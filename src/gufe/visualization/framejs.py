# This code is part of gufe and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
"""Native framejs.io visualizations for gufe objects.

This is the *infrastructure* layer that lets any gufe object render as an
interactive framejs.io widget. It is intentionally thin: the heavy lifting is
done by

* the existing :class:`metaframe_widget.MetaframeWidget` anywidget (renders the
  framejs iframe; works in Jupyter, marimo and VSCode; pushes ``inputs`` live
  over the metapage comm channel), and
* the per-object *canonical visualizations* — small, versioned framejs apps
  **published once to framejs.io** as canonical frames.

On-disk frames + two URL forms (framejs.io, 2026-07)
----------------------------------------------------
The viz *source* lives in the repo as a **canonical framejs frame directory**
under ``viz_assets/<frame>/`` — the on-disk format documented at
https://framejs.io/docs/guide/local-file-io :

* ``code.js``     — the ``js`` hash param (the viz JavaScript)   [required]
* ``og.json``     — the ``og`` hash param (title/description)    [optional]
* ``modules.json``— the ``modules`` hash param (classic scripts) [optional]
* ``inputs.json`` / ``options.json`` / ``definition.json``       [optional]

From that single on-disk source gufe can build **two equivalent URL forms**,
and the runtime chooses between them (see :meth:`VizRef.resolve_url` and the
``GUFE_VIZ_SOURCE`` env var):

* **local** (default) — a full hash-param URL derived directly from the frame
  directory: ``https://framejs.io/#?js=<b64>&og=<b64>…``. Self-contained, never
  expires, always matches the repo source, needs no account. See
  :meth:`VizRef.local_url`.
* **canonical** — a pinned short URL ``https://framejs.io/j/<uuid>`` minted by
  *publishing* the frame directory once (the ``publish-viz`` recipe). Short and
  shareable; updatable without re-releasing gufe. See :meth:`VizRef.canonical_url`.

The default is **auto**: use the local URL when the frame dir is present on
disk, else fall back to the pinned ``uuid``. The per-object data is supplied
*per render* on top of whichever base URL is chosen — the viz JS is never
re-derived per object:

* **Notebook** (``.view()`` / bare-cell auto-display): the widget loads the base
  URL and the serialized object is pushed as live ``inputs`` over the comm
  channel — no URL size limit. See :func:`view_object`.
* **CLI** (browser, no live Python channel): the serialized object is appended to
  the base URL's hash as ``inputs=<b64>``. Appended ``inputs`` **always take
  priority** over anything baked into the frame, so one shared viz renders each
  object's data correctly. See :func:`build_cli_url`.

Round-trip / live-edit the on-disk frames with the ``framejs_frame.py`` helper
and the framejs local server (``just viz-url`` / ``just viz-save`` /
``just viz-edit``) — see ``visualization-demo/scripts/framejs_frame.py``.

Everything here is **non-destructive and optional**. ``import gufe`` works
without ``metaframe-widget`` installed; ``.view()`` then falls back to the
legacy RDKit / py3Dmol renderers (see :func:`legacy_view`). Install the extra
with ``pip install gufe[viz]``.
"""

from __future__ import annotations

import base64
import json
import os
import urllib.error
import urllib.parse
import warnings
from dataclasses import dataclass
from importlib import resources
from typing import TYPE_CHECKING, Any, Callable

if TYPE_CHECKING:
    from metaframe_widget import MetaframeWidget

FRAMEJS_BASE = "https://framejs.io"

# Sentinel for a viz whose framejs.io frame has not been published yet. A VizRef
# still carrying this (and with no env override) raises FramejsUnavailable from
# `canonical_url()`, so `.view()` falls back cleanly and the CLI errors friendly.
UNPUBLISHED = "REPLACE_WITH_PUBLISHED_UUID"

# CLI-only: encoded payloads above this many bytes should be uploaded + referenced
# rather than inlined in the URL hash. Exact tuning is deferred (PLAN §3/§8); for
# now the CLI always hash-appends and warns past the threshold.
CLI_HASH_INLINE_MAX_BYTES = 16 * 1024


class FramejsUnavailable(RuntimeError):
    """Raised when a framejs widget/URL cannot be built (missing dep, no viz,
    unpublished frame, …).

    Callers that want graceful degradation should catch this (plus
    ``ImportError`` / ``OSError``) and fall back to a legacy renderer.
    """


# Equivalent to JavaScript ``btoa(encodeURIComponent(value))`` — the canonical
# framejs hash-param encoder (mirrors ``metaframe_widget.string_to_base64_string``;
# duplicated so this module works even if that symbol moves).
def string_to_base64_string(value: str) -> str:
    """Encode a string the way framejs.io expects in hash params."""
    encoded_uri_component = urllib.parse.quote(value, safe="-_.!~*'()")
    return base64.b64encode(encoded_uri_component.encode("ascii")).decode("ascii")


@dataclass(frozen=True)
class VizRef:
    """A pinned reference to a canonical framejs.io visualization.

    Parameters
    ----------
    id
        Stable short identifier, e.g. ``"network"``. Also the suffix of the
        per-viz env override ``GUFE_VIZ_<ID>_UUID`` (uppercased).
    uuid
        The framejs.io canonical frame id: the viz is published once to
        ``https://framejs.io/j/<uuid>``. Defaults to :data:`UNPUBLISHED` until a
        real frame exists.
    frame
        Name of the on-disk **frame directory** under
        ``gufe/visualization/viz_assets/`` (canonical framejs local-file-io
        format: ``<frame>/code.js`` plus optional ``og.json`` / ``modules.json``
        / ``inputs.json``). This is both the authoring source that gets
        **published** to framejs.io (``publish-viz``) and — unlike the old bare
        ``.js`` asset — the source the runtime reads to build a self-contained
        local URL (:meth:`local_url`).
    version
        Optional immutable ``?v=<sha256>`` pin (deferred — see PLAN §7). When set,
        :meth:`canonical_url` appends it.
    modules
        Classic (non-ESM) script URLs the viz needs (e.g. 3Dmol.js). Baked into
        the frame at publish time via the framejs ``modules`` mechanism; recorded
        here so the publish tooling can pass them. Not part of the runtime URL.
    default_height, default_width
        CSS sizes used when the caller does not override them.
    """

    id: str
    uuid: str = UNPUBLISHED
    frame: str | None = None
    version: str | None = None
    modules: tuple[str, ...] = ()
    default_height: str = "500px"
    default_width: str = "100%"

    def _resolved_uuid(self) -> str:
        """The effective uuid, allowing a per-viz env override.

        ``GUFE_VIZ_<ID>_UUID`` lets a dev/self-hosted publish swap in a freshly
        minted frame without a code change (e.g. in CI or the dev container).
        """
        return os.environ.get(f"GUFE_VIZ_{self.id.upper()}_UUID") or self.uuid

    @property
    def published(self) -> bool:
        uuid = self._resolved_uuid()
        return bool(uuid) and uuid != UNPUBLISHED

    def canonical_url(self) -> str:
        """Return the pinned canonical framejs.io URL for this viz.

        ``https://framejs.io/j/<uuid>`` (plus ``?v=<version>`` when a version is
        pinned). Raises :class:`FramejsUnavailable` if the frame is unpublished.
        """
        if not self.published:
            raise FramejsUnavailable(
                f"canonical viz {self.id!r} has no published framejs.io uuid yet; "
                f"publish it (see the `publish-viz` recipe) and set its uuid, or "
                f"export GUFE_VIZ_{self.id.upper()}_UUID=<uuid>"
            )
        # Canonical route: framejs.io/j/<uuid> (the live, deployed route). On
        # this route, `inputs` appended in the hash (CLI) take priority over the
        # frame's baked-in inputs — the behavior gufe relies on.
        url = f"{FRAMEJS_BASE}/j/{self._resolved_uuid()}"
        return f"{url}?v={self.version}" if self.version else url

    # -- on-disk frame directory (canonical framejs local-file-io format) ---- #

    def _frame_dir(self):
        """The importable frame directory traversable, or raise if not recorded."""
        if not self.frame:
            raise FramejsUnavailable(f"viz {self.id!r} has no on-disk frame recorded")
        return resources.files("gufe.visualization.viz_assets").joinpath(self.frame)

    def has_local(self) -> bool:
        """True if the on-disk frame directory exists and has a ``code.js``."""
        if not self.frame:
            return False
        try:
            return self._frame_dir().joinpath("code.js").is_file()
        except (FileNotFoundError, ModuleNotFoundError, OSError):
            return False

    def js_source(self) -> str:
        """Return the viz JavaScript (``<frame>/code.js``) from package data."""
        try:
            return self._frame_dir().joinpath("code.js").read_text(encoding="utf-8")
        except (FileNotFoundError, ModuleNotFoundError, OSError) as e:
            raise FramejsUnavailable(f"viz frame {self.frame!r} has no readable code.js: {e}") from e

    def _frame_json(self, name: str):
        """Read a JSON sidecar (``og.json``/``modules.json``/…) or return None."""
        try:
            return json.loads(self._frame_dir().joinpath(name).read_text(encoding="utf-8"))
        except (FileNotFoundError, ModuleNotFoundError, OSError):
            return None

    def local_url(self) -> str:
        """Return a self-contained framejs.io URL built from the on-disk frame.

        Encodes ``code.js`` and any JSON sidecars into hash params exactly as the
        framejs runtime expects — ``https://framejs.io/#?js=<b64>&og=<b64>…`` —
        where each value is ``base64(encodeURIComponent(text))`` (raw text for
        ``js``, compact JSON for the rest). This mirrors ``framejs_frame.py
        restore`` and the framejs local server, so the URL round-trips with them.
        Raises :class:`FramejsUnavailable` if the frame dir is missing.
        """
        parts = [f"js={string_to_base64_string(self.js_source())}"]
        # Order mirrors the framejs local server's frameToUrl (hash-query keys).
        for key in ("options", "inputs", "definition", "og", "modules"):
            value = self._frame_json(f"{key}.json")
            if value is None:
                continue
            encoded = string_to_base64_string(json.dumps(value, separators=(",", ":")))
            parts.append(f"{key}={encoded}")
        return f"{FRAMEJS_BASE}/#?" + "&".join(parts)

    def resolve_url(self) -> str:
        """Return the base viz URL, choosing the source per ``GUFE_VIZ_SOURCE``.

        * ``local``     — always :meth:`local_url` (from the on-disk frame dir).
        * ``canonical`` — always :meth:`canonical_url` (the pinned ``/j/<uuid>``).
        * ``auto`` (default) — prefer the on-disk frame when present, else fall
          back to the pinned ``uuid``.

        Raises :class:`FramejsUnavailable` if the requested source is unavailable.
        """
        source = os.environ.get("GUFE_VIZ_SOURCE", "auto").strip().lower()
        if source == "canonical":
            return self.canonical_url()
        if source == "local":
            return self.local_url()
        # auto
        if self.has_local():
            return self.local_url()
        if self.published:
            return self.canonical_url()
        raise FramejsUnavailable(
            f"viz {self.id!r} has neither an on-disk frame ({self.frame!r}) nor a "
            f"published uuid; add its frame dir under viz_assets/ or publish it "
            f"(see the `publish-viz` recipe)."
        )


# Object-type name -> pinned canonical visualization.
# v1 ships LigandNetwork end-to-end. Its `uuid` is the network viz published to
# framejs.io under a claimed account (PLAN §7): https://framejs.io/j/<uuid>
# (== https://framejs.app/j/<uuid>). A dev/CI publish can override it via
# GUFE_VIZ_NETWORK_UUID without a code change.
CANONICAL_VIZ: dict[str, VizRef] = {
    "LigandNetwork": VizRef(id="network", uuid="019f2b55e1f57722af0293acbda78362", frame="network"),
    # Phase 2 (not yet implemented):
    # "ProteinComponent":       VizRef(id="protein",  frame="protein",
    #                                  modules=("https://3dmol.org/build/3Dmol-min.js",)),
    # "SmallMoleculeComponent": VizRef(id="ligand2d", frame="ligand2d"),
    # "LigandAtomMapping":      VizRef(id="mapping2d", frame="mapping2d"),
}


# --------------------------------------------------------------------------- #
# Per-object serializers (payload builders): object -> the viz `inputs` dict.   #
# Each returns a JSON-serializable dict whose keys match what the canonical viz #
# reads in `onInputs(inputs)`.                                                  #
# --------------------------------------------------------------------------- #


def _network_payload(net) -> dict[str, Any]:
    """Serialize a :class:`gufe.LigandNetwork` to the ``network`` viz input.

    The canonical network viz reads ``inputs['network.graphml']`` and parses the
    GraphML itself (accepting a string, Blob, or ArrayBuffer). We push the exact
    ``LigandNetwork.to_graphml()`` output — the same serialization gufe uses
    everywhere — so the viz sees the full network with no lossy re-encoding.
    """
    return {"network.graphml": net.to_graphml()}


# Object-type name -> payload builder.
_PAYLOAD_BUILDERS: dict[str, Callable[[Any], dict[str, Any]]] = {
    "LigandNetwork": _network_payload,
}


def _viz_ref_for(obj) -> VizRef:
    name = type(obj).__name__
    viz = CANONICAL_VIZ.get(name)
    if viz is None:
        raise FramejsUnavailable(f"no canonical framejs viz registered for {name!r}")
    return viz


def _payload_for(obj) -> dict[str, Any]:
    name = type(obj).__name__
    builder = _PAYLOAD_BUILDERS.get(name)
    if builder is None:
        raise FramejsUnavailable(f"no framejs payload serializer registered for {name!r}")
    return builder(obj)


# --------------------------------------------------------------------------- #
# Notebook path — canonical URL + live inputs over the comm channel             #
# --------------------------------------------------------------------------- #


def _viz_widget(viz: VizRef, payload: dict[str, Any], *, width: str, height: str) -> "MetaframeWidget":
    """Build the notebook widget: load the canonical URL, push data as inputs."""
    try:
        from metaframe_widget import MetaframeWidget
    except ImportError as e:
        raise FramejsUnavailable(
            "metaframe-widget is not installed; install the visualization extra "
            "with `pip install gufe[viz]`"
        ) from e

    w = MetaframeWidget(url=viz.resolve_url(), width=width, height=height)
    w.set_inputs(payload)  # pushed live over the comm channel — never in the URL
    return w


def view_object(obj, *, width: str | None = None, height: str | None = None, **opts):
    """Return an interactive framejs widget for ``obj`` (notebook path).

    Falls back to the object's legacy renderer (RDKit / py3Dmol) if framejs is
    unavailable for any reason — see :func:`legacy_view`.
    """
    try:
        viz = _viz_ref_for(obj)
        payload = _payload_for(obj)
        return _viz_widget(
            viz,
            payload,
            width=width or viz.default_width,
            height=height or viz.default_height,
        )
    except (FramejsUnavailable, ImportError, OSError) as e:
        warnings.warn(f"framejs view unavailable ({e}); falling back to legacy renderer.", stacklevel=2)
        return legacy_view(obj, **opts)


def repr_mimebundle(obj, include=None, exclude=None):
    """``_repr_mimebundle_`` implementation for auto-display in notebooks.

    Returns the framejs widget's mimebundle, or ``None`` so the notebook falls
    back to the object's plain ``repr`` when no viz / widget is available.
    """
    try:
        viz = _viz_ref_for(obj)
        payload = _payload_for(obj)
        widget = _viz_widget(viz, payload, width=viz.default_width, height=viz.default_height)
    except (FramejsUnavailable, ImportError, OSError):
        return None
    return widget._repr_mimebundle_(include=include, exclude=exclude)


def legacy_view(obj, **opts):
    """Fall back to the object's pre-framejs renderer, if it has one.

    Returns whatever the legacy renderer produces, or ``None`` if the object has
    no legacy view (e.g. ``LigandNetwork``, whose only interactive view is the
    new framejs one).
    """
    # Each object opts in by providing `_legacy_view`; absence => no fallback.
    fn = getattr(obj, "_legacy_view", None)
    if callable(fn):
        return fn(**opts)
    return None


# --------------------------------------------------------------------------- #
# CLI path — canonical URL + inputs appended in the hash                        #
# --------------------------------------------------------------------------- #


def build_cli_url(obj, *, timeout: float = 15.0) -> str:
    """Return a framejs.io URL that renders ``obj`` (for ``webbrowser.open``).

    The base viz URL is chosen by :meth:`VizRef.resolve_url` (local hash-param
    URL from the on-disk frame by default, else the pinned ``/j/<uuid>``), and the
    per-object ``inputs`` are merged into that URL's hash. Appended ``inputs``
    always take priority over anything baked into the frame, so the shared viz
    renders this object's data.

    ``timeout`` is accepted for API compatibility with the (deferred) large-payload
    upload path; the hash-inline path makes no network call.
    """
    viz = _viz_ref_for(obj)
    payload = _payload_for(obj)
    base = viz.resolve_url()

    encoded = string_to_base64_string(json.dumps(payload))
    if len(encoded) > CLI_HASH_INLINE_MAX_BYTES:
        # TODO(PLAN §3 large-payload): PUT to /api/upload/presign (/f/<sha256>) and
        # reference as a {type:"url"} input instead of inlining. For now we still
        # inline and warn — fine for LigandNetwork; matters for proteins (Phase 2).
        warnings.warn(
            f"framejs inputs are large ({len(encoded)} B encoded); the URL may exceed "
            f"browser limits. Upload+reference transport is not implemented yet.",
            stacklevel=2,
        )
    # The local URL already carries a `#?js=…` hash — merge inputs into it with
    # `&`; the canonical `/j/<uuid>` URL has no hash yet, so start one with `#?`.
    sep = "&" if "#?" in base else "#?"
    return f"{base}{sep}inputs={encoded}"
