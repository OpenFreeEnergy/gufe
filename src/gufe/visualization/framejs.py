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
        Stable short identifier — by convention the snake_case gufe class name,
        e.g. ``"ligand_network"``. Also the suffix of the
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
#
# Naming rule: a viz's `id` (== its `frame` directory under viz_assets/) is the
# snake_case form of the gufe class it renders — `LigandNetwork` ->
# `ligand_network`, `SmallMoleculeComponent` -> `small_molecule_component`, …
# So the on-disk frame is always greppable from the class name and vice versa.
#
# `LigandNetwork`'s `uuid` is the one viz published to framejs.io so far (PLAN §7):
# https://framejs.io/j/<uuid> (== https://framejs.app/j/<uuid>). The rest resolve
# via their on-disk frame (`GUFE_VIZ_SOURCE=auto` -> `local_url()`) until they are
# published; a dev/CI publish can fill each in via GUFE_VIZ_<ID>_UUID with no code
# change.
# No frame declares `modules`: 3Dmol.js is injected on demand by the frames that
# need it (each has a `load3Dmol()` mirroring their `loadRDKit()`), rather than
# being loaded eagerly for every render. That keeps the engine off the critical
# path — `SmallMoleculeComponent.view()` paints its 2D depiction without waiting
# on it, and a `LigandNetwork` only fetches it once an edge is clicked.
CANONICAL_VIZ: dict[str, VizRef] = {
    "LigandNetwork": VizRef(
        id="ligand_network",
        uuid="019f2b55e1f57722af0293acbda78362",
        frame="ligand_network",
    ),
    "AlchemicalNetwork": VizRef(id="alchemical_network", frame="alchemical_network"),
    # Registered on the base so `Transformation` and `NonTransformation` (which are
    # siblings, not parent/child) both resolve to it.
    "TransformationBase": VizRef(id="transformation", frame="transformation"),
    "ChemicalSystem": VizRef(id="chemical_system", frame="chemical_system"),
    "LigandAtomMapping": VizRef(id="ligand_atom_mapping", frame="ligand_atom_mapping"),
    "SmallMoleculeComponent": VizRef(id="small_molecule_component", frame="small_molecule_component"),
    "ProteinComponent": VizRef(id="protein_component", frame="protein_component", default_height="560px"),
    "SolventComponent": VizRef(id="solvent_component", frame="solvent_component", default_height="320px"),
}


# --------------------------------------------------------------------------- #
# Per-object serializers (payload builders): object -> the viz `inputs` dict.   #
# Each returns a JSON-serializable dict whose keys match what the canonical viz #
# reads in `onInputs(inputs)`.                                                  #
# --------------------------------------------------------------------------- #


def _ligand_network_payload(net) -> dict[str, Any]:
    """Serialize a :class:`gufe.LigandNetwork` to the ``ligand_network`` viz input.

    The canonical network viz reads ``inputs['network.graphml']`` and parses the
    GraphML itself (accepting a string, Blob, or ArrayBuffer). We push the exact
    ``LigandNetwork.to_graphml()`` output — the same serialization gufe uses
    everywhere — so the viz sees the full network with no lossy re-encoding.
    """
    return {"network.graphml": net.to_graphml()}


# -- component descriptors (shared by the container vizzes) ----------------- #
#
# ChemicalSystem / Transformation / AlchemicalNetwork all need to describe the
# components they contain. They share one JSON shape, produced here, so the three
# frames can render a component with the same code:
#
#     {"type": <class name>, "name": str, "sdf"|"pdb": str, ...type-specific...}


def _component_descriptor(comp) -> dict[str, Any]:
    """Describe any :class:`gufe.Component` as a plain JSON dict for a viz."""
    kind = type(comp).__name__
    desc: dict[str, Any] = {"type": kind, "name": getattr(comp, "name", "") or ""}
    try:
        if hasattr(comp, "to_sdf"):
            desc["sdf"] = comp.to_sdf()
            desc["smiles"] = comp.smiles
        elif hasattr(comp, "to_pdb_file"):
            desc["pdb"] = _protein_pdb_string(comp)
        elif hasattr(comp, "smiles"):  # SolventComponent & friends
            desc.update(_solvent_fields(comp))
    except Exception as e:  # a bad component must not take down the whole view
        desc["error"] = f"{type(e).__name__}: {e}"
    return desc


def _protein_pdb_string(protein) -> str:
    """Render a :class:`gufe.ProteinComponent` to a PDB string (in memory)."""
    import io

    buf = io.StringIO()
    protein.to_pdb_file(buf)
    return buf.getvalue()


def _json_safe(value: Any) -> Any:
    """Coerce ``value`` into something ``json.dumps`` can handle, stringifying the rest.

    Free-form gufe metadata — ``LigandAtomMapping.annotations`` most of all — can
    hold ``openff.units.Quantity`` and other rich objects. A viz only ever
    *displays* these, so rendering the leftovers with ``str()`` (``"1.2 nanometer"``)
    is both lossless enough and more readable than gufe's ``:custom:`` JSON codec.
    """
    return json.loads(json.dumps(value, default=str))


def _solvent_fields(solvent) -> dict[str, Any]:
    return {
        "smiles": solvent.smiles,
        "positive_ion": solvent.positive_ion,
        "negative_ion": solvent.negative_ion,
        "neutralize": solvent.neutralize,
        "ion_concentration": str(solvent.ion_concentration),
    }


def _small_molecule_component_payload(mol) -> dict[str, Any]:
    """Serialize a :class:`gufe.SmallMoleculeComponent` (2D depiction + 3D)."""
    return {
        "molecule.sdf": mol.to_sdf(),
        "name": mol.name,
        "smiles": mol.smiles,
        "total_charge": mol.total_charge,
    }


def _protein_component_payload(protein) -> dict[str, Any]:
    """Serialize a :class:`gufe.ProteinComponent` as a PDB string for 3Dmol."""
    return {"protein.pdb": _protein_pdb_string(protein), "name": protein.name}


def _solvent_component_payload(solvent) -> dict[str, Any]:
    """Serialize a :class:`gufe.SolventComponent` — a small settings card."""
    return {"solvent": _solvent_fields(solvent)}


def _ligand_atom_mapping_payload(mapping) -> dict[str, Any]:
    """Serialize a :class:`gufe.LigandAtomMapping` for the mapping viewer.

    Both endpoint molecules go over as SDF plus the raw
    ``componentA_to_componentB`` index map, so the viz can colour core / unique
    atoms and draw the correspondence lines itself.
    """
    return {
        "molA.sdf": mapping.componentA.to_sdf(),
        "molB.sdf": mapping.componentB.to_sdf(),
        "nameA": mapping.componentA.name,
        "nameB": mapping.componentB.name,
        "mapping": {str(k): v for k, v in mapping.componentA_to_componentB.items()},
        "annotations": _json_safe(mapping.annotations),
    }


def _chemical_system_payload(system) -> dict[str, Any]:
    """Serialize a :class:`gufe.ChemicalSystem` — its labelled components."""
    return {
        "system": {
            "name": system.name,
            "components": [
                dict(_component_descriptor(comp), label=label) for label, comp in sorted(system.components.items())
            ],
        }
    }


def _transformation_payload(transformation) -> dict[str, Any]:
    """Serialize a :class:`gufe.Transformation` — stateA vs stateB, plus mapping.

    ``NonTransformation`` has the same ``stateA``/``stateB`` properties (both its
    single system), so it renders through this builder unchanged.
    """
    mapping = getattr(transformation, "mapping", None)
    mappings = mapping if isinstance(mapping, list) else ([] if mapping is None else [mapping])
    return {
        "transformation": {
            "name": transformation.name,
            "protocol": type(transformation.protocol).__name__,
            "stateA": _chemical_system_payload(transformation.stateA)["system"],
            "stateB": _chemical_system_payload(transformation.stateB)["system"],
            "mappings": [_ligand_atom_mapping_payload(m) for m in mappings if hasattr(m, "componentA_to_componentB")],
        }
    }


def _alchemical_network_payload(net) -> dict[str, Any]:
    """Serialize a :class:`gufe.AlchemicalNetwork` as a node/edge graph.

    ``AlchemicalNetwork.to_graphml()`` is not implemented, so — unlike
    ``LigandNetwork`` — we build the graph JSON here. Nodes are ``ChemicalSystem``
    summaries keyed by gufe key; edges are ``Transformation`` summaries. Component
    SDF/PDB payloads are *not* inlined (a solvated network would be enormous); the
    viz shows composition and topology.
    """

    def node_id(system):
        return str(system.key)

    return {
        "alchemical_network": {
            "name": net.name,
            "nodes": [
                {
                    "id": node_id(system),
                    "name": system.name,
                    "components": [
                        {"label": label, "type": type(c).__name__, "name": getattr(c, "name", "") or ""}
                        for label, c in sorted(system.components.items())
                    ],
                }
                for system in sorted(net.nodes, key=node_id)
            ],
            "edges": [
                {
                    "id": str(edge.key),
                    "name": edge.name,
                    "source": node_id(edge.stateA),
                    "target": node_id(edge.stateB),
                    "protocol": type(edge.protocol).__name__,
                }
                for edge in sorted(net.edges, key=lambda e: str(e.key))
            ],
        }
    }


# Object-type name -> payload builder. Looked up over the object's MRO (see
# `_registry_lookup`), so subclasses such as `SolvatedPDBComponent` (a
# `ProteinComponent`) and `NonTransformation` (a `Transformation`) inherit the
# viz + serializer of their base class for free.
_PAYLOAD_BUILDERS: dict[str, Callable[[Any], dict[str, Any]]] = {
    "LigandNetwork": _ligand_network_payload,
    "AlchemicalNetwork": _alchemical_network_payload,
    "TransformationBase": _transformation_payload,
    "ChemicalSystem": _chemical_system_payload,
    "LigandAtomMapping": _ligand_atom_mapping_payload,
    "SmallMoleculeComponent": _small_molecule_component_payload,
    "ProteinComponent": _protein_component_payload,
    "SolventComponent": _solvent_component_payload,
}


def _registry_lookup(obj, registry: dict[str, Any]):
    """Look ``obj``'s type up in a registry, walking the MRO for the best match.

    Exact class first, then base classes in MRO order — so a subclass gets its
    parent's viz unless it registers its own. Returns ``None`` on no match.
    """
    for klass in type(obj).__mro__:
        hit = registry.get(klass.__name__)
        if hit is not None:
            return hit
    return None


def _viz_ref_for(obj) -> VizRef:
    viz = _registry_lookup(obj, CANONICAL_VIZ)
    if viz is None:
        raise FramejsUnavailable(f"no canonical framejs viz registered for {type(obj).__name__!r}")
    return viz


def _payload_for(obj) -> dict[str, Any]:
    builder = _registry_lookup(obj, _PAYLOAD_BUILDERS)
    if builder is None:
        raise FramejsUnavailable(f"no framejs payload serializer registered for {type(obj).__name__!r}")
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
