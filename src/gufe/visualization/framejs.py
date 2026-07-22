# This code is part of gufe and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
"""Native framejs.io visualizations for gufe objects.

This is the *infrastructure* layer that lets a gufe object render as an
interactive framejs.io widget. It is intentionally thin: the rendering is done by
:class:`metaframe_widget.MetaframeWidget` (an anywidget that hosts the framejs
iframe in Jupyter / marimo / VSCode and pushes ``inputs`` live over the metapage
comm channel), and the drawing is done by the frames under ``viz_assets/``.

Adding a visualization
----------------------
One entry in :data:`VIZ_REGISTRY` plus one frame directory. The registry maps a
gufe class name to a :class:`VizRef` holding both halves of the contract — which
frame draws it (``frame``) and how the object is serialized for that frame
(``payload``). Lookups walk the MRO, so a subclass inherits its base's viz.

By convention ``frame`` is the snake_case form of the class it renders
(``LigandNetwork`` -> ``ligand_network``), so the frame directory is always
greppable from the class name and vice versa.

Payload keys (the contract with ``onInputs``)
---------------------------------------------
A payload is a **flat** dict. Two kinds of key, and no others:

* **file-shaped values** use a descriptive ``<thing>.<ext>`` key — ``molecule.sdf``,
  ``protein.pdb``, ``network.graphml``, ``molA.sdf``/``molB.sdf``. Frames read these
  through their ``asText()`` helper, which accepts a string, ``Blob``,
  ``ArrayBuffer`` or a dropped file, so these keys double as file-drop targets.
  They are descriptive rather than frame-named because one frame may take several
  (the mapping viewer takes two SDFs).
* **everything else** is a bare snake_case field (``name``, ``smiles``,
  ``total_charge``, ``mapping``, ``annotations``). Where a viz takes a single
  structured object rather than loose fields, that object's key is the **frame
  name** — ``chemical_system``, ``transformation``, ``alchemical_network``,
  ``solvent_component``.

Renaming a payload key is a breaking change for any frame already published to a
``/j/<uuid>`` (that frame is pinned and keeps reading the old key), so treat these
names as the API they are.

Two URL forms
-------------
The viz source lives in the repo as a framejs frame directory under
``viz_assets/<frame>/`` — the on-disk format documented at
https://framejs.io/docs/guide/local-file-io (``code.js`` required; ``og.json`` /
``modules.json`` / ``inputs.json`` / ``options.json`` / ``definition.json``
optional). From that one source gufe builds either of two equivalent URLs, chosen
by ``GUFE_VIZ_SOURCE`` (see :meth:`VizRef.resolve_url`):

* **local** — a self-contained hash-param URL built from the frame directory,
  ``https://framejs.io/#?js=<b64>&og=<b64>…``. Never expires, always matches the
  repo source, needs no account. See :meth:`VizRef.local_url`.
* **canonical** — a pinned short URL ``https://framejs.io/j/<uuid>`` minted by
  publishing the frame directory once (``just publish-viz``). Short and shareable,
  and updatable without re-releasing gufe. See :meth:`VizRef.canonical_url`.

Default is **auto**: local when the frame directory is present, else the pinned
uuid. Either way the per-object data rides on top of the chosen base URL — the
viz JavaScript is never re-derived per object:

* **Notebook** (``.view()`` / bare-cell): the object is pushed as live ``inputs``
  over the comm channel, so there is no URL size limit. See :func:`view_object`.
* **CLI** (a browser, no live Python channel): the object is appended to the base
  URL's hash as ``inputs=<b64>``. Appended ``inputs`` take priority over anything
  baked into the frame. See :func:`build_cli_url`.

Everything here is optional. ``import gufe`` works without ``metaframe-widget``
installed; rendering then falls back to the object's legacy RDKit / py3Dmol
renderer if it has one (see :func:`legacy_view`). Install the extra with
``pip install gufe[viz]``.
"""

from __future__ import annotations

import base64
import json
import os
import urllib.parse
import warnings
from dataclasses import dataclass
from importlib import resources
from typing import TYPE_CHECKING, Any, Callable

if TYPE_CHECKING:
    from metaframe_widget import MetaframeWidget

FRAMEJS_BASE = "https://framejs.io"


class FramejsUnavailable(RuntimeError):
    """Raised when a framejs widget/URL cannot be built (missing dependency, no
    viz registered for the object, no frame on disk and no published uuid, …).

    Callers that want graceful degradation should catch this (plus ``OSError``)
    and fall back to a legacy renderer.
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
    """How one gufe class is visualized: which frame draws it, and how it is fed.

    Parameters
    ----------
    frame
        Name of the frame directory under ``gufe/visualization/viz_assets/``, in
        the framejs local-file-io format (``<frame>/code.js`` plus optional JSON
        sidecars). By convention the snake_case gufe class name. Doubles as the
        viz's short id: it is the suffix of the ``GUFE_VIZ_<FRAME>_UUID`` env
        override, uppercased.
    payload
        Serializer turning the object into the frame's ``inputs`` dict. Its keys
        are the contract with the frame's ``onInputs`` — see the module docstring
        for the naming convention.
    uuid
        The framejs.io canonical frame id, if this viz has been published to
        ``https://framejs.io/j/<uuid>``. ``None`` until it is, in which case only
        the local (on-disk) URL form is available.
    default_height, default_width
        CSS sizes used when the caller does not override them.
    """

    frame: str
    payload: Callable[[Any], dict[str, Any]]
    uuid: str | None = None
    default_height: str = "500px"
    default_width: str = "100%"

    # -- canonical form: the pinned https://framejs.io/j/<uuid> ---------------- #

    def _resolved_uuid(self) -> str | None:
        """The effective uuid, allowing a per-viz env override.

        ``GUFE_VIZ_<FRAME>_UUID`` lets a dev/CI publish swap in a freshly minted
        frame without a code change.
        """
        return os.environ.get(f"GUFE_VIZ_{self.frame.upper()}_UUID") or self.uuid

    @property
    def published(self) -> bool:
        """True if this viz has a canonical framejs.io uuid to point at."""
        return bool(self._resolved_uuid())

    def canonical_url(self) -> str:
        """Return the pinned ``https://framejs.io/j/<uuid>`` URL for this viz.

        Raises :class:`FramejsUnavailable` if the viz has not been published.
        """
        uuid = self._resolved_uuid()
        if not uuid:
            raise FramejsUnavailable(
                f"viz {self.frame!r} has no published framejs.io uuid yet; publish it "
                f"(`just publish-viz`) and set its uuid, or export "
                f"GUFE_VIZ_{self.frame.upper()}_UUID=<uuid>"
            )
        # On the /j/<uuid> route, `inputs` appended in the hash (the CLI path)
        # take priority over the frame's baked-in inputs — what gufe relies on.
        return f"{FRAMEJS_BASE}/j/{uuid}"

    # -- local form: a self-contained URL built from the frame directory ------- #

    def _frame_dir(self):
        """The importable frame directory traversable."""
        return resources.files("gufe.visualization.viz_assets").joinpath(self.frame)

    def has_local(self) -> bool:
        """True if the on-disk frame directory exists and has a ``code.js``."""
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
        Raises :class:`FramejsUnavailable` if the frame directory is missing.
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
        """Return the base viz URL, choosing the form per ``GUFE_VIZ_SOURCE``.

        * ``local``     — always :meth:`local_url` (from the on-disk frame dir).
        * ``canonical`` — always :meth:`canonical_url` (the pinned ``/j/<uuid>``).
        * ``auto`` (default) — the on-disk frame when present, else the uuid.

        Raises :class:`FramejsUnavailable` if the requested form is unavailable.
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
            f"viz {self.frame!r} has neither an on-disk frame directory nor a published "
            f"uuid; add viz_assets/{self.frame}/ or publish it (`just publish-viz`)."
        )


# --------------------------------------------------------------------------- #
# Payload builders: object -> the frame's `inputs` dict.                        #
# Keys follow the convention documented in the module docstring.                #
# --------------------------------------------------------------------------- #


def _ligand_network_payload(net) -> dict[str, Any]:
    """Serialize a :class:`gufe.LigandNetwork` for the ``ligand_network`` frame.

    We push the exact ``LigandNetwork.to_graphml()`` output — the same
    serialization gufe uses everywhere — and the frame parses the GraphML itself,
    so the viz sees the full network with no lossy re-encoding.
    """
    return {"network.graphml": net.to_graphml()}


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
    return {"solvent_component": _solvent_fields(solvent)}


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
    return {"chemical_system": _chemical_system_fields(system)}


def _transformation_payload(transformation) -> dict[str, Any]:
    """Serialize a :class:`gufe.Transformation` — stateA vs stateB, plus mappings.

    ``NonTransformation`` has the same ``stateA``/``stateB`` properties (both its
    single system), so it renders through this builder unchanged.
    """
    mapping = getattr(transformation, "mapping", None)
    mappings = mapping if isinstance(mapping, list) else ([] if mapping is None else [mapping])
    return {
        "transformation": {
            "name": transformation.name,
            "protocol": type(transformation.protocol).__name__,
            "stateA": _chemical_system_fields(transformation.stateA),
            "stateB": _chemical_system_fields(transformation.stateB),
            "mappings": [_ligand_atom_mapping_payload(m) for m in mappings if hasattr(m, "componentA_to_componentB")],
        }
    }


def _alchemical_network_payload(net) -> dict[str, Any]:
    """Serialize a :class:`gufe.AlchemicalNetwork` as a node/edge graph.

    ``AlchemicalNetwork.to_graphml()`` is not implemented, so — unlike
    ``LigandNetwork`` — we build the graph JSON here. Nodes are ``ChemicalSystem``
    summaries keyed by gufe key; edges are ``Transformation`` summaries. Component
    structures are summarized rather than inlined (a solvated network's SDF/PDB
    would be enormous); the viz shows composition and topology.
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
                    "components": [_component_summary(label, c) for label, c in sorted(system.components.items())],
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


# -- shared component shapes ------------------------------------------------- #
#
# ChemicalSystem, Transformation and AlchemicalNetwork all describe the components
# they contain. They share these two shapes so their frames can render a component
# with the same code:
#
#   summary    {"label", "type", "name"}                     — topology only
#   descriptor summary + {"sdf"|"pdb"|solvent fields, ...}   — enough to draw it


def _component_summary(label: str, comp) -> dict[str, Any]:
    """The cheap shape: what a component *is*, with no structure data."""
    return {"label": label, "type": type(comp).__name__, "name": getattr(comp, "name", "") or ""}


def _component_descriptor(label: str, comp) -> dict[str, Any]:
    """The full shape: a summary plus whatever the viz needs to draw the component."""
    desc = _component_summary(label, comp)
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


def _chemical_system_fields(system) -> dict[str, Any]:
    """The inner object shared by the chemical_system and transformation vizzes."""
    return {
        "name": system.name,
        "components": [_component_descriptor(label, comp) for label, comp in sorted(system.components.items())],
    }


def _solvent_fields(solvent) -> dict[str, Any]:
    return {
        "smiles": solvent.smiles,
        "positive_ion": solvent.positive_ion,
        "negative_ion": solvent.negative_ion,
        "neutralize": solvent.neutralize,
        "ion_concentration": str(solvent.ion_concentration),
    }


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


# --------------------------------------------------------------------------- #
# The registry: gufe class name -> its visualization.                           #
# --------------------------------------------------------------------------- #
#
# Looked up over the object's MRO (see `_registry_lookup`), so subclasses such as
# `SolvatedPDBComponent` (a `ProteinComponent`) and `NonTransformation` (a
# `TransformationBase`) inherit their base's viz and serializer for free.
#
# `LigandNetwork` is the one viz published to framejs.io so far; the rest resolve
# via their on-disk frame until they are published, and a dev/CI publish can fill
# each in via GUFE_VIZ_<FRAME>_UUID with no code change.
VIZ_REGISTRY: dict[str, VizRef] = {
    "LigandNetwork": VizRef(
        frame="ligand_network",
        payload=_ligand_network_payload,
        uuid="019f2b55e1f57722af0293acbda78362",
    ),
    "AlchemicalNetwork": VizRef(frame="alchemical_network", payload=_alchemical_network_payload),
    # Registered on the base so `Transformation` and `NonTransformation` (which are
    # siblings, not parent/child) both resolve to it.
    "TransformationBase": VizRef(frame="transformation", payload=_transformation_payload),
    "ChemicalSystem": VizRef(frame="chemical_system", payload=_chemical_system_payload),
    "LigandAtomMapping": VizRef(frame="ligand_atom_mapping", payload=_ligand_atom_mapping_payload),
    "SmallMoleculeComponent": VizRef(frame="small_molecule_component", payload=_small_molecule_component_payload),
    "ProteinComponent": VizRef(
        frame="protein_component", payload=_protein_component_payload, default_height="560px"
    ),
    "SolventComponent": VizRef(
        frame="solvent_component", payload=_solvent_component_payload, default_height="320px"
    ),
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


def _viz_for(obj) -> VizRef:
    """The :class:`VizRef` registered for ``obj``'s class (or a base of it)."""
    viz = _registry_lookup(obj, VIZ_REGISTRY)
    if viz is None:
        raise FramejsUnavailable(f"no framejs viz registered for {type(obj).__name__!r}")
    return viz


# --------------------------------------------------------------------------- #
# Notebook path — base URL + live inputs over the comm channel                  #
# --------------------------------------------------------------------------- #


def _build_widget(obj, *, width: str | None = None, height: str | None = None) -> "MetaframeWidget":
    """Build the notebook widget for ``obj``: load its viz URL, push its data.

    Raises :class:`FramejsUnavailable` if no viz is registered, the widget
    dependency is missing, or no URL form is available.
    """
    viz = _viz_for(obj)
    try:
        from metaframe_widget import MetaframeWidget
    except ImportError as e:
        raise FramejsUnavailable(
            "metaframe-widget is not installed; install the visualization extra "
            "with `pip install gufe[viz]`"
        ) from e

    widget = MetaframeWidget(
        url=viz.resolve_url(),
        width=width or viz.default_width,
        height=height or viz.default_height,
    )
    widget.set_inputs(viz.payload(obj))  # live over the comm channel — never in the URL
    return widget


def view_object(obj, *, width: str | None = None, height: str | None = None, **opts):
    """Return an interactive framejs widget for ``obj`` (notebook path).

    Falls back to the object's legacy renderer (RDKit / py3Dmol) if framejs is
    unavailable for any reason — see :func:`legacy_view`.
    """
    try:
        return _build_widget(obj, width=width, height=height)
    except (FramejsUnavailable, OSError) as e:
        # Explicit call, so be loud about why the interactive view didn't happen.
        warnings.warn(f"framejs view unavailable ({e}); falling back to legacy renderer.", stacklevel=2)
        return legacy_view(obj, **opts)


def repr_mimebundle(obj, include=None, exclude=None):
    """``_repr_mimebundle_`` implementation for auto-display in notebooks.

    Returns the framejs widget's mimebundle; failing that the legacy renderer's,
    so an install without the ``viz`` extra keeps the pre-framejs picture; failing
    that ``None``, so the notebook prints the object's plain ``repr``.

    Unlike :func:`view_object` this is silent — it runs on every bare-cell display
    of the object, and a warning per render would be intolerable.
    """
    try:
        widget = _build_widget(obj)
    except (FramejsUnavailable, OSError):
        return _legacy_mimebundle(obj, include=include, exclude=exclude)
    return widget._repr_mimebundle_(include=include, exclude=exclude)


# --------------------------------------------------------------------------- #
# Legacy fallback — the object's pre-framejs RDKit / py3Dmol renderer           #
# --------------------------------------------------------------------------- #


def legacy_view(obj, **opts):
    """Return the object's pre-framejs rendering, or ``None`` if it has none.

    A class opts in by defining ``_legacy_view()``; absence means no fallback
    exists (``LigandNetwork``, for instance, never had a static renderer).
    """
    fn = getattr(obj, "_legacy_view", None)
    if not callable(fn):
        return None
    try:
        return fn(**opts)
    except Exception as e:  # a broken fallback must not mask the original problem
        warnings.warn(f"legacy renderer for {type(obj).__name__} failed: {e}", stacklevel=2)
        return None


# Rich-display hooks an IPython display object may expose, in preference order.
_LEGACY_REPR_METHODS = (
    ("_repr_html_", "text/html"),
    ("_repr_svg_", "image/svg+xml"),
    ("_repr_png_", "image/png"),
    ("_repr_jpeg_", "image/jpeg"),
)


def _legacy_mimebundle(obj, include=None, exclude=None):
    """Adapt :func:`legacy_view`'s result into a mimebundle, or ``None``.

    The legacy renderers return IPython display objects (``Image``, ``HTML``, …),
    which carry their own ``_repr_*_`` hooks; we read those directly rather than
    booting an IPython formatter, so this works outside IPython too.

    The result is whatever mimebundle shape the display object produces — IPython
    accepts either a ``{mime: data}`` dict or a ``(data, metadata)`` pair, and
    ``IPython.display.Image`` returns the latter.
    """
    legacy = legacy_view(obj)
    if legacy is None:
        return None
    if hasattr(legacy, "_repr_mimebundle_"):
        return legacy._repr_mimebundle_(include=include, exclude=exclude)
    for attr, mime in _LEGACY_REPR_METHODS:
        fn = getattr(legacy, attr, None)
        if callable(fn):
            data = fn()
            if data is not None:
                return {mime: data, "text/plain": repr(obj)}
    return None


# --------------------------------------------------------------------------- #
# CLI path — base URL + inputs appended in the hash                             #
# --------------------------------------------------------------------------- #


def build_cli_url(obj) -> str:
    """Return a framejs.io URL that renders ``obj`` (for ``webbrowser.open``).

    The base viz URL is chosen by :meth:`VizRef.resolve_url` (the local hash-param
    URL built from the on-disk frame by default, else the pinned ``/j/<uuid>``),
    and this object's ``inputs`` are merged into that URL's hash. Appended
    ``inputs`` take priority over anything baked into the frame.
    """
    viz = _viz_for(obj)
    base = viz.resolve_url()
    encoded = string_to_base64_string(json.dumps(viz.payload(obj)))
    # The local URL already carries a `#?js=…` hash — merge inputs into it with
    # `&`; the canonical `/j/<uuid>` URL has no hash yet, so start one with `#?`.
    sep = "&" if "#?" in base else "#?"
    return f"{base}{sep}inputs={encoded}"
