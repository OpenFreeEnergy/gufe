# This code is part of gufe and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
"""Tests for the framejs.io visualization infrastructure (gufe.visualization.framejs).

These cover the *infrastructure* (registry, serializers, canonical-URL building,
display hooks, fallback) and deliberately do not assert on the canonical viz
JavaScript itself. They are fully **hermetic** — no network: the canonical viz
``uuid`` is swapped into the registry by fixture so we never need a
live published frame. Tests that need the optional ``metaframe-widget`` dependency
skip cleanly when it is not installed (it lives behind the ``gufe[viz]`` extra).
"""

import base64
import importlib
import json
import urllib.parse
import warnings
from dataclasses import replace

import pytest
from openff.units import unit
from rdkit import Chem

from gufe import LigandAtomMapping, LigandNetwork, SmallMoleculeComponent
from gufe.visualization import framejs

# A fake published frame id, swapped into the registry so canonical_url() resolves
# without any live framejs.io frame.
FAKE_UUID = "0192f0a1-0000-7000-8000-deadbeef0001"


def mol_from_smiles(smi):
    m = Chem.MolFromSmiles(smi)
    m.Compute2DCoords()
    return m


@pytest.fixture
def simple_network():
    mol1 = SmallMoleculeComponent(mol_from_smiles("CCO"), name="ethanol")
    mol2 = SmallMoleculeComponent(mol_from_smiles("CC"), name="ethane")
    mol3 = SmallMoleculeComponent(mol_from_smiles("CO"), name="methanol")
    edges = [
        LigandAtomMapping(mol1, mol2, {0: 0, 1: 1}, {"score": 0.0, "length": 1.0 * unit.angstrom}),
        LigandAtomMapping(mol2, mol3, {0: 0}, {"score": 1.0}),
        LigandAtomMapping(mol1, mol3, {0: 0, 2: 1}, {"score": 0.5}),
    ]
    return LigandNetwork(edges)


@pytest.fixture
def published_network(monkeypatch):
    """Registered LigandNetwork viz pinned to FAKE_UUID, with no frame directory
    on disk — the one situation where the canonical ``/j/<uuid>`` form is what
    `resolve_url()` returns. Hermetic: nothing is fetched."""
    viz = framejs.VIZ_REGISTRY["LigandNetwork"]
    monkeypatch.setitem(
        framejs.VIZ_REGISTRY,
        "LigandNetwork",
        replace(viz, frame="does_not_exist", uuid=FAKE_UUID),
    )
    return FAKE_UUID


@pytest.fixture
def short_url_network(monkeypatch):
    """Registered LigandNetwork viz pinned to FAKE_UUID, frame dir intact — so
    `build_cli_url(short=True)` has a uuid to build on while the default path
    still resolves locally."""
    viz = framejs.VIZ_REGISTRY["LigandNetwork"]
    monkeypatch.setitem(framejs.VIZ_REGISTRY, "LigandNetwork", replace(viz, uuid=FAKE_UUID))
    return FAKE_UUID


@pytest.fixture
def unpublished_network(monkeypatch):
    """Force the registered LigandNetwork viz to be fully unavailable: no pinned
    uuid AND no on-disk frame, so URL building must raise / fall back."""
    monkeypatch.setitem(framejs.VIZ_REGISTRY, "LigandNetwork", framejs.VizRef(frame="does_not_exist", payload=lambda o: {}))


# --------------------------------------------------------------------------- #
# Encoding                                                                      #
# --------------------------------------------------------------------------- #


def test_string_to_base64_string_roundtrip():
    value = "export function onInputs(i){ /* á & ? = # */ }"
    encoded = framejs.string_to_base64_string(value)
    # btoa(encodeURIComponent(value)) == base64(quote(value))
    decoded = base64.b64decode(encoded).decode("ascii")
    assert "%" in decoded  # percent-encoded
    assert "export function onInputs" in importlib.import_module("urllib.parse").unquote(decoded)


# --------------------------------------------------------------------------- #
# Registry / VizRef                                                             #
# --------------------------------------------------------------------------- #


def test_ligandnetwork_is_registered():
    viz = framejs.VIZ_REGISTRY["LigandNetwork"]
    assert viz.frame == "ligand_network"
    assert viz.has_local()  # the on-disk frame dir viz_assets/ligand_network/ ships with gufe


def test_ligandnetwork_is_published():
    # the registered viz carries the published framejs.io uuid
    viz = framejs.VIZ_REGISTRY["LigandNetwork"]
    assert viz.published is True
    assert viz.canonical_url() == f"https://framejs.io/j/{viz.uuid}"


def test_vizref_without_uuid_is_unpublished():
    viz = framejs.VizRef(frame="x", payload=lambda o: {})
    assert viz.published is False
    with pytest.raises(framejs.FramejsUnavailable):
        viz.canonical_url()


def test_vizref_canonical_url_explicit_uuid():
    viz = framejs.VizRef(frame="x", payload=lambda o: {}, uuid="abc123")
    assert viz.canonical_url() == "https://framejs.io/j/abc123"


def test_vizref_js_source_loads_from_frame_dir():
    viz = framejs.VIZ_REGISTRY["LigandNetwork"]
    src = viz.js_source()  # reads viz_assets/ligand_network/code.js
    assert "onInputs" in src  # the framejs render entrypoint (may be `async`)
    assert len(src) > 100


def test_vizref_missing_frame_raises():
    bad = framejs.VizRef(frame="does_not_exist", payload=lambda o: {}, uuid="abc")
    assert bad.has_local() is False
    with pytest.raises(framejs.FramejsUnavailable):
        bad.js_source()


def test_vizref_no_frame_has_no_local():
    assert framejs.VizRef(frame="does_not_exist", payload=lambda o: {}, uuid="abc").has_local() is False


# --------------------------------------------------------------------------- #
# Local URL (built from the on-disk frame directory) + source selection        #
# --------------------------------------------------------------------------- #


def test_vizref_local_url_encodes_frame():
    viz = framejs.VIZ_REGISTRY["LigandNetwork"]
    url = viz.local_url()
    assert url.startswith("https://framejs.io/#?js=")
    # the js hash param decodes back to the on-disk code.js (btoa(encodeURIComponent(js)))
    js_param = urllib.parse.parse_qs(url.split("#?", 1)[1])["js"][0]
    decoded = urllib.parse.unquote(base64.b64decode(js_param).decode("ascii"))
    assert decoded == viz.js_source()
    assert "onInputs" in decoded
    # og.json sidecar is carried as an og hash param
    assert "&og=" in url
    # No frame ships a modules.json: 3Dmol.js is injected on demand by the frame
    # itself (load3Dmol()) so it stays off the initial-render critical path.
    assert "&modules=" not in url
    assert "load3Dmol" in decoded


def test_resolve_url_prefers_local_over_uuid():
    """The shipped frame wins over the pinned uuid: it always matches the
    installed code and cannot expire."""
    viz = framejs.VIZ_REGISTRY["LigandNetwork"]
    assert viz.published and viz.has_local()
    assert viz.resolve_url() == viz.local_url()


def test_resolve_url_falls_back_to_uuid_without_a_frame_dir(published_network):
    viz = framejs.VIZ_REGISTRY["LigandNetwork"]
    assert not viz.has_local()
    assert viz.resolve_url() == f"https://framejs.io/j/{FAKE_UUID}"


def test_build_cli_url_local_default(simple_network):
    # AUTO default builds a self-contained hash-param URL from the on-disk frame,
    # merging the per-object inputs into the SAME hash with `&`.
    url = framejs.build_cli_url(simple_network)
    assert url.startswith("https://framejs.io/#?js=")
    assert "&inputs=" in url
    encoded = url.split("&inputs=", 1)[1]
    decoded = json.loads(urllib.parse.unquote(base64.b64decode(encoded).decode("ascii")))
    assert "network.graphml" in decoded
    assert "<graphml" in decoded["network.graphml"]


def test_unregistered_type_raises():
    with pytest.raises(framejs.FramejsUnavailable):
        framejs._viz_for(object())
    with pytest.raises(framejs.FramejsUnavailable):
        framejs._viz_for(object())


# --------------------------------------------------------------------------- #
# Serializer (payload builder)                                                  #
# --------------------------------------------------------------------------- #


def test_network_payload_is_graphml(simple_network):
    # the canonical viz reads inputs['network.graphml'] and parses GraphML itself
    payload = framejs._viz_for(simple_network).payload(simple_network)
    assert set(payload) == {"network.graphml"}
    xml = payload["network.graphml"]
    assert isinstance(xml, str)
    assert xml == simple_network.to_graphml()  # the exact gufe serialization
    assert "<graphml" in xml


def test_network_payload_is_json_serializable(simple_network):
    payload = framejs._viz_for(simple_network).payload(simple_network)
    # must round-trip cleanly (it travels as JSON to the browser)
    assert json.loads(json.dumps(payload)) == payload


def test_network_payload_stable(simple_network):
    # reproducible: same network -> identical payload
    assert framejs._viz_for(simple_network).payload(simple_network) == framejs._viz_for(simple_network).payload(simple_network)


# --------------------------------------------------------------------------- #
# CLI URL building — canonical URL + appended inputs                            #
# --------------------------------------------------------------------------- #


def test_build_cli_url_appends_inputs_to_canonical_url(simple_network, published_network):
    url = framejs.build_cli_url(simple_network)
    assert url.startswith(f"https://framejs.io/j/{FAKE_UUID}#?inputs=")
    # the appended inputs decode back to the payload (btoa(encodeURIComponent(json)))
    encoded = url.split("#?inputs=", 1)[1]
    decoded = json.loads(urllib.parse.unquote(base64.b64decode(encoded).decode("ascii")))
    assert "network.graphml" in decoded
    assert "<graphml" in decoded["network.graphml"]
    # the viz js is NOT in the URL (it is hosted at the canonical frame)
    assert "js=" not in url


def test_build_cli_url_unpublished_raises(simple_network, unpublished_network):
    # an unpublished registered viz (no uuid) makes the CLI URL unbuildable
    with pytest.raises(framejs.FramejsUnavailable):
        framejs.build_cli_url(simple_network)


def test_build_cli_url_short_uses_the_uuid(simple_network, short_url_network):
    """`short=True` opts into the pinned /j/<uuid> even though the frame dir is
    present — the only reason to do so being URL size."""
    url = framejs.build_cli_url(simple_network, short=True)
    assert url.startswith(f"https://framejs.io/j/{FAKE_UUID}#?inputs=")
    assert "js=" not in url  # the viz JavaScript is hosted, not inlined


def test_build_cli_url_short_is_dramatically_smaller(simple_network, short_url_network):
    long_url = framejs.build_cli_url(simple_network)
    short = framejs.build_cli_url(simple_network, short=True)
    assert long_url.startswith("https://framejs.io/#?js=")  # default stays local
    assert len(short) * 4 < len(long_url)  # inlining code.js is what makes it big


def test_build_cli_url_short_raises_when_unpublished(simple_network):
    """The shipped LigandAtomMapping viz has no uuid, so there is no short form."""
    edge = next(iter(simple_network.edges))
    assert framejs._viz_for(edge).published is False
    framejs.build_cli_url(edge)  # the default (local) form still works
    with pytest.raises(framejs.FramejsUnavailable):
        framejs.build_cli_url(edge, short=True)


def test_module_reads_no_environment_variables():
    """Behaviour comes from the registry and the caller's arguments, not the env.

    Guards the deliberate removal of GUFE_VIZ_SOURCE / GUFE_VIZ_<FRAME>_UUID: both
    were reachable only in combination, and drifted out of sync with the docs.
    """
    import inspect

    src = inspect.getsource(framejs)
    code = "\n".join(l for l in src.splitlines() if not l.strip().startswith("#"))
    assert "os.environ" not in code
    assert "getenv" not in code


# --------------------------------------------------------------------------- #
# Display hooks (need metaframe-widget)                                         #
# --------------------------------------------------------------------------- #

metaframe_widget = pytest.importorskip("metaframe_widget")


def test_view_returns_widget_with_canonical_url_and_inputs(simple_network, published_network):
    w = simple_network.view()
    assert type(w).__name__ == "MetaframeWidget"
    # notebook path: canonical URL loaded, data pushed as inputs (NOT in the URL)
    assert w.url == f"https://framejs.io/j/{FAKE_UUID}"
    assert "network.graphml" in w.inputs
    assert "inputs=" not in w.url


def test_view_respects_size_overrides(simple_network, published_network):
    w = simple_network.view(height="123px", width="456px")
    assert w.height == "123px"
    assert w.width == "456px"


def test_repr_mimebundle_returns_bundle(simple_network, published_network):
    mb = simple_network._repr_mimebundle_()
    assert mb  # non-empty dict/tuple so the notebook auto-displays the widget


def test_repr_mimebundle_none_for_unregistered():
    # a gufe object with no registered viz must not break bare-cell display
    assert framejs.repr_mimebundle(object()) is None


def test_view_falls_back_when_unpublished(simple_network, unpublished_network):
    # no published uuid -> FramejsUnavailable inside view_object -> legacy (None here)
    with pytest.warns(UserWarning):
        assert simple_network.view() is None


def test_view_falls_back_when_viz_unavailable(simple_network, published_network, monkeypatch):
    # simulate framejs being unavailable -> falls back to legacy_view (None here)
    monkeypatch.setattr(
        framejs,
        "_build_widget",
        lambda *a, **k: (_ for _ in ()).throw(framejs.FramejsUnavailable("boom")),
    )
    with pytest.warns(UserWarning):
        assert simple_network.view() is None


# --------------------------------------------------------------------------- #
# Registry coverage for every visualizable gufe type                            #
# --------------------------------------------------------------------------- #
#
# The rest of this module exercises the infrastructure through LigandNetwork.
# These tests assert the *breadth* of the registry instead: every type that mixes
# in FramejsViewable must have both a viz and a serializer, its frame must ship
# on disk, and its payload must survive the JSON transport both hops use.

# Every gufe class that opts in to a framejs view, and the frame it must resolve to.
EXPECTED_FRAMES = {
    "LigandNetwork": "ligand_network",
    "AlchemicalNetwork": "alchemical_network",
    "TransformationBase": "transformation",
    "ChemicalSystem": "chemical_system",
    "LigandAtomMapping": "ligand_atom_mapping",
    "SmallMoleculeComponent": "small_molecule_component",
    "ProteinComponent": "protein_component",
    "SolventComponent": "solvent_component",
}


@pytest.mark.parametrize("cls_name,frame", sorted(EXPECTED_FRAMES.items()))
def test_every_registered_viz_ships_its_frame(cls_name, frame):
    """Naming rule: a viz's frame dir is the snake_case gufe class name."""
    viz = framejs.VIZ_REGISTRY[cls_name]
    assert viz.frame == frame
    assert viz.has_local(), f"viz_assets/{frame}/code.js is missing from the package"
    assert "onInputs" in viz.js_source()


@pytest.mark.parametrize("cls_name", sorted(EXPECTED_FRAMES))
def test_every_registered_viz_has_a_serializer(cls_name):
    assert callable(framejs.VIZ_REGISTRY[cls_name].payload)


def test_registry_lookup_walks_the_mro():
    """Subclasses inherit their base's viz — that is how SolvatedPDBComponent /
    ProteinMembraneComponent get the protein viz, and NonTransformation the
    transformation viz, without registering anything themselves."""

    class Sub(LigandNetwork):
        pass

    assert framejs._registry_lookup(Sub.__new__(Sub), framejs.VIZ_REGISTRY).frame == "ligand_network"
    assert framejs._registry_lookup(object(), framejs.VIZ_REGISTRY) is None


def test_json_safe_stringifies_units():
    """openff Quantity et al. are display-only, so they go over the wire as str."""
    safe = framejs._json_safe({"score": 0.5, "length": 1.0 * unit.angstrom})
    assert safe["score"] == 0.5
    assert isinstance(safe["length"], str)
    json.dumps(safe)  # must not raise


def test_mapping_payload_is_json_serializable(simple_network):
    """Regression: mapping annotations carry openff Quantity objects, which broke
    both transports (widget set_inputs and the CLI hash) with a TypeError."""
    edge = next(e for e in simple_network.edges if "length" in e.annotations)
    payload = framejs._viz_for(edge).payload(edge)
    json.dumps(payload)  # must not raise
    assert set(payload) == {"molA.sdf", "molB.sdf", "nameA", "nameB", "mapping", "annotations"}
    # mapping keys are stringified for JSON; values stay ints
    assert all(isinstance(k, str) for k in payload["mapping"])


def test_small_molecule_payload(simple_network):
    mol = next(iter(simple_network.nodes))
    payload = framejs._viz_for(mol).payload(mol)
    assert payload["name"] == mol.name
    assert payload["smiles"] == mol.smiles
    assert "V2000" in payload["molecule.sdf"] or "V3000" in payload["molecule.sdf"]
    json.dumps(payload)


def test_solvent_payload():
    from gufe import SolventComponent

    payload = framejs._viz_for(SolventComponent()).payload(SolventComponent())
    assert set(payload["solvent_component"]) == {
        "smiles",
        "positive_ion",
        "negative_ion",
        "neutralize",
        "ion_concentration",
    }
    json.dumps(payload)


def test_chemical_system_payload(simple_network):
    from gufe import ChemicalSystem, SolventComponent

    mol = next(iter(simple_network.nodes))
    system = ChemicalSystem({"ligand": mol, "solvent": SolventComponent()}, name="sys")
    payload = framejs._viz_for(system).payload(system)["chemical_system"]
    assert payload["name"] == "sys"
    assert [c["label"] for c in payload["components"]] == ["ligand", "solvent"]
    # each component is described by type + name, plus its own structural payload
    by_label = {c["label"]: c for c in payload["components"]}
    assert by_label["ligand"]["type"] == "SmallMoleculeComponent"
    assert "sdf" in by_label["ligand"]
    assert by_label["solvent"]["type"] == "SolventComponent"
    json.dumps(payload)


def test_component_descriptor_records_errors_instead_of_raising():
    """A component that cannot be serialized must degrade to an error string, not
    take down the whole ChemicalSystem view."""

    class Exploding(SmallMoleculeComponent):
        def to_sdf(self):
            raise ValueError("nope")

    mol = SmallMoleculeComponent(mol_from_smiles("CCO"), name="boom")
    mol.__class__ = Exploding
    desc = framejs._component_descriptor("ligand", mol)
    assert "nope" in desc["error"]
    json.dumps(desc)


# --------------------------------------------------------------------------- #
# Legacy fallback — the pre-framejs RDKit / py3Dmol renderers                    #
# --------------------------------------------------------------------------- #


def _break_framejs(monkeypatch):
    """Make every framejs render path fail, as if the viz extra were missing."""
    monkeypatch.setattr(
        framejs,
        "_build_widget",
        lambda *a, **k: (_ for _ in ()).throw(framejs.FramejsUnavailable("boom")),
    )


def test_no_viewable_class_defines_ipython_display(simple_network):
    """`_ipython_display_` short-circuits `_repr_mimebundle_` in IPython, so a
    class defining both would show a different thing in a bare cell than through
    `.view()`. Display belongs to FramejsViewable alone."""
    offenders = [
        cls_name
        for cls_name in EXPECTED_FRAMES
        if any("_ipython_display_" in vars(k) for k in _class_for(cls_name).__mro__)
    ]
    assert offenders == []


def _class_for(cls_name):
    import gufe
    import gufe.transformations

    return {
        "LigandNetwork": gufe.LigandNetwork,
        "AlchemicalNetwork": gufe.AlchemicalNetwork,
        "TransformationBase": gufe.transformations.transformation.TransformationBase,
        "ChemicalSystem": gufe.ChemicalSystem,
        "LigandAtomMapping": gufe.LigandAtomMapping,
        "SmallMoleculeComponent": gufe.SmallMoleculeComponent,
        "ProteinComponent": gufe.ProteinComponent,
        "SolventComponent": gufe.SolventComponent,
    }[cls_name]


def test_mapping_legacy_view_returns_an_image(simple_network):
    edge = next(iter(simple_network.edges))
    img = edge._legacy_view()
    assert hasattr(img, "_repr_png_")  # IPython.display.Image
    assert img._repr_png_()


def test_view_falls_back_to_legacy_renderer(simple_network, monkeypatch):
    """LigandAtomMapping has a legacy renderer, so `.view()` must return it."""
    _break_framejs(monkeypatch)
    edge = next(iter(simple_network.edges))
    with pytest.warns(UserWarning):
        result = edge.view()
    assert result is not None
    assert result._repr_png_()


def test_repr_mimebundle_falls_back_to_legacy_renderer(simple_network, monkeypatch):
    """A bare cell must keep the pre-framejs picture when framejs is unavailable."""
    _break_framejs(monkeypatch)
    edge = next(iter(simple_network.edges))
    bundle = edge._repr_mimebundle_()
    assert bundle is not None
    # a mimebundle is either the data dict or a (data, metadata) pair
    data = bundle[0] if isinstance(bundle, tuple) else bundle
    assert "image/png" in data


def test_repr_mimebundle_is_none_without_a_legacy_renderer(simple_network, monkeypatch):
    """LigandNetwork has no legacy renderer, so the notebook gets the plain repr."""
    _break_framejs(monkeypatch)
    assert simple_network._repr_mimebundle_() is None


def test_repr_mimebundle_does_not_warn(simple_network, monkeypatch):
    """It runs on every display of the object — a warning per render is intolerable."""
    _break_framejs(monkeypatch)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert simple_network._repr_mimebundle_() is None
