# This code is part of gufe and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
"""Tests for the framejs.io visualization infrastructure (gufe.visualization.framejs).

These cover the *infrastructure* (registry, serializers, URL building, display
hooks, fallback) and deliberately do not assert on the viz JavaScript itself.
They are fully **hermetic** — no network: the canonical frame ``uuid`` is
monkeypatched in by fixture so we never need a live published frame. Tests that
need the optional ``metaframe-widget`` dependency skip cleanly when it is not
installed (it lives behind the ``gufe[viz]`` extra).

One frame draws every gufe object, so the base URL does not depend on the object
and these tests exercise it through whichever object is convenient. What *is*
per-class is the serializer, and the breadth tests at the bottom pin that down.
"""

import base64
import importlib
import json
import urllib.parse
import warnings

import pytest
from openff.units import unit
from rdkit import Chem

from gufe import LigandAtomMapping, LigandNetwork, SmallMoleculeComponent
from gufe.visualization import framejs

# A fake published frame id, monkeypatched in so canonical_url() resolves
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
def published_frame(monkeypatch):
    """A pinned uuid with no frame directory on disk — the one situation where
    the canonical ``/j/<uuid>`` form is what `resolve_url()` returns. Hermetic:
    nothing is fetched."""
    monkeypatch.setattr(framejs, "FRAME", "does_not_exist")
    monkeypatch.setattr(framejs, "FRAME_UUID", FAKE_UUID)
    return FAKE_UUID


@pytest.fixture
def short_url_frame(monkeypatch):
    """A pinned uuid with the frame dir intact — so `build_cli_url(short=True)`
    has a uuid to build on while the default path still resolves locally."""
    monkeypatch.setattr(framejs, "FRAME_UUID", FAKE_UUID)
    return FAKE_UUID


@pytest.fixture
def unavailable_frame(monkeypatch):
    """Neither form available: no pinned uuid AND no on-disk frame, so URL
    building must raise / fall back."""
    monkeypatch.setattr(framejs, "FRAME", "does_not_exist")
    monkeypatch.setattr(framejs, "FRAME_UUID", None)


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
# The one frame                                                                 #
# --------------------------------------------------------------------------- #


def test_the_frame_ships_on_disk():
    assert framejs.FRAME == "gufe"
    assert framejs.has_local()  # viz_assets/gufe/code.js ships with gufe


def test_js_source_loads_from_the_frame_dir():
    src = framejs.js_source()  # reads viz_assets/gufe/code.js
    assert "onInputs" in src  # the framejs render entrypoint (may be `async`)
    assert len(src) > 100


def test_js_source_raises_without_a_frame_dir(unavailable_frame):
    assert framejs.has_local() is False
    with pytest.raises(framejs.FramejsUnavailable):
        framejs.js_source()


def test_canonical_url_raises_when_unpublished(unavailable_frame):
    with pytest.raises(framejs.FramejsUnavailable):
        framejs.canonical_url()


def test_canonical_url_uses_the_pinned_uuid(short_url_frame):
    assert framejs.canonical_url() == f"https://framejs.io/j/{FAKE_UUID}"


# --------------------------------------------------------------------------- #
# Local URL (built from the on-disk frame directory) + source selection        #
# --------------------------------------------------------------------------- #


def test_local_url_encodes_the_frame():
    url = framejs.local_url()
    assert url.startswith("https://framejs.io/#?js=")
    # the js hash param decodes back to the on-disk code.js (btoa(encodeURIComponent(js)))
    js_param = urllib.parse.parse_qs(url.split("#?", 1)[1])["js"][0]
    decoded = urllib.parse.unquote(base64.b64decode(js_param).decode("ascii"))
    assert decoded == framejs.js_source()
    assert "onInputs" in decoded
    # og.json sidecar is carried as an og hash param
    assert "&og=" in url
    # The frame ships no modules.json: d3, RDKit and 3Dmol are all fetched on
    # demand by the frame itself, so they stay off the initial-render path.
    assert "&modules=" not in url
    for lazy in ("load3Dmol", "loadRDKit", "loadD3"):
        assert lazy in decoded


def test_resolve_url_prefers_local_over_uuid(short_url_frame):
    """The shipped frame wins over the pinned uuid: it always matches the
    installed code and cannot expire."""
    assert framejs.has_local()
    assert framejs.resolve_url() == framejs.local_url()


def test_resolve_url_falls_back_to_uuid_without_a_frame_dir(published_frame):
    assert not framejs.has_local()
    assert framejs.resolve_url() == f"https://framejs.io/j/{FAKE_UUID}"


def test_resolve_url_raises_when_neither_form_exists(unavailable_frame):
    with pytest.raises(framejs.FramejsUnavailable):
        framejs.resolve_url()


def test_one_frame_serves_every_object(simple_network):
    """The base URL does not depend on the object — that is the whole point of
    collapsing the per-class frames into one."""
    mol = next(iter(simple_network.nodes))
    edge = next(iter(simple_network.edges))
    urls = {framejs.resolve_url() for _ in (simple_network, mol, edge)}
    assert len(urls) == 1


def test_build_cli_url_local_default(simple_network):
    # The default builds a self-contained hash-param URL from the on-disk frame,
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
        framejs.payload_for_object(object())


# --------------------------------------------------------------------------- #
# Serializer (payload builder)                                                  #
# --------------------------------------------------------------------------- #


def test_network_payload_is_graphml(simple_network):
    # the viz reads inputs['network.graphml'] and parses GraphML itself
    payload = framejs.payload_for_object(simple_network)
    assert set(payload) == {"network.graphml"}
    xml = payload["network.graphml"]
    assert isinstance(xml, str)
    assert xml == simple_network.to_graphml()  # the exact gufe serialization
    assert "<graphml" in xml


def test_network_payload_is_json_serializable(simple_network):
    payload = framejs.payload_for_object(simple_network)
    # must round-trip cleanly (it travels as JSON to the browser)
    assert json.loads(json.dumps(payload)) == payload


def test_network_payload_stable(simple_network):
    # reproducible: same network -> identical payload
    assert framejs.payload_for_object(simple_network) == framejs.payload_for_object(simple_network)


# --------------------------------------------------------------------------- #
# CLI URL building — canonical URL + appended inputs                            #
# --------------------------------------------------------------------------- #


def test_build_cli_url_appends_inputs_to_canonical_url(simple_network, published_frame):
    url = framejs.build_cli_url(simple_network)
    assert url.startswith(f"https://framejs.io/j/{FAKE_UUID}#?inputs=")
    # the appended inputs decode back to the payload (btoa(encodeURIComponent(json)))
    encoded = url.split("#?inputs=", 1)[1]
    decoded = json.loads(urllib.parse.unquote(base64.b64decode(encoded).decode("ascii")))
    assert "network.graphml" in decoded
    assert "<graphml" in decoded["network.graphml"]
    # the viz js is NOT in the URL (it is hosted at the canonical frame)
    assert "js=" not in url


def test_build_cli_url_unavailable_raises(simple_network, unavailable_frame):
    # no frame on disk and no uuid makes the CLI URL unbuildable
    with pytest.raises(framejs.FramejsUnavailable):
        framejs.build_cli_url(simple_network)


def test_build_cli_url_short_uses_the_uuid(simple_network, short_url_frame):
    """`short=True` opts into the pinned /j/<uuid> even though the frame dir is
    present — the only reason to do so being URL size."""
    url = framejs.build_cli_url(simple_network, short=True)
    assert url.startswith(f"https://framejs.io/j/{FAKE_UUID}#?inputs=")
    assert "js=" not in url  # the viz JavaScript is hosted, not inlined


def test_build_cli_url_short_is_dramatically_smaller(simple_network, short_url_frame):
    long_url = framejs.build_cli_url(simple_network)
    short = framejs.build_cli_url(simple_network, short=True)
    assert long_url.startswith("https://framejs.io/#?js=")  # default stays local
    assert len(short) * 4 < len(long_url)  # inlining code.js is what makes it big


def test_build_cli_url_short_raises_when_unpublished(simple_network):
    """The shipped frame has no uuid yet, so there is no short form."""
    assert framejs.FRAME_UUID is None
    framejs.build_cli_url(simple_network)  # the default (local) form still works
    with pytest.raises(framejs.FramejsUnavailable):
        framejs.build_cli_url(simple_network, short=True)


# --------------------------------------------------------------------------- #
# Standalone HTML — one file, no server, no framejs.io                          #
# --------------------------------------------------------------------------- #


def test_build_html_inlines_the_frame_and_the_payload(simple_network):
    page = framejs.build_html(simple_network)
    assert page.startswith("<!doctype html>")
    # the viz JavaScript is in the file, not linked to framejs.io
    assert "framejs.io" not in page
    assert "export async function onInputs" in page
    assert '<script type="module">' in page
    # ...and so is the object
    inputs = json.loads(page.split('type="application/json">', 1)[1].split("</script>", 1)[0].replace("<\\/", "</"))
    assert "<graphml" in inputs["network.graphml"]
    assert inputs == framejs.payload_for_object(simple_network)


def test_build_html_bootstraps_the_frame(simple_network):
    """The frame is an ES module reading a global `root`, so the page must set
    that up before the module runs and then hand it the payload itself."""
    page = framejs.build_html(simple_network)
    assert 'window.root = document.getElementById("gufe-root");' in page
    # the classic script setting `root` must come before the deferred module
    assert page.index("window.root =") < page.index('<script type="module">')
    assert "await onInputs(JSON.parse(" in page
    assert "onResize()" in page


def test_build_html_titles_the_page_with_the_object_type(simple_network):
    assert "<title>LigandNetwork</title>" in framejs.build_html(simple_network)
    assert "<title>my run</title>" in framejs.build_html(simple_network, title="my run")


def test_build_html_escapes_the_title(simple_network):
    page = framejs.build_html(simple_network, title="<script>alert(1)</script>")
    assert "<title>&lt;script&gt;alert(1)&lt;/script&gt;</title>" in page


def test_build_html_defaults_to_cdn_engines(simple_network):
    """The default file is small and pulls the engines when opened; that is what
    keeps a report emailable."""
    page = framejs.build_html(simple_network)
    # the frame *reads* the pre-seed hook (it is in code.js), but nothing sets it
    assert "window.__gufeEngines = {" not in page
    assert "locateFile" not in page
    assert len(page.encode()) < 1_000_000


def test_build_html_payload_cannot_close_the_script_block():
    """A molecule name (or any payload string) containing `</script>` must not be
    able to break out of the JSON block."""
    from gufe import SolventComponent

    obj = SolventComponent()
    page = framejs.build_html(obj)
    body = page.split('type="application/json">', 1)[1].split("</script>", 1)[0]
    assert "</" not in body  # every one was escaped to `<\/`
    # and it is still valid JSON once unescaped, which is what the page does
    json.loads(body.replace("<\\/", "</"))


def test_script_safe_neutralizes_a_closing_tag():
    assert framejs._script_safe('var s = "</script>";') == 'var s = "<\\/script>";'
    assert framejs._script_safe("var s = '</SCRIPT>';") == "var s = '<\\/SCRIPT>';"
    assert framejs._script_safe("a < b") == "a < b"  # untouched


def test_the_frame_code_needs_no_escaping():
    """Belt and braces: `_script_safe` is applied to code.js, but the frame should
    not contain the sequence in the first place."""
    assert "</script" not in framejs.js_source().lower()


def test_build_html_requires_the_frame_on_disk(simple_network, published_frame):
    """Unlike the URL forms, a pinned uuid is no substitute — the whole point is
    to have the code in the file."""
    with pytest.raises(framejs.FramejsUnavailable):
        framejs.build_html(simple_network)


def test_build_html_unregistered_type_raises():
    with pytest.raises(framejs.FramejsUnavailable):
        framejs.build_html(object())


def test_build_html_self_contained_inlines_the_engines(simple_network, monkeypatch):
    """Hermetic: the downloader is stubbed, so this pins the *wiring* — the
    engines end up inline and pre-seeded through the hook code.js reads."""
    fetched = []

    def fake_fetch(url, **kwargs):
        fetched.append(url)
        return b"/* wasm */" if url.endswith(".wasm") else f"/* {url} */".encode()

    monkeypatch.setattr(framejs, "_fetch_engine", fake_fetch)
    page = framejs.build_html(simple_network, self_contained=True)

    assert sorted(fetched) == sorted(framejs.ENGINE_URLS.values())
    for url in framejs.ENGINE_URLS.values():
        if not url.endswith(".wasm"):
            assert f"/* {url} */" in page
    # the hook the frame's loaders check, with all three engines on it
    assert "window.__gufeEngines = { threeDmol: window.$3Dmol, d3: window.d3, rdkit: rdkitReady };" in page
    # RDKit's wasm has no neighbouring file to find, so it rides as a data: URL
    assert 'locateFile: () => "data:application/wasm;base64,' in page
    assert base64.b64encode(b"/* wasm */").decode() in page
    # engines must be inlined before the module that consumes them
    assert page.index("__gufeEngines") < page.index('<script type="module">')


def test_the_frame_reads_the_engine_preseed_hook():
    """The contract between build_html(self_contained=True) and code.js."""
    src = framejs.js_source()
    assert "__gufeEngines" in src
    for engine in ("threeDmol", "d3", "rdkit"):
        assert f"preseededEngine('{engine}')" in src


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


def test_view_returns_widget_with_canonical_url_and_inputs(simple_network, published_frame):
    w = simple_network.view()
    assert type(w).__name__ == "MetaframeWidget"
    # notebook path: canonical URL loaded, data pushed as inputs (NOT in the URL)
    assert w.url == f"https://framejs.io/j/{FAKE_UUID}"
    assert "network.graphml" in w.inputs
    assert "inputs=" not in w.url


def test_view_respects_size_overrides(simple_network, published_frame):
    w = simple_network.view(height="123px", width="456px")
    assert w.height == "123px"
    assert w.width == "456px"


def test_repr_mimebundle_returns_bundle(simple_network, published_frame):
    mb = simple_network._repr_mimebundle_()
    assert mb  # non-empty dict/tuple so the notebook auto-displays the widget


def test_repr_mimebundle_none_for_unregistered():
    # a gufe object with no registered viz must not break bare-cell display
    assert framejs.repr_mimebundle(object()) is None


def test_view_falls_back_when_frame_unavailable(simple_network, unavailable_frame):
    # no uuid and no frame -> FramejsUnavailable inside view_object -> legacy (None here)
    with pytest.warns(UserWarning):
        assert simple_network.view() is None


def test_view_falls_back_when_viz_unavailable(simple_network, published_frame, monkeypatch):
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
# in FramejsViewable must have a serializer, that serializer's keys must be ones
# the frame dispatches on, and its payload must survive the JSON transport.

# Every gufe class that opts in to a framejs view.
EXPECTED_CLASSES = {
    "LigandNetwork",
    "AlchemicalNetwork",
    "TransformationBase",
    "ChemicalSystem",
    "LigandAtomMapping",
    "SmallMoleculeComponent",
    "ProteinComponent",
    "SolventComponent",
}

# The payload keys the frame's VIEWS table dispatches on — the contract between
# PAYLOAD_REGISTRY and code.js. A serializer must produce exactly one of these.
DISPATCH_KEYS = {
    "network.graphml",
    "molA.sdf",
    "molecule.sdf",
    "protein.pdb",
    "solvent_component",
    "chemical_system",
    "transformation",
    "alchemical_network",
}


def test_registry_covers_exactly_the_viewable_classes():
    assert set(framejs.PAYLOAD_REGISTRY) == EXPECTED_CLASSES


@pytest.mark.parametrize("cls_name", sorted(EXPECTED_CLASSES))
def test_every_registered_class_has_a_serializer(cls_name):
    assert callable(framejs.PAYLOAD_REGISTRY[cls_name])


def test_the_frame_dispatches_on_every_key_the_registry_emits():
    """Each dispatch key the frame knows about is claimed by the frame's VIEWS
    table, so a payload can never arrive at a frame with no view for it."""
    src = framejs.js_source()
    for key in DISPATCH_KEYS:
        assert f"key: '{key}'" in src, f"the frame's VIEWS table has no entry for {key!r}"


def test_registry_lookup_walks_the_mro():
    """Subclasses inherit their base's serializer — that is how SolvatedPDBComponent /
    ProteinMembraneComponent get the protein payload, and NonTransformation the
    transformation payload, without registering anything themselves."""

    class Sub(LigandNetwork):
        pass

    assert framejs._registry_lookup(Sub.__new__(Sub), framejs.PAYLOAD_REGISTRY) is framejs._ligand_network_payload
    assert framejs._registry_lookup(object(), framejs.PAYLOAD_REGISTRY) is None


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
    payload = framejs.payload_for_object(edge)
    json.dumps(payload)  # must not raise
    assert set(payload) == {"molA.sdf", "molB.sdf", "nameA", "nameB", "mapping", "annotations"}
    # mapping keys are stringified for JSON; values stay ints
    assert all(isinstance(k, str) for k in payload["mapping"])


def test_small_molecule_payload(simple_network):
    mol = next(iter(simple_network.nodes))
    payload = framejs.payload_for_object(mol)
    assert payload["name"] == mol.name
    assert payload["smiles"] == mol.smiles
    assert "V2000" in payload["molecule.sdf"] or "V3000" in payload["molecule.sdf"]
    json.dumps(payload)


def test_solvent_payload():
    from gufe import SolventComponent

    payload = framejs.payload_for_object(SolventComponent())
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
    payload = framejs.payload_for_object(system)["chemical_system"]
    assert payload["name"] == "sys"
    assert [c["label"] for c in payload["components"]] == ["ligand", "solvent"]
    # each component is described by type + name, plus its own structural payload
    by_label = {c["label"]: c for c in payload["components"]}
    assert by_label["ligand"]["type"] == "SmallMoleculeComponent"
    assert "sdf" in by_label["ligand"]
    assert by_label["solvent"]["type"] == "SolventComponent"
    json.dumps(payload)


def test_every_payload_claims_exactly_one_dispatch_key(simple_network):
    """The frame picks its view from the first dispatch key present, so a payload
    that claimed two would render the wrong one."""
    from gufe import ChemicalSystem, SolventComponent

    mol = next(iter(simple_network.nodes))
    objects = [
        simple_network,
        next(iter(simple_network.edges)),
        mol,
        SolventComponent(),
        ChemicalSystem({"ligand": mol}, name="sys"),
    ]
    for obj in objects:
        claimed = DISPATCH_KEYS & set(framejs.payload_for_object(obj))
        assert len(claimed) == 1, f"{type(obj).__name__} claims {claimed}"


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


def test_no_viewable_class_defines_ipython_display():
    """`_ipython_display_` short-circuits `_repr_mimebundle_` in IPython, so a
    class defining both would show a different thing in a bare cell than through
    `.view()`. Display belongs to FramejsViewable alone."""
    offenders = [
        cls_name
        for cls_name in EXPECTED_CLASSES
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
