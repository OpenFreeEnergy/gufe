# This code is part of gufe and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
"""Tests for the framejs.io visualization infrastructure (gufe.visualization.framejs).

These cover the *infrastructure* (registry, serializers, canonical-URL building,
display hooks, fallback) and deliberately do not assert on the canonical viz
JavaScript itself. They are fully **hermetic** — no network: the canonical viz
``uuid`` is injected via the ``GUFE_VIZ_<ID>_UUID`` env override so we never need a
live published frame. Tests that need the optional ``metaframe-widget`` dependency
skip cleanly when it is not installed (it lives behind the ``gufe[viz]`` extra).
"""

import base64
import importlib
import json
import urllib.parse

import pytest
from openff.units import unit
from rdkit import Chem

from gufe import LigandAtomMapping, LigandNetwork, SmallMoleculeComponent
from gufe.visualization import framejs

# A fake published frame id, injected via the env override so canonical_url()
# resolves without any live framejs.io frame.
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
    """Registered LigandNetwork viz, pinned to FAKE_UUID and forced to the
    *canonical* source (``/j/<uuid>``) so these tests exercise that branch
    regardless of whether the on-disk frame dir is present."""
    monkeypatch.setenv("GUFE_VIZ_NETWORK_UUID", FAKE_UUID)
    monkeypatch.setenv("GUFE_VIZ_SOURCE", "canonical")
    return FAKE_UUID


@pytest.fixture
def local_network(monkeypatch):
    """Force the registered LigandNetwork viz to build its URL from the on-disk
    frame directory (``GUFE_VIZ_SOURCE=local``) — no network, no uuid needed."""
    monkeypatch.setenv("GUFE_VIZ_SOURCE", "local")


@pytest.fixture
def unpublished_network(monkeypatch):
    """Force the registered LigandNetwork viz to be fully unavailable: no pinned
    uuid AND no on-disk frame, so URL building must raise / fall back."""
    monkeypatch.setitem(framejs.CANONICAL_VIZ, "LigandNetwork", framejs.VizRef(id="network"))


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
    viz = framejs.CANONICAL_VIZ["LigandNetwork"]
    assert viz.id == "network"
    assert viz.frame == "network"
    assert viz.has_local()  # the on-disk frame dir viz_assets/network/ ships with gufe


def test_ligandnetwork_is_published():
    # the registered viz carries the published framejs.io uuid
    viz = framejs.CANONICAL_VIZ["LigandNetwork"]
    assert viz.published is True
    assert viz.canonical_url() == f"https://framejs.io/j/{viz.uuid}"


def test_vizref_without_uuid_is_unpublished():
    viz = framejs.VizRef(id="x")
    assert viz.published is False
    with pytest.raises(framejs.FramejsUnavailable):
        viz.canonical_url()


def test_vizref_canonical_url_from_env(published_network):
    viz = framejs.CANONICAL_VIZ["LigandNetwork"]
    assert viz.published is True
    assert viz.canonical_url() == f"https://framejs.io/j/{FAKE_UUID}"


def test_vizref_canonical_url_explicit_uuid():
    viz = framejs.VizRef(id="x", uuid="abc123")
    assert viz.canonical_url() == "https://framejs.io/j/abc123"


def test_vizref_canonical_url_with_version_pin():
    viz = framejs.VizRef(id="x", uuid="abc123", version="sha256hex")
    assert viz.canonical_url() == "https://framejs.io/j/abc123?v=sha256hex"


def test_vizref_js_source_loads_from_frame_dir():
    viz = framejs.CANONICAL_VIZ["LigandNetwork"]
    src = viz.js_source()  # reads viz_assets/network/code.js
    assert "onInputs" in src  # the framejs render entrypoint (may be `async`)
    assert len(src) > 100


def test_vizref_missing_frame_raises():
    bad = framejs.VizRef(id="nope", uuid="abc", frame="does_not_exist")
    assert bad.has_local() is False
    with pytest.raises(framejs.FramejsUnavailable):
        bad.js_source()


def test_vizref_no_frame_has_no_local():
    assert framejs.VizRef(id="x", uuid="abc").has_local() is False


# --------------------------------------------------------------------------- #
# Local URL (built from the on-disk frame directory) + source selection        #
# --------------------------------------------------------------------------- #


def test_vizref_local_url_encodes_frame():
    viz = framejs.CANONICAL_VIZ["LigandNetwork"]
    url = viz.local_url()
    assert url.startswith("https://framejs.io/#?js=")
    # the js hash param decodes back to the on-disk code.js (btoa(encodeURIComponent(js)))
    js_param = urllib.parse.parse_qs(url.split("#?", 1)[1])["js"][0]
    decoded = urllib.parse.unquote(base64.b64decode(js_param).decode("ascii"))
    assert decoded == viz.js_source()
    assert "onInputs" in decoded
    # og.json sidecar is carried as an og hash param
    assert "&og=" in url
    # modules.json sidecar (3Dmol.js) MUST be carried — the viz is broken without it
    assert "&modules=" in url
    params = urllib.parse.parse_qs(url.split("#?", 1)[1])
    modules = json.loads(urllib.parse.unquote(base64.b64decode(params["modules"][0]).decode("ascii")))
    assert any("3Dmol" in m for m in modules)


def test_resolve_url_auto_prefers_local():
    # default (no GUFE_VIZ_SOURCE): the on-disk frame wins over the pinned uuid
    viz = framejs.CANONICAL_VIZ["LigandNetwork"]
    assert viz.resolve_url() == viz.local_url()


def test_resolve_url_canonical_env(published_network):
    viz = framejs.CANONICAL_VIZ["LigandNetwork"]
    assert viz.resolve_url() == viz.canonical_url()
    assert viz.resolve_url() == f"https://framejs.io/j/{FAKE_UUID}"


def test_resolve_url_local_env(local_network):
    viz = framejs.CANONICAL_VIZ["LigandNetwork"]
    assert viz.resolve_url() == viz.local_url()


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
        framejs._viz_ref_for(object())
    with pytest.raises(framejs.FramejsUnavailable):
        framejs._payload_for(object())


# --------------------------------------------------------------------------- #
# Serializer (payload builder)                                                  #
# --------------------------------------------------------------------------- #


def test_network_payload_is_graphml(simple_network):
    # the canonical viz reads inputs['network.graphml'] and parses GraphML itself
    payload = framejs._payload_for(simple_network)
    assert set(payload) == {"network.graphml"}
    xml = payload["network.graphml"]
    assert isinstance(xml, str)
    assert xml == simple_network.to_graphml()  # the exact gufe serialization
    assert "<graphml" in xml


def test_network_payload_is_json_serializable(simple_network):
    payload = framejs._payload_for(simple_network)
    # must round-trip cleanly (it travels as JSON to the browser)
    assert json.loads(json.dumps(payload)) == payload


def test_network_payload_stable(simple_network):
    # reproducible: same network -> identical payload
    assert framejs._payload_for(simple_network) == framejs._payload_for(simple_network)


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
        "_viz_widget",
        lambda *a, **k: (_ for _ in ()).throw(framejs.FramejsUnavailable("boom")),
    )
    with pytest.warns(UserWarning):
        assert simple_network.view() is None
