# This code is part of gufe and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
"""Tests for the local viz server (gufe.visualization.server) — the CLI path.

Hermetic: the server is bound to an ephemeral loopback port in-process and driven
with ``urllib``. Nothing here reaches framejs.io — the frame URL is built from the
on-disk frame directory, exactly as it is at runtime, and never fetched.

The invariant these guard is the one that makes the terminal path match the
notebook path: the page carries the viz JavaScript and *no object data*, and the
object data comes back separately from ``?inputs=1``.
"""

import functools
import http.server
import json
import threading
import urllib.error
import urllib.request
from pathlib import Path

import pytest
from rdkit import Chem

from gufe import LigandAtomMapping, LigandNetwork, SmallMoleculeComponent
from gufe.visualization import server


def mol_from_smiles(smi):
    m = Chem.MolFromSmiles(smi)
    m.Compute2DCoords()
    return m


@pytest.fixture
def simple_network():
    mol1 = SmallMoleculeComponent(mol_from_smiles("CCO"), name="ethanol")
    mol2 = SmallMoleculeComponent(mol_from_smiles("CC"), name="ethane")
    return LigandNetwork([LigandAtomMapping(mol1, mol2, {0: 0, 1: 1}, {"score": 0.0})])


@pytest.fixture
def data_dir(tmp_path, simple_network):
    """A directory holding one file of each interesting kind."""
    (tmp_path / "network.graphml").write_text(simple_network.to_graphml())
    mol = next(iter(simple_network.nodes))
    (tmp_path / "ethanol.sdf").write_text(mol.to_sdf())
    simple_network.to_json(file=tmp_path / "network.json")
    (tmp_path / "notes.txt").write_text("not a gufe object")
    (tmp_path / "sub").mkdir()
    (tmp_path / "sub" / "ethanol.sdf").write_text(mol.to_sdf())
    return tmp_path


@pytest.fixture
def client(data_dir):
    """A running server rooted at ``data_dir``; returns ``get(path)``."""
    handler = functools.partial(server._VizHandler, root=data_dir.resolve(), quiet=True)
    httpd = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{httpd.server_address[1]}"

    def get(path):
        """Return (status, body-text) — HTTP errors come back, they don't raise."""
        try:
            with urllib.request.urlopen(base + path, timeout=10) as r:
                return r.status, r.read().decode("utf-8")
        except urllib.error.HTTPError as e:
            return e.code, e.read().decode("utf-8")

    try:
        yield get
    finally:
        httpd.shutdown()
        httpd.server_close()


# --------------------------------------------------------------------------- #
# Loading: file extension -> gufe object                                        #
# --------------------------------------------------------------------------- #


def test_load_object_dispatches_on_extension(data_dir):
    assert isinstance(server.load_object(data_dir / "network.graphml"), LigandNetwork)
    assert isinstance(server.load_object(data_dir / "ethanol.sdf"), SmallMoleculeComponent)


def test_load_object_reads_gufe_json(data_dir, simple_network):
    """One `.json` loader covers every tokenizable — the class comes out of the file."""
    loaded = server.load_object(data_dir / "network.json")
    assert isinstance(loaded, LigandNetwork)
    assert loaded == simple_network


def test_load_object_rejects_unknown_extension(data_dir):
    with pytest.raises(server.UnviewableFile, match="unrecognized extension"):
        server.load_object(data_dir / "notes.txt")


def test_is_viewable(data_dir):
    assert server.is_viewable(data_dir / "network.graphml")
    assert not server.is_viewable(data_dir / "notes.txt")


def test_every_loader_extension_is_lowercase_and_dotted():
    assert all(ext.startswith(".") and ext.islower() for ext in server.LOADER_REGISTRY)


# --------------------------------------------------------------------------- #
# The viewer page                                                               #
# --------------------------------------------------------------------------- #


def test_viewer_page_bakes_in_the_frame_url_but_no_data(client, data_dir):
    status, body = client("/network.graphml")
    assert status == 200
    config = json.loads(body.split('type="application/json">', 1)[1].split("</script>", 1)[0])
    # the viz JavaScript is baked in by Python...
    assert config["frame_url"] == server._viz_for(server.load_object(data_dir / "network.graphml")).resolve_url()
    assert config["object_type"] == "LigandNetwork"
    # ...and the object is NOT: it is fetched from the server, like the notebook
    # pushes it over the comm channel.
    assert "inputs=" not in config["frame_url"]
    assert config["inputs_url"] == "/network.graphml?inputs=1"
    assert "<graphml" not in body


def test_viewer_page_is_a_single_metaframe(client):
    """One iframe, whose URL is the one Python baked into the config."""
    _, body = client("/network.graphml")
    assert "metaframes: { mf: { url: cfg.frame_url } }" in body
    assert body.count("metaframes:") == 1


def test_viewer_page_pins_the_same_metapage_as_the_widget(client):
    """The page drives the frame through the runtime the Jupyter widget uses."""
    _, body = client("/ethanol.sdf")
    config = json.loads(body.split('type="application/json">', 1)[1].split("</script>", 1)[0])
    assert config["metapage_module"] == server.METAPAGE_MODULE
    assert "@metapages/metapage@" in server.METAPAGE_MODULE


def test_viewer_url_path_is_the_data_file_path(client):
    """The URL says what it shows, including for a file in a subdirectory."""
    status, body = client("/sub/ethanol.sdf")
    assert status == 200
    config = json.loads(body.split('type="application/json">', 1)[1].split("</script>", 1)[0])
    assert config["path"] == "sub/ethanol.sdf"
    assert config["inputs_url"] == "/sub/ethanol.sdf?inputs=1"


# --------------------------------------------------------------------------- #
# The payload endpoint                                                          #
# --------------------------------------------------------------------------- #


def test_inputs_endpoint_returns_the_widget_payload(client, data_dir):
    status, body = client("/network.graphml?inputs=1")
    assert status == 200
    payload = json.loads(body)
    # byte-identical to what the notebook widget pushes for the same object
    assert payload == server.payload_for(data_dir / "network.graphml")
    assert "<graphml" in payload["network.graphml"]


def test_inputs_are_rebuilt_from_disk_on_every_request(client, data_dir, simple_network):
    """Editing the file and re-fetching re-renders — no restart, no stale link."""
    before = json.loads(client("/ethanol.sdf?inputs=1")[1])
    mol = SmallMoleculeComponent(mol_from_smiles("CCCO"), name="propanol")
    (data_dir / "ethanol.sdf").write_text(mol.to_sdf())
    after = json.loads(client("/ethanol.sdf?inputs=1")[1])
    assert before["name"] == "ethanol"
    assert after["name"] == "propanol"


def test_raw_endpoint_serves_the_file_itself(client, data_dir):
    status, body = client("/network.graphml?raw=1")
    assert status == 200
    assert body == (data_dir / "network.graphml").read_text()


# --------------------------------------------------------------------------- #
# Listings, and refusing what it cannot (or must not) serve                      #
# --------------------------------------------------------------------------- #


def test_directory_listing_links_viewable_files_only(client):
    status, body = client("/")
    assert status == 200
    assert 'href="/network.graphml"' in body
    assert 'href="/sub/"' in body
    assert "notes.txt" in body  # shown...
    assert 'href="/notes.txt"' not in body  # ...but not linked: nothing can view it


def test_unviewable_file_gets_an_explanation_not_a_traceback(client):
    status, body = client("/notes.txt")
    assert status == 415
    assert "unrecognized extension" in body


def test_malformed_file_reports_the_error(client, data_dir):
    (data_dir / "broken.graphml").write_text("<not really graphml>")
    status, body = client("/broken.graphml")
    assert status == 422
    assert "broken.graphml" in body


def test_missing_file_is_404(client):
    assert client("/nope.sdf")[0] == 404


def test_paths_outside_the_root_are_refused(client):
    # normalized before it ever touches the filesystem, so this is a 403/404 and
    # never a read of /etc/passwd
    assert client("/../../etc/passwd")[0] in (403, 404)
    assert client("/%2e%2e/%2e%2e/etc/passwd")[0] in (403, 404)


# --------------------------------------------------------------------------- #
# Root/URL selection                                                            #
# --------------------------------------------------------------------------- #


def test_root_is_the_working_directory_when_the_file_is_under_it(data_dir, monkeypatch):
    monkeypatch.chdir(data_dir)
    root, url_path = server._root_and_url_path(Path("sub/ethanol.sdf"))
    assert root == data_dir.resolve()
    assert url_path == "/sub/ethanol.sdf"


def test_root_is_the_files_own_directory_when_it_is_elsewhere(data_dir, tmp_path_factory, monkeypatch):
    monkeypatch.chdir(tmp_path_factory.mktemp("elsewhere"))
    root, url_path = server._root_and_url_path(data_dir / "network.graphml")
    assert root == data_dir.resolve()
    assert url_path == "/network.graphml"


def test_a_directory_target_serves_its_listing(data_dir, monkeypatch):
    monkeypatch.chdir(data_dir)
    root, url_path = server._root_and_url_path(Path("sub"))
    assert root == data_dir.resolve()
    assert url_path == "/sub"
