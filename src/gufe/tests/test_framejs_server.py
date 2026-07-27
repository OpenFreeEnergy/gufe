# This code is part of gufe and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
"""Tests for the local viz server (gufe.visualization.server) — the CLI path.

Hermetic: the server is bound to an ephemeral loopback port in-process and driven
with ``urllib``. Nothing here reaches framejs.io — the frame URL is built from the
on-disk frame directory, exactly as it is at runtime, and never fetched.

The invariant these guard: the page carries the viz JavaScript and *no object
data*, and the object data comes back separately from ``?inputs=1`` — fetched by
the page, which is same-origin with this server. It deliberately is **not** the
frame that fetches: a frame served from framejs.io is a public origin, and
browsers now refuse those a loopback address without a user permission prompt.
"""

import base64
import functools
import http.server
import json
import threading
import urllib.error
import urllib.parse
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
    # by name, not by iteration order: `nodes` is a frozenset, so `next(iter(…))`
    # picks a different molecule from run to run
    mol = next(m for m in simple_network.nodes if m.name == "ethanol")
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

    def get(path, headers=False):
        """Return (status, body-text) — HTTP errors come back, they don't raise."""
        try:
            with urllib.request.urlopen(base + path, timeout=10) as r:
                return (r.status, dict(r.headers)) if headers else (r.status, r.read().decode("utf-8"))
        except urllib.error.HTTPError as e:
            return (e.code, dict(e.headers)) if headers else (e.code, e.read().decode("utf-8"))

    get.base = base
    try:
        yield get
    finally:
        httpd.shutdown()
        httpd.server_close()


def page_config(body):
    """The JSON config block the page is rendered with."""
    return json.loads(body.split('type="application/json">', 1)[1].split("</script>", 1)[0])


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
    config = page_config(body)
    # the viz JavaScript is baked in by Python...
    assert config["frame_url"] == server.resolve_url()
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
    assert page_config(body)["metapage_module"] == server.METAPAGE_MODULE
    assert "@metapages/metapage@" in server.METAPAGE_MODULE


def test_the_page_fetches_the_data_not_the_frame(client):
    """Load-bearing: a framejs.io frame is a public origin, and browsers refuse
    those a loopback address without a user permission prompt (Chrome's Local
    Network Access). The page is same-origin with this server, so it can."""
    _, body = client("/network.graphml")
    config = page_config(body)
    assert config["inputs_url"].startswith("/")  # relative == same-origin, no address space crossed
    assert "://" not in config["inputs_url"]
    # nothing in what the frame is handed points back at this server
    assert "http://" not in config["frame_url"]
    assert "setInputs" in body  # the data goes over postMessage instead


def test_viewer_url_path_is_the_data_file_path(client):
    """The URL says what it shows, including for a file in a subdirectory."""
    status, body = client("/sub/ethanol.sdf")
    assert status == 200
    config = page_config(body)
    assert config["path"] == "sub/ethanol.sdf"
    assert config["inputs_url"] == "/sub/ethanol.sdf?inputs=1"


# --------------------------------------------------------------------------- #
# Proxy mode — framejs re-hosted here, so the frame can fetch its own data      #
# --------------------------------------------------------------------------- #

UPSTREAM_DOC = b"<!doctype html><title>framejs</title><script>/* app is inline */</script>"


@pytest.fixture
def proxy_client(data_dir, monkeypatch):
    """A proxy-mode server, with framejs.io stubbed out (these stay hermetic)."""
    fetched = []

    def fake_fetch(url, **kwargs):
        fetched.append(url)
        return UPSTREAM_DOC, "text/html; charset=utf-8"

    monkeypatch.setattr(server, "_fetch_upstream", fake_fetch)
    handler = functools.partial(server._VizHandler, root=data_dir.resolve(), quiet=True, proxy=True)
    httpd = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{httpd.server_address[1]}"

    def get(path):
        try:
            with urllib.request.urlopen(base + path, timeout=10) as r:
                return r.status, r.read().decode("utf-8")
        except urllib.error.HTTPError as e:
            return e.code, e.read().decode("utf-8")

    get.base = base
    get.fetched = fetched
    try:
        yield get
    finally:
        httpd.shutdown()
        httpd.server_close()


def iframe_src(body):
    """The one iframe's src, un-escaped."""
    return body.split('<iframe src="', 1)[1].split('"', 1)[0].replace("&amp;", "&")


def test_proxy_page_has_no_script_and_one_iframe(proxy_client):
    """With the frame feeding itself there is nothing left for the page to do."""
    status, body = proxy_client("/network.graphml")
    assert status == 200
    assert body.count("<iframe") == 1
    assert "<script" not in body
    assert "metapage" not in body.lower()


def test_proxy_frame_and_data_share_one_origin(proxy_client):
    """The whole point: no address space is crossed, so Local Network Access
    never applies and Chrome does not block the frame's fetch."""
    _, body = proxy_client("/network.graphml")
    src = iframe_src(body)
    assert src.startswith(proxy_client.base + "/_framejs/#?js=")

    encoded = src.split("inputs=", 1)[1].split("&", 1)[0]
    inputs = json.loads(urllib.parse.unquote(base64.b64decode(encoded).decode("ascii")))
    assert inputs == {"network.graphml": f"{proxy_client.base}/network.graphml?input=network.graphml"}
    # frame origin == data origin, which is the invariant that makes this work
    assert urllib.parse.urlsplit(src).netloc == urllib.parse.urlsplit(inputs["network.graphml"]).netloc


def test_proxy_serves_the_framejs_document_from_here(proxy_client):
    status, body = proxy_client("/_framejs/")
    assert status == 200
    assert body == UPSTREAM_DOC.decode()
    assert proxy_client.fetched == ["https://framejs.io/"]


def test_proxy_forwards_the_apps_own_root_relative_requests(proxy_client):
    """framejs asks for its assets relative to whatever origin it is on — which is
    now us, so they have to reach upstream unchanged."""
    assert proxy_client("/favicon.svg")[0] == 200
    assert proxy_client.fetched == ["https://framejs.io/favicon.svg"]


def test_proxy_refuses_to_install_a_service_worker(proxy_client):
    """A service worker registered here would take scope `/` over the data
    endpoints and outlive both the page and the server — see PROXY_DENY."""
    status, body = proxy_client("/sw.js")
    assert status == 404
    assert "not proxied" in body
    assert proxy_client.fetched == []  # never even asked upstream


def test_proxy_does_not_shadow_the_data_files(proxy_client):
    """A path that exists locally is still ours; only unknown paths go upstream."""
    status, body = proxy_client("/network.graphml?raw=1")
    assert status == 200
    assert "<graphml" in body
    assert proxy_client.fetched == []


def test_unknown_paths_are_404_without_proxy(client):
    assert client("/sw.js")[0] == 404


def test_frame_src_for_proxy_rehosts_both_url_forms():
    from gufe.visualization.framejs import FRAMEJS_BASE

    local = server.frame_src_for_proxy(f"{FRAMEJS_BASE}/#?js=abc", {}, base="http://h:1")
    assert local.startswith("http://h:1/_framejs/#?js=abc&inputs=")
    # the /j/<uuid> form maps just as directly; the catch-all proxy resolves it
    canonical = server.frame_src_for_proxy(f"{FRAMEJS_BASE}/j/deadbeef", {}, base="http://h:1")
    assert canonical.startswith("http://h:1/_framejs/j/deadbeef#?inputs=")


def test_hash_inputs_replaces_only_the_big_values(client):
    """`hash_inputs` is not on the viewer path (see its docstring), but the frames
    do accept URL inputs, so its contract is still worth pinning down."""
    payload = {"molecule.sdf": "V2000\n" * 100, "name": "ethanol", "total_charge": 0}
    rewritten = server.hash_inputs(payload, data_url="http://host/x.sdf")
    assert rewritten["molecule.sdf"] == "http://host/x.sdf?input=molecule.sdf"
    assert rewritten["name"] == "ethanol"  # small scalars are not worth a round-trip
    assert rewritten["total_charge"] == 0


def test_input_endpoint_serves_one_payload_value(client, data_dir):
    """What the frame fetches, for exactly the key the hash pointed it at."""
    status, body = client("/network.graphml?input=network.graphml")
    assert status == 200
    assert body == server.payload_for(data_dir / "network.graphml")["network.graphml"]
    assert "<graphml" in body


def test_input_endpoint_serves_structured_values_as_json(client, data_dir, simple_network):
    edge = next(iter(simple_network.edges))
    (data_dir / "mapping.json").write_text(edge.to_json())
    status, body = client("/mapping.json?input=mapping")
    assert status == 200
    assert json.loads(body) == server.payload_for(data_dir / "mapping.json")["mapping"]


def test_unknown_input_key_is_404(client):
    assert client("/network.graphml?input=nope")[0] == 404


def test_responses_are_readable_by_the_framejs_origin(client):
    """The frame is cross-origin by construction, so CORS is not optional here."""
    for path in ("/network.graphml", "/network.graphml?input=network.graphml"):
        _, headers = client(path, headers=True)
        assert headers["Access-Control-Allow-Origin"] == "*"


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


# --------------------------------------------------------------------------- #
# The `openfe view` integration contract                                        #
# --------------------------------------------------------------------------- #
#
# `openfe view` (openfecli/commands/view.py) is a thin click wrapper: it imports a
# fixed set of names from gufe and calls them. openfe is NOT installed in gufe's
# test env, so these tests do not import openfecli — they instead pin the exact
# gufe-side surface the wrapper depends on, from gufe alone. If one of these
# fails, `openfe view` is broken even though everything else here is green; fix
# the wrapper and gufe together, and update both sides of the contract.
#
# The wrapper's imports, verbatim (keep in sync with openfecli/commands/view.py):
#   from gufe.visualization.server import (
#       DEFAULT_HOST, DEFAULT_PORT, UnviewableFile, load_object, serve)
#   from gufe.visualization.framejs import (
#       FramejsUnavailable, build_cli_url, build_html)

import inspect


def test_openfe_view_imports_resolve():
    """Every name `openfe view` imports from gufe still exists."""
    from gufe.visualization.framejs import FramejsUnavailable, build_cli_url, build_html  # noqa: F401
    from gufe.visualization.server import (  # noqa: F401
        DEFAULT_HOST,
        DEFAULT_PORT,
        UnviewableFile,
        load_object,
        serve,
    )


def test_serve_accepts_every_flag_the_cli_maps_to():
    """The click options must line up with `serve`'s keyword parameters; a rename
    on either side silently drops a flag."""
    params = inspect.signature(server.serve).parameters
    assert "target" in params  # the positional the CLI passes the file/dir as
    for kw in ("host", "port", "open_browser", "quiet", "proxy"):
        assert params[kw].kind == inspect.Parameter.KEYWORD_ONLY, kw


def test_build_cli_url_and_build_html_match_the_cli_calls():
    """`--url-only` calls build_cli_url(obj); `--html [--self-contained]` calls
    build_html(obj, self_contained=...)."""
    from gufe.visualization.framejs import build_cli_url, build_html

    inspect.signature(build_cli_url).bind(object())  # openfe view --url-only
    inspect.signature(build_html).bind(object(), self_contained=True)  # --html --self-contained


def test_serve_reaches_the_wrapper_the_same_way_the_cli_does(simple_network, tmp_path, monkeypatch):
    """Reproduce exactly what `_run()` does for the serve path — load a file, then
    `serve(path, host=, port=, open_browser=, proxy=)` — without importing openfe.
    This is the whole body of `openfe view <file>`."""
    graphml = tmp_path / "net.graphml"
    graphml.write_text(simple_network.to_graphml())

    captured = {}
    monkeypatch.setattr(server, "serve", lambda target, **kw: captured.update(target=str(target), **kw))

    # ---- verbatim from openfecli/commands/view.py `_run`, serve branch ----
    obj = server.load_object(Path(graphml))
    server.serve(graphml, host="0.0.0.0", port=8899, open_browser=False, proxy=True)
    # -----------------------------------------------------------------------

    assert isinstance(obj, LigandNetwork)
    assert captured == {"target": str(graphml), "host": "0.0.0.0", "port": 8899, "open_browser": False, "proxy": True}


def test_the_demo_service_command_drives_the_real_entrypoint(tmp_path, monkeypatch):
    """The `openfe-view` demo service runs `python -m gufe.visualization.server`
    with these args. Drive the real `main()` (mocking only the blocking serve) so
    a change to its argparse is caught — this is the demo's contract, and the
    fallback way to smoke-test the CLI path with no openfe present."""
    target = tmp_path / "data"
    target.mkdir()

    captured = {}
    monkeypatch.setattr(server, "serve", lambda t, **kw: captured.update(target=str(t), **kw))

    # mirrors visualization-demo/docker-compose.yml `openfe-view` command
    rc = server.main(["--host=0.0.0.0", "--port=8899", "--no-browser", "--proxy", str(target)])

    assert rc == 0
    assert captured == {
        "target": str(target),
        "host": "0.0.0.0",
        "port": 8899,
        "open_browser": False,
        "quiet": False,
        "proxy": True,
    }
