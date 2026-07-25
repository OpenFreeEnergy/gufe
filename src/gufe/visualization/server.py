# This code is part of gufe and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
"""A small local web server that shows a gufe object's framejs visualization.

This is the **terminal** counterpart to the notebook path in :mod:`gufe.visualization.framejs`,
and it deliberately works the same way. In a notebook the widget loads the viz
into an iframe and then *pushes the object's data in* over the metapage channel;
here a locally served page does exactly that, in the browser, with no Python
kernel behind it::

    $ python -m gufe.visualization.server ligand_network.graphml
    gufe viz  →  http://localhost:8899/ligand_network.graphml

    open the URL above (Ctrl-C to stop the server)

What the browser gets at that URL is a page rendered by :func:`viewer_html`: a
header, one metaframe, and ~30 lines of JavaScript. The framejs URL of the viz
registered for that file's object type is **baked into the page** by Python; the
page then fetches the object's ``inputs`` from ``?inputs=1`` on its own path and
hands them to the frame through ``renderMetapage`` — the same
``@metapages/metapage`` entrypoint, at the same pinned version, that
``metaframe_widget`` uses in Jupyter.

Why *this page* does the fetching
---------------------------------
It would be tidier to put the inputs straight in the iframe's URL hash and let the
frame fetch its own data from ``?input=<key>`` — the frames support exactly that
(see :func:`hash_inputs`, and ``asText()``/``asObject()`` in every frame, which
fetch a whitespace-free ``http(s)://`` string). **Chrome no longer allows it**
across origins. Chrome 142's Local Network Access blocks a public-origin page —
which the frame, served from framejs.io, is — from reaching a loopback address
without an explicit user permission prompt, and inside a cross-origin iframe the
``local-network-access`` Permissions Policy has to be delegated on top. It fails
as a CORS error: *"Permission was denied for this request to access the
`loopback` address space"*. Firefox has not shipped this and still works. (The
older Private Network Access opt-in header, ``Access-Control-Allow-Private-Network``,
is superseded and ignored.)

This page has no such problem: it is served from ``localhost:8899`` itself, so its
fetch is same-origin — no address space is crossed and nothing prompts — and
handing the result to the frame is a ``postMessage``, not a network request.

Proxy mode (``--proxy``) removes the barrier instead
-----------------------------------------------------
The restriction is about *crossing* address spaces, so it disappears if the frame
is on this origin too. In proxy mode the server re-hosts the framejs.io document
at ``/_framejs/`` (and forwards anything else it does not recognize upstream, for
the app's own ``/sw.js`` and friends), so the iframe, the page and the data are
all ``http://localhost:8899``. Nothing is rewritten on the way through: framejs's
application code is inline in its HTML and its dependencies are absolute
``https://cdn.jsdelivr.net`` imports, which an ``http://localhost`` page — a
secure context, and an *upgrade* rather than mixed content — loads normally.

That buys back the simpler design: inputs go in the iframe's hash as URLs the
frame fetches for itself (:func:`hash_inputs`), the page needs no JavaScript at
all, and the URL stays small however large the object. It is opt-in because it
depends on framejs.io tolerating being re-hosted, which only a browser can
confirm.

Why serve at all, rather than just open a URL
---------------------------------------------
:func:`gufe.visualization.framejs.build_cli_url` puts everything — the viz
JavaScript *and* the object — in the URL's hash. That is perfect for a link you
want to paste somewhere, and it is still available (``--url-only`` in the CLI),
but it is a dead end for anything real: a solvated ``AlchemicalNetwork`` blows
past what a browser will accept in a URL, and the data is frozen at the moment
the link was made. Serving the payload instead means

* **no size limit** — the object travels as a normal HTTP response body;
* **live data** — the payload is rebuilt from disk on every fetch, so editing the
  file and hitting ⟳ in the page re-renders it without restarting anything;
* **one code path with the notebook** — both transports push ``inputs`` into an
  already-loaded frame, so a frame never has to care where it is embedded.

URL layout
----------
The server is rooted at a directory (the working directory when the file lives
under it, otherwise the file's own directory), and **the URL path is the data
file's path** relative to that root — a link that says what it shows::

    http://localhost:8899/setup/ligand_network.graphml       the viewer page
    http://localhost:8899/setup/ligand_network.graphml?input=network.graphml
                                                            one payload value (what the frame fetches)
    http://localhost:8899/setup/ligand_network.graphml?inputs=1   the whole payload (JSON)
    http://localhost:8899/setup/ligand_network.graphml?raw=1      the file itself
    http://localhost:8899/setup/                             a listing of what is viewable

Anything under the root can be viewed from the one server, so a single
``openfe view`` session covers a whole results directory. Nothing outside the
root is reachable.

Which objects can be viewed is :data:`LOADER_REGISTRY` (file extension -> loader)
composed with ``framejs.PAYLOAD_REGISTRY`` (gufe class -> serializer): a file type
is viewable when it loads to an object gufe knows how to feed the frame.
"""

from __future__ import annotations

import functools
import http.server
import json
import mimetypes
import posixpath
import sys
import threading
import urllib.parse
import webbrowser
from pathlib import Path
from typing import Any, Callable

from .framejs import FramejsUnavailable, payload_for_object, resolve_url

__all__ = [
    "DEFAULT_HOST",
    "DEFAULT_PORT",
    "LOADER_REGISTRY",
    "UnviewableFile",
    "frame_src_for_proxy",
    "hash_inputs",
    "load_object",
    "proxy_viewer_html",
    "serve",
    "viewer_html",
]

DEFAULT_HOST = "127.0.0.1"
DEFAULT_PORT = 8899

# The metapage runtime the page uses to host the frame and push it inputs. Pinned
# to the same version `metaframe_widget`'s anywidget ESM loads, so the terminal
# and the notebook drive a frame through identical code. Unused in proxy mode,
# where the frame feeds itself.
METAPAGE_MODULE = "https://cdn.jsdelivr.net/npm/@metapages/metapage@1.10.11/+esm"

# Where the re-hosted framejs.io document lives in proxy mode. A prefix rather
# than the root so it cannot collide with the data files, which keep the plain
# paths; the app's own root-relative requests (/sw.js, /favicon.ico, …) fall
# through to the catch-all proxy, so nothing in its HTML needs rewriting.
FRAMEJS_PROXY_PREFIX = "/_framejs/"

#: Upstream paths the proxy refuses rather than forwards.
#:
#: Only the service worker, and for two reasons: registered at this origin it
#: would take scope ``/`` over the *data* endpoints, whose whole contract is that
#: they are re-read from disk; and a service worker outlives the page and the
#: server, so the next thing to use this port would inherit it. framejs registers
#: it inside a try/catch and runs fine without one — it is an asset cache for
#: framejs.io itself, which is not what this origin is.
PROXY_DENY = frozenset({"/sw.js"})

# Payload values at or under this many bytes ride inline in the URL hash; bigger
# ones are replaced by a URL back to this server. Only reachable through
# `hash_inputs`, which the viewer page does not use — see its docstring.
INLINE_LIMIT = 200


class UnviewableFile(ValueError):
    """Raised for a file this server cannot turn into a viewable gufe object.

    Either the extension is not in :data:`LOADER_REGISTRY`, or it loaded to an
    object with no serializer in ``framejs.PAYLOAD_REGISTRY``.
    """


# --------------------------------------------------------------------------- #
# File extension -> gufe object                                                 #
# --------------------------------------------------------------------------- #


def _load_graphml(path: Path):
    from gufe import LigandNetwork

    return LigandNetwork.from_graphml(path.read_text())


def _load_pdb(path: Path):
    from gufe import ProteinComponent

    return ProteinComponent.from_pdb_file(str(path))


def _load_pdbx(path: Path):
    from gufe import ProteinComponent

    return ProteinComponent.from_pdbx_file(str(path))


def _load_sdf(path: Path):
    from gufe import SmallMoleculeComponent

    return SmallMoleculeComponent.from_sdf_file(str(path))


def _load_json(path: Path):
    """Load any gufe keyed-chain JSON — the ``.to_json()`` of any tokenizable.

    This is the interesting one for openfe users: ``AlchemicalNetwork``,
    ``Transformation`` and ``ChemicalSystem`` all reach disk this way, and the
    concrete class comes back out of the file, so one loader covers them all.
    """
    from gufe.tokenization import GufeTokenizable

    return GufeTokenizable.from_json(file=path)


#: File extension -> loader. The other half of "what can be viewed"; see the
#: module docstring. Keys are lowercase and include the dot.
LOADER_REGISTRY: dict[str, Callable[[Path], Any]] = {
    ".graphml": _load_graphml,
    ".pdb": _load_pdb,
    ".cif": _load_pdbx,
    ".pdbx": _load_pdbx,
    ".sdf": _load_sdf,
    ".mol": _load_sdf,
    ".json": _load_json,
}


def load_object(path: Path):
    """Load the gufe object stored in ``path``, dispatching on its extension.

    Raises :class:`UnviewableFile` if the extension is not registered, and lets
    the loader's own exception through if the file is registered but malformed.
    """
    loader = LOADER_REGISTRY.get(path.suffix.lower())
    if loader is None:
        raise UnviewableFile(
            f"don't know how to view {path.name!r} (unrecognized extension "
            f"{path.suffix!r}); supported: {', '.join(sorted(LOADER_REGISTRY))}"
        )
    return loader(path)


def is_viewable(path: Path) -> bool:
    """True if ``path``'s extension has a loader (cheap — the file is not read)."""
    return path.suffix.lower() in LOADER_REGISTRY


def _viz_url_for(path: Path) -> tuple[Any, str]:
    """Load ``path`` and return ``(object, framejs URL of the viz frame)``.

    The URL carries the viz JavaScript but *no* object data — the page fetches
    that separately. It is the same URL whatever the object is (one frame draws
    them all), but serializing here still makes an object nothing is registered
    for fail loudly rather than render an empty frame. Raises
    :class:`UnviewableFile` if nothing can render it.
    """
    obj = load_object(path)
    try:
        payload_for_object(obj)
        return obj, resolve_url()
    except FramejsUnavailable as e:
        raise UnviewableFile(f"no visualization available for {type(obj).__name__}: {e}") from e


def payload_for(path: Path) -> dict[str, Any]:
    """The frame ``inputs`` for the object in ``path`` — what ``?inputs=1`` returns.

    Identical to what the notebook widget pushes over the comm channel for the
    same object; the frame cannot tell the two apart.
    """
    obj = load_object(path)
    try:
        return payload_for_object(obj)
    except FramejsUnavailable as e:
        raise UnviewableFile(f"no visualization available for {type(obj).__name__}: {e}") from e


# --------------------------------------------------------------------------- #
# Payload -> URL hash, with the big values left here to be fetched              #
# --------------------------------------------------------------------------- #


def _fits_inline(value: Any) -> bool:
    """True if ``value`` is small enough to ride in the URL hash as itself."""
    if value is None or isinstance(value, (bool, int, float)):
        return True
    text = value if isinstance(value, str) else json.dumps(value, separators=(",", ":"))
    return len(text) <= INLINE_LIMIT


def hash_inputs(payload: dict[str, Any], *, data_url: str) -> dict[str, Any]:
    """Rewrite ``payload`` for the URL hash, replacing big values with URLs.

    ``data_url`` is the absolute URL of the data file on this server; each
    replaced value becomes ``<data_url>?input=<key>``, which serves exactly that
    value. The frames accept either form, so this is invisible to them.

    .. warning::
       **The viewer page does not use this**, because a frame served from
       framejs.io can no longer fetch a loopback address — see the module
       docstring on Local Network Access. It is kept for the case the frames were
       taught to handle, and which does work: inputs pointing at data the frame
       is allowed to reach (a public ``https://`` URL, or same-origin data when
       the frame is hosted alongside it).
    """
    return {
        key: value if _fits_inline(value) else f"{data_url}?input={urllib.parse.quote(key)}"
        for key, value in payload.items()
    }


# --------------------------------------------------------------------------- #
# The page: one metaframe, fed by a same-origin fetch                           #
# --------------------------------------------------------------------------- #

# Everything variable rides in a single JSON <script> block, so the page needs no
# templating of its HTML/CSS/JS (and nothing user-controlled is ever interpolated
# into executable text).
_VIEWER_HTML = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>gufe viz</title>
<style>
  html, body { height: 100%; margin: 0; }
  body { display: flex; flex-direction: column; font: 13px/1.4 ui-sans-serif, system-ui, sans-serif;
         color: #222; background: #fff; }
  header { flex: 0 0 auto; display: flex; align-items: center; gap: 10px; padding: 6px 10px;
           border-bottom: 1px solid #e2e2e2; background: #fafafa; }
  header .path { font-family: ui-monospace, monospace; font-weight: 600; }
  header .type { color: #777; }
  header .spacer { flex: 1 1 auto; }
  header button, header a { font: inherit; color: #555; text-decoration: none; cursor: pointer;
                            background: none; border: 1px solid #d5d5d5; border-radius: 4px;
                            padding: 2px 8px; }
  header button:hover, header a:hover { background: #efefef; }
  #frame { flex: 1 1 auto; min-height: 0; }
  #frame iframe { border: 0; display: block; width: 100%; height: 100%; }
  #frame > div, #frame .metapage-metaframe { height: 100%; }
  #error { display: none; padding: 10px 12px; background: #fff4f4; border-bottom: 1px solid #f0c0c0;
           color: #a00; white-space: pre-wrap; font-family: ui-monospace, monospace; }
  @media (prefers-color-scheme: dark) {
    body { color: #ddd; background: #1b1b1b; }
    header { background: #242424; border-bottom-color: #3a3a3a; }
    header button, header a { color: #bbb; border-color: #444; }
    header button:hover, header a:hover { background: #333; }
    #error { background: #2a1a1a; border-bottom-color: #5a2a2a; color: #ff9a9a; }
  }
</style>
</head>
<body>
<header>
  <span class="path"></span>
  <span class="type"></span>
  <span class="spacer"></span>
  <button id="reload" title="re-read the file and re-render, without reloading the frame">&#8635; data</button>
  <a id="raw" href="" download>file</a>
  <a id="up" href="">directory</a>
</header>
<div id="error"></div>
<div id="frame"></div>
<script id="gufe-viz-config" type="application/json">__CONFIG__</script>
<script type="module">
const cfg = JSON.parse(document.getElementById("gufe-viz-config").textContent);
const err = document.getElementById("error");
const fail = (e) => { err.textContent = String(e && e.stack || e); err.style.display = "block"; };

document.title = cfg.path + " — gufe viz";
document.querySelector("header .path").textContent = cfg.path;
document.querySelector("header .type").textContent = cfg.object_type;
document.getElementById("raw").href = cfg.raw_url;
document.getElementById("up").href = cfg.parent_url;

// THIS page fetches the data, not the frame: we are same-origin with the server,
// and the frame (on framejs.io) is not allowed to reach a loopback address at all
// — see the Local Network Access note in the module docstring. Handing the result
// over is a postMessage, which no such policy applies to.
const fetchInputs = () => fetch(cfg.inputs_url, { cache: "no-store" })
  .then(async (r) => r.ok ? r.json() : Promise.reject(new Error(await r.text())));

try {
  const { renderMetapage } = await import(cfg.metapage_module);
  const metapage = await renderMetapage({
    definition: { version: "0.3", metaframes: { mf: { url: cfg.frame_url } } },
    rootDiv: document.getElementById("frame"),
  });
  const push = () => fetchInputs().then((inputs) => metapage.setInputs({ mf: inputs })).catch(fail);
  document.getElementById("reload").addEventListener("click", push);
  await push();
} catch (e) {
  fail(e);
}
</script>
</body>
</html>
"""

# Proxy mode's page: no script at all, because the frame feeds itself. Four
# substitutions, each HTML-escaped text or an attribute value.
_PROXY_VIEWER_HTML = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__PATH__ — gufe viz</title>
<style>
  html, body { height: 100%; margin: 0; }
  body { display: flex; flex-direction: column; font: 13px/1.4 ui-sans-serif, system-ui, sans-serif;
         color: #222; background: #fff; }
  header { flex: 0 0 auto; display: flex; align-items: center; gap: 10px; padding: 6px 10px;
           border-bottom: 1px solid #e2e2e2; background: #fafafa; }
  header .path { font-family: ui-monospace, monospace; font-weight: 600; }
  header .type { color: #777; }
  header .spacer { flex: 1 1 auto; }
  header a { font: inherit; color: #555; text-decoration: none; background: none;
             border: 1px solid #d5d5d5; border-radius: 4px; padding: 2px 8px; }
  header a:hover { background: #efefef; }
  iframe { flex: 1 1 auto; min-height: 0; border: 0; display: block; width: 100%; }
  @media (prefers-color-scheme: dark) {
    body { color: #ddd; background: #1b1b1b; }
    header { background: #242424; border-bottom-color: #3a3a3a; }
    header a { color: #bbb; border-color: #444; }
    header a:hover { background: #333; }
  }
</style>
</head>
<body>
<header>
  <span class="path">__PATH__</span>
  <span class="type">__TYPE__</span>
  <span class="spacer"></span>
  <a href="" title="re-read the file from disk and re-render">&#8635; reload</a>
  <a href="__RAW_URL__" download>file</a>
  <a href="__PARENT_URL__">directory</a>
</header>
<iframe src="__FRAME_SRC__" allow="clipboard-write; fullscreen"></iframe>
</body>
</html>
"""

_LISTING_HTML = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>gufe viz</title>
<style>
  body { font: 14px/1.6 ui-sans-serif, system-ui, sans-serif; margin: 2rem auto; max-width: 46rem;
         padding: 0 1rem; color: #222; background: #fff; }
  h1 { font-size: 1rem; font-family: ui-monospace, monospace; }
  ul { list-style: none; padding: 0; }
  li { font-family: ui-monospace, monospace; }
  a { text-decoration: none; }
  .dim { color: #999; }
  @media (prefers-color-scheme: dark) { body { color: #ddd; background: #1b1b1b; } .dim { color: #777; } }
</style>
</head>
<body>
__BODY__
</body>
</html>
"""


def viewer_html(
    *,
    path: str,
    frame_url: str,
    object_type: str,
    inputs_url: str,
    raw_url: str,
    parent_url: str,
) -> str:
    """Render the viewer page for one data file.

    ``frame_url`` is the framejs URL of the registered viz (JavaScript only, no
    object data); ``inputs_url`` is where the page fetches that data from.
    """
    config = {
        "path": path,
        "object_type": object_type,
        "frame_url": frame_url,
        "inputs_url": inputs_url,
        "raw_url": raw_url,
        "parent_url": parent_url,
        "metapage_module": METAPAGE_MODULE,
    }
    # `</` cannot appear inside a <script> block, whatever the JSON contains.
    blob = json.dumps(config).replace("</", "<\\/")
    return _VIEWER_HTML.replace("__CONFIG__", blob)


def proxy_viewer_html(
    *,
    path: str,
    frame_src: str,
    object_type: str,
    raw_url: str,
    parent_url: str,
) -> str:
    """Render the proxy-mode viewer page: a header and one iframe, no script.

    ``frame_src`` is the whole iframe URL — the re-hosted framejs document plus
    the object's ``inputs`` hash, as built by :func:`frame_src_for_proxy`.
    """
    html = _PROXY_VIEWER_HTML
    for token, value in (
        ("__PATH__", path),
        ("__TYPE__", object_type),
        ("__RAW_URL__", raw_url),
        ("__PARENT_URL__", parent_url),
        ("__FRAME_SRC__", frame_src),
    ):
        html = html.replace(token, _escape(value))
    return html


def frame_src_for_proxy(frame_url: str, inputs: dict[str, Any], *, base: str) -> str:
    """Rewrite a framejs.io viz URL to this server, with ``inputs`` in its hash.

    ``https://framejs.io/#?js=…`` becomes ``<base>/_framejs/#?js=…`` — same
    document, served from here, so the frame shares an origin with the data it
    is about to fetch. The ``/j/<uuid>`` form maps just as directly
    (``<base>/_framejs/j/<uuid>``) because the catch-all proxy resolves the path
    upstream unchanged.

    Appended ``inputs`` take priority over anything baked into the frame — the
    same contract :func:`~gufe.visualization.framejs.build_cli_url` relies on.
    """
    from .framejs import FRAMEJS_BASE, string_to_base64_string

    if not frame_url.startswith(FRAMEJS_BASE):
        raise FramejsUnavailable(f"cannot re-host a non-framejs.io viz URL: {frame_url}")
    rehosted = base + FRAMEJS_PROXY_PREFIX + frame_url[len(FRAMEJS_BASE):].lstrip("/")
    encoded = string_to_base64_string(json.dumps(inputs))
    return f"{rehosted}{'&' if '#?' in rehosted else '#?'}inputs={encoded}"


def _escape(text: str) -> str:
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")


def listing_html(url_path: str, entries: list[tuple[str, str, bool]]) -> str:
    """Render a directory listing: ``entries`` of ``(href, label, viewable)``."""
    items = "\n".join(
        f'  <li><a href="{_escape(href)}">{_escape(label)}</a></li>'
        if viewable
        else f'  <li class="dim">{_escape(label)}</li>'
        for href, label, viewable in entries
    )
    body = f"<h1>{_escape(url_path)}</h1>\n<ul>\n{items or '  <li class=dim>(nothing to view here)</li>'}\n</ul>"
    return _LISTING_HTML.replace("__BODY__", body)


# --------------------------------------------------------------------------- #
# The server                                                                    #
# --------------------------------------------------------------------------- #


class _VizHandler(http.server.BaseHTTPRequestHandler):
    """Serves, under ``root``: viewer pages, their payloads, and the raw files.

    Held to GET only and to paths inside ``root`` — this is a developer tool
    bound to localhost by default, but it still never hands out anything the
    user did not point it at.
    """

    server_version = "gufe-viz"
    protocol_version = "HTTP/1.1"

    def __init__(self, *args, root: Path, quiet: bool = False, proxy: bool = False, **kwargs):
        self.root = root
        self.quiet = quiet
        self.proxy = proxy
        super().__init__(*args, **kwargs)

    # -- plumbing ---------------------------------------------------------- #

    def log_message(self, fmt, *args):  # noqa: D102 - indented, and silenced by --quiet
        if not self.quiet:
            sys.stderr.write(f"  {fmt % args}\n")

    def _send_cors_headers(self) -> None:
        """Let anything that is *allowed* to reach this server read the response.

        The viewer page does not need this (it is same-origin), but a notebook,
        a script or a frame hosted alongside the data might, and the alternative
        would be a same-origin lock on a server that already binds loopback and
        serves only what it was pointed at.
        """
        self.send_header("Access-Control-Allow-Origin", "*")

    def _send(self, body: bytes, content_type: str, status: int = 200) -> None:
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self._send_cors_headers()
        self.end_headers()
        self.wfile.write(body)

    def _send_text(self, text: str, content_type: str, status: int = 200) -> None:
        self._send(text.encode("utf-8"), f"{content_type}; charset=utf-8", status)

    def _send_error_page(self, message: str, status: int) -> None:
        self._send_text(listing_html(f"{status}", [("", message, False)]), "text/html", status)

    def _local_path(self, url_path: str) -> Path | None:
        """Map a URL path to a file under ``root``, or ``None`` if it escapes it."""
        # posixpath.normpath collapses `..` before we ever touch the filesystem;
        # resolve() then catches symlinks that point out of the tree.
        relative = posixpath.normpath(urllib.parse.unquote(url_path)).lstrip("/")
        if relative in (".", ""):
            return self.root
        candidate = (self.root / relative).resolve()
        if candidate != self.root and self.root not in candidate.parents:
            return None
        return candidate

    # -- routing ------------------------------------------------------------ #

    def do_GET(self):  # noqa: N802 - stdlib naming
        split = urllib.parse.urlsplit(self.path)
        query = urllib.parse.parse_qs(split.query)

        if self.proxy and split.path.startswith(FRAMEJS_PROXY_PREFIX):
            return self._serve_proxied("/" + split.path[len(FRAMEJS_PROXY_PREFIX):], split.query)

        path = self._local_path(split.path)

        if path is None:
            return self._send_error_page("outside the served directory", 403)
        if not path.exists():
            # In proxy mode an unknown path is most likely the re-hosted app
            # asking for one of its own root-relative files (/sw.js, favicons).
            if self.proxy:
                return self._serve_proxied(split.path, split.query)
            return self._send_error_page(f"no such file: {split.path}", 404)
        if path.is_dir():
            return self._serve_listing(path, split.path)
        if "raw" in query:
            return self._serve_raw(path)
        if "input" in query:
            return self._serve_input(path, query["input"][0])
        if "inputs" in query:
            return self._serve_inputs(path)
        return self._serve_viewer(path, split.path)

    def do_OPTIONS(self):  # noqa: N802 - stdlib naming
        """Answer CORS preflights, for anything fetching these responses directly.

        Note this is *not* what makes the viewer page work — that fetch is
        same-origin. Nor can it re-enable a fetch from the framejs.io frame:
        Local Network Access is a permission the browser asks the user for, not
        something a response header can grant (``Access-Control-Allow-Private-Network``
        was that header, and it is superseded).
        """
        self.send_response(204)
        self.send_header("Access-Control-Allow-Methods", "GET, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", self.headers.get("Access-Control-Request-Headers", "*"))
        self.send_header("Access-Control-Max-Age", "600")
        self.send_header("Content-Length", "0")
        self._send_cors_headers()
        self.end_headers()

    def _serve_viewer(self, path: Path, url_path: str) -> None:
        try:
            obj, frame_url = _viz_url_for(path)
            payload = payload_for_object(obj) if self.proxy else None
        except UnviewableFile as e:
            return self._send_error_page(str(e), 415)
        except Exception as e:  # a malformed file is the user's, not a server fault
            return self._send_error_page(f"could not read {path.name}: {type(e).__name__}: {e}", 422)

        quoted = urllib.parse.quote(url_path)
        common = dict(
            path=path.relative_to(self.root).as_posix(),
            object_type=type(obj).__name__,
            raw_url=f"{quoted}?raw=1",
            parent_url=posixpath.dirname(quoted.rstrip("/")) or "/",
        )

        if self.proxy:
            # Frame, page and data all on this origin: the frame can fetch its
            # own inputs, so the page carries no script.
            base = self._origin()
            html = proxy_viewer_html(
                frame_src=frame_src_for_proxy(
                    frame_url, hash_inputs(payload, data_url=base + quoted), base=base
                ),
                **common,
            )
        else:
            html = viewer_html(frame_url=frame_url, inputs_url=f"{quoted}?inputs=1", **common)
        self._send_text(html, "text/html")

    def _origin(self) -> str:
        """This server's origin as the client sees it (``http://localhost:8899``).

        Taken from the request's ``Host`` header rather than the bound address,
        so a container that publishes its port yields the *browser's* address —
        the bound ``0.0.0.0`` would be useless to it, and getting this wrong is
        exactly what re-introduces a cross-origin fetch.
        """
        host = self.headers.get("Host") or f"{self.server.server_address[0]}:{self.server.server_address[1]}"
        return f"http://{host}"

    def _serve_proxied(self, upstream_path: str, query: str) -> None:
        """Serve ``https://framejs.io<upstream_path>`` from this origin.

        Byte-for-byte: framejs's application code is inline in its HTML and its
        dependencies are absolute ``https://`` imports, so nothing needs
        rewriting — only the document's *origin* has to change, which is the
        whole point (see the module docstring).
        """
        from .framejs import FRAMEJS_BASE

        if upstream_path in PROXY_DENY:
            return self._send_text(f"{upstream_path} is not proxied (see PROXY_DENY)", "text/plain", 404)
        url = FRAMEJS_BASE + upstream_path + (f"?{query}" if query else "")
        try:
            body, content_type = _fetch_upstream(url)
        except OSError as e:
            return self._send_error_page(f"could not reach {url}: {e}", 502)
        self._send(body, content_type)

    def _serve_input(self, path: Path, key: str) -> None:
        """Serve one payload value — the endpoint a frame fetches a big input from."""
        try:
            payload = payload_for(path)
        except UnviewableFile as e:
            return self._send_text(str(e), "text/plain", 415)
        except Exception as e:
            return self._send_text(f"could not read {path.name}: {type(e).__name__}: {e}", "text/plain", 422)
        if key not in payload:
            return self._send_text(f"no input {key!r} in the payload for {path.name}", "text/plain", 404)
        value = payload[key]
        if isinstance(value, str):
            return self._send_text(value, "text/plain")
        self._send_text(json.dumps(value), "application/json")

    def _serve_inputs(self, path: Path) -> None:
        try:
            payload = payload_for(path)
        except UnviewableFile as e:
            return self._send_text(str(e), "text/plain", 415)
        except Exception as e:
            return self._send_text(f"could not read {path.name}: {type(e).__name__}: {e}", "text/plain", 422)
        self._send_text(json.dumps(payload), "application/json")

    def _serve_raw(self, path: Path) -> None:
        guessed, _ = mimetypes.guess_type(path.name)
        self._send(path.read_bytes(), guessed or "application/octet-stream")

    def _serve_listing(self, path: Path, url_path: str) -> None:
        prefix = url_path if url_path.endswith("/") else url_path + "/"
        entries: list[tuple[str, str, bool]] = []
        if path != self.root:
            entries.append((posixpath.dirname(url_path.rstrip("/")) or "/", "../", True))
        for child in sorted(path.iterdir(), key=lambda p: (p.is_file(), p.name)):
            if child.name.startswith("."):
                continue
            href = urllib.parse.quote(prefix + child.name)
            if child.is_dir():
                entries.append((href + "/", child.name + "/", True))
            elif is_viewable(child):
                entries.append((href, child.name, True))
            else:
                entries.append(("", child.name, False))
        here = path.relative_to(self.root).as_posix()
        self._send_text(listing_html("/" if here in (".", "") else here, entries), "text/html")


#: Proxied upstream responses, keyed by URL. The framejs document is ~80 kB and
#: is re-fetched on every viewer page load, so this is worth having; it is a
#: session cache with no expiry, which is what `--proxy` being a dev-loop feature
#: makes appropriate (restart the server to pick up a framejs.io release).
_UPSTREAM_CACHE: dict[str, tuple[bytes, str]] = {}


def _fetch_upstream(url: str, *, timeout: float = 20.0) -> tuple[bytes, str]:
    """GET ``url``, returning ``(body, content-type)``. Cached for the session."""
    import urllib.request

    hit = _UPSTREAM_CACHE.get(url)
    if hit is not None:
        return hit
    with urllib.request.urlopen(url, timeout=timeout) as response:
        result = (response.read(), response.headers.get("Content-Type", "application/octet-stream"))
    _UPSTREAM_CACHE[url] = result
    return result


def _root_and_url_path(target: Path) -> tuple[Path, str]:
    """Pick what to serve and the URL path of ``target`` within it.

    The working directory when ``target`` is under it — so the URL reads like the
    path the user typed and its sibling files are reachable from the same server
    — otherwise the target's own directory, so an absolute path elsewhere never
    exposes the tree above it.
    """
    target = target.resolve()
    base = target if target.is_dir() else target.parent
    cwd = Path.cwd().resolve()
    root = cwd if (base == cwd or cwd in base.parents) and cwd.parent != cwd else base
    relative = target.relative_to(root).as_posix()
    return root, "/" + urllib.parse.quote(relative if relative != "." else "")


def _bind(host: str, port: int, handler) -> http.server.ThreadingHTTPServer:
    """Bind a threaded server, walking forward from ``port`` if it is taken."""
    for candidate in range(port, port + 20):
        try:
            return http.server.ThreadingHTTPServer((host, candidate), handler)
        except OSError as e:
            if candidate == port + 19 or e.errno not in (48, 98, 10048):  # EADDRINUSE
                raise
            if candidate == port:
                print(f"  port {port} is in use, trying {port + 1}…", file=sys.stderr)
    raise AssertionError("unreachable")


def serve(
    target: str | Path,
    *,
    host: str = DEFAULT_HOST,
    port: int = DEFAULT_PORT,
    open_browser: bool = True,
    quiet: bool = False,
    proxy: bool = False,
) -> None:
    """Serve the visualization of ``target`` and block until Ctrl-C.

    Parameters
    ----------
    target
        A data file to view, or a directory to list. Its directory (or the
        working directory, when ``target`` is under it) becomes the server root:
        every viewable file under that root is reachable from this one server.
    host
        Interface to bind. The default is loopback-only; pass ``0.0.0.0`` to make
        the server reachable from outside a container.
    port
        Preferred port; the next free one is used (and announced) if it is taken.
    open_browser
        Try to open the URL in a browser. A no-op where there is no browser (a
        container, a remote shell) — the URL is always printed either way.
    quiet
        Suppress the per-request log lines.
    proxy
        Re-host the framejs.io document at ``/_framejs/`` so the frame shares this
        origin, and let it fetch its own inputs instead of the page pushing them.
        Removes the Local Network Access barrier described in the module
        docstring, and leaves the page script-free. See there for the caveat.
    """
    target = Path(target)
    if not target.exists():
        raise FileNotFoundError(target)
    root, url_path = _root_and_url_path(target)

    handler = functools.partial(_VizHandler, root=root, quiet=quiet, proxy=proxy)
    httpd = _bind(host, port, handler)
    bound_port = httpd.server_address[1]
    # 0.0.0.0 is not an address a browser can go to; localhost is, and that is
    # also the name that works through `docker run -p`.
    displayed_host = "localhost" if host in ("0.0.0.0", "::", "") else host
    url = f"http://{displayed_host}:{bound_port}{url_path}"

    print(f"gufe viz  →  {url}")
    print(f"  serving {root}{'' if host == DEFAULT_HOST else f' on {host}:{bound_port}'}")
    if proxy:
        print(f"  proxying framejs.io at {FRAMEJS_PROXY_PREFIX} (frame and data share this origin)")
    print("  open the URL above (Ctrl-C to stop the server)", flush=True)

    if open_browser:
        threading.Thread(target=webbrowser.open, args=(url,), daemon=True).start()

    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nstopped.")
    finally:
        httpd.shutdown()
        httpd.server_close()


def main(argv: list[str] | None = None) -> int:
    """``python -m gufe.visualization.server`` — the CLI-free entry point.

    openfe wraps this as ``openfe view``; it lives here so a gufe-only install
    (and the docker demo, which has no openfe) can use it too.
    """
    import argparse

    parser = argparse.ArgumentParser(
        prog="python -m gufe.visualization.server",
        description="Serve an interactive framejs visualization of a gufe data file.",
    )
    parser.add_argument("target", help="data file to view, or a directory to list")
    parser.add_argument("--host", default=DEFAULT_HOST, help=f"interface to bind (default {DEFAULT_HOST})")
    parser.add_argument("--port", type=int, default=DEFAULT_PORT, help=f"port (default {DEFAULT_PORT})")
    parser.add_argument("--no-browser", action="store_true", help="only print the URL")
    parser.add_argument("--quiet", action="store_true", help="no per-request logging")
    parser.add_argument(
        "--proxy",
        action="store_true",
        help="re-host framejs.io on this origin so the frame fetches its own data "
        "(needed for Chrome, whose Local Network Access blocks the cross-origin form)",
    )
    args = parser.parse_args(argv)

    try:
        serve(
            args.target,
            host=args.host,
            port=args.port,
            open_browser=not args.no_browser,
            quiet=args.quiet,
            proxy=args.proxy,
        )
    except (FileNotFoundError, UnviewableFile, OSError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
