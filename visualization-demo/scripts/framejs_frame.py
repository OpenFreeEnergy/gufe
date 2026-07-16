#!/usr/bin/env python3
# This code is part of gufe and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
"""restore/save a framejs.io page <-> a local frame directory (stdlib only).

A framejs.io page is fully described by its URL hash — a self-contained,
clean-room implementation of the framejs *local-file-io on-disk format*
(https://framejs.io/docs/guide/local-file-io). It vendors **no** framejs.io code;
it just reads/writes the documented on-disk layout and the documented hash-param
wire format, exactly as :mod:`gufe.visualization.framejs` does at runtime:

    <dir>/code.js       <->  js       hash param   (raw JavaScript text)
    <dir>/og.json       <->  og       hash param   (JSON)
    <dir>/modules.json  <->  modules  hash param   (JSON: classic-script URLs)
    <dir>/inputs.json / options.json / definition.json           (JSON)

Each hash value is ``base64(encodeURIComponent(text))`` — raw text for ``js``,
compact JSON for the rest.

Usage (host-side, needs only python3):
    python3 scripts/framejs_frame.py restore <dir> [base-url]   # dir -> prints URL
    python3 scripts/framejs_frame.py save    <url> <dir>        # URL -> dir/*
"""

from __future__ import annotations

import base64
import json
import sys
import urllib.parse
from pathlib import Path

BASE_URL_DEFAULT = "https://framejs.io/"

# key -> on-disk filename.  `js` is raw text; the rest are JSON. Order matters
# only for a stable URL (mirrors the framejs local server's param order).
JS_KEY = "js"
JSON_KEYS = ("options", "inputs", "definition", "og", "modules")
FILE = {"js": "code.js", **{k: f"{k}.json" for k in JSON_KEYS}}


def _b64_encode(text: str) -> str:
    """text -> base64(encodeURIComponent(text)) — the framejs hash encoding."""
    quoted = urllib.parse.quote(text, safe="-_.!~*'()")
    return base64.b64encode(quoted.encode("ascii")).decode("ascii")


def _b64_decode(value: str) -> str:
    """base64(encodeURIComponent(text)) -> text."""
    return urllib.parse.unquote(base64.b64decode(value).decode("ascii"))


def _hash_params(url: str) -> dict[str, str]:
    frag = url.split("#", 1)[1] if "#" in url else ""
    frag = frag[1:] if frag.startswith("?") else frag
    params: dict[str, str] = {}
    for kv in frag.split("&"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            params[k] = v
    return params


def restore(directory: str, base: str = BASE_URL_DEFAULT) -> str:
    """Build the framejs.io URL for a frame directory (offline; no server)."""
    d = Path(directory)
    if not d.is_dir():
        raise SystemExit(f"error: no such dir: {directory}")
    parts: list[str] = []
    code = d / FILE[JS_KEY]
    if code.is_file():
        parts.append(f"{JS_KEY}={_b64_encode(code.read_text(encoding='utf-8'))}")
    for key in JSON_KEYS:
        f = d / FILE[key]
        if not f.is_file():
            continue
        value = json.loads(f.read_text(encoding="utf-8"))
        parts.append(f"{key}={_b64_encode(json.dumps(value, separators=(',', ':')))}")
    if not parts:
        raise SystemExit(f"error: no frame files found in {directory}")
    return f"{base.rstrip('/')}/#?" + "&".join(parts)


def save(url: str, directory: str) -> None:
    """Decode a framejs.io URL into a frame directory (captures browser edits)."""
    d = Path(directory)
    d.mkdir(parents=True, exist_ok=True)
    params = _hash_params(url)
    written: list[str] = []
    if JS_KEY in params:
        (d / FILE[JS_KEY]).write_text(_b64_decode(params[JS_KEY]), encoding="utf-8")
        written.append(FILE[JS_KEY])
    for key in JSON_KEYS:
        if key not in params:
            continue
        value = json.loads(_b64_decode(params[key]))
        # pretty, sorted keys for clean diffs
        (d / FILE[key]).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        written.append(FILE[key])
    if not written:
        raise SystemExit(f"error: no framejs hash params found in URL")
    print(f"saved -> {directory}", file=sys.stderr)
    print("\n".join(written), file=sys.stderr)


def main(argv: list[str]) -> None:
    match argv:
        case ["restore", directory]:
            print(restore(directory))
        case ["restore", directory, base]:
            print(restore(directory, base))
        case ["save", url, directory]:
            save(url, directory)
        case _:
            raise SystemExit(
                "usage: framejs_frame.py {restore <dir> [base-url] | save <url> <dir>}"
            )


if __name__ == "__main__":
    main(sys.argv[1:])
