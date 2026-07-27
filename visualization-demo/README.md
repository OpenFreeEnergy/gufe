# gufe framejs visualization — demos

Two self-contained demos of the [framejs.io](https://framejs.io) visualizations
for gufe objects. Both run in Docker with plain `docker compose` — nothing is
installed on your host but Docker, and there is no `just`.

Everything builds one fast native-arm64 gufe image and editable-installs the
checkout one level up, so your host edits to `src/gufe/` — including the frame
itself, `src/gufe/visualization/viz_assets/gufe/code.js` — are live.

Run these from this directory (`visualization-demo/`).

---

## Demo 1 — the notebook implementation (`.view()`)

How a gufe object renders itself in Jupyter: `obj.view()`, or just an object on
the last line of a cell.

```bash
docker compose up notebook          # builds on first run (a few minutes), then cached
```

→ open **<http://localhost:8888/lab>** and run **`demo/framejs_demo.ipynb`**.

The notebook builds one of every visualizable gufe type from gufe's bundled test
data and renders each — small molecule, solvent, protein, atom mapping, chemical
system, transformation, ligand network, alchemical network. It needs only gufe
(no `openfe`, no external data). Edit the frame on your host, restart the kernel,
re-run a cell to see the change.

---

## Demo 2 — the `openfe view` implementation (the local server)

How a gufe **file** opens in a browser from the terminal. The service runs
exactly what `openfe view <file>` runs — the openfe CLI command is a thin wrapper
around `gufe.visualization.server`.

```bash
docker compose up openfe-view       # http://localhost:8899
```

→ open **<http://localhost:8899>**: a browsable listing of gufe's test fixtures.
Click any file to see its visualization.

To view one file (or your own directory) instead of the fixtures, set `VIZ_TARGET`
to a path inside the container (the gufe checkout is mounted at `/src/gufe`):

```bash
VIZ_TARGET=/src/gufe/src/gufe/tests/data/181l.pdb docker compose up openfe-view
```

The service runs with `--proxy`, which re-hosts the framejs viewer on this same
origin so it works in **Chrome** as well as Firefox — see
[`../src/gufe/visualization/README.md`](../src/gufe/visualization/README.md)
("Why Chrome needs `--proxy`") for why that flag exists. It also runs
`--host 0.0.0.0 --no-browser` because it
is in a container: it binds all interfaces so the published port is reachable,
and prints the URL rather than trying to open a browser that isn't there
(`docker compose logs openfe-view` shows it).

---

## Both at once, logs, teardown

```bash
docker compose up -d notebook openfe-view   # both, in the background
docker compose logs -f openfe-view          # follow one service's logs
docker compose down                         # stop everything (keeps the built image)
```

---

## Handy checks (plain docker compose)

`docker compose run --rm` starts a throwaway container without needing a service
up. Prefix the command with **`dev-start.sh`** so it editable-installs gufe from
the mount first (that is the service entrypoint; `run` replaces it otherwise):

```bash
# gufe + the widget import from the editable mount
docker compose run --rm notebook dev-start.sh python -c \
  "import gufe, metaframe_widget as m; print(gufe.__file__); print(m.__file__)"

# the hermetic framejs tests (no network) — keep these green
docker compose run --rm notebook dev-start.sh \
  python -m pytest /src/gufe/src/gufe/tests/test_framejs_visualization.py \
                   /src/gufe/src/gufe/tests/test_framejs_server.py -q

# a shell in the image (gufe already installed)
docker compose run --rm notebook dev-start.sh bash
```

If a service is already up, use `exec` instead — gufe is installed, so no
`dev-start.sh` is needed:

```bash
docker compose exec notebook python -m pytest /src/gufe/src/gufe/tests -k view -q
```

---

## Optional extras

**marimo** — the same `.view()` rendering in marimo, which renders anywidgets
natively where legacy py3Dmol output breaks:

```bash
docker compose --profile marimo up marimo   # http://localhost:2718
```

**Developing the widget itself** — to hack on `metaframe-widget` (the framejs.io
anywidget) instead of using the PyPI build, layer the widget override and point
`FRAMEJS_IO_REPO` at your `metaframe-js` checkout root:

```bash
FRAMEJS_IO_REPO=~/path/to/metaframe-js \
  docker compose -f docker-compose.yml -f docker-compose.widget.yml up notebook
```

**Editing the frame live in the browser** (no rebuild, auto-save to disk) uses
the framejs local server and needs only [`deno`](https://deno.com) — see the
"Editing the visualization" section of
[`../src/gufe/visualization/README.md`](../src/gufe/visualization/README.md).

---

## What's here

| Path | What it is |
|------|------------|
| `docker-compose.yml` | the two demo services (`notebook`, `openfe-view`) + optional `marimo` |
| `docker-compose.widget.yml` | override to develop `metaframe-widget` in-place |
| `Dockerfile`, `dev-start.sh` | the fast gufe image; the entrypoint that editable-installs the mount |
| `demo/framejs_demo.ipynb` | Demo 1's notebook |
| `demo/marimo_demo.py` | the optional marimo notebook |
| `frames/` | framejs frames that are **not** the gufe viz (e.g. the architecture diagram) |
| `scripts/framejs_frame.py` | round-trip a frame directory ↔ a framejs.io URL |
| `scratch/` | throwaway notebooks / experiments |
