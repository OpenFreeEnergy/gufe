# gufe visualization

Interactive [framejs.io](https://framejs.io) visualizations for gufe objects.
`LigandNetwork.view()` (and bare-cell auto-display in Jupyter / marimo / VSCode)
renders an interactive radial network — click an edge to drive a 3D atom-mapping
viewer for its two endpoint molecules. Rendering is powered by the optional
`viz` extra (`pip install gufe[viz]`); without it, `.view()` falls back to the
legacy RDKit / py3Dmol renderers.

- `framejs.py` — builds the self-contained framejs.io URL / `anywidget` from a
  gufe object plus an on-disk frame directory.
- `mapping_visualization.py` — the legacy RDKit atom-mapping renderer.
- `viz_assets/<frame>/` — the on-disk framejs frames that ship with gufe.

## Running the Jupyter notebook demo

The demo stack (`visualization-demo/`) is fully driven by Docker Compose — no
`just` required. It builds a native gufe env, editable-installs this checkout,
and serves a JupyterLab notebook that exercises `LigandNetwork.view()`.

From the repository root:

```bash
cd visualization-demo

# build the image + start JupyterLab in the background
docker compose up -d --build jupyter

# → open http://localhost:8888/lab and run demo/framejs_demo.ipynb
```

The mounted gufe source is editable-installed, so host edits to `src/gufe/` (or
to a `viz_assets/` frame) are live — just restart the notebook kernel.

Optional extras:

```bash
docker compose up -d --build marimo   # marimo demo at http://localhost:2718
docker compose logs -f jupyter        # follow logs
docker compose down                   # stop the stack (keeps built images)
```


## Editing the visualizations

Each visualization is a framejs **frame directory** under `viz_assets/`, in the
canonical [framejs local-file-io](https://framejs.io/docs/guide/local-file-io)
format. For example `viz_assets/network/`:

| file           | contents                                                        |
| -------------- | --------------------------------------------------------------- |
| `code.js`      | the frame's JavaScript (the actual visualization)               |
| `og.json`      | title + description metadata                                    |
| `modules.json` | external classic-script URLs the frame loads (e.g. 3Dmol.js)    |

To edit a frame live in the browser, run the framejs **local server**, which
watches these files and auto-saves your edits back to disk. It needs only
[`deno`](https://docs.deno.com/runtime/getting_started/installation/) — the
canonical command from
[the docs](https://framejs.io/docs/guide/local-file-io#run-it-deno) is:

```bash
deno run --allow-read=$PWD --allow-write=$PWD --allow-net --allow-env \
  https://raw.githubusercontent.com/metapages/framejs.io/main/local-server/server.ts \
  --root "$PWD" --port 4700
```

Point `--root` at this package's `viz_assets/` so the server sees the shipped
frames. From the repository root:

```bash
deno run --allow-read="$PWD/src/gufe/visualization/viz_assets" \
         --allow-write="$PWD/src/gufe/visualization/viz_assets" \
         --allow-net --allow-env \
  https://raw.githubusercontent.com/metapages/framejs.io/main/local-server/server.ts \
  --root "$PWD/src/gufe/visualization/viz_assets" --port 4700
```

Then open <http://localhost:4700>, pick a frame (dirs with a `code.js` show as
`◆ FRAME`), edit, and watch it save (`saved ✓`). gufe's `.view()` builds its URL
straight from this on-disk source, so your edits show up on the next render.

> The `visualization-demo/justfile` wraps this as `just viz-edit` (plus
> `just viz-url` / `just viz-save` to round-trip a frame ↔ a framejs.io URL), but
> the raw `deno` command above is all you need.

