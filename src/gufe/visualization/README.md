# gufe visualization

Interactive [framejs.io](https://framejs.io) visualizations for gufe objects.
Calling `.view()` — or just leaving an object as the last line of a notebook cell
in Jupyter / marimo / VSCode — renders it as an interactive widget. Rendering is
powered by the optional `viz` extra (`pip install gufe[viz]`); without it,
`.view()` falls back to the legacy RDKit / py3Dmol renderers.

> **Architecture, interactively:**
> <https://framejs.app/j/019f8b367fac7e469f26b99bb30e06e4> — both render paths
> and the registry table below, clickable. Source:
> `visualization-demo/frames/architecture/` (`just viz-edit frames` to edit it).

- `framejs.py` — `VIZ_REGISTRY` maps each gufe class to a `VizRef` holding both
  halves of the contract: the frame that draws it and the serializer that turns
  an object into that frame's `inputs`. It then builds the framejs.io URL /
  `anywidget`.
- `../_viewable.py` — the `FramejsViewable` mixin that gives a class `.view()` and
  `_repr_mimebundle_`. Mixing it in is the whole opt-in.
- `mapping_visualization.py` — the legacy RDKit atom-mapping renderer, reached
  through `LigandAtomMapping._legacy_view()` when framejs is unavailable. It is
  the only legacy renderer left.
- `server.py` — the terminal path: a small local web server that shows a **file's**
  visualization in a browser (`python -m gufe.visualization.server <file>`, which
  is what `openfe view <file>` runs). Its `LOADER_REGISTRY` maps a file extension
  to the gufe object it loads to; `VIZ_REGISTRY` takes it from there.
- `viz_assets/<frame>/` — the on-disk framejs frames that ship with gufe.
- `../tests/test_framejs_visualization.py`, `../tests/test_framejs_server.py` —
  hermetic (no network) tests for all of the above, including the invariants
  stated here.

Adding a visualization is one `VIZ_REGISTRY` entry plus one frame directory. A
class must **not** define `_ipython_display_`: IPython checks that hook before
`_repr_mimebundle_` and short-circuits on it, so a bare cell and `.view()` would
disagree. Define `_legacy_view()` instead for the non-framejs fallback.

This package reads **no environment variables**. What renders, and from where, is
determined entirely by `VIZ_REGISTRY` and the caller's arguments.

Payload keys are the contract with each frame's `onInputs`. File-shaped values
use a descriptive `<thing>.<ext>` key (`molecule.sdf`, `protein.pdb`,
`network.graphml`) — frames read these through `asText()`, which also accepts a
dropped file. Everything else is a bare snake_case field, and where a viz takes
one structured object that object's key is the frame name (`chemical_system`,
`transformation`, `alchemical_network`, `solvent_component`). Renaming a key
breaks any frame already published to a `/j/<uuid>`, which is pinned and keeps
reading the old name.

## What renders

A frame directory is named after the gufe class it renders, snake_cased, so the
two are greppable from each other.

| gufe class | frame | shows |
| --- | --- | --- |
| `LigandNetwork` | `ligand_network` | radial network; click an edge to drive a 3D atom-mapping viewer |
| `AlchemicalNetwork` | `alchemical_network` | d3 force graph of ChemicalSystem nodes / Transformation edges |
| `Transformation`, `NonTransformation` | `transformation` | stateA↔stateB component diff + the atom mapping |
| `ChemicalSystem` | `chemical_system` | master/detail over the system's labelled components |
| `LigandAtomMapping` | `ligand_atom_mapping` | the mapping viewer standalone (plain/colored/lines/overlay/2d) |
| `SmallMoleculeComponent` | `small_molecule_component` | 2D depiction + 3D conformer + SMILES/charge |
| `ProteinComponent` | `protein_component` | 3Dmol with representation / colour-scheme switchers |
| `SolventComponent` | `solvent_component` | settings card (it has no coordinates) |

Lookup walks the MRO, so subclasses inherit their base's viz: `SolvatedPDBComponent`
and `ProteinMembraneComponent` get the protein frame, with nothing registered for
them. The transformation entry is registered on `TransformationBase` rather than
on `Transformation`, because `Transformation` and `NonTransformation` are siblings
under it, not parent and child.

To add one: mix `FramejsViewable` into the class, add a `VizRef` and a payload
builder in `framejs.py`, and drop a frame directory under `viz_assets/`.

## Render paths, one core

Every path looks the object up in the same `VIZ_REGISTRY`, gets the same `VizRef`,
builds the same base URL from the same on-disk frame, and serializes the object
with the same `VizRef.payload`. They differ **only** in how that payload reaches
the frame.

| | Notebook | Terminal | Shareable link |
| --- | --- | --- | --- |
| entry point | `.view()` / bare cell → `view_object()` / `repr_mimebundle()` | `server.serve(path)` → `openfe view <file>` | `build_cli_url(obj)` → `openfe view --url-only` |
| carrier | `MetaframeWidget.set_inputs()` — live over the widget comm channel | `fetch()` from the served page → `metapage.setInputs()` | `…&inputs=<base64(json)>` appended to the URL hash |
| size limit | none | none | the URL's |
| result | an iframe inline in the cell, updatable in place | a page at `http://localhost:8899/<the file's path>` | one link you can paste anywhere |

The first two are the same design — load the frame, then push data into it — which
is why a frame does not care which one is hosting it, and why neither has a size
limit. The third bakes the object into the link instead: no server, but the data
is frozen at the moment the link was made, and a solvated `AlchemicalNetwork` will
not fit in a URL. Appended `inputs` take priority over anything baked into the
frame, which is what makes that form work at all.

### Two URL forms

Both are built from the one on-disk frame directory:

- **local** — a self-contained `https://framejs.io/#?js=<b64>&og=<b64>…`. **The
  default everywhere.** It always matches the installed gufe, needs no framejs
  account, and cannot expire.
- **canonical** — the pinned short `https://framejs.io/j/<uuid>`, minted once by
  `just publish-viz`. Its one advantage is size, so it is opt-in exactly where
  size can matter: `build_cli_url(obj, short=True)`. For a `LigandNetwork` that is
  a ~10 kB URL rather than ~140 kB — the difference between a link you can paste
  somewhere and one you cannot. It requires that viz to have been published;
  only `ligand_network` has been so far.

`VizRef.resolve_url()` prefers the on-disk frame and falls back to the pinned
uuid only if no frame directory is present.

### The terminal path in practice

`openfe`'s viewer commands are thin wrappers over `server.serve`:

```bash
openfe view network_setup/ligand_network.graphml   # serve it, open a browser, wait for Ctrl-C
openfe view results/                               # a listing: browse every viewable file
openfe view ligand.sdf --no-browser                # only print the URL (ssh, container)
openfe view ligand.sdf --url-only                  # no server: one self-contained link
openfe view-ligand-network network.graphml         # the older .graphml-only entry point
```

What gets served at `http://localhost:8899/<path/to/the/file>` is a page Python
renders per file: a header, one metaframe, and ~30 lines of JavaScript. The
framejs URL of the registered viz is **baked into the page**; the page then
fetches the object's `inputs` from `?inputs=1` on its own path and hands them to
the frame through the same `@metapages/metapage` entry point, pinned to the same
version, that `metaframe_widget` uses in Jupyter.

The URL path *is* the data file's path, relative to the served root — the
directory you pointed at, or the working directory if the file is under it.
Everything viewable under that root is reachable from the one server, and nothing
outside it is. The payload is rebuilt from disk on every fetch, so regenerating a
file and pressing the page's ⟳ re-renders it without restarting anything.

`LOADER_REGISTRY` is what makes a file type viewable: `.json` (any saved gufe
object), `.graphml`, `.pdb`, `.cif`/`.pdbx`, `.sdf`/`.mol`. A file it cannot load,
or one that loads to an object with no registered viz, produces an explanatory
page (415/422) rather than a traceback.

In a container, bind `--host 0.0.0.0` and publish the port: loopback inside the
container is not reachable from the host, and there is no browser in there to
open. The URL is printed either way.

### When framejs is unavailable

Everything here is optional — `import gufe` works without `metaframe-widget`. On
`FramejsUnavailable` or `OSError` (no `viz` extra, nothing registered, no frame on
disk), both paths degrade: the object falls back to `_legacy_view()` if it defines
one, and otherwise to its plain `repr`. `.view()` warns when it falls back;
auto-display stays silent, because it runs on every display of the object.

## Running the Jupyter notebook demo

The demo stack (`visualization-demo/`) builds a native gufe env, editable-installs
this checkout, and serves a JupyterLab notebook that exercises `.view()` on one of
every visualizable gufe object.

```bash
cd visualization-demo
just dev     # build (first time) + start JupyterLab :8888 and marimo :2718
```

→ open <http://localhost:8888/lab> and run `demo/framejs_demo.ipynb`.

`just dev` also clones the pinned `openfe-demo` locally (git-ignored) so the demo
works from a fresh clone. If you would rather not install `just`, the stack is
plain Docker Compose:

```bash
docker compose up -d --build jupyter marimo   # same thing, minus the openfe-demo clone
docker compose logs -f jupyter                # follow logs
docker compose down                           # stop (keeps built images)
```

The mounted gufe source is editable-installed, so host edits to `src/gufe/` (or
to a `viz_assets/` frame) are live — just restart the notebook kernel.

Useful checks, all from `visualization-demo/`:

| command | does |
| --- | --- |
| `just ci` | build → start → wait → run the framejs tests, from cold |
| `just test-framejs` | just the hermetic framejs tests against the running stack |
| `just doctor` | confirm `gufe` + `metaframe_widget` import from the editable mount |

`just dev` also starts the **terminal** path — the `viz` service, which is
`python -m gufe.visualization.server` over gufe's own test fixtures, bound to
`0.0.0.0:8899` because it is in a container. Open <http://localhost:8899> for a
listing and click any fixture; `just serve <file>` does the same for one file of
your choosing (on :8900, so it can run alongside).

## Editing the visualizations

Each visualization is a framejs **frame directory** under `viz_assets/`, in the
canonical [framejs local-file-io](https://framejs.io/docs/guide/local-file-io)
format. For example `viz_assets/ligand_network/`:

| file           | contents                                                        |
| -------------- | --------------------------------------------------------------- |
| `code.js`      | the frame's JavaScript (the actual visualization)               |
| `og.json`      | title + description metadata                                    |
| `modules.json` | external classic-script URLs loaded eagerly — no frame ships one |

Heavy third-party libraries are **not** listed in `modules.json`, because
everything there blocks the frame's first paint. Instead each frame injects what
it needs on demand — `loadRDKit()` for the 2D depictions and `load3Dmol()` for
the 3D viewers — so a frame renders immediately and only pays for the engines a
given view actually opens (`small_molecule_component` draws its 2D depiction
without waiting on 3Dmol; `ligand_network` fetches 3Dmol only once you click an
edge, and never in `2d` mode).

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
> `just viz-url` / `just viz-save` to round-trip a frame ↔ a framejs.io URL, and
> `just publish-viz` to mint a `/j/<uuid>`), but the raw `deno` command above is
> all you need.

Frames that are **not** gufe vizzes — the architecture diagram linked at the top —
live under `visualization-demo/frames/` instead, since every directory in
`viz_assets/` is package data named after the gufe class it draws. Edit those with
`just viz-edit frames`.
