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

- `framejs.py` — `PAYLOAD_REGISTRY` maps each gufe class to the serializer that
  turns an object into the frame's `inputs`. It then builds the framejs.io URL /
  `anywidget`.
- `../_viewable.py` — the `FramejsViewable` mixin that gives a class `.view()` and
  `_repr_mimebundle_`. Mixing it in is the whole opt-in.
- `mapping_visualization.py` — the legacy RDKit atom-mapping renderer, reached
  through `LigandAtomMapping._legacy_view()` when framejs is unavailable. It is
  the only legacy renderer left.
- `server.py` — the terminal path: a small local web server that shows a **file's**
  visualization in a browser (`python -m gufe.visualization.server <file>`, which
  is what `openfe view <file>` runs). Its `LOADER_REGISTRY` maps a file extension
  to the gufe object it loads to; `PAYLOAD_REGISTRY` takes it from there.
- `viz_assets/gufe/` — the one on-disk framejs frame that ships with gufe.
- `../tests/test_framejs_visualization.py`, `../tests/test_framejs_server.py` —
  hermetic (no network) tests for all of the above, including the invariants
  stated here.

## One frame draws everything

Every visualizable gufe class renders through the **same** frame,
`viz_assets/gufe/`. The frame picks its view from the payload it is handed — a
`network.graphml` key means a ligand network, `protein.pdb` a 3Dmol protein,
`chemical_system` a component browser — so Python never has to know which drawing
code applies, and the base URL is built once and is identical for every object.

That is what makes the shared parts shared: the atom-mapping viewer exists once
and is reused by the `LigandAtomMapping` page, the right pane of the ligand
network, and the bottom of a transformation; the theme, the lazy RDKit / 3Dmol /
d3 loaders, and the SDF / PDB parsers likewise.

Adding a visualization is one `PAYLOAD_REGISTRY` entry plus one branch in the
frame's `VIEWS` table. A class must **not** define `_ipython_display_`: IPython
checks that hook before `_repr_mimebundle_` and short-circuits on it, so a bare
cell and `.view()` would disagree. Define `_legacy_view()` instead for the
non-framejs fallback.

This package reads **no environment variables**. What renders, and from where, is
determined entirely by `PAYLOAD_REGISTRY` and the caller's arguments.

Payload keys are the contract with the frame's `onInputs`, and are what it
dispatches on. File-shaped values use a descriptive `<thing>.<ext>` key
(`molecule.sdf`, `protein.pdb`, `network.graphml`) — the frame reads these through
`asText()`, which also accepts a dropped file. Everything else is a bare
snake_case field, and where a view takes one structured object that object's key
names the thing it describes (`chemical_system`, `transformation`,
`alchemical_network`, `solvent_component`). Renaming a key breaks any frame
already published to a `/j/<uuid>`, which is pinned and keeps reading the old
name.

## What renders

The dispatch key is the payload key the serializer emits; the frame's `VIEWS`
table maps it to the drawing code, first match wins.

| gufe class | dispatch key | shows |
| --- | --- | --- |
| `LigandNetwork` | `network.graphml` | radial network; click an edge to drive a 3D atom-mapping viewer |
| `AlchemicalNetwork` | `alchemical_network` | d3 force graph of ChemicalSystem nodes / Transformation edges |
| `Transformation`, `NonTransformation` | `transformation` | stateA↔stateB component diff + the atom mapping |
| `ChemicalSystem` | `chemical_system` | master/detail over the system's labelled components |
| `LigandAtomMapping` | `molA.sdf` | the mapping viewer standalone (plain/colored/lines/overlay/2d) |
| `SmallMoleculeComponent` | `molecule.sdf` | 2D depiction + 3D conformer + SMILES/charge |
| `ProteinComponent` | `protein.pdb` | 3Dmol with representation / colour-scheme switchers |
| `SolventComponent` | `solvent_component` | settings card (it has no coordinates) |

Lookup walks the MRO, so subclasses inherit their base's serializer:
`SolvatedPDBComponent` and `ProteinMembraneComponent` get the protein payload,
with nothing registered for them. The transformation entry is registered on
`TransformationBase` rather than on `Transformation`, because `Transformation` and
`NonTransformation` are siblings under it, not parent and child.

Because the frame dispatches on the first key it recognizes, a payload must claim
exactly one of them — `test_every_payload_claims_exactly_one_dispatch_key` pins
that down, as does `test_the_frame_dispatches_on_every_key_the_registry_emits`
for the other direction.

To add one: mix `FramejsViewable` into the class, add a payload builder and a
`PAYLOAD_REGISTRY` entry in `framejs.py`, and add a `VIEWS` entry plus its render
function in `viz_assets/gufe/code.js`.

## Render paths, one core

Every path serializes the object with the same `PAYLOAD_REGISTRY` builder and
loads the same base URL, built once from the one on-disk frame. They differ
**only** in how that payload reaches the frame.

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

**Inputs may be content or a URL to content**, and the frame accepts both:
`asText()` / `asObject()` fetch a whitespace-free `http(s)://` string and use
anything else as-is. Read file-shaped inputs through those helpers rather than
touching `inputs[…]` directly and a new view gets that for free. Note the one
place it does *not* help — see below.

### Why the page fetches, and not the frame

It would be tidier for the terminal path to put its inputs straight in the
iframe's hash as URLs and let the frame pull its own data from
`?input=<key>`. **Browsers no longer allow that.** Chrome 142's [Local Network
Access](https://developer.chrome.com/blog/local-network-access) blocks a
public-origin page — which the frame, served from framejs.io, is — from reaching a
loopback address without an explicit user permission prompt, and a cross-origin
iframe additionally needs the `local-network-access` Permissions Policy delegated
to it. It surfaces as a CORS failure:

```
Access to fetch at 'http://localhost:8899/181l.pdb?input=protein.pdb' from origin
'https://framejs.io' has been blocked by CORS policy: Permission was denied for
this request to access the `loopback` address space.
```

No response header fixes it; `Access-Control-Allow-Private-Network` was the header
for the older Private Network Access design and is superseded. The viewer page has
no such problem, because it is served from `localhost:8899` itself — its fetch is
same-origin, so no address space is crossed — and handing the result to the frame
is a `postMessage`, which no such policy applies to.

### `--proxy`: remove the barrier instead of working around it

The restriction is about *crossing* address spaces, so it disappears if the frame
is on this origin too. `openfe view <file> --proxy` re-hosts the framejs.io
document at `/_framejs/` and forwards anything else it does not recognize
upstream, so the page, the iframe and the data are all `http://localhost:8899`:

```
page    http://localhost:8899/181l.pdb
iframe  http://localhost:8899/_framejs/#?js=…&inputs=…
fetch   http://localhost:8899/181l.pdb?input=protein.pdb
```

Nothing is rewritten in transit. framejs's application code is inline in its HTML
and its dependencies are absolute `https://cdn.jsdelivr.net` imports, which an
`http://localhost` page loads normally — that direction is an upgrade, not mixed
content, and localhost is a secure context. Its own root-relative requests
(`/favicon.svg`, …) fall through the catch-all proxy. Data paths always win over
the proxy, so a file you pointed the server at is never shadowed.

With the frame able to feed itself, this mode gets the simpler design back:
inputs are URLs in the iframe's hash (`server.hash_inputs()`), the page contains
**no JavaScript at all**, and the URL stays small however large the object.

Two things it deliberately does not do. `/sw.js` is refused (`PROXY_DENY`): a
service worker registered here would take scope `/` over the data endpoints,
whose contract is that they are re-read from disk, and it would outlive both the
page and the server — the next thing on that port would inherit it. And proxied
responses are cached for the life of the process, so restart the server to pick
up a framejs.io release.

### Two URL forms

Both are built from the one on-disk frame directory:

- **local** — a self-contained `https://framejs.io/#?js=<b64>&og=<b64>…`. **The
  default everywhere.** It always matches the installed gufe, needs no framejs
  account, and cannot expire.
- **canonical** — the pinned short `https://framejs.io/j/<uuid>`, minted once by
  `just publish-viz`. Its one advantage is size, so it is opt-in exactly where
  size can matter: `build_cli_url(obj, short=True)`. Inlining `code.js` is what
  makes the local form big — the difference between a link you can paste
  somewhere and one you cannot. It requires the frame to have been published and
  its uuid pinned as `FRAME_UUID`, which has not happened since the per-class
  frames were merged into one.

`resolve_url()` prefers the on-disk frame and falls back to the pinned uuid only
if no frame directory is present.

### The terminal path in practice

`openfe`'s viewer commands are thin wrappers over `server.serve`:

```bash
openfe view network_setup/ligand_network.graphml   # serve it, open a browser, wait for Ctrl-C
openfe view results/                               # a listing: browse every viewable file
openfe view ligand.sdf --no-browser                # only print the URL (ssh, container)
openfe view ligand.sdf --url-only                  # no server: one self-contained link
openfe view ligand.sdf --proxy                     # serve framejs from here too (see below)
openfe view-ligand-network network.graphml         # the older .graphml-only entry point
```

What gets served at `http://localhost:8899/<path/to/the/file>` is a page Python
renders per file: a header and **one iframe**, with no JavaScript at all. The
iframe's `src` is the registered viz's framejs URL with the object's `inputs` in
its hash — so the page is a thin frame around a link you could paste into a
browser yourself. Each input too big to inline is
`http://localhost:8899/<the file>?input=<key>`, which serves exactly that value
and which the frame fetches for itself.

The URL path *is* the data file's path, relative to the served root — the
directory you pointed at, or the working directory if the file is under it.
Everything viewable under that root is reachable from the one server, and nothing
outside it is. The payload is rebuilt from disk on every request, so regenerating
a file and reloading the page re-renders it without restarting anything.

Because the frame is served from framejs.io and fetches from localhost, every
response carries `Access-Control-Allow-Origin: *` and preflights are answered,
including Chrome's Private Network Access check for an https page reaching a
local address.

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

The visualization is one framejs **frame directory**, `viz_assets/gufe/`, in the
canonical [framejs local-file-io](https://framejs.io/docs/guide/local-file-io)
format:

| file           | contents                                                        |
| -------------- | --------------------------------------------------------------- |
| `code.js`      | the frame's JavaScript (every view, plus the dispatcher)        |
| `og.json`      | title + description metadata                                    |
| `modules.json` | external classic-script URLs loaded eagerly — the frame ships none |

Heavy third-party libraries are **not** listed in `modules.json`, because
everything there blocks the frame's first paint. Instead the frame injects what
it needs on demand — `loadRDKit()` for the 2D depictions, `load3Dmol()` for the
3D viewers, `loadD3()` for the two graph views — so it renders immediately and
only pays for the engines a given view actually opens (the small-molecule view
draws its 2D depiction without waiting on 3Dmol; the ligand network fetches 3Dmol
only once you click an edge, and never in `2d` mode; a solvent card fetches
neither 3Dmol nor d3).

`code.js` is laid out top to bottom as: theme, engine loaders, input helpers, DOM
helpers, molecule / mapping / protein helpers, the shared atom-mapping viewer,
then one `view*()` function per gufe type, then the `VIEWS` dispatch table. A
view function is handed a fresh container and the `inputs`, and returns an
optional `{ onResize, cleanup }` handle; the dispatcher tears the previous view
down before mounting the next, so views never have to know about each other.

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

Frames that are **not** the gufe viz — the architecture diagram linked at the top —
live under `visualization-demo/frames/` instead, since `viz_assets/` holds exactly
the one frame gufe ships as package data. Edit those with `just viz-edit frames`.
