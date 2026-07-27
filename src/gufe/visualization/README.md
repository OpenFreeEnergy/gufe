# gufe visualization

Interactive [framejs.io](https://framejs.io) visualizations for gufe objects,
across three surfaces from one frame: a notebook widget (`.view()`), a local
web server (`openfe view <file>`), and a standalone `.html` file. Rendering is
powered by the optional `viz` extra (`pip install gufe[viz]`); without it
`.view()` falls back to the legacy RDKit / py3Dmol renderers.

> **Architecture, interactively:**
> <https://framejs.app/j/019f8b367fac7e469f26b99bb30e06e4> — the render paths and
> the registry table, clickable. Source: `visualization-demo/frames/architecture/`.

The three top sections are, in order: **[Using it](#using-it)**,
**[Editing the visualization](#editing-the-visualization)** (the frame code), and
**[Developing](#developing)** (the demo stack + tests). **[How it works](#how-it-works)**
at the end is reference.

---

## Using it

Every gufe type that mixes in `FramejsViewable` renders the same way on all three
surfaces. Nothing here needs an account or network beyond loading framejs.io.

### In a notebook — `.view()` / bare cell

```python
import gufe
net = gufe.LigandNetwork.from_graphml(open("ligand_network.graphml").read())

net                       # a bare object on the last line renders itself
net.view()                # the same widget, explicitly
net.view(height="600px")  # …with a size override
```

The object pushes its data into the frame live over the widget comm channel, so
there is no URL size limit and the view updates in place. Works in Jupyter,
marimo and VSCode. See **[Demo 1](#developing)**.

### From the terminal — `openfe view`

```bash
openfe view network_setup/ligand_network.graphml   # serve it, open a browser, wait for Ctrl-C
openfe view results/                               # a directory: browse every viewable file
openfe view ligand.sdf --no-browser                # only print the URL (ssh, container)
openfe view ligand.sdf --proxy                     # serve framejs from here too — needed for Chrome
openfe view-ligand-network network.graphml         # the older .graphml-only entry point
```

`openfe view` starts a small local web server (`gufe.visualization.server`), prints
its URL, opens a browser, and stays up until **Ctrl-C**. The URL path *is* the
data file's path, so a link says what it shows; every viewable file under the
served root is reachable from the one server.

Which file types are viewable is `LOADER_REGISTRY`: `.json` (any saved gufe
object), `.graphml`, `.pdb`, `.cif`/`.pdbx`, `.sdf`/`.mol`. A file it cannot load,
or one with no registered viz, gets an explanatory page rather than a traceback.

> **Chrome needs `--proxy`.** By default the framejs frame is served from
> framejs.io, and Chrome 142+ blocks that public-origin frame from fetching your
> data off `localhost` (Local Network Access). `--proxy` re-hosts framejs on the
> server's own origin so frame and data match, which sidesteps it; Firefox works
> either way. Details in [How it works](#how-it-works). In a container, also pass
> `--host 0.0.0.0` and publish the port.

### As a shareable link or a standalone file

```bash
openfe view ligand.sdf --url-only                        # one self-contained framejs.io URL
openfe view ligand.sdf --html report.html                # one self-contained .html file
openfe view ligand.sdf --html report.html --self-contained   # …that needs no network at all
```

`--url-only` bakes the object into a single framejs.io link — paste it anywhere,
but the data is frozen and large objects overflow the URL. `--html` writes a
`.html` you can email; it fetches the RDKit / 3Dmol / d3 engines from a CDN when
opened (~130 kB), and `--self-contained` inlines those too (~10 MB, no network).

---

## Editing the visualization

All the drawing is one framejs **frame**, `viz_assets/gufe/`, in the canonical
[framejs local-file-io](https://framejs.io/docs/guide/local-file-io) format:

| file           | contents                                                        |
| -------------- | --------------------------------------------------------------- |
| `code.js`      | every view, plus the dispatcher — the actual visualization      |
| `og.json`      | title + description metadata                                    |
| `modules.json` | eager external scripts — the frame ships none (see below)       |

`code.js` reads top to bottom as: theme → engine loaders → input helpers → DOM
helpers → molecule / mapping / protein helpers → the shared atom-mapping viewer →
one `view*()` function per gufe type → the `VIEWS` dispatch table. A `view*()`
function is handed a fresh container and the `inputs`, and returns an optional
`{ onResize, cleanup }` handle; the dispatcher tears the previous view down before
mounting the next, so views never have to know about each other.

Heavy libraries are **not** in `modules.json`, because everything there blocks the
first paint. The frame injects what a view needs on demand — `loadRDKit()` for 2D
depictions, `load3Dmol()` for 3D, `loadD3()` for the graph views — so it paints
immediately and only pays for the engines a view actually opens.

### Live-editing in the browser

Run the framejs **local server**: it watches these files and auto-saves your edits
back to disk. It needs only [`deno`](https://docs.deno.com/runtime/getting_started/installation/).
From the gufe repo root:

```bash
deno run --allow-read="$PWD/src/gufe/visualization/viz_assets" \
         --allow-write="$PWD/src/gufe/visualization/viz_assets" \
         --allow-net --allow-env \
  https://raw.githubusercontent.com/metapages/framejs.io/main/local-server/server.ts \
  --root "$PWD/src/gufe/visualization/viz_assets" --port 4700
```

Open <http://localhost:4700>, pick the frame (dirs with a `code.js` show as
`◆ gufe`), edit, watch it save (`saved ✓`). `.view()` and `openfe view` build
their URL straight from this on-disk source, so edits show up on the next render
(restart the notebook kernel; reload the served page).

`scripts/framejs_frame.py` (under `visualization-demo/`) round-trips a frame
directory ↔ a framejs.io URL if you'd rather edit on framejs.io and pull the
result back down:

```bash
python3 visualization-demo/scripts/framejs_frame.py restore src/gufe/visualization/viz_assets/gufe   # dir → URL
python3 visualization-demo/scripts/framejs_frame.py save '<framejs.io URL>' src/gufe/visualization/viz_assets/gufe  # URL → dir
```

### Adding a view for a new gufe type

Three edits — no new frame:

1. `code.js`: write a `view*()` function and add a `VIEWS` entry keyed on the
   payload key that identifies the type.
2. `framejs.py`: write a payload builder and register it in `PAYLOAD_REGISTRY`,
   emitting that dispatch key.
3. the gufe class: mix in `FramejsViewable`. Do **not** define
   `_ipython_display_` (IPython short-circuits on it, so a bare cell and `.view()`
   would disagree); define `_legacy_view()` for the non-framejs fallback instead.

`test_every_payload_claims_exactly_one_dispatch_key` and
`test_the_frame_dispatches_on_every_key_the_registry_emits` keep the two sides in
sync.

---

## Developing

The demo stack under [`visualization-demo/`](../../../visualization-demo/) runs
everything in Docker with plain `docker compose` — no `just`, nothing on the host
but Docker. It builds a fast native-arm64 gufe env and editable-installs this
checkout, so host edits (including to `viz_assets/gufe/code.js`) are live. Its
[README](../../../visualization-demo/README.md) is the full guide; the two demos:

```bash
cd visualization-demo

docker compose up notebook       # Demo 1: .view() in JupyterLab   → http://localhost:8888/lab
docker compose up openfe-view    # Demo 2: `openfe view` server     → http://localhost:8899
```

Demo 1 opens `demo/framejs_demo.ipynb`, which renders one of every visualizable
gufe type from the bundled test data. Demo 2 runs exactly what `openfe view` runs
(with `--proxy`, so it works in Chrome).

Run the hermetic tests (no network) in a throwaway container — prefix
`dev-start.sh` so it editable-installs gufe first:

```bash
docker compose run --rm notebook dev-start.sh \
  python -m pytest /src/gufe/src/gufe/tests/test_framejs_visualization.py \
                   /src/gufe/src/gufe/tests/test_framejs_server.py -q
```

- `test_framejs_visualization.py` — the registry, serializers, URL/HTML building,
  display hooks, and fallback.
- `test_framejs_server.py` — the `openfe view` server: routing, the payload
  endpoints, `--proxy`, and path-traversal refusal.

The Python side is `framejs.py` (`PAYLOAD_REGISTRY` + URL/HTML builders),
`server.py` (the local server), `../_viewable.py` (the `FramejsViewable` mixin),
and `mapping_visualization.py` (the one remaining legacy renderer).

---

## How it works

### One frame draws everything

Every visualizable gufe class renders through the **same** frame,
`viz_assets/gufe/`. The frame picks its view from the payload it is handed — a
`network.graphml` key means a ligand network, `protein.pdb` a 3Dmol protein,
`chemical_system` a component browser — so Python never has to know which drawing
code applies, and the base URL is built once and is identical for every object.
That is what lets the shared parts be shared: the atom-mapping viewer exists once
and drives the `LigandAtomMapping` page, the right pane of the ligand network, and
the bottom of a transformation; likewise the theme, the lazy engine loaders, and
the SDF / PDB parsers.

The dispatch key is the payload key the serializer emits; the frame's `VIEWS`
table maps it to the drawing code, first match wins.

| gufe class | dispatch key | shows |
| --- | --- | --- |
| `LigandNetwork` | `network.graphml` | radial network; click an edge → 3D atom-mapping viewer |
| `AlchemicalNetwork` | `alchemical_network` | d3 force graph of ChemicalSystem nodes / Transformation edges |
| `Transformation`, `NonTransformation` | `transformation` | stateA↔stateB component diff + the atom mapping |
| `ChemicalSystem` | `chemical_system` | master/detail over the labelled components |
| `LigandAtomMapping` | `molA.sdf` | the mapping viewer standalone (plain/colored/lines/overlay/2d) |
| `SmallMoleculeComponent` | `molecule.sdf` | 2D depiction + 3D conformer + SMILES/charge |
| `ProteinComponent` | `protein.pdb` | 3Dmol with representation / colour-scheme switchers |
| `SolventComponent` | `solvent_component` | settings card (it has no coordinates) |

Lookup walks the MRO, so subclasses inherit their base's serializer:
`SolvatedPDBComponent` and `ProteinMembraneComponent` get the protein payload with
nothing registered for them. The transformation entry is on `TransformationBase`,
because `Transformation` and `NonTransformation` are siblings under it. A payload
must claim exactly one dispatch key — the two tests named above pin that down.

This package reads **no environment variables**: what renders, and from where, is
determined entirely by `PAYLOAD_REGISTRY` and the caller's arguments.

### Payload keys are the contract

Payload keys are what the frame's `onInputs` dispatches on. File-shaped values use
a descriptive `<thing>.<ext>` key (`molecule.sdf`, `protein.pdb`,
`network.graphml`) — the frame reads these through `asText()`, which also accepts
a dropped file. Everything else is a bare snake_case field, and where a view takes
one structured object that object's key names what it describes
(`chemical_system`, `transformation`, `alchemical_network`, `solvent_component`).
Renaming a key breaks any frame already published to a `/j/<uuid>`, which is
pinned and keeps reading the old name.

### Render paths, one core

Every path serializes the object with the same `PAYLOAD_REGISTRY` builder and
loads the same base URL, built once from the one on-disk frame. They differ
**only** in how that payload reaches the frame.

| | Notebook | Terminal | Shareable link | Standalone file |
| --- | --- | --- | --- | --- |
| entry point | `.view()` / bare cell | `openfe view <file>` | `openfe view --url-only` | `openfe view --html FILE` |
| carrier | `set_inputs()` over the comm channel | `fetch()` from the served page → `setInputs()` | `#?inputs=<base64>` in the URL hash | a JSON `<script>` block in the page |
| size limit | none | none | the URL's | none |
| needs | a kernel | a running server | framejs.io | nothing (`--self-contained`: not even a network) |
| result | an iframe in the cell, updatable | a page at `localhost:8899/<the file>` | one link you can paste | one `.html` you can email |

The first two are the same design — load the frame, then push data into it — which
is why a frame doesn't care which is hosting it, and why neither has a size limit.
The link and file bake the object in instead; appended `inputs` take priority over
anything baked into the frame, which is what makes those forms work. `build_html()`
is the only one that loads no framejs.io URL at all: it inlines `code.js` and the
payload directly, with a classic-`<script>` bootstrap that sets up the `root`
global the frame reads before the deferred module runs.

### Why Chrome needs `--proxy`

By default the terminal path serves a page that carries the framejs.io frame URL
and fetches the object's data from `?inputs=1` on its own origin, then pushes it
in over `postMessage`. It would be tidier to skip the page's fetch and let the
frame pull its own data from a `?input=<key>` URL in its hash — the frame supports
exactly that (`asText()`/`asObject()` fetch a whitespace-free `http(s)://` string
and use anything else as-is). **Chrome no longer allows it across origins.**
Chrome 142's [Local Network Access](https://developer.chrome.com/blog/local-network-access)
blocks a public-origin page — which the framejs.io frame is — from reaching a
loopback address without a user permission prompt:

```
Access to fetch at 'http://localhost:8899/181l.pdb?input=protein.pdb' from origin
'https://framejs.io' has been blocked by CORS policy: Permission was denied for
this request to access the `loopback` address space.
```

No response header fixes it (`Access-Control-Allow-Private-Network` was the older
Private Network Access opt-in and is superseded). Firefox has not shipped this and
works either way.

`--proxy` removes the barrier instead of working around it: the server re-hosts
the framejs.io document at `/_framejs/` and forwards anything it doesn't recognize
upstream, so page, iframe and data are all `http://localhost:8899`. Nothing is
rewritten in transit — framejs's app code is inline and its dependencies are
absolute `https://cdn.jsdelivr.net` imports, which an `http://localhost` page (a
secure context, an upgrade not mixed content) loads normally. With one origin the
frame can feed itself, so the served page needs no JavaScript at all. Two
deliberate limits: `/sw.js` is refused (`PROXY_DENY`) so a service worker can't
take scope `/` over the data endpoints or outlive the server, and proxied
responses are cached for the life of the process (restart to pick up a framejs.io
release).

### Two URL forms, and fallback

The base viz URL is built from the on-disk frame in one of two forms, both
carrying the same `code.js`:

- **local** (default everywhere) — a self-contained `https://framejs.io/#?js=<b64>…`.
  Always matches the installed gufe, needs no account, cannot expire.
- **canonical** — the pinned short `https://framejs.io/j/<uuid>`, opt-in via
  `build_cli_url(obj, short=True)` where URL size matters. Requires the frame to
  have been published and its uuid pinned as `FRAME_UUID`, which has not happened
  since the per-class frames were merged into one. `resolve_url()` prefers the
  on-disk frame and only falls back to the uuid if no frame directory is present.

Everything here is optional — `import gufe` works without `metaframe-widget`. On
`FramejsUnavailable`/`OSError`, both paths degrade: the object falls back to
`_legacy_view()` if it defines one, else its plain `repr`. `.view()` warns on
fallback; auto-display stays silent, since it runs on every display.
