# gufe visualization

Interactive [framejs.io](https://framejs.io) visualizations for gufe objects.
Calling `.view()` — or just leaving an object as the last line of a notebook cell
in Jupyter / marimo / VSCode — renders it as an interactive widget. Rendering is
powered by the optional `viz` extra (`pip install gufe[viz]`); without it,
`.view()` falls back to the legacy RDKit / py3Dmol renderers.

- `framejs.py` — the registry: maps each gufe class to its frame (`CANONICAL_VIZ`)
  and to the serializer that turns an object into that frame's `inputs`
  (`_PAYLOAD_BUILDERS`), then builds the self-contained framejs.io URL /
  `anywidget`.
- `../_viewable.py` — the `FramejsViewable` mixin that gives a class `.view()` and
  `_repr_mimebundle_`. Mixing it in is the whole opt-in.
- `mapping_visualization.py` — the legacy RDKit atom-mapping renderer.
- `viz_assets/<frame>/` — the on-disk framejs frames that ship with gufe.

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
and `ProteinMembraneComponent` get the protein frame, `NonTransformation` gets the
transformation frame, with nothing registered for them.

To add one: mix `FramejsViewable` into the class, add a `VizRef` and a payload
builder in `framejs.py`, and drop a frame directory under `viz_assets/`.

## Running the Jupyter notebook demo

The demo stack (`visualization-demo/`) is fully driven by Docker Compose — no
`just` required. It builds a native gufe env, editable-installs this checkout,
and serves a JupyterLab notebook that exercises `.view()` on one of every
visualizable gufe object.

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
format. For example `viz_assets/ligand_network/`:

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

