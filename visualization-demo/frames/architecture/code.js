// ============================================================================
// gufe framejs visualization — ARCHITECTURE
//
// A documentation frame, not a data viz: it draws the two render paths that
// turn a gufe object into an interactive framejs.io view, and the shared
// registry core both of them route through.
//
//   Notebook   .view() / bare cell  → MetaframeWidget → live `inputs` over comm
//   CLI        openfe view <file>   → a URL with `inputs=<b64>` in its hash
//
// Click any stage for its detail + the file it lives in. The REGISTRY tab is
// the same table as visualization/README.md, so the two are checkable against
// each other.
//
// This frame ships under visualization-demo/frames/ (dev/doc material), NOT
// under gufe/visualization/viz_assets/ — that directory is package data, and
// every directory in it is named after the gufe class it renders.
//   just viz-edit frames        → live-edit this frame
// ============================================================================

// ============================================================================
// THEME — follows the viz_assets house style; `auto` tracks the host page
// ============================================================================
var DARK_MODE = 'auto'; // true | false | 'auto'

var THEMES = {
  dark: {
    appBg:        '#1a1a2e',
    panelBg:      '#16213e',
    cardBg:       '#16213e',
    cardBorder:   '#2a4a7f',
    cardBorderOn: '#7ecfff',
    titleColor:   '#7ecfff',
    textPrimary:  '#e6f3ff',
    textMuted:    '#9bb8d6',
    textMuted2:   '#7a96b8',
    chipBg:       '#0f3460',
    chipBorder:   '#2a4a7f',
    tabBg:        '#0f3460',
    tabOnBg:      '#7ecfff',
    tabOnFg:      '#0b1a2e',
    connector:    '#2a4a7f',
    codeBg:       '#0b1a2e',
    kinds: {
      user:     '#c792ea',
      gufe:     '#7ecfff',
      framejs:  '#82e0aa',
      core:     '#ffcc66',
      browser:  '#ff79c6',
      fallback: '#93a3b8'
    }
  },
  light: {
    appBg:        '#ffffff',
    panelBg:      '#f8fafc',
    cardBg:       '#ffffff',
    cardBorder:   '#e2e8f0',
    cardBorderOn: '#0369a1',
    titleColor:   '#0369a1',
    textPrimary:  '#1e293b',
    textMuted:    '#475569',
    textMuted2:   '#64748b',
    chipBg:       '#f1f5f9',
    chipBorder:   '#cbd5e1',
    tabBg:        '#f1f5f9',
    tabOnBg:      '#0369a1',
    tabOnFg:      '#ffffff',
    connector:    '#cbd5e1',
    codeBg:       '#f1f5f9',
    kinds: {
      user:     '#7c3aed',
      gufe:     '#0369a1',
      framejs:  '#0f766e',
      core:     '#b45309',
      browser:  '#be185d',
      fallback: '#64748b'
    }
  }
};

function pickTheme() {
  if (DARK_MODE === true) return THEMES.dark;
  if (DARK_MODE === false) return THEMES.light;
  const mq = window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)');
  return mq && mq.matches ? THEMES.dark : THEMES.light;
}

var T = pickTheme();

const MONO = "ui-monospace,'SF Mono',Menlo,Consolas,monospace";
const SANS = "'Inter',system-ui,-apple-system,sans-serif";

// Width below which the flow and the detail panel stack instead of sitting
// side by side. Everything else is fluid.
const NARROW_PX = 700;

// ============================================================================
// THE ARCHITECTURE — the whole content of this frame is these three tables
// ============================================================================
//
// kind → the colour band on a stage card:
//   user     something the user holds or types
//   gufe     core gufe (the mixin, the object)
//   framejs  gufe's framejs layer (visualization/framejs.py)
//   core     the shared registry core — identical on both paths
//   browser  what finally renders
//   fallback the degraded path taken when framejs is unavailable

const PATHS = {
  notebook: {
    label: 'Notebook',
    icon: '📓',
    tagline: 'Jupyter · marimo · VSCode — the object pushes live inputs over the widget comm channel.',
    steps: [
      {
        id: 'nb-obj',
        kind: 'user',
        title: 'a gufe object',
        sub: 'LigandNetwork · AlchemicalNetwork · Transformation · ChemicalSystem · LigandAtomMapping · SmallMoleculeComponent · ProteinComponent · SolventComponent',
        file: 'src/gufe/*.py',
        detail: 'Any of the eight classes that mix in FramejsViewable. Subclasses ' +
                'need nothing of their own: SolvatedPDBComponent renders as a ' +
                'ProteinComponent and NonTransformation as a Transformation, ' +
                'because the registry lookup walks the MRO.'
      },
      {
        id: 'nb-mixin',
        kind: 'gufe',
        title: 'FramejsViewable',
        sub: '.view()  ·  _repr_mimebundle_()',
        file: 'src/gufe/_viewable.py',
        arrow: 'bare cell, or an explicit .view()',
        detail: 'Mixing this in is the entire opt-in. Both entry points land in ' +
                'the same place, so a bare cell and .view() can never disagree — ' +
                'which is also why no viewable class may define ' +
                '_ipython_display_: IPython checks that hook first and would ' +
                'short-circuit past _repr_mimebundle_. A test enforces it.'
      },
      {
        id: 'nb-entry',
        kind: 'framejs',
        title: 'view_object()  /  repr_mimebundle()',
        sub: 'the two notebook entry points into the framejs layer',
        file: 'visualization/framejs.py',
        detail: 'Both call _build_widget(). They differ only in how loudly they ' +
                'fail: view_object() was asked for explicitly, so it warns when ' +
                'it has to fall back; repr_mimebundle() runs on every display of ' +
                'the object, so it is silent.'
      },
      {
        id: 'nb-registry',
        kind: 'core',
        shared: true,
        title: 'VIZ_REGISTRY  →  VizRef',
        sub: '_viz_for(obj) — class name lookup, walking the MRO',
        file: 'visualization/framejs.py',
        detail: 'The shared core, identical on both paths. One dict maps a gufe ' +
                'class name to a VizRef holding both halves of the contract: ' +
                'which frame draws the object, and how the object is serialized ' +
                'for that frame. No match → FramejsUnavailable.'
      },
      {
        id: 'nb-url',
        kind: 'core',
        shared: true,
        title: 'VizRef.resolve_url()',
        sub: 'local: framejs.io/#?js=<b64>&og=<b64>… built from viz_assets/<frame>/',
        file: 'visualization/viz_assets/<frame>/code.js',
        detail: 'The base URL — the viz JavaScript, and nothing about this ' +
                'particular object. It is built by base64-encoding the on-disk ' +
                'frame directory (code.js plus its JSON sidecars), so it always ' +
                'matches the installed gufe, needs no framejs account and cannot ' +
                'expire. The pinned /j/<uuid> form is used only if no frame dir ' +
                'is on disk.'
      },
      {
        id: 'nb-payload',
        kind: 'core',
        shared: true,
        title: 'VizRef.payload(obj)',
        sub: 'the object → a flat dict of inputs',
        file: 'visualization/framejs.py',
        detail: 'The other half of the contract: a flat dict whose keys are read ' +
                'by the frame\'s onInputs. File-shaped values get a descriptive ' +
                '<thing>.<ext> key (molecule.sdf, protein.pdb, network.graphml) ' +
                'and double as file-drop targets; everything else is a bare ' +
                'snake_case field. Renaming a key breaks any already-published ' +
                'frame, so the keys are API.'
      },
      {
        id: 'nb-widget',
        kind: 'framejs',
        title: 'MetaframeWidget(url, width, height)',
        sub: '.set_inputs(payload)  — anywidget, from the gufe[viz] extra',
        file: 'metaframe-widget (PyPI)',
        detail: 'The anywidget that hosts the framejs iframe in the notebook. ' +
                'gufe hands it the base URL once and then pushes the payload as ' +
                'live inputs over the comm channel — never through the URL, so ' +
                'this path has no size limit at all.'
      },
      {
        id: 'nb-render',
        kind: 'browser',
        title: 'the frame renders, inline in the cell',
        sub: 'inputs arrive live; re-pushing updates the same iframe',
        file: 'framejs.io',
        detail: 'Because inputs travel over the comm channel rather than the ' +
                'URL, the widget stays live: pushing a new payload re-renders ' +
                'the frame in place.'
      }
    ],
    fallback: {
      id: 'nb-fallback',
      kind: 'fallback',
      title: 'legacy_view()  →  _legacy_view()',
      sub: 'RDKit / py3Dmol — else the plain repr',
      file: 'visualization/framejs.py · visualization/mapping_visualization.py',
      when: 'FramejsUnavailable / OSError — no gufe[viz] extra, nothing registered, no frame on disk',
      detail: 'Nothing here is required for gufe to import or run. If the framejs ' +
              'path cannot be taken, the object falls back to its pre-framejs ' +
              'renderer if it defines _legacy_view() (only LigandAtomMapping ' +
              'does, via mapping_visualization.py), and otherwise to the plain ' +
              'repr. .view() warns; auto-display stays silent.'
    }
  },

  cli: {
    label: 'CLI',
    icon: '💻',
    tagline: 'No live Python channel — so the object rides along inside the URL itself.',
    steps: [
      {
        id: 'cli-file',
        kind: 'user',
        title: 'a file on disk',
        sub: '.graphml · .pdb · .cif / .pdbx · .sdf / .mol',
        file: 'network_setup/ligand_network.graphml',
        detail: 'The CLI starts from a file rather than a live object, which is ' +
                'the only real difference between the two paths.'
      },
      {
        id: 'cli-cmd',
        kind: 'user',
        title: 'openfe view <file>',
        sub: '--no-browser to print the URL instead of opening it',
        file: 'openfecli/commands/view.py',
        arrow: 'typed at a terminal',
        detail: 'A generic viewer command: _load_object() dispatches on the file ' +
                'extension to the right gufe loader. openfe view-ligand-network ' +
                'is the older, .graphml-only entry point onto the same machinery.'
      },
      {
        id: 'cli-obj',
        kind: 'gufe',
        title: 'the loaded gufe object',
        sub: 'LigandNetwork.from_graphml() · ProteinComponent.from_pdb_file() · …',
        file: 'openfecli/commands/view.py',
        detail: 'From here on the CLI holds exactly what a notebook cell holds — ' +
                'so everything below is shared with the notebook path.'
      },
      {
        id: 'cli-build',
        kind: 'framejs',
        title: 'build_cli_url(obj, short=False)',
        sub: 'the one CLI entry point into the framejs layer',
        file: 'visualization/framejs.py',
        detail: 'Returns a URL that renders this object, for webbrowser.open(). ' +
                'Raises FramejsUnavailable if the object has no registered viz, ' +
                'which the CLI reports as a friendly "not viewable yet".'
      },
      {
        id: 'cli-registry',
        kind: 'core',
        shared: true,
        title: 'VIZ_REGISTRY  →  VizRef',
        sub: '_viz_for(obj) — the same lookup the notebook path makes',
        file: 'visualization/framejs.py',
        detail: 'The shared core. There is one registry, one serializer per ' +
                'class and one frame per class; the two paths differ only in how ' +
                'the payload is delivered.'
      },
      {
        id: 'cli-url',
        kind: 'core',
        shared: true,
        title: 'the base URL',
        sub: 'default: local #?js=<b64>…   ·   short=True: the pinned /j/<uuid>',
        file: 'visualization/framejs.py',
        detail: 'Size is the only reason to choose the pinned form, and it is ' +
                'the one place size can matter — a self-contained LigandNetwork ' +
                'URL inlines the frame\'s JavaScript at ~140 kB, against ~10 kB ' +
                'on /j/<uuid>: the difference between a link you can paste ' +
                'somewhere and one you cannot. short=True requires that viz to ' +
                'have been published (just publish-viz); only ligand_network has.'
      },
      {
        id: 'cli-inputs',
        kind: 'core',
        shared: true,
        title: '…&inputs=<base64(json(payload))>',
        sub: 'the same VizRef.payload(obj), appended to the hash',
        file: 'visualization/framejs.py',
        detail: 'The same payload the notebook pushes over the comm channel, ' +
                'base64-encoded into the hash instead. Appended inputs take ' +
                'priority over anything baked into the frame. Merged with & onto ' +
                'the local form\'s existing #?js=…, or as a fresh #? on /j/<uuid>.'
      },
      {
        id: 'cli-open',
        kind: 'browser',
        title: 'webbrowser.open(url)',
        sub: 'the whole visualization travels in the link',
        file: 'openfecli/commands/view.py',
        detail: 'The URL is self-contained by default: it carries both the viz ' +
                'and the data, needs nothing from framejs.io\'s frame store, and ' +
                'cannot go stale against the installed gufe.'
      }
    ]
  }
};

// The registry, mirroring the table in visualization/README.md.
const REGISTRY = [
  ['LigandNetwork',          'ligand_network',           'radial network; click an edge to drive a 3D atom-mapping viewer', true],
  ['AlchemicalNetwork',      'alchemical_network',       'd3 force graph of ChemicalSystem nodes / Transformation edges', false],
  ['TransformationBase',     'transformation',           'stateA↔stateB component diff + the atom mapping', false],
  ['ChemicalSystem',         'chemical_system',          'master/detail over the system\'s labelled components', false],
  ['LigandAtomMapping',      'ligand_atom_mapping',      'the mapping viewer standalone (plain/colored/lines/overlay/2d)', false],
  ['SmallMoleculeComponent', 'small_molecule_component', '2D depiction + 3D conformer + SMILES/charge', false],
  ['ProteinComponent',       'protein_component',        '3Dmol with representation / colour-scheme switchers', false],
  ['SolventComponent',       'solvent_component',        'settings card (it has no coordinates)', false]
];

const REGISTRY_NOTE =
  'Registered on the class name and looked up over the MRO, so subclasses ' +
  'inherit: NonTransformation resolves to TransformationBase, and ' +
  'SolvatedPDBComponent / ProteinMembraneComponent to ProteinComponent. A ' +
  'frame directory is named after the gufe class it draws, snake_cased, so the ' +
  'two are greppable from each other. Only ligand_network has been published ' +
  'to a /j/<uuid>; the rest resolve from disk.';

// ============================================================================
// DOM SCAFFOLD
// ============================================================================
var mount = (typeof root !== 'undefined' && root) ? root : document.body;

mount.innerHTML = '';
mount.style.background = T.appBg;
mount.style.color = T.textPrimary;
mount.style.fontFamily = SANS;

const app = document.createElement('div');
app.style.cssText =
  'width:100%;height:100%;box-sizing:border-box;display:flex;flex-direction:column;' +
  'gap:10px;padding:14px 16px;overflow:hidden;background:' + T.appBg + ';';
mount.appendChild(app);

// ─── Header: title + path tabs ───
const header = document.createElement('div');
header.style.cssText = 'display:flex;flex-wrap:wrap;align-items:baseline;gap:10px 14px;flex-shrink:0;';
app.appendChild(header);

const titleEl = document.createElement('div');
titleEl.textContent = 'gufe → framejs';
titleEl.style.cssText =
  'font-weight:700;font-size:17px;color:' + T.titleColor + ';letter-spacing:.02em;';
header.appendChild(titleEl);

const subtitleEl = document.createElement('div');
subtitleEl.textContent = 'how an object becomes an interactive view';
subtitleEl.style.cssText = 'font-size:12px;color:' + T.textMuted2 + ';flex:1;min-width:120px;';
header.appendChild(subtitleEl);

const tabRow = document.createElement('div');
tabRow.style.cssText = 'display:flex;gap:6px;flex-shrink:0;';
header.appendChild(tabRow);

// ─── Tagline for the active tab ───
const taglineEl = document.createElement('div');
taglineEl.style.cssText =
  'font-size:12px;color:' + T.textMuted + ';flex-shrink:0;line-height:1.45;';
app.appendChild(taglineEl);

// ─── Body: the flow, and the detail panel ───
const body = document.createElement('div');
body.style.cssText = 'flex:1;min-height:0;display:flex;gap:14px;';
app.appendChild(body);

const flowPane = document.createElement('div');
flowPane.style.cssText =
  'flex:1;min-width:0;min-height:0;overflow-y:auto;overflow-x:hidden;' +
  'display:flex;flex-direction:column;align-items:stretch;padding-right:4px;';
body.appendChild(flowPane);

const detailPane = document.createElement('div');
detailPane.style.cssText =
  'width:300px;flex-shrink:0;min-height:0;overflow-y:auto;box-sizing:border-box;' +
  'padding:12px 14px;border-radius:10px;background:' + T.panelBg + ';' +
  'border:1px solid ' + T.cardBorder + ';';
body.appendChild(detailPane);

// ============================================================================
// SMALL BUILDERS
// ============================================================================

function makeTab(key, label, icon) {
  const b = document.createElement('button');
  b.textContent = icon + '  ' + label;
  b.dataset.key = key;
  b.style.cssText =
    'border:1px solid ' + T.chipBorder + ';border-radius:999px;padding:5px 14px;' +
    'font-family:' + SANS + ';font-size:12px;font-weight:600;cursor:pointer;' +
    'background:' + T.tabBg + ';color:' + T.textMuted + ';white-space:nowrap;';
  b.addEventListener('click', function() { selectTab(key); });
  tabRow.appendChild(b);
  return b;
}

function paintTabs() {
  Array.prototype.forEach.call(tabRow.children, function(b) {
    const on = b.dataset.key === activeTab;
    b.style.background = on ? T.tabOnBg : T.tabBg;
    b.style.color = on ? T.tabOnFg : T.textMuted;
    b.style.borderColor = on ? T.tabOnBg : T.chipBorder;
  });
}

function makeConnector(label) {
  const wrap = document.createElement('div');
  wrap.style.cssText =
    'display:flex;align-items:center;justify-content:center;gap:8px;' +
    'padding:3px 0;flex-shrink:0;';
  const line = document.createElement('div');
  line.textContent = '↓';
  line.style.cssText = 'color:' + T.connector + ';font-size:15px;line-height:1;';
  wrap.appendChild(line);
  if (label) {
    const l = document.createElement('div');
    l.textContent = label;
    l.style.cssText = 'font-size:10.5px;color:' + T.textMuted2 + ';font-style:italic;';
    wrap.appendChild(l);
  }
  flowPane.appendChild(wrap);
  return wrap;
}

function makeStageCard(step) {
  const card = document.createElement('div');
  card.dataset.id = step.id;
  card.style.cssText =
    'box-sizing:border-box;display:flex;gap:10px;align-items:flex-start;cursor:pointer;' +
    'padding:9px 12px;border-radius:9px;flex-shrink:0;' +
    'background:' + T.cardBg + ';border:1px solid ' + T.cardBorder + ';' +
    'border-left:4px solid ' + T.kinds[step.kind] + ';' +
    (step.dashed ? 'border-style:dashed;border-left-style:solid;' : '');

  const textCol = document.createElement('div');
  textCol.style.cssText = 'flex:1;min-width:0;display:flex;flex-direction:column;gap:2px;';
  card.appendChild(textCol);

  const t = document.createElement('div');
  t.textContent = step.title;
  t.style.cssText =
    'font-family:' + MONO + ';font-size:12.5px;font-weight:600;' +
    'color:' + T.textPrimary + ';overflow-wrap:anywhere;';
  textCol.appendChild(t);

  const s = document.createElement('div');
  s.textContent = step.sub;
  s.style.cssText =
    'font-size:11px;color:' + T.textMuted + ';line-height:1.4;overflow-wrap:anywhere;';
  textCol.appendChild(s);

  if (step.shared) {
    const badge = document.createElement('div');
    badge.textContent = 'shared core';
    badge.style.cssText =
      'align-self:flex-start;flex-shrink:0;font-size:9px;font-weight:700;' +
      'text-transform:uppercase;letter-spacing:.07em;padding:2px 7px;border-radius:999px;' +
      'background:' + T.chipBg + ';border:1px solid ' + T.kinds.core + ';color:' + T.kinds.core + ';';
    card.appendChild(badge);
  }

  card.addEventListener('click', function() { selectStep(step); });
  flowPane.appendChild(card);
  return card;
}

function paintCards() {
  Array.prototype.forEach.call(flowPane.querySelectorAll('[data-id]'), function(c) {
    const on = c.dataset.id === (activeStep && activeStep.id);
    // Only the three non-left sides — the left border is the kind colour band,
    // and setting `borderColor` wholesale would wipe it.
    const color = on ? T.cardBorderOn : T.cardBorder;
    c.style.borderTopColor = color;
    c.style.borderRightColor = color;
    c.style.borderBottomColor = color;
    c.style.boxShadow = on ? '0 0 0 1px ' + T.cardBorderOn : 'none';
  });
}

function codeSpan(text) {
  return '<span style="font-family:' + MONO + ';font-size:11px;background:' + T.codeBg +
         ';padding:1px 5px;border-radius:4px;overflow-wrap:anywhere;">' + esc(text) + '</span>';
}

function esc(s) {
  return String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
}

// ============================================================================
// PANES
// ============================================================================

let activeTab = 'notebook';
let activeStep = null;

function renderDetailPlaceholder() {
  const path = PATHS[activeTab];
  detailPane.innerHTML =
    '<div style="font-size:10px;font-weight:700;text-transform:uppercase;letter-spacing:.07em;' +
    'color:' + T.textMuted2 + ';margin-bottom:8px;">' + esc(path ? path.label + ' path' : 'registry') + '</div>' +
    '<div style="font-size:12px;line-height:1.55;color:' + T.textMuted + ';">' +
    (path
      ? 'Click any stage to see what it does and where it lives.<br><br>' +
        'The four <span style="color:' + T.kinds.core + ';font-weight:600;">shared core</span> ' +
        'stages are byte-identical between the two paths — one registry, one ' +
        'serializer per class, one frame per class. The paths differ only in ' +
        'how the payload reaches the frame: live over a comm channel, or ' +
        'base64ed into the URL.'
      : esc(REGISTRY_NOTE)) +
    '</div>';
}

function selectStep(step) {
  activeStep = step;
  const rows = [];
  rows.push(
    '<div style="font-family:' + MONO + ';font-size:13px;font-weight:700;color:' + T.titleColor +
    ';overflow-wrap:anywhere;margin-bottom:6px;">' + esc(step.title) + '</div>'
  );
  if (step.when) {
    rows.push(
      '<div style="font-size:11px;color:' + T.kinds.fallback + ';font-style:italic;' +
      'line-height:1.45;margin-bottom:8px;">taken when: ' + esc(step.when) + '</div>'
    );
  }
  rows.push(
    '<div style="font-size:12px;line-height:1.6;color:' + T.textPrimary + ';">' +
    esc(step.detail) + '</div>'
  );
  rows.push(
    '<div style="margin-top:12px;padding-top:10px;border-top:1px solid ' + T.cardBorder + ';">' +
    '<div style="font-size:9.5px;font-weight:700;text-transform:uppercase;letter-spacing:.07em;' +
    'color:' + T.textMuted2 + ';margin-bottom:5px;">lives in</div>' +
    codeSpan(step.file) + '</div>'
  );
  detailPane.innerHTML = rows.join('');
  paintCards();
}

function renderPath(key) {
  const path = PATHS[key];
  flowPane.innerHTML = '';
  path.steps.forEach(function(step, i) {
    if (i > 0) makeConnector(step.arrow || '');
    makeStageCard(step);
  });
  if (path.fallback) {
    const wrap = document.createElement('div');
    wrap.style.cssText =
      'display:flex;align-items:center;justify-content:center;gap:8px;padding:8px 0 3px;flex-shrink:0;';
    const l = document.createElement('div');
    l.textContent = '⤷  or, at any stage above';
    l.style.cssText = 'font-size:10.5px;color:' + T.textMuted2 + ';font-style:italic;';
    wrap.appendChild(l);
    flowPane.appendChild(wrap);
    makeStageCard(Object.assign({ dashed: true }, path.fallback));
  }
  paintCards();
}

function renderRegistry() {
  flowPane.innerHTML = '';
  const table = document.createElement('div');
  table.style.cssText = 'display:flex;flex-direction:column;gap:6px;';
  flowPane.appendChild(table);

  REGISTRY.forEach(function(row) {
    const cls = row[0], frame = row[1], shows = row[2], published = row[3];
    const card = document.createElement('div');
    card.dataset.id = 'reg-' + frame;
    card.style.cssText =
      'box-sizing:border-box;padding:9px 12px;border-radius:9px;cursor:pointer;flex-shrink:0;' +
      'background:' + T.cardBg + ';border:1px solid ' + T.cardBorder + ';' +
      'border-left:4px solid ' + T.kinds.gufe + ';';
    card.innerHTML =
      '<div style="display:flex;flex-wrap:wrap;align-items:baseline;gap:8px;">' +
        '<span style="font-family:' + MONO + ';font-size:12.5px;font-weight:600;color:' +
          T.textPrimary + ';">' + esc(cls) + '</span>' +
        '<span style="color:' + T.textMuted2 + ';font-size:11px;">→</span>' +
        '<span style="font-family:' + MONO + ';font-size:11.5px;color:' + T.kinds.framejs +
          ';">viz_assets/' + esc(frame) + '/</span>' +
        (published
          ? '<span style="font-size:9px;font-weight:700;text-transform:uppercase;letter-spacing:.06em;' +
            'padding:2px 7px;border-radius:999px;background:' + T.chipBg + ';border:1px solid ' +
            T.kinds.core + ';color:' + T.kinds.core + ';">published</span>'
          : '') +
      '</div>' +
      '<div style="font-size:11px;color:' + T.textMuted + ';margin-top:3px;line-height:1.4;">' +
        esc(shows) + '</div>';
    card.addEventListener('click', function() {
      selectStep({
        id: 'reg-' + frame,
        title: cls,
        sub: shows,
        detail: shows.charAt(0).toUpperCase() + shows.slice(1) + '. ' +
                (published
                  ? 'This is the one viz published to a canonical framejs.io ' +
                    '/j/<uuid>, so build_cli_url(obj, short=True) works for it.'
                  : 'Not published, so it always resolves from its on-disk frame ' +
                    'directory — which is the default everywhere anyway.'),
        file: 'src/gufe/visualization/viz_assets/' + frame + '/code.js'
      });
    });
    table.appendChild(card);
  });

  const note = document.createElement('div');
  note.textContent = REGISTRY_NOTE;
  note.style.cssText =
    'font-size:11px;color:' + T.textMuted2 + ';line-height:1.55;font-style:italic;' +
    'margin-top:10px;padding-top:10px;border-top:1px solid ' + T.cardBorder + ';';
  table.appendChild(note);

  paintCards();
}

function selectTab(key) {
  activeTab = key;
  activeStep = null;
  paintTabs();
  if (key === 'registry') {
    taglineEl.textContent =
      'One entry, one frame directory, one serializer — the whole of "what renders".';
    renderRegistry();
  } else {
    taglineEl.textContent = PATHS[key].tagline;
    renderPath(key);
  }
  renderDetailPlaceholder();
  flowPane.scrollTop = 0;
}

// ============================================================================
// RESPONSIVE — stack the detail panel under the flow on narrow screens
// ============================================================================

function applyLayout() {
  const w = mount.getBoundingClientRect().width || window.innerWidth;
  const narrow = w < NARROW_PX;
  body.style.flexDirection = narrow ? 'column' : 'row';
  body.style.overflowY = narrow ? 'auto' : 'hidden';
  detailPane.style.width = narrow ? '100%' : '300px';
  detailPane.style.overflowY = narrow ? 'visible' : 'auto';
  flowPane.style.overflowY = narrow ? 'visible' : 'auto';
  flowPane.style.flex = narrow ? '0 0 auto' : '1';
}

// ============================================================================
// Module entry points
// ============================================================================

// This frame is its own content — it needs no data. `inputs.highlight` is
// honoured so a link can deep-link straight to one path or stage, e.g.
// #?inputs=<b64 of {"highlight":"cli"}>.
export function onInputs(inputs) {
  try {
    const h = inputs && inputs.highlight;
    if (!h) return;
    if (PATHS[h] || h === 'registry') { selectTab(h); return; }
    for (const key of Object.keys(PATHS)) {
      const hit = PATHS[key].steps.filter(function(s) { return s.id === h; })[0];
      if (hit) { selectTab(key); selectStep(hit); return; }
    }
  } catch (err) {
    console.warn('[architecture] onInputs failed:', err);
  }
}

export function onResize() {
  applyLayout();
}

export function cleanup() {
  // No timers, listeners on window, or viewers to tear down.
}

// ── boot ──
makeTab('notebook', PATHS.notebook.label, PATHS.notebook.icon);
makeTab('cli', PATHS.cli.label, PATHS.cli.icon);
makeTab('registry', 'Registry', '🗂');
applyLayout();
selectTab('notebook');
