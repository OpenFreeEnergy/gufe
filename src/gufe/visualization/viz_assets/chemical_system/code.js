// ============================================================================
// Chemical system browser — component list (left) + detail pane (right).
//
// Renders a gufe.ChemicalSystem: a labelled bag of Components. Each component
// is one card in the left list; selecting a card renders it on the right
// according to its type — SDF small molecules and PDB proteins get a 3Dmol
// viewer, solvents get a settings card.
//
// Payload:  { "system": { "name": <str>, "components": [ <descriptor>, ... ] } }
// A descriptor is { label, type, name } plus ONE of:
//   { sdf, smiles }                                            small molecules
//   { pdb }                                                    proteins
//   { smiles, positive_ion, negative_ion, neutralize, ion_concentration }
//   { error }                                       serialization failed
// ============================================================================

var ThreeDmol = $3Dmol;

// ============================================================================
// THEME — flip DARK_MODE to switch the entire app
// ============================================================================
var DARK_MODE = false;

var THEMES = {
  dark: {
    appBg:         '#1a1a2e',
    splitBorder:   '#2a4a7f',
    toolbarBg:     '#16213e',
    toolbarBorder: '#2a4a7f',
    titleColor:    '#7ecfff',
    textPrimary:   '#e6f3ff',
    textMuted:     '#9bb8d6',
    textMuted2:    '#7a96b8',
    cardBg:        '#16213e',
    cardBgHover:   '#1a2b4f',
    cardBgActive:  '#0f3460',
    cardBorder:    '#2a4a7f',
    cardBorderActive: '#7ecfff',
    codeBg:        '#0f1c33',
    viewerBg:      '0x1a1a2e',
    btnBg:         '#0f3460',
    btnBgHover:    '#1a4a8a',
    btnBgActive:   '#2a6ab5',
    btnFg:         '#7ecfff',
    btnBorder:     '#2a4a7f',
    warnBg:        '#3b1d1d',
    warnFg:        '#ffb4b4',
    warnBorder:    '#7f2a2a',
    loadingFg:     '#888',
    errorFg:       '#ff8080'
  },
  light: {
    appBg:         '#ffffff',
    splitBorder:   '#e2e8f0',
    toolbarBg:     '#f8fafc',
    toolbarBorder: '#e2e8f0',
    titleColor:    '#0369a1',
    textPrimary:   '#1e293b',
    textMuted:     '#475569',
    textMuted2:    '#64748b',
    cardBg:        '#ffffff',
    cardBgHover:   '#f1f5f9',
    cardBgActive:  '#e0f2fe',
    cardBorder:    '#e2e8f0',
    cardBorderActive: '#0369a1',
    codeBg:        '#f1f5f9',
    viewerBg:      '0xffffff',
    btnBg:         '#e1e8f2',
    btnBgHover:    '#c8d4e8',
    btnBgActive:   '#7ab0e5',
    btnFg:         '#1a4a8a',
    btnBorder:     '#a8bcd6',
    warnBg:        '#fee2e2',
    warnFg:        '#991b1b',
    warnBorder:    '#fecaca',
    loadingFg:     '#888',
    errorFg:       '#c33'
  }
};

var T = DARK_MODE ? THEMES.dark : THEMES.light;

// ============================================================================
// CONFIG
// ============================================================================
const CONFIG = {
  listWidth: 240,
  mol:   { stickRadius: 0.15, sphereScale: 0.25 },
  water: { stickRadius: 0.05, sphereScale: 0.18 },
  hetero:{ stickRadius: 0.20, sphereScale: 0.28 }
};

// Residue names treated as water (hidden by default in the protein view).
const WATER_RESN = ['HOH', 'WAT', 'SOL', 'TIP3'];

// Stable badge colour per component type. Unknown types fall back to a hash
// over the type name, so the same type always gets the same colour.
const TYPE_COLORS = {
  SmallMoleculeComponent: '#7c3aed',
  ProteinComponent:       '#0f766e',
  SolvatedPDBComponent:   '#0369a1',
  SolventComponent:       '#0284c7',
  ProteinMembraneComponent: '#b45309'
};
const FALLBACK_COLORS = ['#be123c', '#a16207', '#4d7c0f', '#1d4ed8', '#86198f', '#0f766e'];

function typeColor(type) {
  if (TYPE_COLORS[type]) return TYPE_COLORS[type];
  let h = 0;
  for (let i = 0; i < type.length; i++) h = (h * 31 + type.charCodeAt(i)) >>> 0;
  return FALLBACK_COLORS[h % FALLBACK_COLORS.length];
}

// Short, human word for a type — used in the composition summary line.
function typeWord(type) {
  if (type === 'SmallMoleculeComponent') return 'ligand';
  if (type === 'SolventComponent') return 'solvent';
  if (type.indexOf('Protein') !== -1) return 'protein';
  return type.replace(/Component$/, '').toLowerCase() || type;
}

// ============================================================================
// HELPERS
// ============================================================================

// Input values may be string | Blob | ArrayBuffer | TypedArray.
async function asText(v) {
  if (v == null) return null;
  if (v instanceof Blob) return await v.text();
  if (v instanceof ArrayBuffer) return new TextDecoder().decode(v);
  if (ArrayBuffer.isView(v)) return new TextDecoder().decode(v);
  return typeof v === 'string' ? v : String(v);
}

// Object-valued inputs may arrive as a JSON string.
async function asObject(v) {
  if (v == null) return null;
  if (typeof v === 'object' && !(v instanceof Blob) && !(v instanceof ArrayBuffer) && !ArrayBuffer.isView(v)) {
    return v;
  }
  const text = await asText(v);
  if (!text) return null;
  return JSON.parse(text);
}

function fmt(n) {
  return n.toLocaleString('en-US');
}

/**
 * Single pass over the PDB text for the chain / residue / atom readout.
 * Only the first MODEL is counted (NMR ensembles would otherwise multiply
 * every number). Residues are keyed by chain + sequence number + insertion
 * code + name, which is what makes them unique in a PDB.
 */
function parsePdbStats(pdbText) {
  const chains = new Set();
  const residues = new Set();
  let atoms = 0, hetatms = 0, waters = 0;

  const lines = pdbText.split(/\r?\n/);
  for (let i = 0; i < lines.length; i++) {
    const line = lines[i];
    const rec = line.slice(0, 6);
    if (rec === 'ENDMDL') break;                 // first model only
    const isAtom = (rec === 'ATOM  ');
    const isHet  = (rec === 'HETATM');
    if (!isAtom && !isHet) continue;

    atoms++;
    if (isHet) hetatms++;

    const resName = line.slice(17, 20).trim();
    const chainId = line.slice(21, 22).trim() || '_';
    const resSeq  = line.slice(22, 26).trim();
    const iCode   = line.slice(26, 27).trim();

    if (WATER_RESN.indexOf(resName) !== -1) waters++;

    chains.add(chainId);
    residues.add(chainId + '|' + resSeq + iCode + '|' + resName);
  }

  return { chains: chains.size, residues: residues.size, atoms, hetatms, waters };
}

function el(tag, css, text) {
  const n = document.createElement(tag);
  if (css) n.style.cssText = css;
  if (text != null) n.textContent = text;
  return n;
}

// ============================================================================
// DOM SCAFFOLD: header strip + master/detail split
// ============================================================================
var mount = (typeof root !== 'undefined' && root) ? root : document.body;

mount.innerHTML = '';
mount.style.background = T.appBg;
mount.style.fontFamily = "'Inter',system-ui,sans-serif";

const appWrap = el('div',
  'width:100%;height:100%;display:flex;flex-direction:column;overflow:hidden;background:' + T.appBg + ';');
mount.appendChild(appWrap);

// ─── Header strip ───
const header = el('div',
  'display:flex;align-items:baseline;gap:14px;padding:9px 16px;background:' + T.toolbarBg +
  ';border-bottom:1px solid ' + T.toolbarBorder + ';flex-shrink:0;flex-wrap:wrap;');
appWrap.appendChild(header);

const sysNameEl = el('span',
  'font-weight:700;font-size:15px;color:' + T.titleColor + ';letter-spacing:.02em;', 'Chemical System');
header.appendChild(sysNameEl);

const summaryEl = el('span', 'font-size:12px;color:' + T.textMuted2 + ';');
header.appendChild(summaryEl);

// ─── Split ───
const splitWrap = el('div', 'flex:1;display:flex;flex-direction:row;min-height:0;overflow:hidden;');
appWrap.appendChild(splitWrap);

const listPane = el('div',
  'width:' + CONFIG.listWidth + 'px;flex:0 0 ' + CONFIG.listWidth + 'px;height:100%;overflow-y:auto;' +
  'padding:10px;box-sizing:border-box;display:flex;flex-direction:column;gap:8px;background:' + T.appBg + ';');
splitWrap.appendChild(listPane);

splitWrap.appendChild(el('div', 'width:1px;background:' + T.splitBorder + ';flex-shrink:0;'));

const detailPane = el('div',
  'flex:1;min-width:0;min-height:0;position:relative;display:flex;flex-direction:column;background:' + T.appBg + ';');
splitWrap.appendChild(detailPane);

// ============================================================================
// STATE
// ============================================================================
let components = [];
let selectedIdx = -1;
let cardEls = [];
let viewer = null;          // active 3Dmol viewer, if the detail pane has one
let lastPayloadJson = null; // avoid rebuilding on an identical re-delivery

function destroyViewer() {
  if (!viewer) return;
  try { viewer.clear(); } catch (e) {}
  viewer = null;
}

// ============================================================================
// LEFT: component cards
// ============================================================================
function cardCss(active) {
  return 'display:flex;flex-direction:column;gap:4px;padding:8px 10px;border-radius:8px;cursor:pointer;' +
    'border:1px solid ' + (active ? T.cardBorderActive : T.cardBorder) + ';' +
    'background:' + (active ? T.cardBgActive : T.cardBg) + ';' +
    'flex-shrink:0;transition:background .12s;';
}

function refreshCardStyles() {
  cardEls.forEach((c, i) => { c.style.cssText = cardCss(i === selectedIdx); });
}

function buildList() {
  listPane.innerHTML = '';
  cardEls = [];

  if (components.length === 0) {
    listPane.appendChild(el('div',
      'font-size:12px;color:' + T.textMuted + ';padding:6px;', 'No components.'));
    return;
  }

  components.forEach((comp, i) => {
    const card = el('div', cardCss(false));

    card.appendChild(el('div',
      'font-size:13px;font-weight:700;color:' + T.textPrimary + ';word-break:break-word;',
      comp.label || '(unlabelled)'));

    const badgeRow = el('div', 'display:flex;align-items:center;gap:6px;flex-wrap:wrap;');
    const type = comp.type || 'Component';
    badgeRow.appendChild(el('span',
      'font-size:10px;font-weight:700;letter-spacing:.02em;padding:2px 6px;border-radius:4px;' +
      'color:#fff;background:' + typeColor(type) + ';white-space:nowrap;',
      type.replace(/Component$/, '')));
    if (comp.error) {
      badgeRow.appendChild(el('span', 'font-size:11px;color:' + T.errorFg + ';', '⚠'));
    }
    card.appendChild(badgeRow);

    if (comp.name) {
      card.appendChild(el('div',
        'font-size:11px;color:' + T.textMuted + ';word-break:break-word;', comp.name));
    }

    card.onmouseover = () => { if (i !== selectedIdx) card.style.background = T.cardBgHover; };
    card.onmouseout  = () => { if (i !== selectedIdx) card.style.background = T.cardBg; };
    card.onclick     = () => select(i);

    listPane.appendChild(card);
    cardEls.push(card);
  });
}

// ============================================================================
// RIGHT: detail renderers
// ============================================================================

function detailHeader(comp) {
  const wrap = el('div',
    'display:flex;align-items:baseline;gap:10px;padding:9px 16px;flex-wrap:wrap;flex-shrink:0;' +
    'border-bottom:1px solid ' + T.toolbarBorder + ';background:' + T.toolbarBg + ';');
  wrap.appendChild(el('span',
    'font-size:14px;font-weight:700;color:' + T.textPrimary + ';', comp.label || '(unlabelled)'));
  const type = comp.type || 'Component';
  wrap.appendChild(el('span',
    'font-size:10px;font-weight:700;padding:2px 6px;border-radius:4px;color:#fff;background:' + typeColor(type) + ';',
    type));
  if (comp.name) {
    wrap.appendChild(el('span', 'font-size:12px;color:' + T.textMuted + ';', comp.name));
  }
  return wrap;
}

function warnBanner(message) {
  return el('div',
    'margin:12px 16px;padding:8px 12px;border-radius:6px;font-size:12px;white-space:pre-wrap;' +
    'background:' + T.warnBg + ';color:' + T.warnFg + ';border:1px solid ' + T.warnBorder + ';',
    '⚠ ' + message);
}

function fieldRow(label, value) {
  const row = el('div', 'display:flex;gap:10px;align-items:baseline;padding:5px 0;');
  row.appendChild(el('div',
    'flex:0 0 130px;font-size:11px;color:' + T.textMuted2 + ';text-transform:uppercase;letter-spacing:.04em;',
    label));
  row.appendChild(el('div',
    'flex:1;min-width:0;font-size:13px;color:' + T.textPrimary + ';word-break:break-all;' +
    "font-family:ui-monospace,SFMono-Regular,Menlo,monospace;",
    value == null || value === '' ? '—' : String(value)));
  return row;
}

/** A one-line SMILES strip shown under a molecule viewer. */
function smilesStrip(smiles) {
  const wrap = el('div',
    'flex-shrink:0;padding:8px 16px;border-top:1px solid ' + T.toolbarBorder +
    ';background:' + T.toolbarBg + ';display:flex;gap:8px;align-items:baseline;');
  wrap.appendChild(el('span',
    'font-size:11px;color:' + T.textMuted2 + ';text-transform:uppercase;letter-spacing:.04em;', 'SMILES'));
  wrap.appendChild(el('div',
    'flex:1;min-width:0;font-size:12px;color:' + T.textPrimary + ';overflow-x:auto;white-space:nowrap;' +
    "font-family:ui-monospace,SFMono-Regular,Menlo,monospace;",
    smiles || '—'));
  return wrap;
}

/** A stats strip (chains / residues / atoms) shown under the protein viewer. */
function statsStrip(stats) {
  const text =
    fmt(stats.chains) + ' chains · ' +
    fmt(stats.residues) + ' residues · ' +
    fmt(stats.atoms) + ' atoms · ' +
    fmt(stats.hetatms) + ' HETATM' +
    (stats.waters > 0 ? ' (' + fmt(stats.waters) + ' water, hidden)' : '');
  return el('div',
    'flex-shrink:0;padding:8px 16px;border-top:1px solid ' + T.toolbarBorder +
    ';background:' + T.toolbarBg + ';font-size:12px;color:' + T.textMuted2 + ';',
    text);
}

/** Shared viewer host: a flexible box the 3Dmol canvas is created inside. */
function viewerHost() {
  const wrap = el('div', 'flex:1;position:relative;min-height:0;min-width:0;');
  const container = el('div', 'position:absolute;inset:0;');
  wrap.appendChild(container);
  return { wrap, container };
}

function renderSmallMolecule(comp, sdf) {
  const host = viewerHost();
  detailPane.appendChild(host.wrap);
  detailPane.appendChild(smilesStrip(comp.smiles));

  try {
    viewer = ThreeDmol.createViewer(host.container, { backgroundColor: T.viewerBg });
    viewer.addModel(sdf, 'sdf');
    viewer.setStyle({}, {
      stick:  { radius: CONFIG.mol.stickRadius, colorscheme: 'Jmol' },
      sphere: { scale: CONFIG.mol.sphereScale,  colorscheme: 'Jmol' }
    });
    viewer.zoomTo();
    viewer.render();
  } catch (e) {
    destroyViewer();
    host.container.appendChild(warnBanner('Failed to render molecule: ' + e.message));
  }
}

function renderProtein(comp, pdb) {
  let stats = null;
  try { stats = parsePdbStats(pdb); } catch (e) { stats = null; }

  const host = viewerHost();
  detailPane.appendChild(host.wrap);
  if (stats) detailPane.appendChild(statsStrip(stats));

  try {
    viewer = ThreeDmol.createViewer(host.container, { backgroundColor: T.viewerBg });
    viewer.addModel(pdb, 'pdb');
    // Order matters: each setStyle replaces the style of the atoms it matches,
    // so the water rule (last) wins over the blanket hetero rule.
    viewer.setStyle({}, {});
    viewer.setStyle({ hetflag: false }, { cartoon: { colorscheme: 'chain' } });
    viewer.setStyle({ hetflag: true }, {
      stick:  { radius: CONFIG.hetero.stickRadius, colorscheme: 'Jmol' },
      sphere: { scale: CONFIG.hetero.sphereScale,  colorscheme: 'Jmol' }
    });
    viewer.setStyle({ resn: WATER_RESN }, {});   // waters hidden by default
    viewer.zoomTo();
    viewer.render();
  } catch (e) {
    destroyViewer();
    host.container.appendChild(warnBanner('Failed to render structure: ' + e.message));
  }
}

function renderSolvent(comp) {
  const scroll = el('div', 'flex:1;overflow-y:auto;padding:16px;min-height:0;');
  const card = el('div',
    'max-width:520px;padding:14px 16px;border-radius:10px;border:1px solid ' + T.cardBorder +
    ';background:' + T.cardBg + ';');

  card.appendChild(el('div',
    'font-size:12px;font-weight:700;color:' + T.titleColor + ';margin-bottom:6px;letter-spacing:.03em;',
    'SOLVENT SETTINGS'));

  card.appendChild(fieldRow('SMILES', comp.smiles));
  card.appendChild(fieldRow('Positive ion', comp.positive_ion));
  card.appendChild(fieldRow('Negative ion', comp.negative_ion));
  card.appendChild(fieldRow('Concentration', comp.ion_concentration));
  card.appendChild(fieldRow('Neutralize', comp.neutralize == null ? null : (comp.neutralize ? 'yes' : 'no')));

  scroll.appendChild(card);
  detailPane.appendChild(scroll);
}

/** Anything with no sdf/pdb/solvent fields — show whatever scalars we got. */
function renderGeneric(comp) {
  const scroll = el('div', 'flex:1;overflow-y:auto;padding:16px;min-height:0;');
  const card = el('div',
    'max-width:520px;padding:14px 16px;border-radius:10px;border:1px solid ' + T.cardBorder +
    ';background:' + T.cardBg + ';');
  card.appendChild(el('div',
    'font-size:12px;color:' + T.textMuted + ';margin-bottom:8px;',
    'No structural data was serialized for this component.'));
  Object.keys(comp).forEach(k => {
    if (k === 'label' || k === 'type' || k === 'error') return;
    const v = comp[k];
    if (v == null || typeof v === 'object') return;
    card.appendChild(fieldRow(k.replace(/_/g, ' '), v));
  });
  scroll.appendChild(card);
  detailPane.appendChild(scroll);
}

function renderDetail(comp) {
  destroyViewer();
  detailPane.innerHTML = '';

  if (!comp) {
    detailPane.appendChild(el('div',
      'flex:1;display:flex;align-items:center;justify-content:center;font-size:13px;color:' + T.textMuted + ';',
      'Select a component.'));
    return;
  }

  detailPane.appendChild(detailHeader(comp));

  if (comp.error) {
    detailPane.appendChild(warnBanner(comp.error));
    // Fall through: some fields may still have made it across.
  }

  if (comp.sdf) {
    renderSmallMolecule(comp, comp.sdf);
  } else if (comp.pdb) {
    renderProtein(comp, comp.pdb);
  } else if (comp.smiles != null || comp.ion_concentration != null || comp.neutralize != null) {
    renderSolvent(comp);
  } else if (!comp.error) {
    renderGeneric(comp);
  } else {
    detailPane.appendChild(el('div', 'flex:1;'));   // banner only
  }
}

function select(i) {
  if (i < 0 || i >= components.length) return;
  selectedIdx = i;
  refreshCardStyles();
  renderDetail(components[i]);
}

// ============================================================================
// HEADER SUMMARY
// ============================================================================
function compositionSummary(comps) {
  if (comps.length === 0) return 'no components';
  const counts = new Map();
  comps.forEach(c => {
    const w = typeWord(c.type || 'Component');
    counts.set(w, (counts.get(w) || 0) + 1);
  });
  const parts = [];
  counts.forEach((n, w) => parts.push(n + ' ' + w + (n > 1 ? 's' : '')));
  return comps.length + ' component' + (comps.length > 1 ? 's' : '') + ' — ' + parts.join(', ');
}

// ============================================================================
// INIT
// ============================================================================
function init(system) {
  destroyViewer();

  if (!system || typeof system !== 'object') {
    components = [];
    selectedIdx = -1;
    sysNameEl.textContent = 'Chemical System';
    summaryEl.textContent = 'waiting for data…';
    buildList();
    renderDetail(null);
    return;
  }

  components = Array.isArray(system.components) ? system.components.filter(c => c && typeof c === 'object') : [];
  sysNameEl.textContent = system.name || 'Chemical System';
  summaryEl.textContent = compositionSummary(components);

  buildList();
  selectedIdx = -1;
  if (components.length > 0) {
    select(0);                 // auto-select the first component
  } else {
    renderDetail(null);
  }
}

// Keep the canvas in step with the pane (the frame can be resized by its host
// without onResize firing).
const detailRO = new ResizeObserver(() => {
  if (viewer) { viewer.resize(); viewer.render(); }
});
detailRO.observe(detailPane);

init(null);

// ============================================================================
// Module entry points
// ============================================================================
export async function onInputs(inputs) {
  if (!inputs) return;
  try {
    const system = await asObject(inputs['system']);
    const json = JSON.stringify(system);
    if (json === lastPayloadJson) return;   // identical re-delivery
    lastPayloadJson = json;
    init(system);
  } catch (e) {
    destroyViewer();
    detailPane.innerHTML = '';
    detailPane.appendChild(warnBanner('Failed to read system payload: ' + e.message));
    if (typeof logStderr === 'function') logStderr('chemical system onInputs failed: ' + e.message);
  }
}

export function onResize() {
  if (viewer) { viewer.resize(); viewer.render(); }
}

export function cleanup() {
  try { detailRO.disconnect(); } catch (e) {}
  destroyViewer();
}
