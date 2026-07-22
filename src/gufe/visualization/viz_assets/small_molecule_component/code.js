// ============================================================================
// SmallMoleculeComponent viewer.
//
// Header strip with the molecule name, then a 50/50 split:
//   LEFT  — 2D depiction from the SDF molblock, via RDKit-js.
//   RIGHT — 3Dmol.js viewer with a stick / ball+stick / sphere switcher and a
//           spin toggle.
// A compact info bar along the bottom shows name, SMILES, total charge and the
// atom/bond counts read off the SDF counts line.
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
    labelBg:       '#16213e',
    labelFg:       '#7ecfff',
    switcherBg:    'rgba(22,33,62,0.9)',
    btnBg:         '#0f3460',
    btnBgHover:    '#1a4a8a',
    btnBgActive:   '#2a6ab5',
    btnFg:         '#7ecfff',
    btnBorder:     '#2a4a7f',
    viewerBg:      '0x1a1a2e',
    canvas2DBg:    '#fff',
    loadingFg:     '#888',
    errorFg:       '#c33'
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
    labelBg:       '#f0f4fa',
    labelFg:       '#1a4a8a',
    switcherBg:    'rgba(240,244,250,0.95)',
    btnBg:         '#e1e8f2',
    btnBgHover:    '#c8d4e8',
    btnBgActive:   '#7ab0e5',
    btnFg:         '#1a4a8a',
    btnBorder:     '#a8bcd6',
    viewerBg:      '0xffffff',
    canvas2DBg:    '#fff',
    loadingFg:     '#888',
    errorFg:       '#c33'
  }
};

var T = DARK_MODE ? THEMES.dark : THEMES.light;

// 3D style presets offered by the switcher.
const STYLE_MODES = [
  { id: 'stick',  label: 'Stick',      title: 'Sticks only' },
  { id: 'ball',   label: 'Ball+Stick', title: 'Ball and stick' },
  { id: 'sphere', label: 'Sphere',     title: 'Space-filling spheres' }
];

const STYLE_SPECS = {
  stick:  { stick: { radius: 0.15, colorscheme: 'Jmol' } },
  ball:   { stick: { radius: 0.12, colorscheme: 'Jmol' },
            sphere: { scale: 0.28, colorscheme: 'Jmol' } },
  sphere: { sphere: { scale: 1.0, colorscheme: 'Jmol' } }
};

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

// RDKit loader (singleton) — copied from ligand_network/code.js
let rdkitPromise = null;
function loadRDKit() {
  if (rdkitPromise) return rdkitPromise;
  rdkitPromise = new Promise(function(resolve, reject) {
    if (window.RDKit) { resolve(window.RDKit); return; }
    const script = document.createElement('script');
    script.src = 'https://unpkg.com/@rdkit/rdkit/dist/RDKit_minimal.js';
    script.onload = function() {
      window.initRDKitModule().then(function(rdkit) {
        window.RDKit = rdkit;
        resolve(rdkit);
      }).catch(reject);
    };
    script.onerror = function() { reject(new Error('Failed to load RDKit')); };
    document.head.appendChild(script);
  });
  return rdkitPromise;
}

// 2D depiction SVG — adapted from ligand_network/code.js to take a molblock
// string directly (the payload already hands us one).
function depictMoleculeSVG(RDKit, molblock, size) {
  let rdmol = null;
  try {
    rdmol = RDKit.get_mol(molblock, JSON.stringify({ removeHs: true }));
    if (!rdmol) return null;
    try { rdmol.set_new_coords(true); } catch (e) {}
    const svg = rdmol.get_svg(size, size);
    return svg || null;
  } catch (e) {
    console.warn('[smc] depictMoleculeSVG threw -', e.message);
    return null;
  } finally {
    if (rdmol) { try { rdmol.delete(); } catch (e) {} }
  }
}

// The V2000 counts line is the 4th line: aaabbb... (3 chars each).
function parseCounts(sdf) {
  try {
    const lines = sdf.split(/\r?\n/);
    if (lines.length < 4) return null;
    const counts = lines[3];
    const nAtoms = parseInt(counts.slice(0, 3), 10);
    const nBonds = parseInt(counts.slice(3, 6), 10);
    if (isNaN(nAtoms) || isNaN(nBonds)) return null;
    return { atoms: nAtoms, bonds: nBonds };
  } catch (e) {
    return null;
  }
}

// 3Dmol wants the "$$$$" record terminator; to_sdf() usually includes it.
function ensureSDFTerminator(sdf) {
  const DOLLAR = String.fromCharCode(36);
  const term = DOLLAR + DOLLAR + DOLLAR + DOLLAR;
  return sdf.indexOf(term) >= 0 ? sdf : sdf + String.fromCharCode(10) + term;
}

// ============================================================================
// DOM SCAFFOLD: header / split / info bar
// ============================================================================
const mount = (typeof root !== 'undefined' && root) ? root : document.body;
mount.innerHTML = '';
mount.style.background = T.appBg;
mount.style.fontFamily = "'Inter',system-ui,sans-serif";
mount.style.color = T.textPrimary;

const shell = document.createElement('div');
shell.style.cssText = 'width:100%;height:100%;display:flex;flex-direction:column;overflow:hidden;background:' + T.appBg + ';';
mount.appendChild(shell);

// ─── Header ───
const header = document.createElement('div');
header.style.cssText = 'flex-shrink:0;padding:10px 16px;background:' + T.toolbarBg + ';border-bottom:1px solid ' + T.toolbarBorder + ';display:flex;align-items:baseline;gap:10px;';
shell.appendChild(header);

const titleEl = document.createElement('span');
titleEl.textContent = '—';
titleEl.style.cssText = 'font-weight:700;font-size:15px;color:' + T.titleColor + ';letter-spacing:.02em;';
header.appendChild(titleEl);

const subtitleEl = document.createElement('span');
subtitleEl.textContent = 'SmallMoleculeComponent';
subtitleEl.style.cssText = 'font-size:11px;color:' + T.textMuted2 + ';text-transform:uppercase;letter-spacing:.08em;';
header.appendChild(subtitleEl);

// ─── Split ───
const splitWrap = document.createElement('div');
splitWrap.style.cssText = 'flex:1;display:flex;flex-direction:row;overflow:hidden;min-height:0;';
shell.appendChild(splitWrap);

const leftPane = document.createElement('div');
leftPane.style.cssText = 'flex:1 1 50%;width:50%;height:100%;display:flex;flex-direction:column;min-width:0;min-height:0;';
splitWrap.appendChild(leftPane);

const divider = document.createElement('div');
divider.style.cssText = 'width:1px;background:' + T.splitBorder + ';flex-shrink:0;';
splitWrap.appendChild(divider);

const rightPane = document.createElement('div');
rightPane.style.cssText = 'flex:1 1 50%;width:50%;height:100%;display:flex;flex-direction:column;min-width:0;min-height:0;position:relative;';
splitWrap.appendChild(rightPane);

// LEFT: 2D depiction
const depictLabel = document.createElement('div');
depictLabel.textContent = '2D';
depictLabel.style.cssText = 'flex-shrink:0;padding:4px 10px;font-size:12px;font-weight:bold;color:' + T.labelFg + ';background:' + T.labelBg + ';';
leftPane.appendChild(depictLabel);

const depictBox = document.createElement('div');
depictBox.style.cssText = 'flex:1;min-height:0;display:flex;align-items:center;justify-content:center;background:' + T.canvas2DBg + ';overflow:hidden;padding:8px;';
leftPane.appendChild(depictBox);

// RIGHT: 3D viewer
const viewerLabel = document.createElement('div');
viewerLabel.textContent = '3D';
viewerLabel.style.cssText = 'flex-shrink:0;padding:4px 10px;font-size:12px;font-weight:bold;color:' + T.labelFg + ';background:' + T.labelBg + ';';
rightPane.appendChild(viewerLabel);

const viewerBox = document.createElement('div');
viewerBox.style.cssText = 'flex:1;min-height:0;position:relative;';
rightPane.appendChild(viewerBox);

const switcher = document.createElement('div');
switcher.style.cssText = 'position:absolute;bottom:10px;right:10px;display:flex;gap:4px;background:' + T.switcherBg + ';padding:4px;border-radius:6px;z-index:10;box-shadow:0 2px 8px rgba(0,0,0,0.25);';
rightPane.appendChild(switcher);

const BTN_CSS = 'background:' + T.btnBg + ';color:' + T.btnFg + ';border:1px solid ' + T.btnBorder +
  ';padding:4px 8px;font-size:11px;font-weight:bold;border-radius:3px;cursor:pointer;font-family:inherit;';

let currentStyle = 'stick';
let spinning = false;

STYLE_MODES.forEach(m => {
  const btn = document.createElement('button');
  btn.textContent = m.label;
  btn.title = m.title;
  btn.dataset.style = m.id;
  btn.style.cssText = BTN_CSS;
  btn.onmouseover = () => { btn.style.background = T.btnBgHover; };
  btn.onmouseout  = () => { btn.style.background = (currentStyle === m.id) ? T.btnBgActive : T.btnBg; };
  btn.onclick = () => setStyle(m.id);
  switcher.appendChild(btn);
});

const spinBtn = document.createElement('button');
spinBtn.textContent = 'Spin';
spinBtn.title = 'Toggle continuous rotation';
spinBtn.style.cssText = BTN_CSS + 'margin-left:4px;';
spinBtn.onmouseover = () => { spinBtn.style.background = T.btnBgHover; };
spinBtn.onmouseout  = () => { spinBtn.style.background = spinning ? T.btnBgActive : T.btnBg; };
spinBtn.onclick = () => setSpin(!spinning);
switcher.appendChild(spinBtn);

// ─── Info bar ───
const infoBar = document.createElement('div');
infoBar.style.cssText = 'flex-shrink:0;display:flex;flex-wrap:wrap;align-items:baseline;gap:6px 20px;padding:8px 16px;background:' + T.toolbarBg + ';border-top:1px solid ' + T.toolbarBorder + ';font-size:12px;color:' + T.textPrimary + ';';
shell.appendChild(infoBar);

// Builds one "LABEL value" cell. `mono` renders the value as selectable
// monospace text (used for SMILES).
function infoCell(label, value, mono) {
  const cell = document.createElement('div');
  cell.style.cssText = 'display:flex;align-items:baseline;gap:6px;min-width:0;';

  const k = document.createElement('span');
  k.textContent = label;
  k.style.cssText = 'font-size:10px;font-weight:700;letter-spacing:.08em;text-transform:uppercase;color:' + T.textMuted2 + ';flex-shrink:0;';
  cell.appendChild(k);

  const v = document.createElement('span');
  v.textContent = value;
  v.title = value;
  v.style.cssText = 'color:' + T.textPrimary + ';user-select:text;-webkit-user-select:text;cursor:text;' +
    (mono ? "font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:11px;overflow-wrap:anywhere;" : '');
  cell.appendChild(v);

  infoBar.appendChild(cell);
  return cell;
}

function renderInfoBar(name, smiles, charge, counts) {
  infoBar.innerHTML = '';
  infoCell('Name', name || '—', false);
  infoCell('SMILES', smiles || '—', true);
  infoCell('Charge', (charge === null || charge === undefined) ? '—' : String(charge), false);
  infoCell('Atoms', counts ? String(counts.atoms) : '—', false);
  infoCell('Bonds', counts ? String(counts.bonds) : '—', false);
}

// ============================================================================
// STATUS / PLACEHOLDER RENDERING
// ============================================================================
function showMessage(el, text, color) {
  el.innerHTML = '';
  const d = document.createElement('div');
  d.textContent = text;
  d.style.cssText = 'color:' + color + ';font-size:12px;padding:10px;text-align:center;';
  el.appendChild(d);
}

function showPlaceholder(text) {
  destroyViewer();
  showMessage(depictBox, text, T.loadingFg);
  showMessage(viewerBox, text, T.loadingFg);
  titleEl.textContent = '—';
  renderInfoBar(null, null, null, null);
}

// Error banner, mirroring the ligand_network init() pattern.
function showErrorBanner(message) {
  const warn = document.createElement('div');
  warn.textContent = '⚠ ' + message;
  warn.style.cssText = 'position:absolute;top:8px;left:50%;transform:translateX(-50%);background:#fee2e2;color:#991b1b;border:1px solid #fecaca;padding:6px 14px;border-radius:6px;font-size:12px;z-index:20;max-width:90%;';
  rightPane.appendChild(warn);
}

// ============================================================================
// 3D VIEWER
// ============================================================================
let viewer = null;

function destroyViewer() {
  if (viewer) {
    try { viewer.spin(false); } catch (e) {}
    try { viewer.clear(); } catch (e) {}
    viewer = null;
  }
  spinning = false;
  viewerBox.innerHTML = '';
}

function updateSwitcherHighlight() {
  switcher.querySelectorAll('button[data-style]').forEach(b => {
    b.style.background = (b.dataset.style === currentStyle) ? T.btnBgActive : T.btnBg;
  });
  spinBtn.style.background = spinning ? T.btnBgActive : T.btnBg;
}

function setStyle(styleId) {
  currentStyle = styleId;
  updateSwitcherHighlight();
  if (!viewer) return;
  viewer.setStyle({}, STYLE_SPECS[styleId] || STYLE_SPECS.stick);
  viewer.render();
}

function setSpin(on) {
  spinning = !!on;
  updateSwitcherHighlight();
  if (!viewer) return;
  try { viewer.spin(spinning ? 'y' : false); } catch (e) {}
}

function render3D(sdf) {
  destroyViewer();
  try {
    viewer = ThreeDmol.createViewer(viewerBox, { backgroundColor: T.viewerBg });
    viewer.addModel(ensureSDFTerminator(sdf), 'sdf');
    viewer.setStyle({}, STYLE_SPECS[currentStyle] || STYLE_SPECS.stick);
    viewer.zoomTo();
    viewer.render();
    updateSwitcherHighlight();
  } catch (e) {
    viewer = null;
    showMessage(viewerBox, '3D render failed: ' + e.message, T.errorFg);
  }
}

// ============================================================================
// 2D DEPICTION
// ============================================================================
const DEPICT_SIZE = 400;

// Bumped on every payload so a slow RDKit load can't overwrite a newer render.
let depictToken = 0;

function render2D(sdf) {
  showMessage(depictBox, 'Loading 2D depiction…', T.loadingFg);
  const token = ++depictToken;

  loadRDKit().then(RDKit => {
    if (token !== depictToken) return;   // a newer payload superseded this one
    const svg = depictMoleculeSVG(RDKit, sdf, DEPICT_SIZE);
    if (!svg) {
      showMessage(depictBox, 'Failed to parse molecule for 2D depiction', T.errorFg);
      return;
    }
    depictBox.innerHTML = svg;
    const svgEl = depictBox.querySelector('svg');
    if (svgEl) {
      // Let the depiction scale to fill the pane rather than sit at 400px.
      svgEl.removeAttribute('width');
      svgEl.removeAttribute('height');
      if (!svgEl.getAttribute('viewBox')) {
        svgEl.setAttribute('viewBox', '0 0 ' + DEPICT_SIZE + ' ' + DEPICT_SIZE);
      }
      svgEl.setAttribute('preserveAspectRatio', 'xMidYMid meet');
      svgEl.style.width = '100%';
      svgEl.style.height = '100%';
      svgEl.style.maxWidth = '100%';
      svgEl.style.maxHeight = '100%';
    }
  }).catch(err => {
    if (token !== depictToken) return;
    showMessage(depictBox, 'RDKit failed to load: ' + err.message, T.errorFg);
  });
}

// ============================================================================
// INIT
// ============================================================================
function init(sdf, name, smiles, charge) {
  titleEl.textContent = name || 'Unnamed molecule';

  let counts = null;
  try {
    counts = parseCounts(sdf);
  } catch (e) {
    showErrorBanner('SDF counts line unreadable: ' + e.message);
  }
  renderInfoBar(name, smiles, charge, counts);

  render2D(sdf);
  render3D(sdf);
}

// Keep the 3D canvas sized to its pane.
const viewerRO = new ResizeObserver(() => {
  if (viewer) { try { viewer.resize(); viewer.render(); } catch (e) {} }
});
viewerRO.observe(viewerBox);

showPlaceholder('Waiting for molecule…');

// ============================================================================
// Module entry points
// ============================================================================
export async function onInputs(inputs) {
  if (!inputs) { showPlaceholder('No molecule provided'); return; }

  let sdf;
  try {
    sdf = await asText(inputs['molecule.sdf']);
  } catch (e) {
    showPlaceholder('Could not read molecule.sdf: ' + e.message);
    return;
  }
  if (!sdf || !sdf.trim()) { showPlaceholder('No molecule provided'); return; }

  const name = (typeof inputs.name === 'string') ? inputs.name : (inputs.name == null ? null : String(inputs.name));
  const smiles = (typeof inputs.smiles === 'string') ? inputs.smiles : (inputs.smiles == null ? null : String(inputs.smiles));

  let charge = inputs.total_charge;
  if (typeof charge === 'string') {
    const parsed = parseInt(charge, 10);
    charge = isNaN(parsed) ? charge : parsed;
  }

  try {
    init(sdf, name, smiles, charge);
  } catch (e) {
    if (typeof logStderr === 'function') logStderr('smc render failed: ' + e.message);
    showErrorBanner('Render error: ' + e.message);
  }
}

export function onResize() {
  if (viewer) { try { viewer.resize(); viewer.render(); } catch (e) {} }
}

export function cleanup() {
  try { viewerRO.disconnect(); } catch (e) {}
  destroyViewer();
}
