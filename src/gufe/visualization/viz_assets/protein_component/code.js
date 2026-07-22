// ============================================================================
// Protein component viewer — toolbar + full-pane 3Dmol.js viewer.
//
// Renders a single gufe.ProteinComponent (or a subclass such as
// SolvatedPDBComponent / ProteinMembraneComponent, so the PDB may also carry
// waters, ions and lipids).
//
// Payload:  { "protein.pdb": <PDB string>, "name": <str> }
//
// Everything is driven through 3Dmol *selections* — never per-atom DOM work —
// so tens of thousands of atoms stay responsive.
// ============================================================================

// 3Dmol.js is fetched on demand rather than eagerly through modules.json, so
// the toolbar and stats readout paint before the 3D engine arrives.
var ThreeDmol = null;
var threeDmolPromise = null;

function load3Dmol() {
  if (threeDmolPromise) return threeDmolPromise;
  threeDmolPromise = new Promise(function(resolve, reject) {
    if (window.$3Dmol) { ThreeDmol = window.$3Dmol; resolve(ThreeDmol); return; }
    const script = document.createElement('script');
    script.src = 'https://3dmol.org/build/3Dmol-min.js';
    script.onload = function() {
      ThreeDmol = window.$3Dmol;
      if (ThreeDmol) resolve(ThreeDmol);
      else reject(new Error('3Dmol.js loaded but $3Dmol is undefined'));
    };
    script.onerror = function() { reject(new Error('Failed to load 3Dmol.js')); };
    document.head.appendChild(script);
  });
  return threeDmolPromise;
}

// ============================================================================
// THEME — flip DARK_MODE to switch the entire app
// ============================================================================
var DARK_MODE = false;

var THEMES = {
  dark: {
    appBg:         '#1a1a2e',
    toolbarBg:     '#16213e',
    toolbarBorder: '#2a4a7f',
    titleColor:    '#7ecfff',
    textPrimary:   '#e6f3ff',
    textMuted:     '#9bb8d6',
    textMuted2:    '#7a96b8',
    selectBg:      '#0f3460',
    selectBorder:  '#2a4a7f',
    viewerBg:      '0x1a1a2e',
    btnBg:         '#0f3460',
    btnBgHover:    '#1a4a8a',
    btnBgActive:   '#2a6ab5',
    btnFg:         '#7ecfff',
    btnBorder:     '#2a4a7f',
    surfaceColor:  '0xbfd8f0',
    loadingFg:     '#888',
    errorFg:       '#c33'
  },
  light: {
    appBg:         '#ffffff',
    toolbarBg:     '#f8fafc',
    toolbarBorder: '#e2e8f0',
    titleColor:    '#0369a1',
    textPrimary:   '#1e293b',
    textMuted:     '#475569',
    textMuted2:    '#64748b',
    selectBg:      '#ffffff',
    selectBorder:  '#cbd5e1',
    viewerBg:      '0xffffff',
    btnBg:         '#e1e8f2',
    btnBgHover:    '#c8d4e8',
    btnBgActive:   '#7ab0e5',
    btnFg:         '#1a4a8a',
    btnBorder:     '#a8bcd6',
    surfaceColor:  '0x9db8d6',
    loadingFg:     '#888',
    errorFg:       '#c33'
  }
};

var T = DARK_MODE ? THEMES.dark : THEMES.light;

// ============================================================================
// CONFIG
// ============================================================================
const CONFIG = {
  stick:   { radius: 0.15 },
  sphere:  { scale: 0.30 },
  hetero:  { stickRadius: 0.20, sphereScale: 0.28 },
  water:   { stickRadius: 0.05, sphereScale: 0.18 },
  surface: { opacity: 0.85 },
  // Above this many atoms a surface is slow enough to be worth warning about.
  surfaceAtomWarn: 40000
};

// Residue names treated as water.
const WATER_RESN = ['HOH', 'WAT', 'SOL', 'TIP3'];

const REPS = [
  { id: 'cartoon', label: 'Cartoon', title: 'Ribbon / cartoon backbone' },
  { id: 'surface', label: 'Surface', title: 'Molecular (VDW) surface' },
  { id: 'stick',   label: 'Stick',   title: 'All-atom sticks' },
  { id: 'sphere',  label: 'Sphere',  title: 'Space-filling spheres' }
];

const COLOR_SCHEMES = [
  { id: 'chain',    label: 'Chain' },
  { id: 'spectrum', label: 'Spectrum' },
  { id: 'ss',       label: 'Secondary structure' },
  { id: 'element',  label: 'Element' }
];

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

/**
 * Single pass over the PDB text for the toolbar stats readout.
 *
 * Only the first MODEL is counted (NMR ensembles would otherwise multiply every
 * number). Residues are keyed by chain + sequence number + insertion code +
 * name, which is what makes them unique in a PDB.
 */
function parsePdbStats(pdbText) {
  const chains = new Set();
  const residues = new Set();
  let atoms = 0;
  let hetatms = 0;
  let waters = 0;
  let resiMin = Infinity;
  let resiMax = -Infinity;

  const lines = pdbText.split(/\r?\n/);
  for (let i = 0; i < lines.length; i++) {
    const line = lines[i];
    const rec = line.slice(0, 6);
    if (rec === 'ENDMDL') break;                     // first model only
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

    const resi = parseInt(resSeq, 10);
    if (!isNaN(resi)) {
      if (resi < resiMin) resiMin = resi;
      if (resi > resiMax) resiMax = resi;
    }
  }

  return {
    chains: chains.size,
    residues: residues.size,
    atoms,
    hetatms,
    waters,
    // Non-water HETATM records — what the "hetero/ligands" toggle governs.
    heteroNonWater: hetatms - waters,
    resiMin: resiMin === Infinity ? 0 : resiMin,
    resiMax: resiMax === -Infinity ? 0 : resiMax
  };
}

function fmt(n) {
  return n.toLocaleString('en-US');
}

// ============================================================================
// DOM SCAFFOLD
// ============================================================================
var mount = (typeof root !== 'undefined' && root) ? root : document.body;

mount.innerHTML = '';
mount.style.background = T.appBg;
mount.style.fontFamily = "'Inter',system-ui,sans-serif";

const appWrap = document.createElement('div');
appWrap.style.cssText = 'width:100%;height:100%;display:flex;flex-direction:column;overflow:hidden;background:' + T.appBg + ';';
mount.appendChild(appWrap);

// ─── Toolbar ───
const toolbar = document.createElement('div');
toolbar.style.cssText =
  'display:flex;align-items:center;gap:14px;padding:8px 14px;background:' + T.toolbarBg +
  ';border-bottom:1px solid ' + T.toolbarBorder + ';flex-shrink:0;flex-wrap:wrap;font-size:12px;color:' + T.textPrimary + ';';
appWrap.appendChild(toolbar);

const nameEl = document.createElement('span');
nameEl.style.cssText = 'font-weight:700;font-size:14px;color:' + T.titleColor + ';letter-spacing:.02em;';
nameEl.textContent = 'Protein';
toolbar.appendChild(nameEl);

const BTN_CSS =
  'background:' + T.btnBg + ';color:' + T.btnFg + ';border:1px solid ' + T.btnBorder +
  ';padding:4px 9px;font-size:11px;font-weight:bold;border-radius:3px;cursor:pointer;font-family:inherit;';

function makeGroupLabel(text) {
  const el = document.createElement('span');
  el.textContent = text;
  el.style.cssText = 'font-size:11px;color:' + T.textMuted + ';';
  return el;
}

// Representation switcher
toolbar.appendChild(makeGroupLabel('Style:'));

const repGroup = document.createElement('div');
repGroup.style.cssText = 'display:flex;gap:4px;';
toolbar.appendChild(repGroup);

let currentRep = 'cartoon';

REPS.forEach(r => {
  const btn = document.createElement('button');
  btn.textContent = r.label;
  btn.title = r.title;
  btn.dataset.rep = r.id;
  btn.style.cssText = BTN_CSS;
  btn.onmouseover = () => { btn.style.background = T.btnBgHover; };
  btn.onmouseout  = () => { btn.style.background = (currentRep === r.id) ? T.btnBgActive : T.btnBg; };
  btn.onclick = () => { currentRep = r.id; refreshToggleStyles(); applyStyles(); };
  repGroup.appendChild(btn);
});

// Colour-scheme switcher
toolbar.appendChild(makeGroupLabel('Color:'));

const colorSelect = document.createElement('select');
colorSelect.style.cssText =
  'background:' + T.selectBg + ';color:' + T.textPrimary + ';border:1px solid ' + T.selectBorder +
  ';border-radius:6px;padding:4px 8px;font-size:12px;cursor:pointer;font-family:inherit;';
COLOR_SCHEMES.forEach(c => {
  const o = document.createElement('option');
  o.value = c.id;
  o.textContent = c.label;
  o.style.background = T.selectBg;
  o.style.color = T.textPrimary;
  colorSelect.appendChild(o);
});
colorSelect.value = 'chain';
colorSelect.addEventListener('change', () => applyStyles());
toolbar.appendChild(colorSelect);

// Toggles
const toggleGroup = document.createElement('div');
toggleGroup.style.cssText = 'display:flex;gap:4px;';
toolbar.appendChild(toggleGroup);

const toggleState = { waters: false, hetero: true, spin: false };

function makeToggle(key, label, title, onChange) {
  const btn = document.createElement('button');
  btn.textContent = label;
  btn.title = title;
  btn.dataset.toggle = key;
  btn.style.cssText = BTN_CSS;
  btn.onmouseover = () => { btn.style.background = T.btnBgHover; };
  btn.onmouseout  = () => { btn.style.background = toggleState[key] ? T.btnBgActive : T.btnBg; };
  btn.onclick = () => {
    toggleState[key] = !toggleState[key];
    refreshToggleStyles();
    onChange();
  };
  toggleGroup.appendChild(btn);
  return btn;
}

makeToggle('waters', 'Waters', 'Show water molecules', () => applyStyles());
makeToggle('hetero', 'Hetero', 'Show hetero atoms / ligands / ions / lipids', () => applyStyles());
makeToggle('spin',   'Spin',   'Rotate the view continuously', () => {
  if (viewer) { viewer.spin(toggleState.spin); }
});

function refreshToggleStyles() {
  repGroup.querySelectorAll('button').forEach(b => {
    b.style.background = (b.dataset.rep === currentRep) ? T.btnBgActive : T.btnBg;
  });
  toggleGroup.querySelectorAll('button').forEach(b => {
    b.style.background = toggleState[b.dataset.toggle] ? T.btnBgActive : T.btnBg;
  });
}

// Stats readout (pushed to the right)
const statsEl = document.createElement('span');
statsEl.style.cssText = 'margin-left:auto;font-size:11px;color:' + T.textMuted2 + ';white-space:nowrap;';
toolbar.appendChild(statsEl);

// ─── Viewer pane ───
const viewerWrap = document.createElement('div');
viewerWrap.style.cssText = 'flex:1;position:relative;min-height:0;min-width:0;';
appWrap.appendChild(viewerWrap);

const viewerContainer = document.createElement('div');
viewerContainer.style.cssText = 'position:absolute;inset:0;';
viewerWrap.appendChild(viewerContainer);

// Status / placeholder / error overlay
const statusEl = document.createElement('div');
statusEl.style.cssText =
  'position:absolute;top:12px;left:50%;transform:translateX(-50%);padding:6px 14px;border-radius:6px;' +
  'font-size:12px;z-index:20;display:none;pointer-events:none;';
viewerWrap.appendChild(statusEl);

function showStatus(msg, kind) {
  statusEl.textContent = msg;
  statusEl.style.display = 'block';
  if (kind === 'error') {
    statusEl.style.background = '#fee2e2';
    statusEl.style.color = '#991b1b';
    statusEl.style.border = '1px solid #fecaca';
  } else {
    statusEl.style.background = T.toolbarBg;
    statusEl.style.color = T.textMuted;
    statusEl.style.border = '1px solid ' + T.toolbarBorder;
  }
}

function hideStatus() {
  statusEl.style.display = 'none';
}

refreshToggleStyles();

// ============================================================================
// VIEWER STATE
// ============================================================================
let viewer = null;
let stats = null;
let hasModel = false;
let lastPdb = null;
let initToken = 0;

function ensureViewer() {
  if (viewer) return viewer;
  viewer = ThreeDmol.createViewer(viewerContainer, { backgroundColor: T.viewerBg });
  return viewer;
}

// Selection specs. Waters are matched by residue name so that PDBs which write
// TIP3 as ATOM (rather than HETATM) are still caught.
const SEL_POLYMER = { hetflag: false };
const SEL_HETERO  = { hetflag: true };
const SEL_WATER   = { resn: WATER_RESN };

/** Colour arguments for the currently-selected scheme, valid for any rep. */
function colorArgs() {
  switch (colorSelect.value) {
    case 'chain':
      return { colorscheme: 'chain' };
    case 'spectrum':
      // A residue-index gradient works for every representation, unlike
      // 3Dmol's cartoon-only `color:'spectrum'`.
      if (stats && stats.resiMax > stats.resiMin) {
        return { colorscheme: { prop: 'resi', gradient: 'roygb', min: stats.resiMin, max: stats.resiMax } };
      }
      return { colorscheme: 'Jmol' };
    case 'ss':
      return { colorscheme: 'ssJmol' };
    default:
      return { colorscheme: 'Jmol' };
  }
}

/** Style object for the polymer, for the currently-selected representation. */
function polymerStyle() {
  const c = colorArgs();
  switch (currentRep) {
    case 'stick':
      return { stick: Object.assign({ radius: CONFIG.stick.radius }, c) };
    case 'sphere':
      return { sphere: Object.assign({ scale: CONFIG.sphere.scale }, c) };
    case 'surface':
      // The surface itself is added separately; leave the atoms bare so the
      // shell is not cluttered from the inside.
      return {};
    default:
      return { cartoon: Object.assign({}, c) };
  }
}

/**
 * Re-apply every style from scratch.
 *
 * The four setStyle calls are ordered deliberately: each one *replaces* the
 * style of the atoms it matches, so waters (matched last) win over the blanket
 * hetero rule they would otherwise fall under.
 */
function applyStyles() {
  if (!viewer || !hasModel) return;

  try {
    viewer.removeAllSurfaces();
  } catch (e) { /* no surfaces yet */ }

  viewer.setStyle({}, {});                                    // 1. hide everything
  viewer.setStyle(SEL_POLYMER, polymerStyle());               // 2. polymer

  viewer.setStyle(SEL_HETERO, toggleState.hetero ? {          // 3. all hetero
    stick:  { radius: CONFIG.hetero.stickRadius, colorscheme: 'Jmol' },
    sphere: { scale: CONFIG.hetero.sphereScale,  colorscheme: 'Jmol' }
  } : {});

  viewer.setStyle(SEL_WATER, toggleState.waters ? {           // 4. waters override
    stick:  { radius: CONFIG.water.stickRadius, colorscheme: 'Jmol' },
    sphere: { scale: CONFIG.water.sphereScale,  colorscheme: 'Jmol' }
  } : {});

  if (currentRep === 'surface') {
    addSurfaceAsync();
  } else {
    hideStatus();
    viewer.render();
  }
}

/**
 * Surface generation is the one genuinely expensive operation here, so it gets
 * a status message and is deferred a frame to let that message paint.
 */
function addSurfaceAsync() {
  const heavy = stats && stats.atoms > CONFIG.surfaceAtomWarn;
  showStatus(heavy ? 'Computing surface (large structure, this may take a while)…' : 'Computing surface…');
  viewer.render();

  setTimeout(() => {
    if (!viewer || !hasModel || currentRep !== 'surface') return;
    try {
      const surfStyle = Object.assign({ opacity: CONFIG.surface.opacity }, colorArgs());
      const result = viewer.addSurface(ThreeDmol.SurfaceType.VDW, surfStyle, SEL_POLYMER);
      // 3Dmol v2 returns a promise; v1 returns a surface id.
      Promise.resolve(result).then(() => {
        hideStatus();
        if (viewer) viewer.render();
      }).catch(err => {
        showStatus('Surface failed: ' + err.message, 'error');
      });
    } catch (err) {
      showStatus('Surface failed: ' + err.message, 'error');
    }
  }, 30);
}

function updateStatsReadout() {
  if (!stats) { statsEl.textContent = ''; return; }
  statsEl.textContent =
    fmt(stats.chains) + ' chains · ' +
    fmt(stats.residues) + ' residues · ' +
    fmt(stats.atoms) + ' atoms · ' +
    fmt(stats.hetatms) + ' HETATM' +
    (stats.waters > 0 ? ' (' + fmt(stats.waters) + ' water)' : '');
}

// ============================================================================
// INIT
// ============================================================================
function init(pdbText, name) {
  nameEl.textContent = name || 'Protein';

  if (!pdbText || !pdbText.trim()) {
    hasModel = false;
    stats = null;
    updateStatsReadout();
    showStatus('No protein data — waiting for a PDB payload.');
    return;
  }

  let parseError = null;
  try {
    stats = parsePdbStats(pdbText);
  } catch (e) {
    stats = null;
    parseError = '⚠ PDB parse error: ' + e.message;
    showStatus(parseError, 'error');
  }
  updateStatsReadout();

  // Bumped per payload so a slow 3Dmol load can't overwrite a newer structure.
  const token = ++initToken;
  if (!parseError) showStatus('Loading 3D viewer…');

  load3Dmol().then(() => {
    if (token !== initToken) return;   // a newer payload superseded this one
    const v = ensureViewer();
    v.clear();
    try { v.removeAllSurfaces(); } catch (e) { /* none yet */ }
    v.addModel(pdbText, 'pdb');
    hasModel = true;

    applyStyles();
    v.zoomTo();
    v.spin(toggleState.spin);
    v.render();
    // applyStyles() already cleared the "Loading…" status (or replaced it with
    // the surface-computing message), so nothing to do here.
  }).catch(e => {
    if (token !== initToken) return;
    hasModel = false;
    showStatus('⚠ Failed to render structure: ' + e.message, 'error');
    if (typeof logStderr === 'function') logStderr('protein render failed: ' + e.message);
  });
}

// Keep the canvas in step with the pane (the frame can be resized by its host
// without onResize firing).
const viewerRO = new ResizeObserver(() => {
  if (viewer) { viewer.resize(); viewer.render(); }
});
viewerRO.observe(viewerWrap);

showStatus('Waiting for protein data…');

// ============================================================================
// Module entry points
// ============================================================================
export async function onInputs(inputs) {
  if (!inputs) return;
  try {
    const pdb = await asText(inputs['protein.pdb']);
    const name = await asText(inputs['name']);
    if (pdb === lastPdb && hasModel) {
      // Same structure re-delivered — just refresh the label.
      nameEl.textContent = name || 'Protein';
      return;
    }
    lastPdb = pdb;
    init(pdb, name);
  } catch (e) {
    showStatus('⚠ ' + e.message, 'error');
    if (typeof logStderr === 'function') logStderr('protein onInputs failed: ' + e.message);
  }
}

export function onResize() {
  if (viewer) { viewer.resize(); viewer.render(); }
}

export function cleanup() {
  try { viewerRO.disconnect(); } catch (e) {}
  if (viewer) {
    try { viewer.spin(false); } catch (e) {}
    try { viewer.clear(); } catch (e) {}
  }
}
