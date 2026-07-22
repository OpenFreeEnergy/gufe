// ============================================================================
// Ligand network (left) + atom-mapping 3D viewer (right) — split 50/50.
//
// Clicking an edge on the left immediately drives the right panel: it sends
// the two endpoint molecules and the edge's mapping data to the viewer,
// which renders the selected pair in whatever view mode is active.
// ============================================================================

import * as d3 from "https://cdn.jsdelivr.net/npm/d3@7/+esm";

// 3Dmol.js is fetched on demand rather than eagerly through modules.json, so
// the network graph paints immediately — the 3D engine is only needed once an
// edge is clicked, and the 2D mode never needs it at all.
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
    // Shared / app chrome
    appBg:           '#1a1a2e',
    splitBorder:     '#2a4a7f',
    toolbarBg:       '#16213e',
    toolbarBorder:   '#2a4a7f',
    titleColor:      '#7ecfff',
    textPrimary:     '#e6f3ff',
    textMuted:       '#9bb8d6',
    textMuted2:      '#7a96b8',
    selectBg:        '#0f3460',
    selectBorder:    '#2a4a7f',
    tooltipBg:       '#16213e',
    tooltipBorder:   '#2a4a7f',

    // Network panel
    netNodeFill:     '#16213e',
    netNodeStroke:   '#2a4a7f',
    netNodeLabel:    '#cfe6ff',
    netInitials:     '#7ecfff',
    netBondLine:     '#cfe6ff',
    netEdgeRamp:     ['#3a4a6a', '#7ecfff'],   // dim slate → bright cyan
    netEdgeLabel:    '#cfe6ff',
    netLabelBg:      '#16213e',
    netHaloColor:    '#ff79c6',                // pink halo for selection
    netCanvasBg:     '#1a1a2e',

    // Viewer panel
    viewerBg:        '0x1a1a2e',
    labelBg:         '#16213e',
    labelFg:         '#7ecfff',
    switcherBg:      'rgba(22,33,62,0.9)',
    btnBg:           '#0f3460',
    btnBgHover:      '#1a4a8a',
    btnBgActive:     '#2a6ab5',
    btnFg:           '#7ecfff',
    btnBorder:       '#2a4a7f',
    colorCore:       '0xaaaaaa',
    colorUniqueA:    '0xff4d4d',
    colorUniqueB:    '0x4dff88',
    linesMolA:       '0xff8888',
    linesMolB:       '0x88ffaa',
    linesDash:       '0xffee55',
    overlayMolA:     '0xff6666',
    overlayMolB:     '0x66ff99',
    canvas2DBg:      '#fff',
    rgbCore:         [0.7, 0.7, 0.7],
    rgbUniqueA:      [1.0, 0.3, 0.3],
    rgbUniqueB:      [0.3, 1.0, 0.5],
    loadingFg:       '#888',
    errorFg:         '#c33'
  },
  light: {
    appBg:           '#ffffff',
    splitBorder:     '#e2e8f0',
    toolbarBg:       '#f8fafc',
    toolbarBorder:   '#e2e8f0',
    titleColor:      '#0369a1',
    textPrimary:     '#1e293b',
    textMuted:       '#475569',
    textMuted2:      '#64748b',
    selectBg:        '#ffffff',
    selectBorder:    '#cbd5e1',
    tooltipBg:       '#ffffff',
    tooltipBorder:   '#cbd5e1',

    netNodeFill:     '#ffffff',
    netNodeStroke:   '#ffffff',
    netNodeLabel:    '#334155',
    netInitials:     '#0369a1',
    netBondLine:     '#1e293b',
    netEdgeRamp:     ['#cbd5e1', '#0f766e'],
    netEdgeLabel:    '#334155',
    netLabelBg:      '#ffffff',
    netHaloColor:    '#fbcfe8',
    netCanvasBg:     '#ffffff',

    viewerBg:        '0xffffff',
    labelBg:         '#f0f4fa',
    labelFg:         '#1a4a8a',
    switcherBg:      'rgba(240,244,250,0.95)',
    btnBg:           '#e1e8f2',
    btnBgHover:      '#c8d4e8',
    btnBgActive:     '#7ab0e5',
    btnFg:           '#1a4a8a',
    btnBorder:       '#a8bcd6',
    colorCore:       '0x888888',
    colorUniqueA:    '0xd62828',
    colorUniqueB:    '0x2a9d4a',
    linesMolA:       '0xd62828',
    linesMolB:       '0x2a9d4a',
    linesDash:       '0xd9a300',
    overlayMolA:     '0xd62828',
    overlayMolB:     '0x2a9d4a',
    canvas2DBg:      '#fff',
    rgbCore:         [0.55, 0.55, 0.55],
    rgbUniqueA:      [0.84, 0.16, 0.16],
    rgbUniqueB:      [0.16, 0.62, 0.29],
    loadingFg:       '#888',
    errorFg:         '#c33'
  }
};

var T = DARK_MODE ? THEMES.dark : THEMES.light;

// ============================================================================
// NETWORK-SPECIFIC CONFIG (geometry / forces / sizes — palette comes from T)
// ============================================================================
const CONFIG = {
  node: {
    radius:           38,
    strokeWidth:      1.5,
    labelFontSize:    11,
    labelFontWeight:  '600',
    labelMaxChars:    14,
    initialsSize:     18,
    depictionPadding: 4,
  },
  edge: {
    minWidth:        1.5,
    maxWidth:        6.5,
    opacity:         0.9,
    labelFontSize:   10,
    labelBgPadding:  3,
    labelBgOpacity:  0.92,
    arrowWidthPx:    8,
    arrowHeightPx:   8,
    arrowRefX:       46,
    selectionHaloPadding: 4,
    selectionHaloOpacity: 0.95,
  },
  force: {
    linkBaseDistance:    18,
    linkScoreBonus:      10,
    linkStrength:        0.5,
    chargeStrength:     -2500,
    chargeDistanceMin:   20,
    chargeDistanceMax:   5000,
    centerStrength:      0.08,
    collisionPadding:    12,
    collisionIterations: 4,
    driftX:              0.04,
    driftY:              0.04,
    tickMultiplier:      2,
  },
  radial:   { ringStep: 0.18, ringOffset: 40 },
  circular: { radiusFraction: 0.36 },
};

// Element map (atomic number → symbol)
const ELEMENT_MAP = {
  1: 'H', 6: 'C', 7: 'N', 8: 'O', 9: 'F',
  15: 'P', 16: 'S', 17: 'Cl', 35: 'Br', 53: 'I'
};

// ============================================================================
// SHARED HELPERS — molecule parsing, MOL block, RDKit
// ============================================================================

const NL = String.fromCharCode(10);
const DOLLAR = String.fromCharCode(36);

function parseNpyCoords(confStr) {
  const len = confStr.length;
  const buf = new ArrayBuffer(len);
  const vi = new Uint8Array(buf);
  for (let i = 0; i < len; i++) vi[i] = confStr.charCodeAt(i);
  const headerLen = new DataView(buf).getUint16(8, true);
  const headerStr = String.fromCharCode.apply(null, new Uint8Array(buf, 10, headerLen));
  const shapeMatch = headerStr.match(/'shape':[ ]*[(]([0-9]+),[ ]*([0-9]+)[)]/);
  const nAtoms = parseInt(shapeMatch[1]);
  const dataOffset = 10 + headerLen;
  const dv = new DataView(buf, dataOffset);
  const coords = [];
  for (let i = 0; i < nAtoms; i++) {
    coords.push([
      dv.getFloat64(i * 24, true),
      dv.getFloat64(i * 24 + 8, true),
      dv.getFloat64(i * 24 + 16, true)
    ]);
  }
  return coords;
}

function alignFlat(coords) {
  const n = coords.length;
  let cx = 0, cy = 0, cz = 0;
  for (let i = 0; i < n; i++) { cx += coords[i][0]; cy += coords[i][1]; cz += coords[i][2]; }
  cx /= n; cy /= n; cz /= n;
  const pts = coords.map(c => [c[0]-cx, c[1]-cy, c[2]-cz]);
  const cov = [[0,0,0],[0,0,0],[0,0,0]];
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < 3; j++) {
      for (let k = j; k < 3; k++) {
        cov[j][k] += pts[i][j] * pts[i][k];
      }
    }
  }
  cov[1][0] = cov[0][1]; cov[2][0] = cov[0][2]; cov[2][1] = cov[1][2];
  let a = cov.map(r => r.slice());
  let v = [[1,0,0],[0,1,0],[0,0,1]];
  for (let iter = 0; iter < 100; iter++) {
    let p = 0, q = 1, mx = Math.abs(a[0][1]);
    if (Math.abs(a[0][2]) > mx) { p = 0; q = 2; mx = Math.abs(a[0][2]); }
    if (Math.abs(a[1][2]) > mx) { p = 1; q = 2; }
    if (Math.abs(a[p][q]) < 1e-12) break;
    const theta = 0.5 * Math.atan2(2 * a[p][q], a[q][q] - a[p][p]);
    const cs = Math.cos(theta), sn = Math.sin(theta);
    const na = a.map(r => r.slice());
    na[p][p] = cs*cs*a[p][p] - 2*sn*cs*a[p][q] + sn*sn*a[q][q];
    na[q][q] = sn*sn*a[p][p] + 2*sn*cs*a[p][q] + cs*cs*a[q][q];
    na[p][q] = 0; na[q][p] = 0;
    for (let i = 0; i < 3; i++) {
      if (i !== p && i !== q) {
        na[i][p] = cs*a[i][p] - sn*a[i][q]; na[p][i] = na[i][p];
        na[i][q] = sn*a[i][p] + cs*a[i][q]; na[q][i] = na[i][q];
      }
    }
    a = na;
    const nv = v.map(r => r.slice());
    for (let i = 0; i < 3; i++) {
      nv[i][p] = cs*v[i][p] - sn*v[i][q];
      nv[i][q] = sn*v[i][p] + cs*v[i][q];
    }
    v = nv;
  }
  const eigs = [{val: a[0][0], col: 0}, {val: a[1][1], col: 1}, {val: a[2][2], col: 2}];
  eigs.sort((x, y) => y.val - x.val);
  const rot = [];
  for (let i = 0; i < 3; i++) {
    rot.push([v[i][eigs[0].col], v[i][eigs[1].col], v[i][eigs[2].col]]);
  }
  const result = [];
  for (let i = 0; i < n; i++) {
    const px = pts[i][0], py = pts[i][1], pz = pts[i][2];
    result.push([
      rot[0][0]*px + rot[1][0]*py + rot[2][0]*pz,
      rot[0][1]*px + rot[1][1]*py + rot[2][1]*pz,
      rot[0][2]*px + rot[1][2]*py + rot[2][2]*pz
    ]);
  }
  return result;
}

function parseMolecule(nodeEl) {
  const dataEls = nodeEl.getElementsByTagName('data');
  let dataEl;
  for (let i = 0; i < dataEls.length; i++) {
    if (dataEls[i].getAttribute('key') === 'd0') { dataEl = dataEls[i]; break; }
  }
  if (!dataEl) throw new Error('no d0 data element');
  const json = JSON.parse(dataEl.textContent);
  const symbols = json.atoms.map(a => ELEMENT_MAP[a[0]] || 'X');
  const bonds = json.bonds.map(b => [b[0], b[1], b[2]]);
  const rawCoords = parseNpyCoords(json.conformer[0]);
  const coords = alignFlat(rawCoords);
  const name = (json.molprops && json.molprops['ofe-name']) || nodeEl.getAttribute('id');
  return { name, symbols, bonds, coords };
}

// MOL block (V2000) without the SDF terminator. The viewer wraps this with
// $$$$ for 3Dmol; the network passes it straight to RDKit.
function buildMolBlock(mol) {
  const nAtoms = mol.symbols.length;
  const nBonds = mol.bonds.length;
  const lines = [];
  lines.push(mol.name || '');
  lines.push('  Generated');
  lines.push('');
  lines.push(String(nAtoms).padStart(3) + String(nBonds).padStart(3) + '  0  0  0  0  0  0  0  0999 V2000');
  for (let i = 0; i < nAtoms; i++) {
    const c = mol.coords[i];
    lines.push(
      c[0].toFixed(4).padStart(10) +
      c[1].toFixed(4).padStart(10) +
      c[2].toFixed(4).padStart(10) +
      ' ' + mol.symbols[i].padEnd(3) +
      ' 0  0  0  0  0  0  0  0  0  0  0  0'
    );
  }
  for (let j = 0; j < nBonds; j++) {
    const b = mol.bonds[j];
    const sdfType = b[2] === 12 ? 4 : b[2];
    lines.push(String(b[0] + 1).padStart(3) + String(b[1] + 1).padStart(3) + String(sdfType).padStart(3) + '  0  0  0  0');
  }
  lines.push('M  END');
  return lines.join(NL);
}

// SDF = MOL block + the "$$$$" record terminator that 3Dmol expects.
function buildSDF(mol) {
  return buildMolBlock(mol) + NL + DOLLAR + DOLLAR + DOLLAR + DOLLAR;
}

// RDKit loader (singleton)
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

// 2D depiction SVG (no highlights — used by network nodes)
function depictMoleculeSVG(RDKit, mol, size) {
  let rdmol = null;
  try {
    const molblock = buildMolBlock(mol);
    rdmol = RDKit.get_mol(molblock, JSON.stringify({ removeHs: true }));
    if (!rdmol) return null;
    try { rdmol.set_new_coords(true); } catch (e) {}
    const svg = rdmol.get_svg(size, size);
    return svg || null;
  } catch (e) {
    console.warn('[network] depictMoleculeSVG threw for', mol.name, '-', e.message);
    return null;
  } finally {
    if (rdmol) { try { rdmol.delete(); } catch (e) {} }
  }
}

// ============================================================================
// MAPPING DECODE (used by viewer panel for 3D-Map / Pairs / 2D modes)
// ============================================================================

// Accepts many shapes and tries to find a usable atom-mapping inside.
function parseMappingValue(val) {
  if (val == null) return null;
  if (Array.isArray(val)) {
    const m = new Map();
    for (let i = 0; i < val.length; i++) {
      const pair = val[i];
      if (Array.isArray(pair) && pair.length >= 2 &&
          typeof pair[0] === 'number' && typeof pair[1] === 'number') {
        m.set(pair[0], pair[1]);
      }
    }
    return m.size > 0 ? m : null;
  }
  if (typeof val === 'object') {
    const keys = Object.keys(val);
    const looksLikeIntDict = keys.length > 0 && keys.every(k =>
      /^-?\d+$/.test(k) && typeof val[k] === 'number'
    );
    if (looksLikeIntDict) {
      const dm = new Map();
      keys.forEach(k => dm.set(parseInt(k, 10), val[k]));
      return dm.size > 0 ? dm : null;
    }
    for (let ki = 0; ki < keys.length; ki++) {
      let inner = val[keys[ki]];
      if (typeof inner === 'string') {
        try { inner = JSON.parse(inner); } catch (e) { continue; }
      }
      const got = parseMappingValue(inner);
      if (got) return got;
    }
  }
  return null;
}

// ============================================================================
// GraphML parser
// Produces both:
//   - graphData.nodes / edges / molecules  (network panel)
//   - graphData.edgeMappings: per-edge { aToB, bToA } maps  (viewer panel)
//   - graphData.edgeRawData:  per-edge raw decoded data object  (passed to viewer)
// Plus a name||name → Map fallback table (used by the viewer's resolveMapping).
// ============================================================================
function parseGraphML(xmlStr) {
  const parser = new DOMParser();
  const doc = parser.parseFromString(xmlStr, 'text/xml');
  const parseErr = doc.getElementsByTagName('parsererror');
  if (parseErr.length > 0) {
    throw new Error(parseErr[0].textContent);
  }

  // Key registry
  const keyDefs = {};
  const keyEls = doc.getElementsByTagName('key');
  for (let i = 0; i < keyEls.length; i++) {
    const id = keyEls[i].getAttribute('id');
    const name = keyEls[i].getAttribute('attr.name') || id;
    keyDefs[id] = name;
  }

  // Nodes
  const molecules = {};   // id -> mol
  const nodeName  = {};   // id -> ofe-name
  const nodes     = [];
  const nodeEls   = doc.getElementsByTagName('node');
  for (let i = 0; i < nodeEls.length; i++) {
    const id = nodeEls[i].getAttribute('id');
    try {
      const mol = parseMolecule(nodeEls[i]);
      molecules[id] = mol;
      nodeName[id]  = mol.name;
      nodes.push({ id, name: mol.name });
    } catch (e) {
      if (typeof logStderr === 'function') logStderr('Failed to parse node ' + id + ': ' + e.message);
    }
  }

  // Edges
  const edges          = [];
  const edgeMappings   = [];   // parallel array, per edge: { aToB, bToA } | null
  const edgeRawData    = [];   // parallel array, per edge: raw decoded object | null
  const nameMap        = {};   // "nameA||nameB" -> Map<aIdx,bIdx>  (viewer fallback)
  const edgeEls = doc.getElementsByTagName('edge');

  for (let j = 0; j < edgeEls.length; j++) {
    const src = edgeEls[j].getAttribute('source');
    const tgt = edgeEls[j].getAttribute('target');
    if (!molecules[src] || !molecules[tgt]) continue;

    let score = 0.5;
    let mapping = null;
    let rawData = null;

    const dEls = edgeEls[j].getElementsByTagName('data');
    for (let k = 0; k < dEls.length; k++) {
      const keyId   = dEls[k].getAttribute('key');
      const keyName = (keyDefs[keyId] || keyId || '').toLowerCase();
      const text    = dEls[k].textContent.trim();
      if (!text) continue;

      // Score
      if (keyName.includes('score') || keyName.includes('lomap') || keyName === 'weight') {
        const v = parseFloat(text);
        if (!isNaN(v)) score = v;
      }

      // Try to interpret as JSON for mapping extraction
      try {
        const parsed = JSON.parse(text);
        if (rawData == null) rawData = parsed;
        if (mapping == null) {
          const m = parseMappingValue(parsed);
          if (m && m.size > 0) mapping = m;
        }
      } catch (e) { /* not JSON, that's fine */ }
    }

    edges.push({ source: src, target: tgt, score });

    if (mapping) {
      const bToA = new Map();
      mapping.forEach((b, a) => bToA.set(b, a));
      edgeMappings.push({ aToB: mapping, bToA });

      const nameA = nodeName[src], nameB = nodeName[tgt];
      if (nameA && nameB) {
        nameMap[nameA + '||' + nameB] = mapping;
        nameMap[nameB + '||' + nameA] = bToA;
      }
    } else {
      edgeMappings.push(null);
    }
    edgeRawData.push(rawData);
  }

  return { nodes, edges, molecules, edgeMappings, edgeRawData, nameMap };
}

// ============================================================================
// DOM SCAFFOLD: split-pane layout
// ============================================================================
root.innerHTML = '';
root.style.background = T.appBg;
root.style.fontFamily = "'Inter',system-ui,sans-serif";

const splitWrap = document.createElement('div');
splitWrap.style.cssText = 'width:100%;height:100%;display:flex;flex-direction:row;overflow:hidden;';
root.appendChild(splitWrap);

const leftPane = document.createElement('div');
leftPane.style.cssText = 'flex:1 1 50%;width:50%;height:100%;display:flex;flex-direction:column;min-width:0;min-height:0;background:' + T.appBg + ';';
splitWrap.appendChild(leftPane);

const divider = document.createElement('div');
divider.style.cssText = 'width:1px;background:' + T.splitBorder + ';flex-shrink:0;';
splitWrap.appendChild(divider);

const rightPane = document.createElement('div');
rightPane.style.cssText = 'flex:1 1 50%;width:50%;height:100%;display:flex;flex-direction:column;min-width:0;min-height:0;background:' + T.appBg + ';position:relative;';
splitWrap.appendChild(rightPane);

// ============================================================================
// VIEWER PANEL (right) — atom-mapping 3D viewer
// All state and DOM live inside `rightPane`.
// ============================================================================

const viewerWrapper = document.createElement('div');
viewerWrapper.style.cssText = 'position:relative;width:100%;height:100%;display:flex;flex-direction:column;';
rightPane.appendChild(viewerWrapper);

const viewerArea = document.createElement('div');
viewerArea.style.cssText = 'flex:1;display:flex;flex-direction:column;min-height:0;';
viewerWrapper.appendChild(viewerArea);

// View-mode switcher
const MODES = [
  { id: 'plain',   label: '3D',      title: 'Plain 3D view' },
  { id: 'colored', label: '3D-Map',  title: 'Color-coded by mapping' },
  { id: 'lines',   label: 'Pairs',   title: 'Dashed lines between mapped atoms' },
  { id: 'overlay', label: 'Overlay', title: 'Both molecules superimposed' },
  { id: '2d',      label: '2D',      title: '2D depictions with highlights' }
];

let currentMode = 'plain';
let currentMols = null;
let explicitMapping = null;
let viewerNameMap = {};   // populated from graphData.nameMap when graph loads

const switcher = document.createElement('div');
switcher.style.cssText = 'position:absolute;bottom:10px;right:10px;display:flex;gap:4px;background:' + T.switcherBg + ';padding:4px;border-radius:6px;z-index:10;box-shadow:0 2px 8px rgba(0,0,0,0.4);';
MODES.forEach(m => {
  const btn = document.createElement('button');
  btn.textContent = m.label;
  btn.title = m.title;
  btn.dataset.mode = m.id;
  btn.style.cssText = 'background:' + T.btnBg + ';color:' + T.btnFg + ';border:1px solid ' + T.btnBorder + ';padding:4px 8px;font-size:11px;font-weight:bold;border-radius:3px;cursor:pointer;font-family:inherit;';
  btn.onmouseover = () => { btn.style.background = T.btnBgHover; };
  btn.onmouseout  = () => { btn.style.background = (currentMode === m.id) ? T.btnBgActive : T.btnBg; };
  btn.onclick = () => setMode(m.id);
  switcher.appendChild(btn);
});
viewerWrapper.appendChild(switcher);

function updateSwitcherHighlight() {
  switcher.querySelectorAll('button').forEach(b => {
    b.style.background = (b.dataset.mode === currentMode) ? T.btnBgActive : T.btnBg;
  });
}

function setMode(mode) {
  currentMode = mode;
  updateSwitcherHighlight();
  if (currentMols) renderMode();
}

// Viewer-box management
let boxes = [];
let syncing = false;
let syncRunning = false;
const SYNC_VIEWS = true;

function makeViewerBox(labelText) {
  const box = document.createElement('div');
  box.style.cssText = 'flex:1;display:flex;flex-direction:column;position:relative;min-height:0;';
  const label = document.createElement('div');
  label.style.cssText = 'padding:4px 10px;font-size:13px;font-weight:bold;color:' + T.labelFg + ';background:' + T.labelBg + ';';
  label.textContent = labelText;
  box.appendChild(label);
  const container = document.createElement('div');
  container.style.cssText = 'flex:1;position:relative;';
  box.appendChild(container);
  viewerArea.appendChild(box);
  return { label, container, box, viewer: null };
}

function clearAllBoxes() {
  boxes.forEach(b => {
    if (b.viewer) { try { b.viewer.clear(); } catch (e) {} }
  });
  viewerArea.innerHTML = '';
  boxes = [];
}

function startSyncLoop() {
  if (syncRunning || !SYNC_VIEWS || boxes.length < 2) return;
  syncRunning = true;
  const lastViews = boxes.map(() => '');
  function loop() {
    if (!syncing && boxes.length >= 2) {
      for (let i = 0; i < boxes.length; i++) {
        if (!boxes[i].viewer) continue;
        const cur = JSON.stringify(boxes[i].viewer.getView());
        if (cur !== lastViews[i]) {
          syncing = true;
          for (let j = 0; j < boxes.length; j++) {
            if (j !== i && boxes[j].viewer) {
              boxes[j].viewer.setView(boxes[i].viewer.getView());
              boxes[j].viewer.render();
            }
            lastViews[j] = cur;
          }
          syncing = false;
          break;
        }
      }
    }
    requestAnimationFrame(loop);
  }
  loop();
}

function resolveMapping() {
  if (explicitMapping && explicitMapping.size > 0) return explicitMapping;
  if (currentMols && currentMols.length >= 2) {
    return viewerNameMap[currentMols[0].name + '||' + currentMols[1].name] || null;
  }
  return null;
}

// ─── MODE: plain ───
function renderPlain(mols) {
  clearAllBoxes();
  for (let i = 0; i < mols.length; i++) {
    const b = makeViewerBox(mols[i].name);
    b.viewer = ThreeDmol.createViewer(b.container, { backgroundColor: T.viewerBg });
    b.viewer.addModel(buildSDF(mols[i]), 'sdf');
    b.viewer.setStyle({}, {
      stick:  { radius: 0.15, colorscheme: 'Jmol' },
      sphere: { scale: 0.25, colorscheme: 'Jmol' }
    });
    b.viewer.zoomTo();
    b.viewer.render();
    boxes.push(b);
  }
  startSyncLoop();
}

// ─── MODE: colored ───
function renderColored(mols) {
  clearAllBoxes();
  if (mols.length < 2) { renderPlain(mols); return; }

  let mapping = resolveMapping();
  if (!mapping) {
    console.warn('[viewer] No mapping available for', mols[0].name, '->', mols[1].name);
    renderPlain(mols);
    return;
  }

  const nA0 = mols[0].coords.length;
  const nB0 = mols[1].coords.length;
  let maxK = -1, maxV = -1;
  mapping.forEach((v, k) => {
    if (k > maxK) maxK = k;
    if (v > maxV) maxV = v;
  });
  if ((maxK >= nA0 || maxV >= nB0) && maxK < nB0 && maxV < nA0) {
    const flipped0 = new Map();
    mapping.forEach((v, k) => flipped0.set(v, k));
    mapping = flipped0;
  }

  const coreA = new Set(mapping.keys());
  const coreB = new Set(mapping.values());

  for (let i = 0; i < mols.length; i++) {
    const b = makeViewerBox(mols[i].name);
    b.viewer = ThreeDmol.createViewer(b.container, { backgroundColor: T.viewerBg });
    b.viewer.addModel(buildSDF(mols[i]), 'sdf');

    b.viewer.setStyle({}, {
      stick:  { radius: 0.15, color: T.colorCore },
      sphere: { scale: 0.25, color: T.colorCore }
    });

    const coreSet = (i === 0) ? coreA : coreB;
    const uniqueColor = (i === 0) ? T.colorUniqueA : T.colorUniqueB;
    const nAtoms = mols[i].symbols.length;

    for (let ai = 0; ai < nAtoms; ai++) {
      if (!coreSet.has(ai)) {
        b.viewer.addStyle({ serial: ai + 1 }, {
          stick:  { radius: 0.18, color: uniqueColor },
          sphere: { scale: 0.32, color: uniqueColor }
        });
      }
    }

    b.viewer.zoomTo();
    b.viewer.render();
    boxes.push(b);
  }
  startSyncLoop();
}

// ─── Kabsch alignment helpers (used by Pairs mode) ───
function mat3Mul(a, b) {
  const r = new Array(9);
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      r[i*3+j] = a[i*3]*b[j] + a[i*3+1]*b[3+j] + a[i*3+2]*b[6+j];
    }
  }
  return r;
}
function mat3Transpose(m) {
  return [m[0], m[3], m[6], m[1], m[4], m[7], m[2], m[5], m[8]];
}
function mat3Det(m) {
  return m[0]*(m[4]*m[8] - m[5]*m[7])
       - m[1]*(m[3]*m[8] - m[5]*m[6])
       + m[2]*(m[3]*m[7] - m[4]*m[6]);
}
function jacobiSym3(M) {
  const a = M.slice();
  const v = [1,0,0, 0,1,0, 0,0,1];
  for (let sweep = 0; sweep < 50; sweep++) {
    const off = Math.abs(a[1]) + Math.abs(a[2]) + Math.abs(a[5]);
    if (off < 1e-12) break;
    const pairs = [[0,1],[0,2],[1,2]];
    for (let pi = 0; pi < 3; pi++) {
      const p = pairs[pi][0], q = pairs[pi][1];
      const apq = a[p*3+q];
      if (Math.abs(apq) < 1e-14) continue;
      const app = a[p*3+p], aqq = a[q*3+q];
      const theta = (aqq - app) / (2*apq);
      let t;
      if (Math.abs(theta) > 1e10) {
        t = 1/(2*theta);
      } else {
        const sign = theta >= 0 ? 1 : -1;
        t = sign / (Math.abs(theta) + Math.sqrt(theta*theta + 1));
      }
      const c = 1/Math.sqrt(1 + t*t);
      const s = t*c;
      a[p*3+p] = app - t*apq;
      a[q*3+q] = aqq + t*apq;
      a[p*3+q] = 0; a[q*3+p] = 0;
      for (let r = 0; r < 3; r++) {
        if (r !== p && r !== q) {
          const arp = a[r*3+p], arq = a[r*3+q];
          a[r*3+p] = c*arp - s*arq; a[p*3+r] = a[r*3+p];
          a[r*3+q] = s*arp + c*arq; a[q*3+r] = a[r*3+q];
        }
      }
      for (let k = 0; k < 3; k++) {
        const vkp = v[k*3+p], vkq = v[k*3+q];
        v[k*3+p] = c*vkp - s*vkq;
        v[k*3+q] = s*vkp + c*vkq;
      }
    }
  }
  return { values: [a[0], a[4], a[8]], vectors: v };
}

function kabsch(P, Q) {
  const n = Math.min(P.length, Q.length);
  if (n < 1) return null;
  const cP = [0,0,0], cQ = [0,0,0];
  for (let i = 0; i < n; i++) {
    cP[0] += P[i][0]; cP[1] += P[i][1]; cP[2] += P[i][2];
    cQ[0] += Q[i][0]; cQ[1] += Q[i][1]; cQ[2] += Q[i][2];
  }
  cP[0]/=n; cP[1]/=n; cP[2]/=n;
  cQ[0]/=n; cQ[1]/=n; cQ[2]/=n;
  if (n < 3) {
    return { R: [1,0,0, 0,1,0, 0,0,1], t: [cP[0]-cQ[0], cP[1]-cQ[1], cP[2]-cQ[2]] };
  }
  const H = [0,0,0, 0,0,0, 0,0,0];
  for (let k = 0; k < n; k++) {
    const px = P[k][0]-cP[0], py = P[k][1]-cP[1], pz = P[k][2]-cP[2];
    const qx = Q[k][0]-cQ[0], qy = Q[k][1]-cQ[1], qz = Q[k][2]-cQ[2];
    H[0] += px*qx; H[1] += px*qy; H[2] += px*qz;
    H[3] += py*qx; H[4] += py*qy; H[5] += py*qz;
    H[6] += pz*qx; H[7] += pz*qy; H[8] += pz*qz;
  }
  const Ht = mat3Transpose(H);
  const HtH = mat3Mul(Ht, H);
  const HHt = mat3Mul(H, Ht);
  let eV = jacobiSym3(HtH);
  let eU = jacobiSym3(HHt);
  function sortEig(e) {
    const idx = [0,1,2].sort((a, b) => e.values[b] - e.values[a]);
    const sortedVec = new Array(9);
    for (let c = 0; c < 3; c++) {
      const src = idx[c];
      sortedVec[c]   = e.vectors[src];
      sortedVec[3+c] = e.vectors[3+src];
      sortedVec[6+c] = e.vectors[6+src];
    }
    return {
      values: [e.values[idx[0]], e.values[idx[1]], e.values[idx[2]]],
      vectors: sortedVec
    };
  }
  eV = sortEig(eV);
  eU = sortEig(eU);
  const V = eV.vectors;
  const U = eU.vectors;
  for (let col = 0; col < 3; col++) {
    const vx = V[col], vy = V[3+col], vz = V[6+col];
    const hx = H[0]*vx + H[1]*vy + H[2]*vz;
    const hy = H[3]*vx + H[4]*vy + H[5]*vz;
    const hz = H[6]*vx + H[7]*vy + H[8]*vz;
    const ux = U[col], uy = U[3+col], uz = U[6+col];
    const dot = hx*ux + hy*uy + hz*uz;
    if (dot < 0) {
      U[col]   = -ux;
      U[3+col] = -uy;
      U[6+col] = -uz;
    }
  }
  const Vt = mat3Transpose(V);
  let R = mat3Mul(U, Vt);
  if (mat3Det(R) < 0) {
    U[2] = -U[2];
    U[5] = -U[5];
    U[8] = -U[8];
    R = mat3Mul(U, Vt);
  }
  const rcQx = R[0]*cQ[0] + R[1]*cQ[1] + R[2]*cQ[2];
  const rcQy = R[3]*cQ[0] + R[4]*cQ[1] + R[5]*cQ[2];
  const rcQz = R[6]*cQ[0] + R[7]*cQ[1] + R[8]*cQ[2];
  return { R, t: [cP[0]-rcQx, cP[1]-rcQy, cP[2]-rcQz] };
}

function applyRT(coord, R, t) {
  const x = coord[0], y = coord[1], z = coord[2];
  return [
    R[0]*x + R[1]*y + R[2]*z + t[0],
    R[3]*x + R[4]*y + R[5]*z + t[1],
    R[6]*x + R[7]*y + R[8]*z + t[2]
  ];
}

// ─── MODE: lines ───
function renderLines(mols) {
  clearAllBoxes();
  if (mols.length < 2) { renderPlain(mols); return; }

  let mapping = resolveMapping();
  if (!mapping) {
    console.warn('[viewer] No mapping available; falling back to plain');
    renderPlain(mols);
    return;
  }

  const mol0 = mols[0], mol1 = mols[1];
  const nA = mol0.coords.length;
  const nB = mol1.coords.length;

  let maxKey = -1, maxVal = -1;
  mapping.forEach((v, k) => {
    if (k > maxKey) maxKey = k;
    if (v > maxVal) maxVal = v;
  });
  if (maxKey >= nA || maxVal >= nB) {
    if (maxKey < nB && maxVal < nA) {
      const flipped = new Map();
      mapping.forEach((v, k) => flipped.set(v, k));
      mapping = flipped;
    }
  }

  const b = makeViewerBox(mol0.name + ' \u2194 ' + mol1.name + '  (' + mapping.size + ' mapped pairs)');
  b.viewer = ThreeDmol.createViewer(b.container, { backgroundColor: T.viewerBg });

  const P = [], Q = [];
  mapping.forEach((bIdx, aIdx) => {
    const pa = mol0.coords[aIdx];
    const pb = mol1.coords[bIdx];
    if (!pa || !pb) return;
    P.push(pa); Q.push(pb);
  });

  const rt = kabsch(P, Q);
  let alignedCoords;
  if (rt) {
    alignedCoords = mol1.coords.map(c => applyRT(c, rt.R, rt.t));
  } else {
    alignedCoords = mol1.coords.map(c => c.slice());
  }

  function extents(coords, n) {
    const mn = [Infinity, Infinity, Infinity];
    const mx = [-Infinity, -Infinity, -Infinity];
    for (let i = 0; i < n; i++) {
      for (let k = 0; k < 3; k++) {
        const v = coords[i][k];
        if (v < mn[k]) mn[k] = v;
        if (v > mx[k]) mx[k] = v;
      }
    }
    return { min: mn, max: mx, span: [mx[0]-mn[0], mx[1]-mn[1], mx[2]-mn[2]] };
  }
  const ext0 = extents(mol0.coords,    nA);
  const ext1 = extents(alignedCoords,  nB);

  let liftAxis = 0;
  if (ext0.span[1] < ext0.span[liftAxis]) liftAxis = 1;
  if (ext0.span[2] < ext0.span[liftAxis]) liftAxis = 2;
  const maxSpan = Math.max(ext0.span[0], ext0.span[1], ext0.span[2]);

  const GAP = 2.5;
  const clearance = (ext0.max[liftAxis] - ext1.min[liftAxis]) + GAP;
  const minLift = 0.6 * maxSpan + GAP;
  const lift = Math.max(clearance, minLift);

  const mol1Shifted = {
    name: mol1.name,
    symbols: mol1.symbols,
    bonds: mol1.bonds,
    coords: alignedCoords.map(c => {
      const out = [c[0], c[1], c[2]];
      out[liftAxis] += lift;
      return out;
    })
  };

  b.viewer.addModel(buildSDF(mol0), 'sdf');
  b.viewer.addModel(buildSDF(mol1Shifted), 'sdf');

  b.viewer.setStyle({ model: 0 }, {
    stick:  { radius: 0.15, color: T.linesMolA },
    sphere: { scale: 0.22, color: T.linesMolA }
  });
  b.viewer.setStyle({ model: 1 }, {
    stick:  { radius: 0.15, color: T.linesMolB },
    sphere: { scale: 0.22, color: T.linesMolB }
  });

  mapping.forEach((bIdx, aIdx) => {
    const pa = mol0.coords[aIdx];
    const pb = mol1Shifted.coords[bIdx];
    if (!pa || !pb) return;
    b.viewer.addCylinder({
      start:  { x: pa[0], y: pa[1], z: pa[2] },
      end:    { x: pb[0], y: pb[1], z: pb[2] },
      radius: 0.04,
      dashed: true,
      fromCap: 'round',
      toCap:   'round',
      color:  T.linesDash
    });
  });

  b.viewer.zoomTo();
  if (liftAxis === 2) {
    b.viewer.rotate(90, 'x');
  } else if (liftAxis === 0) {
    b.viewer.rotate(-90, 'z');
  }
  b.viewer.render();
  boxes.push(b);
}

// ─── MODE: overlay ───
function renderOverlay(mols) {
  clearAllBoxes();
  if (mols.length < 2) { renderPlain(mols); return; }

  const b = makeViewerBox(mols[0].name + ' + ' + mols[1].name + '  (overlay)');
  b.viewer = ThreeDmol.createViewer(b.container, { backgroundColor: T.viewerBg });
  b.viewer.addModel(buildSDF(mols[0]), 'sdf');
  b.viewer.addModel(buildSDF(mols[1]), 'sdf');
  b.viewer.setStyle({ model: 0 }, {
    stick:  { radius: 0.15, color: T.overlayMolA, opacity: 0.7 },
    sphere: { scale: 0.22, color: T.overlayMolA, opacity: 0.7 }
  });
  b.viewer.setStyle({ model: 1 }, {
    stick:  { radius: 0.15, color: T.overlayMolB, opacity: 0.7 },
    sphere: { scale: 0.22, color: T.overlayMolB, opacity: 0.7 }
  });
  b.viewer.zoomTo();
  b.viewer.render();
  boxes.push(b);
}

// ─── MODE: 2d ───
function molblockOnly(mol) {
  return buildMolBlock(mol);
}

function render2D(mols) {
  clearAllBoxes();
  const pendingBoxes = mols.map(m => {
    const box = document.createElement('div');
    box.style.cssText = 'flex:1;display:flex;flex-direction:column;position:relative;min-height:0;';
    const label = document.createElement('div');
    label.style.cssText = 'padding:4px 10px;font-size:13px;font-weight:bold;color:' + T.labelFg + ';background:' + T.labelBg + ';';
    label.textContent = m.name;
    box.appendChild(label);
    const container = document.createElement('div');
    container.style.cssText = 'flex:1;display:flex;align-items:center;justify-content:center;background:' + T.canvas2DBg + ';overflow:hidden;';
    container.innerHTML = '<div style="color:' + T.loadingFg + ';font-size:12px;">Loading 2D depiction\u2026</div>';
    box.appendChild(container);
    viewerArea.appendChild(box);
    return { container, mol: m };
  });

  loadRDKit().then(RDKit => {
    try {
      const mapping = (mols.length >= 2) ? resolveMapping() : null;
      const coreA = mapping ? new Set(mapping.keys())   : null;
      const coreB = mapping ? new Set(mapping.values()) : null;

      pendingBoxes.forEach((pb, idx) => {
        const rdmol = RDKit.get_mol(molblockOnly(pb.mol));
        if (!rdmol) {
          pb.container.innerHTML = '<div style="color:' + T.errorFg + ';">Failed to parse</div>';
          return;
        }
        try { rdmol.set_new_coords(true); } catch (e) {}

        let drawOpts = {};
        if (mapping && idx < 2) {
          const coreSet = (idx === 0) ? coreA : coreB;
          const nAtoms = pb.mol.symbols.length;
          const highlightAtoms = [];
          const atomColors = {};
          for (let ai = 0; ai < nAtoms; ai++) {
            highlightAtoms.push(ai);
            atomColors[ai] = coreSet.has(ai) ? T.rgbCore : (idx === 0 ? T.rgbUniqueA : T.rgbUniqueB);
          }
          drawOpts = {
            atoms: highlightAtoms,
            highlightAtomColors: atomColors,
            highlightAtomRadii: highlightAtoms.reduce((o, a) => { o[a] = 0.3; return o; }, {})
          };
        }

        let svg;
        try {
          if (mapping && idx < 2) {
            svg = rdmol.get_svg_with_highlights(JSON.stringify(drawOpts));
          } else {
            svg = rdmol.get_svg(300, 300);
          }
        } catch (e) {
          svg = rdmol.get_svg(300, 300);
        }
        rdmol.delete();

        pb.container.innerHTML = svg;
        const svgEl = pb.container.querySelector('svg');
        if (svgEl) {
          svgEl.style.maxWidth = '100%';
          svgEl.style.maxHeight = '100%';
          svgEl.setAttribute('preserveAspectRatio', 'xMidYMid meet');
        }
      });
    } catch (err) {
      pendingBoxes.forEach(pb => {
        pb.container.innerHTML = '<div style="color:' + T.errorFg + ';font-size:12px;padding:10px;">RDKit failed: ' + err.message + '</div>';
      });
    }
  }).catch(err => {
    pendingBoxes.forEach(pb => {
      pb.container.innerHTML = '<div style="color:' + T.errorFg + ';font-size:12px;padding:10px;">RDKit failed to load: ' + err.message + '</div>';
    });
  });
}

// Bumped per render so a slow 3Dmol load can't paint over a newer edge/mode.
let modeToken = 0;

function renderMode() {
  if (!currentMols) return;
  const mols = currentMols;
  const token = ++modeToken;

  // The 2D depiction is pure RDKit — it must never wait on the 3D engine.
  if (currentMode === '2d') { render2D(mols); return; }

  clearAllBoxes();
  const loading = document.createElement('div');
  loading.style.cssText = 'flex:1;display:flex;align-items:center;justify-content:center;' +
    'text-align:center;padding:24px;font-size:13px;color:' + T.textMuted2 + ';';
  loading.textContent = 'Loading 3D viewer…';
  viewerArea.appendChild(loading);

  load3Dmol().then(() => {
    if (token !== modeToken) return;   // a newer edge/mode superseded this one
    clearAllBoxes();                   // drop the loading message
    switch (currentMode) {
      case 'colored': renderColored(mols); break;
      case 'lines':   renderLines(mols);   break;
      case 'overlay': renderOverlay(mols); break;
      default:        renderPlain(mols);
    }
  }).catch(e => {
    if (token !== modeToken) return;
    clearAllBoxes();
    loading.style.color = T.errorFg;
    loading.textContent = '3D viewer unavailable: ' + e.message;
    viewerArea.appendChild(loading);
  });
}

// ResizeObserver for the right pane
const viewerRO = new ResizeObserver(() => {
  for (let i = 0; i < boxes.length; i++) {
    if (boxes[i].viewer) { boxes[i].viewer.resize(); boxes[i].viewer.render(); }
  }
});
viewerRO.observe(rightPane);

// Public entrypoint to drive the viewer panel from the network panel.
// `edgePayload` mirrors what the original viewer expected as `inputs.edge`.
function showInViewer(mols, edgePayload) {
  if (!mols || mols.length === 0) return;
  currentMols = mols;
  explicitMapping = null;

  if (edgePayload) {
    let edgeData = edgePayload.data;
    if (typeof edgeData === 'string') {
      try { edgeData = JSON.parse(edgeData); } catch (e) { edgeData = null; }
    }
    if (edgeData != null) {
      explicitMapping = parseMappingValue(edgeData);
    }
    if (!explicitMapping) {
      explicitMapping = parseMappingValue(edgePayload);
    }
  }

  updateSwitcherHighlight();
  renderMode();
}

// ============================================================================
// NETWORK PANEL (left)
// ============================================================================

const netWrap = document.createElement('div');
netWrap.style.cssText = 'width:100%;height:100%;display:flex;flex-direction:column;background:' + T.netCanvasBg + ';color:' + T.textPrimary + ";font-family:'Inter',system-ui,sans-serif;overflow:hidden;";
leftPane.appendChild(netWrap);

const toolbar = document.createElement('div');
toolbar.style.cssText = 'display:flex;align-items:center;gap:12px;padding:10px 16px;background:' + T.toolbarBg + ';border-top:1px solid ' + T.toolbarBorder + ';flex-shrink:0;flex-wrap:wrap;';

// const titleEl = document.createElement('span');
// titleEl.textContent = 'Ligand Network';
// titleEl.style.cssText = 'font-weight:700;font-size:15px;color:' + T.titleColor + ';letter-spacing:.03em;margin-right:8px;';
// toolbar.appendChild(titleEl);

const layoutLabel = document.createElement('label');
layoutLabel.textContent = 'Layout:';
layoutLabel.style.cssText = 'font-size:12px;color:' + T.textMuted + ';margin-left:auto;';
toolbar.appendChild(layoutLabel);

const layoutSelect = document.createElement('select');
layoutSelect.style.cssText = 'background:' + T.selectBg + ';color:' + T.textPrimary + ';border:1px solid ' + T.selectBorder + ';border-radius:6px;padding:4px 8px;font-size:12px;cursor:pointer;';
['Force-directed', 'Circular', 'Radial Tree'].forEach(l => {
  const o = document.createElement('option');
  o.value = l; o.textContent = l;
  o.style.background = T.selectBg;
  o.style.color = T.textPrimary;
  layoutSelect.appendChild(o);
});
toolbar.appendChild(layoutSelect);

const legendWrap = document.createElement('div');
legendWrap.style.cssText = 'display:flex;align-items:center;gap:6px;font-size:11px;color:' + T.textMuted + ';';
legendWrap.innerHTML =
  '<span>LOMAP score:</span>' +
  '<span style="display:flex;align-items:center;gap:4px;">' +
    '<span style="width:32px;height:4px;background:linear-gradient(to right,' + T.netEdgeRamp.join(',') + ');border-radius:2px;display:inline-block;"></span>' +
    '<span>0 → 1</span>' +
  '</span>';
toolbar.appendChild(legendWrap);

const svgWrap = document.createElement('div');
svgWrap.style.cssText = 'flex:1;position:relative;overflow:hidden;min-height:0;background:' + T.netCanvasBg + ';';
netWrap.appendChild(svgWrap);

netWrap.appendChild(toolbar);   // toolbar at the bottom

const tooltip = document.createElement('div');
tooltip.style.cssText = 'position:absolute;pointer-events:none;opacity:0;background:' + T.tooltipBg + ';border:1px solid ' + T.tooltipBorder + ';border-radius:8px;padding:10px 14px;font-size:12px;color:' + T.textPrimary + ';max-width:220px;box-shadow:0 4px 16px rgba(15,23,42,0.12);transition:opacity .15s;z-index:10;';
svgWrap.appendChild(tooltip);

const scoreColor = d3.scaleSequential()
  .domain([0, 1])
  .interpolator(d3.interpolateRgbBasis(T.netEdgeRamp));

let graphData = null;
let simulation = null;
let svgEl = null;
let activeEdgeIdx = -1;

const rdkitReady = loadRDKit().catch(err => {
  if (typeof logStderr === 'function') logStderr('RDKit load failed: ' + err.message);
  return null;
});

// Drive the right (viewer) panel from the currently-active edge.
function selectEdge(edgeIdx) {
  if (!graphData || edgeIdx < 0 || edgeIdx >= graphData.edges.length) return;
  const edge    = graphData.edges[edgeIdx];
  const sourceMol = graphData.molecules[edge.source];
  const targetMol = graphData.molecules[edge.target];
  if (!sourceMol || !targetMol) return;

  const rawData = graphData.edgeRawData[edgeIdx];
  showInViewer([sourceMol, targetMol], {
    source: edge.source,
    target: edge.target,
    data: rawData,
    sourceNode: sourceMol,
    targetNode: targetMol
  });
}

// ─── Render network ───
function renderNetwork(data, layout) {
  if (svgEl) svgEl.remove();
  if (simulation) { simulation.stop(); simulation = null; }

  const W  = svgWrap.getBoundingClientRect().width  || 800;
  const H  = svgWrap.getBoundingClientRect().height || 600;
  const CX = W / 2, CY = H / 2;
  const NR = CONFIG.node.radius;

  svgEl = d3.select(svgWrap).append('svg')
    .attr('width', W).attr('height', H).style('display', 'block');

  const gRoot = svgEl.append('g');
  svgEl.call(d3.zoom().scaleExtent([0.15, 5]).on('zoom', e => gRoot.attr('transform', e.transform)));

  const defs = svgEl.append('defs');

  function addArrowMarker(id, fill) {
    defs.append('marker')
      .attr('id', id)
      .attr('viewBox', '0 -5 10 10')
      .attr('refX', CONFIG.edge.arrowRefX)
      .attr('refY', 0)
      .attr('markerUnits', 'userSpaceOnUse')
      .attr('markerWidth',  CONFIG.edge.arrowWidthPx)
      .attr('markerHeight', CONFIG.edge.arrowHeightPx)
      .attr('orient', 'auto')
      .append('path').attr('d', 'M0,-5L10,0L0,5').attr('fill', fill);
  }

  const colorToMarkerId = new Map();
  function markerIdForColor(color) {
    if (colorToMarkerId.has(color)) return colorToMarkerId.get(color);
    const id = 'arrow-' + color.replace(/[^a-zA-Z0-9]/g, '');
    colorToMarkerId.set(color, id);
    addArrowMarker(id, color);
    return id;
  }

  defs.append('clipPath')
    .attr('id', 'node-clip')
    .append('circle')
    .attr('r', NR - CONFIG.node.depictionPadding);

  const nodes    = data.nodes.map(n => ({ ...n }));
  const nodeById = Object.fromEntries(nodes.map(n => [n.id, n]));
  const edges    = data.edges.map((e, i) => ({
    ...e,
    _idx: i,
    _color: scoreColor(e.score),
    source: nodeById[e.source] || e.source,
    target: nodeById[e.target] || e.target,
  }));

  edges.forEach(e => markerIdForColor(e._color));

  // Layouts
  if (layout === 'Circular') {
    const r = Math.min(W, H) * CONFIG.circular.radiusFraction;
    nodes.forEach((n, i) => {
      const a = (2 * Math.PI * i / nodes.length) - Math.PI / 2;
      n.x = CX + r * Math.cos(a); n.y = CY + r * Math.sin(a);
      n.fx = n.x; n.fy = n.y;
    });
  } else if (layout === 'Radial Tree') {
    const adj = Object.fromEntries(nodes.map(n => [n.id, []]));
    edges.forEach(e => {
      const sid = e.source.id || e.source;
      const tid = e.target.id || e.target;
      adj[sid]?.push(tid); adj[tid]?.push(sid);
    });
    const visited = new Set([nodes[0].id]);
    const levels  = [[nodes[0].id]];
    while (true) {
      const next = [];
      levels[levels.length - 1].forEach(id =>
        (adj[id] || []).forEach(nb => { if (!visited.has(nb)) { visited.add(nb); next.push(nb); } })
      );
      if (!next.length) break;
      levels.push(next);
    }
    const rStep = Math.min(W, H) * CONFIG.radial.ringStep;
    levels.forEach((level, li) => {
      const r = li === 0 ? 0 : li * rStep + CONFIG.radial.ringOffset;
      level.forEach((id, i) => {
        const a  = (2 * Math.PI * i / level.length) - Math.PI / 2;
        const nd = nodeById[id];
        if (nd) { nd.x = CX + r * Math.cos(a); nd.y = CY + r * Math.sin(a); nd.fx = nd.x; nd.fy = nd.y; }
      });
    });
    nodes.filter(n => !visited.has(n.id)).forEach((n, i) => { n.x = 40 + i * 80; n.y = 40; n.fx = n.x; n.fy = n.y; });
  }

  const ticked = () => {
    link
      .attr('x1', d => d.source.x).attr('y1', d => d.source.y)
      .attr('x2', d => d.target.x).attr('y2', d => d.target.y);
    linkHalo
      .attr('x1', d => d.source.x).attr('y1', d => d.source.y)
      .attr('x2', d => d.target.x).attr('y2', d => d.target.y);
    linkHit
      .attr('x1', d => d.source.x).attr('y1', d => d.source.y)
      .attr('x2', d => d.target.x).attr('y2', d => d.target.y);
    linkLabelG.attr('transform', d =>
      `translate(${(d.source.x + d.target.x) / 2}, ${(d.source.y + d.target.y) / 2 - 6})`
    );
    nodeG.attr('transform', d => `translate(${d.x ?? CX},${d.y ?? CY})`);
  };

  function refreshEdgeStyles() {
    linkHalo.attr('opacity', d =>
      d._idx === activeEdgeIdx ? CONFIG.edge.selectionHaloOpacity : 0
    );
  }

  // Hit layer
  const linkHit = gRoot.append('g').selectAll('line').data(edges).join('line')
    .attr('stroke', 'transparent')
    .attr('stroke-width', 14)
    .style('cursor', 'pointer');

  // Halo layer
  const linkHalo = gRoot.append('g').selectAll('line').data(edges).join('line')
    .attr('stroke', T.netHaloColor)
    .attr('stroke-width', d =>
      CONFIG.edge.minWidth + d.score * (CONFIG.edge.maxWidth - CONFIG.edge.minWidth)
      + CONFIG.edge.selectionHaloPadding * 2)
    .attr('stroke-linecap', 'round')
    .attr('opacity', 0)
    .style('pointer-events', 'none');

  // Visible edges
  const link = gRoot.append('g').selectAll('line').data(edges).join('line')
    .attr('stroke', d => d._color)
    .attr('stroke-width', d =>
      CONFIG.edge.minWidth + d.score * (CONFIG.edge.maxWidth - CONFIG.edge.minWidth))
    .attr('stroke-opacity', CONFIG.edge.opacity)
    .attr('marker-end', d => `url(#${markerIdForColor(d._color)})`)
    .style('pointer-events', 'none');

  // Edge labels
  const linkLabelG = gRoot.append('g').selectAll('g').data(edges).join('g')
    .attr('pointer-events', 'none');

  linkLabelG.append('rect')
    .attr('fill', T.netLabelBg)
    .attr('opacity', CONFIG.edge.labelBgOpacity)
    .attr('rx', 3).attr('ry', 3);

  const linkLabel = linkLabelG.append('text')
    .attr('text-anchor', 'middle')
    .attr('dominant-baseline', 'middle')
    .attr('font-size', CONFIG.edge.labelFontSize + 'px')
    .attr('font-weight', '600')
    .attr('fill', T.netEdgeLabel)
    .text(d => d.score.toFixed(2));

  const pad = CONFIG.edge.labelBgPadding;
  linkLabel.each(function() {
    const bb = this.getBBox();
    const bg = this.previousSibling;
    if (bg && bg.setAttribute) {
      bg.setAttribute('x',      bb.x - pad);
      bg.setAttribute('y',      bb.y - pad);
      bg.setAttribute('width',  bb.width  + 2 * pad);
      bg.setAttribute('height', bb.height + 2 * pad);
    }
  });

  linkHit
    .on('mouseenter', (ev, d) => {
      tooltip.style.opacity = '1';
      const sName = d.source.name || d.source.id || d.source;
      const tName = d.target.name || d.target.id || d.target;
      tooltip.innerHTML =
        '<div style="font-weight:700;margin-bottom:4px;color:' + T.titleColor + ';">Transformation</div>' +
        '<div style="color:' + T.textPrimary + ';">' + sName + ' → ' + tName + '</div>' +
        '<div style="margin-top:6px;color:' + T.titleColor + ';">LOMAP score: <b>' + d.score.toFixed(3) + '</b></div>' +
        '<div style="margin-top:4px;font-size:10px;color:' + T.textMuted2 + ';">Click to select</div>';
    })
    .on('mousemove', ev => { tooltip.style.left = (ev.offsetX + 14) + 'px'; tooltip.style.top = (ev.offsetY - 10) + 'px'; })
    .on('mouseleave', () => { tooltip.style.opacity = '0'; })
    .on('click', (ev, d) => {
      ev.stopPropagation();
      tooltip.style.opacity = '0';
      activeEdgeIdx = d._idx;
      refreshEdgeStyles();
      selectEdge(activeEdgeIdx);
    });

  // Nodes
  const nodeG = gRoot.append('g').selectAll('g').data(nodes).join('g')
    .style('cursor', 'grab')
    .call(d3.drag()
      .on('start', (ev, d) => { d.fx = d.x; d.fy = d.y; })
      .on('drag',  (ev, d) => { d.fx = ev.x; d.fy = ev.y; simulation?.tick(); ticked(); })
      .on('end',   (ev, d) => { d.fx = ev.x; d.fy = ev.y; })
    )
    .on('mouseenter', (ev, d) => {
      tooltip.style.opacity = '1';
      tooltip.innerHTML =
        '<div style="font-weight:700;color:' + T.titleColor + ';">' + d.name + '</div>' +
        '<div style="color:' + T.textMuted2 + ';font-size:10px;margin-top:3px;">id: ' + d.id + '</div>';
    })
    .on('mousemove', ev => { tooltip.style.left = (ev.offsetX + 14) + 'px'; tooltip.style.top = (ev.offsetY - 10) + 'px'; })
    .on('mouseleave', () => { tooltip.style.opacity = '0'; })
    .on('click', ev => ev.stopPropagation());

  nodeG.append('circle')
    .attr('r',            NR)
    .attr('fill',         T.netNodeFill)
    .attr('stroke',       T.netNodeStroke)
    .attr('stroke-width', CONFIG.node.strokeWidth);

  const depictionGroups = nodeG.append('g').attr('class', 'node-depiction');

  const placeholderText = nodeG.append('text')
    .attr('class', 'node-placeholder')
    .attr('text-anchor', 'middle').attr('dominant-baseline', 'middle')
    .attr('font-size', CONFIG.node.initialsSize + 'px')
    .attr('font-weight', '700')
    .attr('fill', T.netInitials)
    .attr('pointer-events', 'none')
    .text(d => (d.name || '?').slice(0, 2).toUpperCase());

  rdkitReady.then(RDKit => {
    if (!RDKit) return;

    const RDKIT_SIZE = 200;
    const innerR = NR - CONFIG.node.depictionPadding;
    const scale = (innerR * 2) / RDKIT_SIZE;
    const parser = new DOMParser();
    let injected = 0;

    depictionGroups.each(function(d) {
      const mol = data.molecules[d.id];
      if (!mol) return;
      const svgStr = depictMoleculeSVG(RDKit, mol, RDKIT_SIZE);
      if (!svgStr) return;
      const doc = parser.parseFromString(svgStr, 'image/svg+xml');
      const rdSvg = doc.documentElement;
      if (!rdSvg || rdSvg.nodeName.toLowerCase() === 'parsererror') return;

      const g = this;
      const half = RDKIT_SIZE / 2;
      g.setAttribute('transform', `translate(${-scale * half}, ${-scale * half}) scale(${scale})`);

      const kids = Array.from(rdSvg.childNodes);
      let appended = 0;
      for (const k of kids) {
        if (k.nodeType !== 1) continue;
        const tag = k.nodeName.toLowerCase();
        if (tag === 'defs' || tag === 'metadata' || tag === 'title') continue;
        if (tag === 'rect') {
          const fill = (k.getAttribute('fill') || '').toLowerCase();
          if (fill === '#ffffff' || fill === 'white' || fill === 'rgb(255,255,255)') continue;
        }
        const imported = document.importNode(k, true);
        g.appendChild(imported);
        appended++;
      }
      if (appended > 0) injected++;
    });

    if (injected > 0) placeholderText.style('display', 'none');
  }).catch(err => {
    if (typeof logStderr === 'function') logStderr('depiction render failed: ' + err.message);
  });

  nodeG.append('text')
    .attr('text-anchor', 'middle')
    .attr('y',           NR + 14)
    .attr('font-size',   CONFIG.node.labelFontSize + 'px')
    .attr('fill',        T.netNodeLabel)
    .attr('font-weight', CONFIG.node.labelFontWeight)
    .attr('pointer-events', 'none')
    .text(d => d.name.length > CONFIG.node.labelMaxChars
      ? d.name.slice(0, CONFIG.node.labelMaxChars - 1) + '…'
      : d.name);

  if (layout === 'Force-directed') {
    const fc = CONFIG.force;
    simulation = d3.forceSimulation(nodes)
      .force('link',      d3.forceLink(edges).id(d => d.id)
                            .distance(d => fc.linkBaseDistance + (1 - d.score) * fc.linkScoreBonus)
                            .strength(fc.linkStrength))
      .force('charge',    d3.forceManyBody()
                            .strength(fc.chargeStrength)
                            .distanceMin(fc.chargeDistanceMin)
                            .distanceMax(fc.chargeDistanceMax))
      .force('center',    d3.forceCenter(CX, CY).strength(fc.centerStrength))
      .force('collision', d3.forceCollide(NR + fc.collisionPadding).strength(1).iterations(fc.collisionIterations))
      .force('x',         d3.forceX(CX).strength(fc.driftX))
      .force('y',         d3.forceY(CY).strength(fc.driftY))
      .stop();

    const n = Math.ceil(Math.log(simulation.alphaMin()) / Math.log(1 - simulation.alphaDecay()));
    for (let i = 0; i < n * fc.tickMultiplier; i++) simulation.tick();
  }

  ticked();
  refreshEdgeStyles();

  // Auto-select first edge so the viewer panel is populated immediately.
  if (activeEdgeIdx < 0 && edges.length > 0) {
    activeEdgeIdx = 0;
    refreshEdgeStyles();
  }
  selectEdge(activeEdgeIdx);
}

// ─── Init ───
function init(xmlStr) {
  try {
    graphData = parseGraphML(xmlStr);
    viewerNameMap = graphData.nameMap || {};
  } catch (e) {
    if (typeof logStderr === 'function') logStderr('GraphML parse failed: ' + e.message);
    const warn = document.createElement('div');
    warn.textContent = '⚠ GraphML parse error: ' + e.message;
    warn.style.cssText = 'position:absolute;top:60px;left:50%;transform:translateX(-50%);background:#fee2e2;color:#991b1b;border:1px solid #fecaca;padding:6px 14px;border-radius:6px;font-size:12px;z-index:20;';
    svgWrap.appendChild(warn);
    return;
  }
  activeEdgeIdx = -1;
  renderNetwork(graphData, layoutSelect.value);
}

layoutSelect.addEventListener('change', () => {
  if (graphData) renderNetwork(graphData, layoutSelect.value);
});

// ============================================================================
// Module entry points
// ============================================================================
export async function onInputs(inputs) {
  let xml = inputs && inputs['network.graphml'];
  if (!xml) return;
  if (xml instanceof Blob) xml = await xml.text();
  else if (xml instanceof ArrayBuffer) xml = new TextDecoder().decode(xml);
  else if (ArrayBuffer.isView(xml))    xml = new TextDecoder().decode(xml);
  if (typeof xml !== 'string') xml = String(xml);
  init(xml);
}

export function onResize() {
  if (graphData) renderNetwork(graphData, layoutSelect.value);
  for (let i = 0; i < boxes.length; i++) {
    if (boxes[i].viewer) { boxes[i].viewer.resize(); boxes[i].viewer.render(); }
  }
}

export function cleanup() {
  if (simulation) simulation.stop();
}