// ============================================================================
// gufe — one frame for every visualizable gufe object.
//
// `onInputs` dispatches on which payload keys arrived and mounts the matching
// view (see VIEWS at the bottom). Everything above that is shared: one theme,
// one set of input/SDF/PDB helpers, one lazily-loaded RDKit / 3Dmol / d3, one
// atom-mapping viewer used by three different views.
//
// Payload key -> view
//   network.graphml     LigandNetwork      radial graph + mapping viewer
//   molA.sdf, molB.sdf  LigandAtomMapping  the mapping viewer, standalone
//   molecule.sdf        SmallMoleculeComponent
//   protein.pdb         ProteinComponent
//   solvent_component   SolventComponent
//   chemical_system     ChemicalSystem
//   transformation      Transformation / NonTransformation
//   alchemical_network  AlchemicalNetwork
// ============================================================================

// ============================================================================
// THEME — flip DARK_MODE to switch the entire app
// ============================================================================
var DARK_MODE = false;

var THEMES = {
  dark: {
    appBg:         '#1a1a2e',
    panelBg:       '#16213e',
    cardBg:        '#16213e',
    cardBgHover:   '#1a2b4f',
    cardBgActive:  '#0f3460',
    cardBorder:    '#2a4a7f',
    cardBorderActive: '#7ecfff',
    splitBorder:   '#2a4a7f',
    toolbarBg:     '#16213e',
    toolbarBorder: '#2a4a7f',
    tooltipBg:     '#16213e',
    tooltipBorder: '#2a4a7f',
    titleColor:    '#7ecfff',
    textPrimary:   '#e6f3ff',
    textMuted:     '#9bb8d6',
    textMuted2:    '#7a96b8',
    selectBg:      '#0f3460',
    selectBorder:  '#2a4a7f',
    labelBg:       '#16213e',
    labelFg:       '#7ecfff',
    badgeBg:       '#0f3460',
    badgeFg:       '#7ecfff',
    switcherBg:    'rgba(22,33,62,0.9)',
    btnBg:         '#0f3460',
    btnBgHover:    '#1a4a8a',
    btnBgActive:   '#2a6ab5',
    btnFg:         '#7ecfff',
    btnBorder:     '#2a4a7f',

    okBg:          '#14532d',
    okFg:          '#86efac',
    okBorder:      '#166534',
    warnBg:        '#3b1d1d',
    warnFg:        '#ffb4b4',
    warnBorder:    '#7f2a2a',
    loadingFg:     '#888',
    errorFg:       '#ff8080',

    // 3D viewers (3Dmol wants 0x-prefixed colours)
    viewerBg:      '0x1a1a2e',
    canvas2DBg:    '#ffffff',
    colorCore:     '0xaaaaaa',
    colorUniqueA:  '0xff4d4d',
    colorUniqueB:  '0x4dff88',
    linesMolA:     '0xff8888',
    linesMolB:     '0x88ffaa',
    linesDash:     '0xffee55',
    overlayMolA:   '0xff6666',
    overlayMolB:   '0x66ff99',
    rgbCore:       [0.7, 0.7, 0.7],
    rgbUniqueA:    [1.0, 0.3, 0.3],
    rgbUniqueB:    [0.3, 1.0, 0.5],
    chipCore:      '#aaaaaa',
    chipUniqueA:   '#ff4d4d',
    chipUniqueB:   '#4dff88',

    // Diff palette (transformation)
    diffUnchanged: '#64748b',
    diffChanged:   '#d9a300',
    diffAdded:     '#2a9d4a',
    diffRemoved:   '#d62828',

    // Graphs
    netCanvasBg:   '#1a1a2e',
    netNodeFill:   '#16213e',
    netNodeStroke: '#2a4a7f',
    netNodeLabel:  '#cfe6ff',
    netInitials:   '#7ecfff',
    netEdgeRamp:   ['#3a4a6a', '#7ecfff'],   // dim slate -> bright cyan
    netEdgeLine:   '#5f7ea8',
    netEdgeLabel:  '#cfe6ff',
    netLabelBg:    '#16213e',
    netHaloColor:  '#ff79c6',

    // Solvent schematic
    boxFill:       '#0f2a4a',
    boxStroke:     '#2a4a7f'
  },
  light: {
    appBg:         '#ffffff',
    panelBg:       '#f8fafc',
    cardBg:        '#ffffff',
    cardBgHover:   '#f1f5f9',
    cardBgActive:  '#e0f2fe',
    cardBorder:    '#e2e8f0',
    cardBorderActive: '#0369a1',
    splitBorder:   '#e2e8f0',
    toolbarBg:     '#f8fafc',
    toolbarBorder: '#e2e8f0',
    tooltipBg:     '#ffffff',
    tooltipBorder: '#cbd5e1',
    titleColor:    '#0369a1',
    textPrimary:   '#1e293b',
    textMuted:     '#475569',
    textMuted2:    '#64748b',
    selectBg:      '#ffffff',
    selectBorder:  '#cbd5e1',
    labelBg:       '#f0f4fa',
    labelFg:       '#1a4a8a',
    badgeBg:       '#e1e8f2',
    badgeFg:       '#1a4a8a',
    switcherBg:    'rgba(240,244,250,0.95)',
    btnBg:         '#e1e8f2',
    btnBgHover:    '#c8d4e8',
    btnBgActive:   '#7ab0e5',
    btnFg:         '#1a4a8a',
    btnBorder:     '#a8bcd6',

    okBg:          '#dcfce7',
    okFg:          '#166534',
    okBorder:      '#bbf7d0',
    warnBg:        '#fee2e2',
    warnFg:        '#991b1b',
    warnBorder:    '#fecaca',
    loadingFg:     '#888',
    errorFg:       '#c33',

    viewerBg:      '0xffffff',
    canvas2DBg:    '#ffffff',
    colorCore:     '0x888888',
    colorUniqueA:  '0xd62828',
    colorUniqueB:  '0x2a9d4a',
    linesMolA:     '0xd62828',
    linesMolB:     '0x2a9d4a',
    linesDash:     '0xd9a300',
    overlayMolA:   '0xd62828',
    overlayMolB:   '0x2a9d4a',
    rgbCore:       [0.55, 0.55, 0.55],
    rgbUniqueA:    [0.84, 0.16, 0.16],
    rgbUniqueB:    [0.16, 0.62, 0.29],
    chipCore:      '#888888',
    chipUniqueA:   '#d62828',
    chipUniqueB:   '#2a9d4a',

    diffUnchanged: '#94a3b8',
    diffChanged:   '#b45309',
    diffAdded:     '#2a9d4a',
    diffRemoved:   '#d62828',

    netCanvasBg:   '#ffffff',
    netNodeFill:   '#ffffff',
    netNodeStroke: '#ffffff',
    netNodeLabel:  '#334155',
    netInitials:   '#0369a1',
    netEdgeRamp:   ['#cbd5e1', '#0f766e'],
    netEdgeLine:   '#94a3b8',
    netEdgeLabel:  '#334155',
    netLabelBg:    '#ffffff',
    netHaloColor:  '#fbcfe8',

    boxFill:       '#eff6ff',
    boxStroke:     '#bfdbfe'
  }
};

const T = DARK_MODE ? THEMES.dark : THEMES.light;

// ============================================================================
// LAZY ENGINE LOADERS
//
// Nothing heavy is listed in modules.json, because everything there blocks the
// frame's first paint. Each engine is fetched the first time a view actually
// needs it: a solvent card never pays for 3Dmol, a 2D depiction never waits on
// it, and only the two graph views pull d3.
//
// A host may instead **pre-seed** an engine through `window.__gufeEngines`, in
// which case nothing is fetched at all. That is how gufe's standalone HTML
// export (`openfe view --html --self-contained`) produces a page that renders
// with no network: it inlines the engines and hands them over here. A seeded
// value may be the module itself or a promise for it.
// ============================================================================

function preseededEngine(name) {
  const seeded = globalThis.__gufeEngines && globalThis.__gufeEngines[name];
  return seeded ? Promise.resolve(seeded) : null;
}

var ThreeDmol = null;
let threeDmolPromise = null;

function load3Dmol() {
  if (threeDmolPromise) return threeDmolPromise;
  const seeded = preseededEngine('threeDmol');
  if (seeded) {
    threeDmolPromise = seeded.then(m => (ThreeDmol = m || window.$3Dmol));
    return threeDmolPromise;
  }
  threeDmolPromise = new Promise((resolve, reject) => {
    if (window.$3Dmol) { ThreeDmol = window.$3Dmol; resolve(ThreeDmol); return; }
    const script = document.createElement('script');
    script.src = 'https://3dmol.org/build/3Dmol-min.js';
    script.onload = () => {
      ThreeDmol = window.$3Dmol;
      if (ThreeDmol) resolve(ThreeDmol);
      else reject(new Error('3Dmol.js loaded but $3Dmol is undefined'));
    };
    script.onerror = () => reject(new Error('Failed to load 3Dmol.js'));
    document.head.appendChild(script);
  });
  return threeDmolPromise;
}

let rdkitPromise = null;

function loadRDKit() {
  if (rdkitPromise) return rdkitPromise;
  const seeded = preseededEngine('rdkit');
  if (seeded) {
    rdkitPromise = seeded.then(m => (window.RDKit = m));
    return rdkitPromise;
  }
  rdkitPromise = new Promise((resolve, reject) => {
    if (window.RDKit) { resolve(window.RDKit); return; }
    const script = document.createElement('script');
    script.src = 'https://unpkg.com/@rdkit/rdkit/dist/RDKit_minimal.js';
    script.onload = () => {
      window.initRDKitModule().then(rdkit => { window.RDKit = rdkit; resolve(rdkit); }).catch(reject);
    };
    script.onerror = () => reject(new Error('Failed to load RDKit'));
    document.head.appendChild(script);
  });
  return rdkitPromise;
}

let d3Promise = null;

function loadD3() {
  if (!d3Promise) {
    d3Promise = preseededEngine('d3') || import('https://cdn.jsdelivr.net/npm/d3@7/+esm');
  }
  return d3Promise;
}

// ============================================================================
// INPUT NORMALIZATION
//
// An input may arrive as an http(s) URL to fetch rather than as the content
// itself: that is how gufe's `openfe view` server hands over data (an SDF, a
// PDB, a GraphML, a JSON blob) that would not fit in the URL hash. A one-token
// string with no whitespace is a URL; real file content never is.
// ============================================================================

const isFetchableUrl = (v) => typeof v === 'string' && /^https?:\/\/[^\s]+$/.test(v.trim());

async function asText(v) {
  if (v == null) return null;
  if (isFetchableUrl(v)) {
    const r = await fetch(v.trim());
    if (!r.ok) throw new Error(`fetch ${v.trim()} failed: ${r.status} ${r.statusText}`);
    return await r.text();
  }
  if (v instanceof Blob) return await v.text();
  if (v instanceof ArrayBuffer) return new TextDecoder().decode(v);
  if (ArrayBuffer.isView(v)) return new TextDecoder().decode(v);
  return typeof v === 'string' ? v : String(v);
}

/** Object-valued inputs may arrive as an object, a JSON string, or a URL to one. */
async function asObject(v) {
  if (v == null) return null;
  if (typeof v === 'object' && !(v instanceof Blob) && !(v instanceof ArrayBuffer) && !ArrayBuffer.isView(v)) {
    return v;
  }
  const text = await asText(v);
  if (text == null || text === '') return null;
  return JSON.parse(text);
}

// ============================================================================
// SMALL DOM HELPERS
// ============================================================================

function el(tag, css, text) {
  const n = document.createElement(tag);
  if (css) n.style.cssText = css;
  if (text != null) n.textContent = text;
  return n;
}

function esc(s) {
  return String(s == null ? '' : s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
}

/**
 * Describe a thrown value. The engines do not always reject with an `Error` —
 * 3Dmol hands back a bare object when WebGL is unavailable, which turned every
 * such failure into the message "undefined".
 */
function errText(e) {
  if (e == null) return 'unknown error';
  if (errText(e)) return errText(e);
  const text = String(e);
  return text === '[object Object]' ? (e.name || 'unknown error') : text;
}

const fmt = (n) => n.toLocaleString('en-US');
const EM_DASH = '—';

const BTN_CSS =
  'background:' + T.btnBg + ';color:' + T.btnFg + ';border:1px solid ' + T.btnBorder +
  ';padding:4px 9px;font-size:11px;font-weight:bold;border-radius:3px;cursor:pointer;font-family:inherit;';

const SELECT_CSS =
  'background:' + T.selectBg + ';color:' + T.textPrimary + ';border:1px solid ' + T.selectBorder +
  ';border-radius:6px;padding:4px 8px;font-size:12px;cursor:pointer;font-family:inherit;';

/**
 * A row of mutually-exclusive buttons. `onPick(id)` fires on click; the active
 * button is highlighted. Returns the container, with `.setActive(id)` on it.
 */
function buttonGroup(items, active, onPick) {
  const group = el('div', 'display:flex;gap:4px;');
  const buttons = items.map(item => {
    const btn = el('button', BTN_CSS, item.label);
    btn.title = item.title || item.label;
    btn.onmouseover = () => { btn.style.background = T.btnBgHover; };
    btn.onmouseout = () => { btn.style.background = (active === item.id) ? T.btnBgActive : T.btnBg; };
    btn.onclick = () => { group.setActive(item.id); onPick(item.id); };
    group.appendChild(btn);
    return { id: item.id, btn };
  });
  group.setActive = (id) => {
    active = id;
    buttons.forEach(b => { b.btn.style.background = (b.id === active) ? T.btnBgActive : T.btnBg; });
  };
  group.setActive(active);
  return group;
}

/** "label <b>value</b>" with an optional colour dot — the stats readouts. */
function statChip(label, value, dotColor) {
  const chip = el('span', 'display:inline-flex;align-items:center;gap:5px;white-space:nowrap;');
  if (dotColor) {
    chip.appendChild(el('span',
      'width:8px;height:8px;border-radius:50%;background:' + dotColor + ';display:inline-block;'));
  }
  const txt = el('span');
  txt.innerHTML = esc(label) + ' <b style="color:' + T.textPrimary + ';">' + esc(value) + '</b>';
  chip.appendChild(txt);
  return chip;
}

function warnBanner(message) {
  return el('div',
    'margin:12px 16px;padding:8px 12px;border-radius:6px;font-size:12px;white-space:pre-wrap;' +
    'background:' + T.warnBg + ';color:' + T.warnFg + ';border:1px solid ' + T.warnBorder + ';',
    '⚠ ' + message);
}

/** A floating banner pinned to the top of a positioned host. */
function floatingWarning(host, message) {
  const warn = el('div', '', '⚠ ' + message);
  warn.style.cssText =
    'position:absolute;top:10px;left:50%;transform:translateX(-50%);max-width:90%;z-index:20;' +
    'padding:6px 14px;border-radius:6px;font-size:12px;' +
    'background:' + T.warnBg + ';color:' + T.warnFg + ';border:1px solid ' + T.warnBorder + ';';
  host.appendChild(warn);
  return warn;
}

/** A centred message filling its container — the placeholder / error state. */
function centredMessage(text, isError) {
  return el('div',
    'flex:1;display:flex;align-items:center;justify-content:center;text-align:center;padding:24px;' +
    'font-size:13px;color:' + (isError ? T.errorFg : T.textMuted2) + ';',
    text);
}

/** The standard header strip: bold title, muted subtitle, right-aligned stats. */
function headerStrip(title, subtitle) {
  const bar = el('div',
    'display:flex;align-items:baseline;gap:12px;flex-wrap:wrap;padding:9px 14px;flex-shrink:0;' +
    'background:' + T.toolbarBg + ';border-bottom:1px solid ' + T.toolbarBorder + ';');
  bar.titleEl = el('span', 'font-weight:700;font-size:15px;color:' + T.titleColor + ';letter-spacing:.02em;', title);
  bar.subtitleEl = el('span', 'font-size:12px;color:' + T.textMuted2 + ';', subtitle || '');
  bar.statsEl = el('div',
    'display:flex;align-items:center;gap:12px;flex-wrap:wrap;margin-left:auto;font-size:11px;color:' + T.textMuted + ';');
  bar.appendChild(bar.titleEl);
  bar.appendChild(bar.subtitleEl);
  bar.appendChild(bar.statsEl);
  return bar;
}

/** A 3D viewer host: an absolutely-filled container inside a flexible box. */
function viewerHost() {
  const wrap = el('div', 'flex:1;position:relative;min-height:0;min-width:0;');
  const container = el('div', 'position:absolute;inset:0;');
  wrap.appendChild(container);
  return { wrap, container };
}

// ============================================================================
// MOLECULE HELPERS — SDF in, internal shape out, MOL block back out
//
// Internal molecule shape:
//   { name, symbols: ['C', …], bonds: [[i, j, order], …], coords: [[x,y,z], …] }
// ============================================================================

const NL = String.fromCharCode(10);
const DOLLAR = String.fromCharCode(36);
const SDF_TERMINATOR = DOLLAR + DOLLAR + DOLLAR + DOLLAR;

/**
 * Parse a V2000 SDF/MOL record. Only the counts line + atom block + bond block
 * are read; property and data blocks are ignored. `bonds` indices are 0-based.
 */
function parseSDF(sdfText, fallbackName) {
  if (!sdfText || !sdfText.trim()) throw new Error('empty SDF');
  const lines = sdfText.replace(/\r/g, '').split(NL);
  if (lines.length < 4) throw new Error('SDF too short');

  const countsLine = lines[3];
  if (countsLine.indexOf('V3000') !== -1) throw new Error('V3000 molfiles are not supported');
  const nAtoms = parseInt(countsLine.substring(0, 3), 10);
  const nBonds = parseInt(countsLine.substring(3, 6), 10);
  if (!isFinite(nAtoms) || nAtoms <= 0) throw new Error('bad counts line: ' + countsLine);

  const coords = [];
  const symbols = [];
  for (let i = 0; i < nAtoms; i++) {
    const ln = lines[4 + i];
    if (ln == null) throw new Error('truncated atom block');
    coords.push([
      parseFloat(ln.substring(0, 10)) || 0,
      parseFloat(ln.substring(10, 20)) || 0,
      parseFloat(ln.substring(20, 30)) || 0
    ]);
    symbols.push(ln.substring(31, 34).trim() || 'X');
  }

  const bonds = [];
  for (let j = 0; j < (isFinite(nBonds) ? nBonds : 0); j++) {
    const ln = lines[4 + nAtoms + j];
    if (ln == null) break;
    const a = parseInt(ln.substring(0, 3), 10);
    const b = parseInt(ln.substring(3, 6), 10);
    const order = parseInt(ln.substring(6, 9), 10);
    if (!isFinite(a) || !isFinite(b)) continue;
    bonds.push([a - 1, b - 1, isFinite(order) ? order : 1]);
  }

  const title = (lines[0] || '').trim();
  return { name: title || fallbackName || 'molecule', symbols, bonds, coords };
}

/** MOL block (V2000) rebuilt from the internal shape — what RDKit is fed. */
function buildMolBlock(mol) {
  const nAtoms = mol.symbols.length;
  const nBonds = mol.bonds.length;
  const lines = [
    mol.name || '',
    '  Generated',
    '',
    String(nAtoms).padStart(3) + String(nBonds).padStart(3) + '  0  0  0  0  0  0  0  0999 V2000'
  ];
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

/** SDF = MOL block + the "$$$$" record terminator that 3Dmol expects. */
const buildSDF = (mol) => buildMolBlock(mol) + NL + SDF_TERMINATOR;

/** to_sdf() usually already terminates its record; make sure of it for 3Dmol. */
const ensureSDFTerminator = (sdf) =>
  sdf.indexOf(SDF_TERMINATOR) >= 0 ? sdf : sdf + NL + SDF_TERMINATOR;

/** The V2000 counts line is the 4th: aaabbb… (3 chars each). */
function parseCounts(sdf) {
  const lines = String(sdf).split(/\r?\n/);
  if (lines.length < 4) return null;
  const atoms = parseInt(lines[3].slice(0, 3), 10);
  const bonds = parseInt(lines[3].slice(3, 6), 10);
  return (isNaN(atoms) || isNaN(bonds)) ? null : { atoms, bonds };
}

/** 2D depiction SVG from a molblock or a SMILES string. */
function depictSVG(RDKit, source, size) {
  let rdmol = null;
  try {
    rdmol = RDKit.get_mol(source, JSON.stringify({ removeHs: true }));
    if (!rdmol) return null;
    try { rdmol.set_new_coords(true); } catch (e) { /* SMILES have no coords */ }
    return rdmol.get_svg(size, size) || null;
  } catch (e) {
    console.warn('[gufe] depictSVG threw -', errText(e));
    return null;
  } finally {
    if (rdmol) { try { rdmol.delete(); } catch (e) { /* already freed */ } }
  }
}

/** Drop an RDKit SVG into `box`, scaled to fill it. */
function placeDepiction(box, svg, size) {
  box.innerHTML = svg;
  const svgEl = box.querySelector('svg');
  if (!svgEl) return;
  svgEl.removeAttribute('width');
  svgEl.removeAttribute('height');
  if (!svgEl.getAttribute('viewBox')) svgEl.setAttribute('viewBox', '0 0 ' + size + ' ' + size);
  svgEl.setAttribute('preserveAspectRatio', 'xMidYMid meet');
  svgEl.style.cssText = 'width:100%;height:100%;max-width:100%;max-height:100%;';
}

// ============================================================================
// ATOM MAPPING
// ============================================================================

/** Pull an atom mapping out of whatever shape it arrived in. */
function parseMappingValue(val) {
  if (val == null) return null;
  if (val instanceof Map) return val.size > 0 ? val : null;
  if (Array.isArray(val)) {
    const m = new Map();
    val.forEach(pair => {
      if (Array.isArray(pair) && typeof pair[0] === 'number' && typeof pair[1] === 'number') {
        m.set(pair[0], pair[1]);
      }
    });
    return m.size > 0 ? m : null;
  }
  if (typeof val !== 'object') return null;

  const keys = Object.keys(val);
  if (keys.length && keys.every(k => /^-?\d+$/.test(k) && typeof val[k] === 'number')) {
    const m = new Map();
    keys.forEach(k => m.set(parseInt(k, 10), val[k]));
    return m.size > 0 ? m : null;
  }
  // Not a mapping itself — look one level in (GraphML edge data nests it).
  for (const k of keys) {
    let inner = val[k];
    if (typeof inner === 'string') {
      try { inner = JSON.parse(inner); } catch (e) { continue; }
    }
    const got = parseMappingValue(inner);
    if (got) return got;
  }
  return null;
}

/**
 * Normalize a mapping to `Map<aIdx, bIdx>` for a specific molecule pair,
 * flipping it if the indices only make sense the other way round.
 */
function toMapping(value, nA, nB) {
  const m = parseMappingValue(value);
  if (!m) return null;
  let maxK = -1, maxV = -1;
  m.forEach((v, k) => { if (k > maxK) maxK = k; if (v > maxV) maxV = v; });
  if ((maxK >= nA || maxV >= nB) && maxK < nB && maxV < nA) {
    const flipped = new Map();
    m.forEach((v, k) => flipped.set(v, k));
    return flipped;
  }
  return m;
}

// ─── Kabsch superposition (Pairs mode) ─────────────────────────────────────

function mat3Mul(a, b) {
  const r = new Array(9);
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      r[i * 3 + j] = a[i * 3] * b[j] + a[i * 3 + 1] * b[3 + j] + a[i * 3 + 2] * b[6 + j];
    }
  }
  return r;
}

const mat3Transpose = (m) => [m[0], m[3], m[6], m[1], m[4], m[7], m[2], m[5], m[8]];

const mat3Det = (m) =>
  m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) + m[2] * (m[3] * m[7] - m[4] * m[6]);

/** Jacobi eigen-decomposition of a symmetric 3×3, columns sorted descending. */
function jacobiSym3(M) {
  const a = M.slice();
  const v = [1, 0, 0, 0, 1, 0, 0, 0, 1];
  for (let sweep = 0; sweep < 50; sweep++) {
    if (Math.abs(a[1]) + Math.abs(a[2]) + Math.abs(a[5]) < 1e-12) break;
    for (const [p, q] of [[0, 1], [0, 2], [1, 2]]) {
      const apq = a[p * 3 + q];
      if (Math.abs(apq) < 1e-14) continue;
      const app = a[p * 3 + p], aqq = a[q * 3 + q];
      const theta = (aqq - app) / (2 * apq);
      const t = Math.abs(theta) > 1e10
        ? 1 / (2 * theta)
        : (theta >= 0 ? 1 : -1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
      const c = 1 / Math.sqrt(1 + t * t);
      const s = t * c;
      a[p * 3 + p] = app - t * apq;
      a[q * 3 + q] = aqq + t * apq;
      a[p * 3 + q] = 0; a[q * 3 + p] = 0;
      for (let r = 0; r < 3; r++) {
        if (r === p || r === q) continue;
        const arp = a[r * 3 + p], arq = a[r * 3 + q];
        a[r * 3 + p] = c * arp - s * arq; a[p * 3 + r] = a[r * 3 + p];
        a[r * 3 + q] = s * arp + c * arq; a[q * 3 + r] = a[r * 3 + q];
      }
      for (let k = 0; k < 3; k++) {
        const vkp = v[k * 3 + p], vkq = v[k * 3 + q];
        v[k * 3 + p] = c * vkp - s * vkq;
        v[k * 3 + q] = s * vkp + c * vkq;
      }
    }
  }
  const order = [0, 1, 2].sort((x, y) => a[y * 4] - a[x * 4]);
  const vectors = new Array(9);
  order.forEach((src, col) => {
    vectors[col] = v[src];
    vectors[3 + col] = v[3 + src];
    vectors[6 + col] = v[6 + src];
  });
  return { vectors };
}

/** Rotation + translation putting Q onto P (paired, same length). */
function kabsch(P, Q) {
  const n = Math.min(P.length, Q.length);
  if (n < 1) return null;
  const cP = [0, 0, 0], cQ = [0, 0, 0];
  for (let i = 0; i < n; i++) {
    for (let k = 0; k < 3; k++) { cP[k] += P[i][k]; cQ[k] += Q[i][k]; }
  }
  for (let k = 0; k < 3; k++) { cP[k] /= n; cQ[k] /= n; }
  if (n < 3) {
    return { R: [1, 0, 0, 0, 1, 0, 0, 0, 1], t: [cP[0] - cQ[0], cP[1] - cQ[1], cP[2] - cQ[2]] };
  }

  const H = [0, 0, 0, 0, 0, 0, 0, 0, 0];
  for (let k = 0; k < n; k++) {
    const p = [P[k][0] - cP[0], P[k][1] - cP[1], P[k][2] - cP[2]];
    const q = [Q[k][0] - cQ[0], Q[k][1] - cQ[1], Q[k][2] - cQ[2]];
    for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) H[i * 3 + j] += p[i] * q[j];
  }

  const Ht = mat3Transpose(H);
  const V = jacobiSym3(mat3Mul(Ht, H)).vectors;
  const U = jacobiSym3(mat3Mul(H, Ht)).vectors;
  // Fix each U column's sign against H·v so U and V describe the same rotation.
  for (let col = 0; col < 3; col++) {
    const vx = V[col], vy = V[3 + col], vz = V[6 + col];
    const hx = H[0] * vx + H[1] * vy + H[2] * vz;
    const hy = H[3] * vx + H[4] * vy + H[5] * vz;
    const hz = H[6] * vx + H[7] * vy + H[8] * vz;
    if (hx * U[col] + hy * U[3 + col] + hz * U[6 + col] < 0) {
      U[col] = -U[col]; U[3 + col] = -U[3 + col]; U[6 + col] = -U[6 + col];
    }
  }

  const Vt = mat3Transpose(V);
  let R = mat3Mul(U, Vt);
  if (mat3Det(R) < 0) {   // reflection, not a rotation — flip the last column
    U[2] = -U[2]; U[5] = -U[5]; U[8] = -U[8];
    R = mat3Mul(U, Vt);
  }
  return {
    R,
    t: [
      cP[0] - (R[0] * cQ[0] + R[1] * cQ[1] + R[2] * cQ[2]),
      cP[1] - (R[3] * cQ[0] + R[4] * cQ[1] + R[5] * cQ[2]),
      cP[2] - (R[6] * cQ[0] + R[7] * cQ[1] + R[8] * cQ[2])
    ]
  };
}

const applyRT = (c, R, t) => [
  R[0] * c[0] + R[1] * c[1] + R[2] * c[2] + t[0],
  R[3] * c[0] + R[4] * c[1] + R[5] * c[2] + t[1],
  R[6] * c[0] + R[7] * c[1] + R[8] * c[2] + t[2]
];

// ============================================================================
// THE MAPPING VIEWER
//
// Used by three views: the LigandAtomMapping page, the right pane of the
// LigandNetwork graph, and the bottom of a Transformation. Owns its 3Dmol
// viewers, its mode switcher and its view-sync loop.
// ============================================================================

const MAPPING_MODES = [
  { id: 'plain',   label: '3D',      title: 'Plain 3D view' },
  { id: 'colored', label: '3D-Map',  title: 'Color-coded by mapping' },
  { id: 'lines',   label: 'Pairs',   title: 'Dashed lines between mapped atoms' },
  { id: 'overlay', label: 'Overlay', title: 'Both molecules superimposed' },
  { id: '2d',      label: '2D',      title: '2D depictions with highlights' }
];

function createMappingViewer(host, initialMode) {
  host.style.position = 'relative';
  const area = el('div', 'flex:1;display:flex;flex-direction:column;min-height:0;');
  host.appendChild(area);

  let mode = initialMode || 'colored';
  let mols = null;         // [molA, molB]
  let mapping = null;      // Map<aIdx, bIdx> | null
  let boxes = [];
  let token = 0;           // guards async engine loads against a newer render
  let syncing = false, syncRunning = false, stopped = false;

  const switcher = el('div',
    'position:absolute;bottom:10px;right:10px;display:flex;gap:4px;padding:4px;border-radius:6px;z-index:10;' +
    'background:' + T.switcherBg + ';box-shadow:0 2px 8px rgba(0,0,0,0.25);');
  const modeButtons = buttonGroup(MAPPING_MODES, mode, id => { mode = id; render(); });
  switcher.appendChild(modeButtons);
  host.appendChild(switcher);

  function makeBox(labelText, is2D) {
    const box = el('div', 'flex:1;display:flex;flex-direction:column;position:relative;min-width:0;min-height:0;');
    box.appendChild(el('div',
      'padding:4px 10px;font-size:12px;font-weight:bold;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;' +
      'color:' + T.labelFg + ';background:' + T.labelBg + ';', labelText));
    const container = el('div', is2D
      ? 'flex:1;display:flex;align-items:center;justify-content:center;overflow:hidden;background:' + T.canvas2DBg + ';'
      : 'flex:1;position:relative;');
    box.appendChild(container);
    area.appendChild(box);
    const entry = { container, viewer: null };
    boxes.push(entry);
    return entry;
  }

  function clearBoxes() {
    boxes.forEach(b => { if (b.viewer) { try { b.viewer.clear(); } catch (e) { /* already gone */ } } });
    area.innerHTML = '';
    boxes = [];
  }

  function message(text, isError) {
    clearBoxes();
    area.appendChild(centredMessage(text, isError));
  }

  /** Keep the two 3D panes pointing the same way as the user rotates either. */
  function startSyncLoop() {
    if (syncRunning || boxes.length < 2) return;
    syncRunning = true;
    const lastViews = boxes.map(() => '');
    (function loop() {
      if (stopped) { syncRunning = false; return; }
      if (!syncing && boxes.length >= 2) {
        for (let i = 0; i < boxes.length; i++) {
          if (!boxes[i].viewer) continue;
          const cur = JSON.stringify(boxes[i].viewer.getView());
          if (cur === lastViews[i]) continue;
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
      requestAnimationFrame(loop);
    })();
  }

  function newViewer(container) {
    return ThreeDmol.createViewer(container, { backgroundColor: T.viewerBg });
  }

  // ─── modes ───

  function renderPlain() {
    clearBoxes();
    mols.forEach(mol => {
      const b = makeBox(mol.name);
      b.viewer = newViewer(b.container);
      b.viewer.addModel(buildSDF(mol), 'sdf');
      b.viewer.setStyle({}, {
        stick: { radius: 0.15, colorscheme: 'Jmol' },
        sphere: { scale: 0.25, colorscheme: 'Jmol' }
      });
      b.viewer.zoomTo();
      b.viewer.render();
    });
    startSyncLoop();
  }

  function renderColored() {
    if (!mapping) { renderPlain(); return; }
    clearBoxes();
    const cores = [new Set(mapping.keys()), new Set(mapping.values())];
    mols.forEach((mol, i) => {
      const b = makeBox(mol.name);
      b.viewer = newViewer(b.container);
      b.viewer.addModel(buildSDF(mol), 'sdf');
      b.viewer.setStyle({}, {
        stick: { radius: 0.15, color: T.colorCore },
        sphere: { scale: 0.25, color: T.colorCore }
      });
      const uniqueColor = i === 0 ? T.colorUniqueA : T.colorUniqueB;
      for (let ai = 0; ai < mol.symbols.length; ai++) {
        if (cores[i].has(ai)) continue;
        b.viewer.addStyle({ serial: ai + 1 }, {
          stick: { radius: 0.18, color: uniqueColor },
          sphere: { scale: 0.32, color: uniqueColor }
        });
      }
      b.viewer.zoomTo();
      b.viewer.render();
    });
    startSyncLoop();
  }

  /**
   * Superimpose B onto A by the mapped atoms, then lift B clear along the
   * molecule's thinnest axis and draw a dashed cylinder per mapped pair.
   */
  function renderLines() {
    if (!mapping) { renderPlain(); return; }
    clearBoxes();
    const [molA, molB] = mols;

    const P = [], Q = [];
    mapping.forEach((bIdx, aIdx) => {
      if (molA.coords[aIdx] && molB.coords[bIdx]) { P.push(molA.coords[aIdx]); Q.push(molB.coords[bIdx]); }
    });
    const rt = kabsch(P, Q);
    const aligned = rt ? molB.coords.map(c => applyRT(c, rt.R, rt.t)) : molB.coords.map(c => c.slice());

    const extents = (coords) => {
      const mn = [Infinity, Infinity, Infinity], mx = [-Infinity, -Infinity, -Infinity];
      coords.forEach(c => {
        for (let k = 0; k < 3; k++) { if (c[k] < mn[k]) mn[k] = c[k]; if (c[k] > mx[k]) mx[k] = c[k]; }
      });
      return { min: mn, max: mx, span: [mx[0] - mn[0], mx[1] - mn[1], mx[2] - mn[2]] };
    };
    const extA = extents(molA.coords), extB = extents(aligned);

    let liftAxis = 0;
    if (extA.span[1] < extA.span[liftAxis]) liftAxis = 1;
    if (extA.span[2] < extA.span[liftAxis]) liftAxis = 2;
    const GAP = 2.5;
    const lift = Math.max(
      (extA.max[liftAxis] - extB.min[liftAxis]) + GAP,
      0.6 * Math.max(...extA.span) + GAP
    );
    const shifted = {
      name: molB.name, symbols: molB.symbols, bonds: molB.bonds,
      coords: aligned.map(c => { const out = c.slice(); out[liftAxis] += lift; return out; })
    };

    const b = makeBox(molA.name + ' ↔ ' + molB.name + '  (' + mapping.size + ' mapped pairs)');
    b.viewer = newViewer(b.container);
    b.viewer.addModel(buildSDF(molA), 'sdf');
    b.viewer.addModel(buildSDF(shifted), 'sdf');
    b.viewer.setStyle({ model: 0 }, {
      stick: { radius: 0.15, color: T.linesMolA }, sphere: { scale: 0.22, color: T.linesMolA }
    });
    b.viewer.setStyle({ model: 1 }, {
      stick: { radius: 0.15, color: T.linesMolB }, sphere: { scale: 0.22, color: T.linesMolB }
    });
    mapping.forEach((bIdx, aIdx) => {
      const pa = molA.coords[aIdx], pb = shifted.coords[bIdx];
      if (!pa || !pb) return;
      b.viewer.addCylinder({
        start: { x: pa[0], y: pa[1], z: pa[2] },
        end: { x: pb[0], y: pb[1], z: pb[2] },
        radius: 0.04, dashed: true, fromCap: 'round', toCap: 'round', color: T.linesDash
      });
    });
    b.viewer.zoomTo();
    // Turn the camera so the lift axis runs up the screen.
    if (liftAxis === 2) b.viewer.rotate(90, 'x');
    else if (liftAxis === 0) b.viewer.rotate(-90, 'z');
    b.viewer.render();
  }

  function renderOverlay() {
    clearBoxes();
    const b = makeBox(mols[0].name + ' + ' + mols[1].name + '  (overlay)');
    b.viewer = newViewer(b.container);
    mols.forEach(mol => b.viewer.addModel(buildSDF(mol), 'sdf'));
    [T.overlayMolA, T.overlayMolB].forEach((color, i) => {
      b.viewer.setStyle({ model: i }, {
        stick: { radius: 0.15, color, opacity: 0.7 },
        sphere: { scale: 0.22, color, opacity: 0.7 }
      });
    });
    b.viewer.zoomTo();
    b.viewer.render();
  }

  function render2D() {
    clearBoxes();
    const pending = mols.map(mol => {
      const b = makeBox(mol.name, true);
      b.container.innerHTML = '<div style="color:' + T.loadingFg + ';font-size:12px;">Loading 2D depiction…</div>';
      return { container: b.container, mol };
    });
    const myToken = token;

    loadRDKit().then(RDKit => {
      if (myToken !== token) return;
      const cores = mapping ? [new Set(mapping.keys()), new Set(mapping.values())] : null;
      pending.forEach((p, idx) => {
        const rdmol = RDKit.get_mol(buildMolBlock(p.mol));
        if (!rdmol) {
          p.container.innerHTML = '<div style="color:' + T.errorFg + ';">Failed to parse</div>';
          return;
        }
        try { rdmol.set_new_coords(true); } catch (e) { /* keep the given coords */ }
        let svg;
        try {
          if (cores) {
            const atoms = p.mol.symbols.map((_, ai) => ai);
            svg = rdmol.get_svg_with_highlights(JSON.stringify({
              atoms,
              highlightAtomColors: Object.fromEntries(atoms.map(ai =>
                [ai, cores[idx].has(ai) ? T.rgbCore : (idx === 0 ? T.rgbUniqueA : T.rgbUniqueB)])),
              highlightAtomRadii: Object.fromEntries(atoms.map(ai => [ai, 0.3]))
            }));
          } else {
            svg = rdmol.get_svg(300, 300);
          }
        } catch (e) {
          svg = rdmol.get_svg(300, 300);
        }
        rdmol.delete();
        p.container.innerHTML = svg;
        const svgEl = p.container.querySelector('svg');
        if (svgEl) {
          svgEl.style.maxWidth = '100%';
          svgEl.style.maxHeight = '100%';
          svgEl.setAttribute('preserveAspectRatio', 'xMidYMid meet');
        }
      });
    }).catch(err => {
      if (myToken !== token) return;
      pending.forEach(p => {
        p.container.innerHTML =
          '<div style="color:' + T.errorFg + ';font-size:12px;padding:10px;">RDKit failed to load: ' + esc(errText(err)) + '</div>';
      });
    });
  }

  function render() {
    if (!mols) return;
    const myToken = ++token;

    // The 2D depiction is pure RDKit — it must never wait on the 3D engine.
    if (mode === '2d') { render2D(); return; }

    message('Loading 3D viewer…');
    load3Dmol().then(() => {
      if (myToken !== token) return;   // a newer mode or molecule pair won
      clearBoxes();
      if (mols.length < 2) renderPlain();
      else if (mode === 'colored') renderColored();
      else if (mode === 'lines') renderLines();
      else if (mode === 'overlay') renderOverlay();
      else renderPlain();
    }).catch(e => {
      if (myToken !== token) return;
      message('3D viewer unavailable: ' + errText(e), true);
    });
  }

  return {
    /** Show a molecule pair. `mappingValue` may be an object, array or Map. */
    show(molA, molB, mappingValue) {
      mols = molB ? [molA, molB] : [molA];
      mapping = molB ? toMapping(mappingValue, molA.coords.length, molB.coords.length) : null;
      render();
      return mapping;
    },
    message,
    resize() {
      boxes.forEach(b => { if (b.viewer) { b.viewer.resize(); b.viewer.render(); } });
    },
    destroy() {
      stopped = true;
      token++;
      clearBoxes();
    }
  };
}

// ============================================================================
// PROTEIN RENDERING — shared by the ProteinComponent view and the
// ChemicalSystem detail pane. Everything goes through 3Dmol *selections*, never
// per-atom DOM work, so tens of thousands of atoms stay responsive.
// ============================================================================

const WATER_RESN = ['HOH', 'WAT', 'SOL', 'TIP3'];
const SEL_POLYMER = { hetflag: false };
const SEL_HETERO = { hetflag: true };
const SEL_WATER = { resn: WATER_RESN };

const PROTEIN_CONFIG = {
  stick: { radius: 0.15 },
  sphere: { scale: 0.30 },
  hetero: { stickRadius: 0.20, sphereScale: 0.28 },
  water: { stickRadius: 0.05, sphereScale: 0.18 },
  surfaceOpacity: 0.85,
  // Above this many atoms a surface is slow enough to be worth warning about.
  surfaceAtomWarn: 40000
};

/**
 * Single pass over the PDB text for the stats readout. Only the first MODEL is
 * counted (NMR ensembles would otherwise multiply every number). Residues are
 * keyed by chain + sequence number + insertion code + name, which is what makes
 * them unique in a PDB.
 */
function parsePdbStats(pdbText) {
  const chains = new Set(), residues = new Set();
  let atoms = 0, hetatms = 0, waters = 0, resiMin = Infinity, resiMax = -Infinity;

  for (const line of pdbText.split(/\r?\n/)) {
    const rec = line.slice(0, 6);
    if (rec === 'ENDMDL') break;
    if (rec !== 'ATOM  ' && rec !== 'HETATM') continue;

    atoms++;
    if (rec === 'HETATM') hetatms++;

    const resName = line.slice(17, 20).trim();
    const chainId = line.slice(21, 22).trim() || '_';
    const resSeq = line.slice(22, 26).trim();
    const iCode = line.slice(26, 27).trim();

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
    chains: chains.size, residues: residues.size, atoms, hetatms, waters,
    // Non-water HETATM records — what the "hetero/ligands" toggle governs.
    heteroNonWater: hetatms - waters,
    resiMin: resiMin === Infinity ? 0 : resiMin,
    resiMax: resiMax === -Infinity ? 0 : resiMax
  };
}

function proteinStatsText(stats) {
  return fmt(stats.chains) + ' chains · ' + fmt(stats.residues) + ' residues · ' +
    fmt(stats.atoms) + ' atoms · ' + fmt(stats.hetatms) + ' HETATM' +
    (stats.waters > 0 ? ' (' + fmt(stats.waters) + ' water)' : '');
}

/** Colour arguments for a scheme, valid for any representation. */
function proteinColorArgs(scheme, stats) {
  if (scheme === 'chain') return { colorscheme: 'chain' };
  if (scheme === 'ss') return { colorscheme: 'ssJmol' };
  if (scheme === 'spectrum' && stats && stats.resiMax > stats.resiMin) {
    // A residue-index gradient works for every representation, unlike 3Dmol's
    // cartoon-only `color: 'spectrum'`.
    return { colorscheme: { prop: 'resi', gradient: 'roygb', min: stats.resiMin, max: stats.resiMax } };
  }
  return { colorscheme: 'Jmol' };
}

/**
 * Re-apply every protein style from scratch.
 *
 * The setStyle calls are ordered deliberately: each one *replaces* the style of
 * the atoms it matches, so waters (matched last) win over the blanket hetero
 * rule they would otherwise fall under.
 *
 * `onStatus(message, kind)` is optional — a surface is the one genuinely
 * expensive operation here, and is computed off a timeout so its message paints.
 */
function applyProteinStyles(viewer, opts, stats, onStatus) {
  const status = onStatus || (() => {});
  const color = proteinColorArgs(opts.color, stats);

  try { viewer.removeAllSurfaces(); } catch (e) { /* none yet */ }

  viewer.setStyle({}, {});
  viewer.setStyle(SEL_POLYMER,
    opts.rep === 'stick' ? { stick: { radius: PROTEIN_CONFIG.stick.radius, ...color } }
      : opts.rep === 'sphere' ? { sphere: { scale: PROTEIN_CONFIG.sphere.scale, ...color } }
        // For 'surface' the shell is added separately; leave the atoms bare so
        // it is not cluttered from the inside.
        : opts.rep === 'surface' ? {}
          : { cartoon: { ...color } });

  viewer.setStyle(SEL_HETERO, opts.hetero ? {
    stick: { radius: PROTEIN_CONFIG.hetero.stickRadius, colorscheme: 'Jmol' },
    sphere: { scale: PROTEIN_CONFIG.hetero.sphereScale, colorscheme: 'Jmol' }
  } : {});

  viewer.setStyle(SEL_WATER, opts.waters ? {
    stick: { radius: PROTEIN_CONFIG.water.stickRadius, colorscheme: 'Jmol' },
    sphere: { scale: PROTEIN_CONFIG.water.sphereScale, colorscheme: 'Jmol' }
  } : {});

  if (opts.rep !== 'surface') {
    status(null);
    viewer.render();
    return;
  }

  status((stats && stats.atoms > PROTEIN_CONFIG.surfaceAtomWarn)
    ? 'Computing surface (large structure, this may take a while)…'
    : 'Computing surface…');
  viewer.render();
  setTimeout(() => {
    try {
      // 3Dmol v2 returns a promise; v1 returns a surface id.
      Promise.resolve(viewer.addSurface(
        ThreeDmol.SurfaceType.VDW,
        { opacity: PROTEIN_CONFIG.surfaceOpacity, ...color },
        SEL_POLYMER
      )).then(() => { status(null); viewer.render(); })
        .catch(err => status('Surface failed: ' + errText(err), 'error'));
    } catch (err) {
      status('Surface failed: ' + errText(err), 'error');
    }
  }, 30);
}

// ============================================================================
// COMPONENT DESCRIPTORS
//
// ChemicalSystem, Transformation and AlchemicalNetwork all describe the
// components they contain with the same shape: { label, type, name } plus one
// of { sdf, smiles } / { pdb } / solvent fields / { error }.
// ============================================================================

const TYPE_COLORS = {
  SmallMoleculeComponent: '#7c3aed',
  ProteinComponent: '#0f766e',
  SolvatedPDBComponent: '#0369a1',
  SolventComponent: '#0284c7',
  ProteinMembraneComponent: '#b45309'
};
const FALLBACK_COLORS = ['#be123c', '#a16207', '#4d7c0f', '#1d4ed8', '#86198f', '#0f766e'];

/** Stable badge colour per component type — unknown types hash to a fallback. */
function typeColor(type) {
  if (TYPE_COLORS[type]) return TYPE_COLORS[type];
  let h = 0;
  for (let i = 0; i < type.length; i++) h = (h * 31 + type.charCodeAt(i)) >>> 0;
  return FALLBACK_COLORS[h % FALLBACK_COLORS.length];
}

/** "SmallMoleculeComponent" -> "SmallMolecule"; keeps signatures readable. */
const shortType = (t) => {
  const s = String(t || '');
  return s.endsWith('Component') && s.length > 9 ? s.slice(0, -9) : s;
};

/** Short human word for a type — used in composition summaries. */
function typeWord(type) {
  if (type === 'SmallMoleculeComponent') return 'ligand';
  if (type === 'SolventComponent') return 'solvent';
  if (type.indexOf('Protein') !== -1) return 'protein';
  return shortType(type).toLowerCase() || type;
}

function typeBadge(type, short) {
  return el('span',
    'font-size:10px;font-weight:700;letter-spacing:.02em;padding:2px 6px;border-radius:4px;' +
    'color:#fff;white-space:nowrap;background:' + typeColor(type) + ';',
    short ? shortType(type) : type);
}

function fieldRow(label, value) {
  const row = el('div', 'display:flex;gap:10px;align-items:baseline;padding:5px 0;');
  row.appendChild(el('div',
    'flex:0 0 130px;font-size:11px;text-transform:uppercase;letter-spacing:.04em;color:' + T.textMuted2 + ';',
    label));
  row.appendChild(el('div',
    'flex:1;min-width:0;font-size:13px;word-break:break-all;color:' + T.textPrimary +
    ";font-family:ui-monospace,SFMono-Regular,Menlo,monospace;",
    value == null || value === '' ? EM_DASH : String(value)));
  return row;
}

/** The solvent settings card — a SolventComponent has no coordinates. */
function solventCard(fields) {
  const card = el('div',
    'max-width:520px;padding:14px 16px;border-radius:10px;border:1px solid ' + T.cardBorder +
    ';background:' + T.cardBg + ';');
  card.appendChild(el('div',
    'font-size:12px;font-weight:700;margin-bottom:6px;letter-spacing:.03em;color:' + T.titleColor + ';',
    'SOLVENT SETTINGS'));
  card.appendChild(fieldRow('SMILES', fields.smiles));
  card.appendChild(fieldRow('Positive ion', fields.positive_ion));
  card.appendChild(fieldRow('Negative ion', fields.negative_ion));
  card.appendChild(fieldRow('Concentration', fields.ion_concentration));
  card.appendChild(fieldRow('Neutralize', fields.neutralize == null ? null : (fields.neutralize ? 'yes' : 'no')));
  return card;
}

// ============================================================================
// VIEW: SmallMoleculeComponent — 2D depiction | 3D viewer, with an info bar
// ============================================================================

const SMALL_MOL_STYLES = [
  { id: 'stick', label: 'Stick', title: 'Sticks only' },
  { id: 'ball', label: 'Ball+Stick', title: 'Ball and stick' },
  { id: 'sphere', label: 'Sphere', title: 'Space-filling spheres' }
];

const SMALL_MOL_SPECS = {
  stick: { stick: { radius: 0.15, colorscheme: 'Jmol' } },
  ball: { stick: { radius: 0.12, colorscheme: 'Jmol' }, sphere: { scale: 0.28, colorscheme: 'Jmol' } },
  sphere: { sphere: { scale: 1.0, colorscheme: 'Jmol' } }
};

const DEPICT_SIZE = 400;

async function viewSmallMolecule(host, inputs) {
  const sdf = await asText(inputs['molecule.sdf']);
  const name = await asText(inputs.name);
  const smiles = await asText(inputs.smiles);
  const charge = inputs.total_charge;

  const header = headerStrip(name || 'Unnamed molecule', 'SmallMoleculeComponent');
  host.appendChild(header);

  const split = el('div', 'flex:1;display:flex;flex-direction:row;overflow:hidden;min-height:0;');
  host.appendChild(split);

  const left = el('div', 'flex:1 1 50%;min-width:0;display:flex;flex-direction:column;');
  const right = el('div', 'flex:1 1 50%;min-width:0;display:flex;flex-direction:column;position:relative;');
  split.appendChild(left);
  split.appendChild(el('div', 'width:1px;flex-shrink:0;background:' + T.splitBorder + ';'));
  split.appendChild(right);

  const paneLabel = (text) => el('div',
    'flex-shrink:0;padding:4px 10px;font-size:12px;font-weight:bold;color:' + T.labelFg +
    ';background:' + T.labelBg + ';', text);

  left.appendChild(paneLabel('2D'));
  const depictBox = el('div',
    'flex:1;min-height:0;display:flex;align-items:center;justify-content:center;overflow:hidden;padding:8px;' +
    'background:' + T.canvas2DBg + ';');
  left.appendChild(depictBox);

  right.appendChild(paneLabel('3D'));
  const host3D = viewerHost();
  right.appendChild(host3D.wrap);

  // ─── info bar ───
  const infoBar = el('div',
    'flex-shrink:0;display:flex;flex-wrap:wrap;align-items:baseline;gap:6px 20px;padding:8px 16px;font-size:12px;' +
    'background:' + T.toolbarBg + ';border-top:1px solid ' + T.toolbarBorder + ';color:' + T.textPrimary + ';');
  host.appendChild(infoBar);

  const counts = sdf ? parseCounts(sdf) : null;
  [
    ['Name', name || EM_DASH, false],
    ['SMILES', smiles || EM_DASH, true],
    ['Charge', charge == null ? EM_DASH : String(charge), false],
    ['Atoms', counts ? String(counts.atoms) : EM_DASH, false],
    ['Bonds', counts ? String(counts.bonds) : EM_DASH, false]
  ].forEach(([label, value, mono]) => {
    const cell = el('div', 'display:flex;align-items:baseline;gap:6px;min-width:0;');
    cell.appendChild(el('span',
      'font-size:10px;font-weight:700;letter-spacing:.08em;text-transform:uppercase;flex-shrink:0;color:' +
      T.textMuted2 + ';', label));
    const v = el('span',
      'user-select:text;cursor:text;color:' + T.textPrimary +
      (mono ? ';font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:11px;overflow-wrap:anywhere;' : ''),
      value);
    v.title = value;
    cell.appendChild(v);
    infoBar.appendChild(cell);
  });

  if (!sdf || !sdf.trim()) {
    depictBox.appendChild(centredMessage('No molecule provided'));
    host3D.container.appendChild(centredMessage('No molecule provided'));
    return {};
  }

  // ─── 2D ───
  depictBox.appendChild(centredMessage('Loading 2D depiction…'));
  loadRDKit().then(RDKit => {
    const svg = depictSVG(RDKit, sdf, DEPICT_SIZE);
    if (svg) placeDepiction(depictBox, svg, DEPICT_SIZE);
    else { depictBox.innerHTML = ''; depictBox.appendChild(centredMessage('Failed to parse molecule', true)); }
  }).catch(err => {
    depictBox.innerHTML = '';
    depictBox.appendChild(centredMessage('RDKit failed to load: ' + errText(err), true));
  });

  // ─── 3D ───
  let viewer = null;
  let style = 'stick';
  let spinning = false;

  const switcher = el('div',
    'position:absolute;bottom:10px;right:10px;display:flex;gap:4px;padding:4px;border-radius:6px;z-index:10;' +
    'background:' + T.switcherBg + ';box-shadow:0 2px 8px rgba(0,0,0,0.25);');
  switcher.appendChild(buttonGroup(SMALL_MOL_STYLES, style, id => {
    style = id;
    if (viewer) { viewer.setStyle({}, SMALL_MOL_SPECS[id]); viewer.render(); }
  }));
  const spinBtn = el('button', BTN_CSS + 'margin-left:4px;', 'Spin');
  spinBtn.title = 'Toggle continuous rotation';
  spinBtn.onclick = () => {
    spinning = !spinning;
    spinBtn.style.background = spinning ? T.btnBgActive : T.btnBg;
    if (viewer) { try { viewer.spin(spinning ? 'y' : false); } catch (e) { /* v1 quirk */ } }
  };
  switcher.appendChild(spinBtn);
  right.appendChild(switcher);

  host3D.container.appendChild(centredMessage('Loading 3D viewer…'));
  load3Dmol().then(() => {
    host3D.container.innerHTML = '';
    viewer = ThreeDmol.createViewer(host3D.container, { backgroundColor: T.viewerBg });
    viewer.addModel(ensureSDFTerminator(sdf), 'sdf');
    viewer.setStyle({}, SMALL_MOL_SPECS[style]);
    viewer.zoomTo();
    viewer.render();
  }).catch(e => {
    host3D.container.innerHTML = '';
    host3D.container.appendChild(centredMessage('3D render failed: ' + errText(e), true));
  });

  return {
    onResize() { if (viewer) { viewer.resize(); viewer.render(); } },
    cleanup() {
      if (!viewer) return;
      try { viewer.spin(false); } catch (e) { /* v1 quirk */ }
      try { viewer.clear(); } catch (e) { /* already gone */ }
    }
  };
}

// ============================================================================
// VIEW: ProteinComponent — toolbar + full-pane 3Dmol viewer
// ============================================================================

const PROTEIN_REPS = [
  { id: 'cartoon', label: 'Cartoon', title: 'Ribbon / cartoon backbone' },
  { id: 'surface', label: 'Surface', title: 'Molecular (VDW) surface' },
  { id: 'stick', label: 'Stick', title: 'All-atom sticks' },
  { id: 'sphere', label: 'Sphere', title: 'Space-filling spheres' }
];

const PROTEIN_COLOR_SCHEMES = [
  { id: 'chain', label: 'Chain' },
  { id: 'spectrum', label: 'Spectrum' },
  { id: 'ss', label: 'Secondary structure' },
  { id: 'element', label: 'Element' }
];

async function viewProtein(host, inputs) {
  const pdb = await asText(inputs['protein.pdb']);
  const name = await asText(inputs.name);

  const opts = { rep: 'cartoon', color: 'chain', waters: false, hetero: true };
  let viewer = null;
  let stats = null;

  // ─── toolbar ───
  const toolbar = el('div',
    'display:flex;align-items:center;gap:14px;flex-wrap:wrap;padding:8px 14px;flex-shrink:0;font-size:12px;' +
    'background:' + T.toolbarBg + ';border-bottom:1px solid ' + T.toolbarBorder + ';color:' + T.textPrimary + ';');
  host.appendChild(toolbar);

  toolbar.appendChild(el('span',
    'font-weight:700;font-size:14px;letter-spacing:.02em;color:' + T.titleColor + ';', name || 'Protein'));

  const groupLabel = (text) => el('span', 'font-size:11px;color:' + T.textMuted + ';', text);

  toolbar.appendChild(groupLabel('Style:'));
  toolbar.appendChild(buttonGroup(PROTEIN_REPS, opts.rep, id => { opts.rep = id; restyle(); }));

  toolbar.appendChild(groupLabel('Color:'));
  const colorSelect = el('select', SELECT_CSS);
  PROTEIN_COLOR_SCHEMES.forEach(c => {
    const o = el('option', '', c.label);
    o.value = c.id;
    colorSelect.appendChild(o);
  });
  colorSelect.value = opts.color;
  colorSelect.addEventListener('change', () => { opts.color = colorSelect.value; restyle(); });
  toolbar.appendChild(colorSelect);

  const toggles = el('div', 'display:flex;gap:4px;');
  toolbar.appendChild(toggles);
  [
    ['waters', 'Waters', 'Show water molecules', restyle],
    ['hetero', 'Hetero', 'Show hetero atoms / ligands / ions / lipids', restyle],
    ['spin', 'Spin', 'Rotate the view continuously', () => { if (viewer) viewer.spin(opts.spin); }]
  ].forEach(([key, label, title, onChange]) => {
    const btn = el('button', BTN_CSS, label);
    btn.title = title;
    btn.style.background = opts[key] ? T.btnBgActive : T.btnBg;
    btn.onclick = () => {
      opts[key] = !opts[key];
      btn.style.background = opts[key] ? T.btnBgActive : T.btnBg;
      onChange();
    };
    toggles.appendChild(btn);
  });

  const statsEl = el('span',
    'margin-left:auto;font-size:11px;white-space:nowrap;color:' + T.textMuted2 + ';');
  toolbar.appendChild(statsEl);

  // ─── viewer + status overlay ───
  const pane = viewerHost();
  host.appendChild(pane.wrap);

  const statusEl = el('div',
    'position:absolute;top:12px;left:50%;transform:translateX(-50%);padding:6px 14px;border-radius:6px;' +
    'font-size:12px;z-index:20;display:none;pointer-events:none;');
  pane.wrap.appendChild(statusEl);

  function showStatus(msg, kind) {
    if (msg == null) { statusEl.style.display = 'none'; return; }
    statusEl.textContent = msg;
    statusEl.style.display = 'block';
    const isError = kind === 'error';
    statusEl.style.background = isError ? T.warnBg : T.toolbarBg;
    statusEl.style.color = isError ? T.warnFg : T.textMuted;
    statusEl.style.border = '1px solid ' + (isError ? T.warnBorder : T.toolbarBorder);
  }

  function restyle() {
    if (viewer) applyProteinStyles(viewer, opts, stats, showStatus);
  }

  if (!pdb || !pdb.trim()) {
    showStatus('No protein data — waiting for a PDB payload.');
    return {};
  }

  try {
    stats = parsePdbStats(pdb);
    statsEl.textContent = proteinStatsText(stats);
  } catch (e) {
    showStatus('⚠ PDB parse error: ' + errText(e), 'error');
  }

  showStatus('Loading 3D viewer…');
  load3Dmol().then(() => {
    viewer = ThreeDmol.createViewer(pane.container, { backgroundColor: T.viewerBg });
    viewer.addModel(pdb, 'pdb');
    // applyProteinStyles clears the "Loading…" status (or replaces it with the
    // surface-computing message), so there is nothing to hide here.
    applyProteinStyles(viewer, opts, stats, showStatus);
    viewer.zoomTo();
    viewer.spin(!!opts.spin);
    viewer.render();
  }).catch(e => {
    showStatus('⚠ Failed to render structure: ' + errText(e), 'error');
  });

  return {
    onResize() { if (viewer) { viewer.resize(); viewer.render(); } },
    cleanup() {
      if (!viewer) return;
      try { viewer.spin(false); } catch (e) { /* v1 quirk */ }
      try { viewer.clear(); } catch (e) { /* already gone */ }
    }
  };
}

// ============================================================================
// VIEW: SolventComponent — a compact settings card
//
// A SolventComponent has no coordinates; it is a *specification* of a solvent
// box. So: a depiction of the solvent molecule, a row of value chips, and a
// purely illustrative schematic whose ion dots scale with the concentration.
// ============================================================================

const SCHEMATIC = {
  width: 320, height: 120, padding: 14, dotRadius: 4,
  // molar -> dot count. 0.15 M ≈ 9 ion pairs, capped so it stays legible.
  dotsPerMolar: 60, maxDots: 26
};

/** "0.15 molar" / "0.15 mol/L" / 0.15 -> 0.15; null if unparseable. */
function parseConcentration(value) {
  if (value == null) return null;
  if (typeof value === 'number') return isFinite(value) ? value : null;
  const m = String(value).match(/-?\d+(\.\d+)?([eE][-+]?\d+)?/);
  const n = m ? parseFloat(m[0]) : NaN;
  return isFinite(n) ? n : null;
}

/** Deterministic PRNG so the ion scatter is stable across re-renders. */
function makeRng(seed) {
  let s = seed >>> 0 || 1;
  return () => {
    s ^= s << 13; s >>>= 0;
    s ^= s >> 17;
    s ^= s << 5; s >>>= 0;
    return s / 4294967296;
  };
}

function solventSchematicSVG(solvent) {
  const { width: W, height: H, padding: pad } = SCHEMATIC;
  const conc = parseConcentration(solvent.ion_concentration);
  const hasIons = solvent.positive_ion != null || solvent.negative_ion != null;

  const nPairs = (hasIons && conc > 0)
    ? Math.min(SCHEMATIC.maxDots, Math.round(conc * SCHEMATIC.dotsPerMolar))
    : 0;
  const dots = [];
  for (let i = 0; i < nPairs; i++) {
    if (solvent.positive_ion != null) dots.push(true);
    if (solvent.negative_ion != null) dots.push(false);
  }

  const rng = makeRng(1337);
  const parts = [
    '<svg viewBox="0 0 ' + W + ' ' + H + '" preserveAspectRatio="xMidYMid meet" ' +
    'style="width:100%;height:100%;display:block;">',
    '<rect x="1" y="1" width="' + (W - 2) + '" height="' + (H - 2) + '" rx="12" ' +
    'fill="' + T.boxFill + '" stroke="' + T.boxStroke + '" stroke-width="1.5"/>'
  ];
  dots.forEach(isCation => {
    const cx = (pad + rng() * (W - 2 * pad)).toFixed(1);
    const cy = (pad + rng() * (H - 2 * pad));
    const color = isCation ? T.chipUniqueA : T.chipUniqueB;
    parts.push('<circle cx="' + cx + '" cy="' + cy.toFixed(1) + '" r="' + SCHEMATIC.dotRadius +
      '" fill="' + color + '" fill-opacity="0.85"/>');
    parts.push('<text x="' + cx + '" y="' + (cy + 2.6).toFixed(1) + '" text-anchor="middle" ' +
      'font-size="7" font-weight="700" fill="#ffffff" pointer-events="none">' +
      (isCation ? '+' : '−') + '</text>');
  });
  if (dots.length === 0) {
    parts.push('<text x="' + (W / 2) + '" y="' + (H / 2 + 4) + '" text-anchor="middle" ' +
      'font-size="11" fill="' + T.textMuted2 + '">no ions</text>');
  }
  parts.push('</svg>');
  return parts.join('');
}

async function viewSolvent(host, inputs) {
  const solvent = await asObject(inputs.solvent_component);

  const card = el('div',
    'flex:1;min-height:0;box-sizing:border-box;padding:14px 18px;overflow:auto;' +
    'display:flex;flex-direction:column;gap:12px;background:' + T.appBg + ';');
  host.appendChild(card);

  if (!solvent || typeof solvent !== 'object') {
    card.appendChild(centredMessage('No solvent specification available.'));
    return {};
  }

  // ─── title row: depiction + SMILES ───
  const titleRow = el('div', 'display:flex;align-items:center;gap:14px;flex-shrink:0;');
  card.appendChild(titleRow);

  const depictionBox = el('div',
    'width:64px;height:64px;flex-shrink:0;border-radius:10px;overflow:hidden;' +
    'display:flex;align-items:center;justify-content:center;' +
    'border:1px solid ' + T.cardBorder + ';background:' + T.canvas2DBg + ';');
  titleRow.appendChild(depictionBox);

  const titleText = el('div', 'display:flex;flex-direction:column;gap:2px;min-width:0;');
  titleText.appendChild(el('div',
    'font-weight:700;font-size:17px;letter-spacing:.03em;color:' + T.titleColor + ';', 'Solvent'));
  const smiles = solvent.smiles == null ? null : String(solvent.smiles);
  const smilesEl = el('div',
    'font-size:12px;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;color:' + T.textMuted +
    ";font-family:ui-monospace,'SF Mono',Menlo,monospace;",
    'SMILES  ' + (smiles || EM_DASH));
  smilesEl.title = smiles || '';
  titleText.appendChild(smilesEl);
  titleRow.appendChild(titleText);

  if (smiles) {
    depictionBox.innerHTML = '<div style="color:' + T.loadingFg + ';font-size:10px;">…</div>';
    loadRDKit().then(RDKit => {
      const svg = depictSVG(RDKit, smiles, 96);
      if (svg) placeDepiction(depictionBox, svg, 96);
      else depictionBox.innerHTML = '<div style="color:' + T.errorFg + ';font-size:10px;">?</div>';
    }).catch(() => {
      depictionBox.innerHTML = '<div style="color:' + T.errorFg + ';font-size:9px;">RDKit</div>';
    });
  } else {
    depictionBox.innerHTML = '<div style="color:' + T.textMuted2 + ';font-size:14px;">' + EM_DASH + '</div>';
  }

  // ─── chips ───
  const chipRow = el('div', 'display:flex;flex-wrap:wrap;gap:8px;flex-shrink:0;');
  card.appendChild(chipRow);

  const addChip = (label, valueNode) => {
    const chip = el('div',
      'display:flex;flex-direction:column;gap:2px;padding:6px 12px;border-radius:8px;min-width:76px;' +
      'background:' + T.cardBg + ';border:1px solid ' + T.cardBorder + ';');
    chip.appendChild(el('div',
      'font-size:10px;text-transform:uppercase;letter-spacing:.06em;color:' + T.textMuted2 + ';', label));
    chip.appendChild(valueNode);
    chipRow.appendChild(chip);
  };
  const valueText = (v) => el('div',
    'font-size:14px;font-weight:600;color:' + T.textPrimary + ';',
    v == null || String(v).trim() === '' ? EM_DASH : String(v));

  addChip('Positive ion', valueText(solvent.positive_ion));
  addChip('Negative ion', valueText(solvent.negative_ion));
  addChip('Ion concentration', valueText(solvent.ion_concentration));

  const pill = el('span', '', solvent.neutralize ? 'yes' : 'no');
  pill.style.cssText =
    'display:inline-block;align-self:flex-start;padding:1px 10px;border-radius:999px;font-size:12px;font-weight:700;' +
    (solvent.neutralize
      ? 'background:' + T.okBg + ';color:' + T.okFg + ';border:1px solid ' + T.okBorder + ';'
      : 'background:' + T.warnBg + ';color:' + T.warnFg + ';border:1px solid ' + T.warnBorder + ';');
  addChip('Neutralize', pill);

  // ─── schematic ───
  const schematicWrap = el('div', 'flex:1;min-height:0;display:flex;flex-direction:column;gap:4px;');
  card.appendChild(schematicWrap);
  const svgWrap = el('div', 'flex:1;min-height:0;');
  svgWrap.innerHTML = solventSchematicSVG(solvent);
  schematicWrap.appendChild(svgWrap);

  const legend = [];
  if (solvent.positive_ion != null) legend.push(solvent.positive_ion + ' ● cation');
  if (solvent.negative_ion != null) legend.push(solvent.negative_ion + ' ● anion');
  schematicWrap.appendChild(el('div',
    'font-size:10px;font-style:italic;color:' + T.textMuted2 + ';',
    'Illustrative only — not a real solvent configuration. Ion dots are drawn in proportion to the ' +
    'requested concentration' + (legend.length ? ' (' + legend.join(', ') + ').' : '.')));

  // The card is pure CSS flex plus a viewBox SVG, so it reflows on its own.
  return {};
}

// ============================================================================
// VIEW: ChemicalSystem — component list (left) + detail pane (right)
// ============================================================================

async function viewChemicalSystem(host, inputs) {
  const system = await asObject(inputs.chemical_system);
  const components = (system && Array.isArray(system.components))
    ? system.components.filter(c => c && typeof c === 'object')
    : [];

  const summary = () => {
    if (!components.length) return 'no components';
    const counts = new Map();
    components.forEach(c => {
      const w = typeWord(c.type || 'Component');
      counts.set(w, (counts.get(w) || 0) + 1);
    });
    const parts = [];
    counts.forEach((n, w) => parts.push(n + ' ' + w + (n > 1 ? 's' : '')));
    return components.length + ' component' + (components.length > 1 ? 's' : '') + ' — ' + parts.join(', ');
  };

  host.appendChild(headerStrip((system && system.name) || 'Chemical System', summary()));

  const split = el('div', 'flex:1;display:flex;flex-direction:row;min-height:0;overflow:hidden;');
  host.appendChild(split);

  const listPane = el('div',
    'width:240px;flex:0 0 240px;overflow-y:auto;padding:10px;box-sizing:border-box;' +
    'display:flex;flex-direction:column;gap:8px;background:' + T.appBg + ';');
  split.appendChild(listPane);
  split.appendChild(el('div', 'width:1px;flex-shrink:0;background:' + T.splitBorder + ';'));
  const detailPane = el('div',
    'flex:1;min-width:0;min-height:0;display:flex;flex-direction:column;background:' + T.appBg + ';');
  split.appendChild(detailPane);

  let viewer = null;
  let detailToken = 0;      // bumped per selection; guards async 3Dmol loads
  let selected = -1;

  const destroyViewer = () => {
    if (!viewer) return;
    try { viewer.clear(); } catch (e) { /* already gone */ }
    viewer = null;
  };

  const cardCss = (active) =>
    'display:flex;flex-direction:column;gap:4px;padding:8px 10px;border-radius:8px;cursor:pointer;flex-shrink:0;' +
    'transition:background .12s;border:1px solid ' + (active ? T.cardBorderActive : T.cardBorder) +
    ';background:' + (active ? T.cardBgActive : T.cardBg) + ';';

  const cards = components.map((comp, i) => {
    const card = el('div', cardCss(false));
    card.appendChild(el('div',
      'font-size:13px;font-weight:700;word-break:break-word;color:' + T.textPrimary + ';',
      comp.label || '(unlabelled)'));

    const badgeRow = el('div', 'display:flex;align-items:center;gap:6px;flex-wrap:wrap;');
    badgeRow.appendChild(typeBadge(comp.type || 'Component', true));
    if (comp.error) badgeRow.appendChild(el('span', 'font-size:11px;color:' + T.errorFg + ';', '⚠'));
    card.appendChild(badgeRow);

    if (comp.name) {
      card.appendChild(el('div',
        'font-size:11px;word-break:break-word;color:' + T.textMuted + ';', comp.name));
    }
    card.onmouseover = () => { if (i !== selected) card.style.background = T.cardBgHover; };
    card.onmouseout = () => { if (i !== selected) card.style.background = T.cardBg; };
    card.onclick = () => select(i);
    listPane.appendChild(card);
    return card;
  });

  if (!components.length) {
    listPane.appendChild(el('div', 'font-size:12px;padding:6px;color:' + T.textMuted + ';', 'No components.'));
  }

  /** A strip pinned under the detail viewer — SMILES or PDB stats. */
  const detailStrip = (label, value) => {
    const wrap = el('div',
      'flex-shrink:0;display:flex;gap:8px;align-items:baseline;padding:8px 16px;font-size:12px;' +
      'background:' + T.toolbarBg + ';border-top:1px solid ' + T.toolbarBorder + ';');
    if (label) {
      wrap.appendChild(el('span',
        'font-size:11px;text-transform:uppercase;letter-spacing:.04em;color:' + T.textMuted2 + ';', label));
    }
    wrap.appendChild(el('div',
      'flex:1;min-width:0;overflow-x:auto;white-space:nowrap;color:' + (label ? T.textPrimary : T.textMuted2) +
      (label ? ";font-family:ui-monospace,SFMono-Regular,Menlo,monospace;" : ''),
      value));
    return wrap;
  };

  /** Mount a 3Dmol viewer for the current selection, guarded by detailToken. */
  function mount3D(styleFn) {
    const pane = viewerHost();
    detailPane.appendChild(pane.wrap);
    const myToken = detailToken;
    load3Dmol().then(() => {
      if (myToken !== detailToken) return;   // a different component was picked
      viewer = ThreeDmol.createViewer(pane.container, { backgroundColor: T.viewerBg });
      styleFn(viewer);
      viewer.zoomTo();
      viewer.render();
    }).catch(e => {
      if (myToken !== detailToken) return;
      destroyViewer();
      pane.container.appendChild(warnBanner('Failed to render: ' + errText(e)));
    });
    return pane;
  }

  function renderDetail(comp) {
    destroyViewer();
    detailPane.innerHTML = '';
    detailToken++;

    if (!comp) {
      detailPane.appendChild(centredMessage('Select a component.'));
      return;
    }

    const header = el('div',
      'display:flex;align-items:baseline;gap:10px;padding:9px 16px;flex-wrap:wrap;flex-shrink:0;' +
      'background:' + T.toolbarBg + ';border-bottom:1px solid ' + T.toolbarBorder + ';');
    header.appendChild(el('span',
      'font-size:14px;font-weight:700;color:' + T.textPrimary + ';', comp.label || '(unlabelled)'));
    header.appendChild(typeBadge(comp.type || 'Component', false));
    if (comp.name) header.appendChild(el('span', 'font-size:12px;color:' + T.textMuted + ';', comp.name));
    detailPane.appendChild(header);

    if (comp.error) detailPane.appendChild(warnBanner(comp.error));

    if (comp.sdf) {
      mount3D(v => {
        v.addModel(ensureSDFTerminator(comp.sdf), 'sdf');
        v.setStyle({}, {
          stick: { radius: 0.15, colorscheme: 'Jmol' },
          sphere: { scale: 0.25, colorscheme: 'Jmol' }
        });
      });
      detailPane.appendChild(detailStrip('SMILES', comp.smiles || EM_DASH));
    } else if (comp.pdb) {
      let stats = null;
      try { stats = parsePdbStats(comp.pdb); } catch (e) { /* show it without stats */ }
      mount3D(v => {
        v.addModel(comp.pdb, 'pdb');
        applyProteinStyles(v, { rep: 'cartoon', color: 'chain', waters: false, hetero: true }, stats);
      });
      if (stats) detailPane.appendChild(detailStrip(null, proteinStatsText(stats) + ' — waters hidden'));
    } else if (comp.smiles != null || comp.ion_concentration != null || comp.neutralize != null) {
      const scroll = el('div', 'flex:1;overflow-y:auto;padding:16px;min-height:0;');
      scroll.appendChild(solventCard(comp));
      detailPane.appendChild(scroll);
    } else if (!comp.error) {
      // Nothing structural made it across — show whatever scalars did.
      const scroll = el('div', 'flex:1;overflow-y:auto;padding:16px;min-height:0;');
      const card = el('div',
        'max-width:520px;padding:14px 16px;border-radius:10px;border:1px solid ' + T.cardBorder +
        ';background:' + T.cardBg + ';');
      card.appendChild(el('div', 'font-size:12px;margin-bottom:8px;color:' + T.textMuted + ';',
        'No structural data was serialized for this component.'));
      Object.keys(comp).forEach(k => {
        if (k === 'label' || k === 'type' || k === 'error') return;
        if (comp[k] == null || typeof comp[k] === 'object') return;
        card.appendChild(fieldRow(k.replace(/_/g, ' '), comp[k]));
      });
      scroll.appendChild(card);
      detailPane.appendChild(scroll);
    } else {
      detailPane.appendChild(el('div', 'flex:1;'));   // the error banner alone
    }
  }

  function select(i) {
    if (i < 0 || i >= components.length) return;
    selected = i;
    cards.forEach((c, j) => { c.style.cssText = cardCss(j === selected); });
    renderDetail(components[i]);
  }

  if (components.length) select(0);
  else renderDetail(null);

  return {
    onResize() { if (viewer) { viewer.resize(); viewer.render(); } },
    cleanup: destroyViewer
  };
}

// ============================================================================
// VIEW: Transformation — component diff on top, mapping viewer below
// ============================================================================

const DIFF_STATUS = {
  unchanged: { label: 'unchanged', color: () => T.diffUnchanged },
  changed: { label: 'changed', color: () => T.diffChanged },
  added: { label: 'added', color: () => T.diffAdded },
  removed: { label: 'removed', color: () => T.diffRemoved }
};

function diffStatus(a, b) {
  if (a && !b) return 'removed';
  if (b && !a) return 'added';
  if (a.type !== b.type) return 'changed';
  if ((a.name || '') !== (b.name || '')) return 'changed';
  return 'unchanged';
}

/** One-line summary of a component descriptor: what it actually is. */
function componentDetail(c) {
  if (!c) return '';
  if (c.error) return 'error: ' + c.error;
  if (c.smiles && c.sdf) return c.smiles;
  if (c.pdb) return c.pdb.split(NL).length + ' PDB lines';
  if (c.smiles) {
    const bits = [c.smiles];
    if (c.ion_concentration) bits.push(c.ion_concentration);
    if (c.positive_ion && c.negative_ion) bits.push(c.positive_ion + '/' + c.negative_ion);
    return bits.join('  ·  ');
  }
  return '';
}

/** One side of a diff row. `null` renders as an explicit empty slot. */
function componentCell(c, side) {
  const cell = el('div',
    'flex:1 1 0;min-width:0;padding:7px 10px;border-radius:6px;border:1px solid ' + T.cardBorder +
    ';background:' + T.cardBg + ';');
  if (!c) {
    cell.style.borderStyle = 'dashed';
    cell.style.opacity = '0.5';
    cell.appendChild(el('span', 'font-size:11px;color:' + T.textMuted2 + ';',
      '— absent in state' + side + ' —'));
    return cell;
  }
  const ellipsis = 'overflow:hidden;text-overflow:ellipsis;white-space:nowrap;';
  const nameRow = el('div',
    'font-size:12px;font-weight:600;color:' + T.textPrimary + ';' + ellipsis, c.name || '(unnamed)');
  nameRow.title = c.name || '';
  cell.appendChild(nameRow);
  cell.appendChild(el('div', 'font-size:10px;margin-top:2px;color:' + T.textMuted2 + ';', c.type || ''));

  const detail = componentDetail(c);
  if (detail) {
    const row = el('div',
      'font-size:10px;font-family:ui-monospace,monospace;margin-top:3px;color:' + T.textMuted + ';' + ellipsis,
      detail);
    row.title = detail;
    cell.appendChild(row);
  }
  return cell;
}

function renderDiff(section, stateA, stateB) {
  const compsA = {}, compsB = {};
  (stateA.components || []).forEach(c => { compsA[c.label] = c; });
  (stateB.components || []).forEach(c => { compsB[c.label] = c; });
  const labels = Array.from(new Set(Object.keys(compsA).concat(Object.keys(compsB)))).sort();

  const heads = el('div',
    'display:flex;gap:8px;align-items:center;margin-bottom:6px;font-size:11px;font-weight:700;' +
    'letter-spacing:.04em;text-transform:uppercase;color:' + T.textMuted + ';');
  heads.appendChild(el('div', 'flex:0 0 120px;', 'label'));
  [stateA.name || 'state A', stateB.name || 'state B'].forEach(n => {
    heads.appendChild(el('div',
      'flex:1 1 0;min-width:0;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;', n));
  });
  section.appendChild(heads);

  const legend = el('div',
    'display:flex;align-items:center;gap:14px;flex-wrap:wrap;margin-bottom:10px;font-size:11px;color:' +
    T.textMuted + ';');
  Object.keys(DIFF_STATUS).forEach(k => {
    const item = el('span', 'display:inline-flex;align-items:center;gap:5px;');
    item.appendChild(el('span',
      'width:9px;height:9px;border-radius:2px;display:inline-block;background:' + DIFF_STATUS[k].color() + ';'));
    item.appendChild(el('span', '', DIFF_STATUS[k].label));
    legend.appendChild(item);
  });
  section.appendChild(legend);

  if (!labels.length) {
    section.appendChild(el('div', 'font-size:12px;padding:8px 0;color:' + T.textMuted2 + ';',
      'Neither chemical system declares any components.'));
    return;
  }

  labels.forEach(label => {
    const a = compsA[label] || null;
    const b = compsB[label] || null;
    const status = diffStatus(a, b);
    const color = DIFF_STATUS[status].color();

    const row = el('div', 'display:flex;gap:8px;align-items:stretch;margin-bottom:6px;');
    const labelCell = el('div',
      'flex:0 0 120px;min-width:0;display:flex;flex-direction:column;justify-content:center;gap:3px;' +
      'padding-left:8px;border-left:3px solid ' + color + ';');
    labelCell.appendChild(el('div',
      'font-size:12px;font-weight:600;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;color:' +
      T.textPrimary + ';', label));
    labelCell.appendChild(el('div',
      'font-size:10px;font-weight:600;color:' + color + ';', DIFF_STATUS[status].label));
    row.appendChild(labelCell);
    row.appendChild(componentCell(a, 'A'));
    row.appendChild(componentCell(b, 'B'));
    section.appendChild(row);
  });
}

async function viewTransformation(host, inputs) {
  const tf = await asObject(inputs.transformation);
  if (!tf || typeof tf !== 'object') {
    host.appendChild(centredMessage('Payload is missing the `transformation` object.'));
    return {};
  }

  const stateA = tf.stateA || { name: 'state A', components: [] };
  const stateB = tf.stateB || { name: 'state B', components: [] };

  const header = headerStrip(
    (stateA.name || 'state A') + '  →  ' + (stateB.name || 'state B'),
    tf.name || '(unnamed transformation)');
  header.statsEl.appendChild(el('span',
    'font-size:11px;font-weight:600;padding:3px 9px;border-radius:10px;white-space:nowrap;' +
    'background:' + T.badgeBg + ';color:' + T.badgeFg + ';',
    tf.protocol || 'unknown protocol'));
  host.appendChild(header);

  const diffSection = el('div',
    'flex:0 0 auto;max-height:55%;overflow:auto;padding:12px 14px;border-bottom:1px solid ' + T.splitBorder + ';');
  host.appendChild(diffSection);
  renderDiff(diffSection, stateA, stateB);

  const mappingSection = el('div',
    'flex:1 1 auto;min-height:0;display:flex;flex-direction:column;position:relative;');
  host.appendChild(mappingSection);

  const mappingBar = el('div',
    'display:flex;align-items:center;gap:12px;flex-wrap:wrap;padding:6px 14px;flex-shrink:0;font-size:11px;' +
    'background:' + T.labelBg + ';color:' + T.textMuted + ';');
  mappingSection.appendChild(mappingBar);

  const viewerBox = el('div', 'flex:1;min-height:0;display:flex;flex-direction:column;');
  mappingSection.appendChild(viewerBox);

  const mappings = Array.isArray(tf.mappings) ? tf.mappings : [];
  const viewer = createMappingViewer(viewerBox, 'colored');

  if (!mappings.length) {
    mappingBar.textContent = 'atom mapping';
    viewer.message('No atom mapping on this transformation — nothing to superimpose.');
    return { onResize: viewer.resize, cleanup: viewer.destroy };
  }

  const stats = el('div', 'display:flex;align-items:center;gap:12px;flex-wrap:wrap;margin-left:auto;');

  function show(payload) {
    let molA, molB;
    try {
      molA = parseSDF(payload['molA.sdf'], payload.nameA || 'molecule A');
      molB = parseSDF(payload['molB.sdf'], payload.nameB || 'molecule B');
    } catch (e) {
      viewer.message('SDF parse error: ' + errText(e), true);
      return;
    }
    if (payload.nameA) molA.name = payload.nameA;
    if (payload.nameB) molB.name = payload.nameB;

    const mapping = viewer.show(molA, molB, payload.mapping);
    const coreA = mapping ? new Set(mapping.keys()).size : 0;
    const coreB = mapping ? new Set(mapping.values()).size : 0;
    stats.innerHTML = '';
    stats.appendChild(statChip('mapped', (mapping ? mapping.size : 0) + ' atoms', T.chipCore));
    stats.appendChild(statChip('unique A', molA.symbols.length - coreA, T.chipUniqueA));
    stats.appendChild(statChip('unique B', molB.symbols.length - coreB, T.chipUniqueB));
  }

  if (mappings.length > 1) {
    const select = el('select', SELECT_CSS + 'padding:3px 8px;font-size:11px;');
    mappings.forEach((m, i) => {
      const o = el('option', '', 'mapping ' + (i + 1) + ': ' + (m.nameA || 'A') + ' → ' + (m.nameB || 'B'));
      o.value = String(i);
      select.appendChild(o);
    });
    select.addEventListener('change', () => show(mappings[parseInt(select.value, 10)]));
    mappingBar.appendChild(select);
  }
  mappingBar.appendChild(stats);
  show(mappings[0]);

  return { onResize: viewer.resize, cleanup: viewer.destroy };
}

// ============================================================================
// VIEW: LigandNetwork — radial graph (left) driving a mapping viewer (right)
//
// The payload is the exact `LigandNetwork.to_graphml()` output; the molecules
// and the per-edge atom mappings are parsed straight out of it.
// ============================================================================

const ELEMENT_MAP = { 1: 'H', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 15: 'P', 16: 'S', 17: 'Cl', 35: 'Br', 53: 'I' };

const NETWORK_CONFIG = {
  node: { radius: 38, strokeWidth: 1.5, labelFontSize: 11, labelMaxChars: 14, initialsSize: 18, depictionPadding: 4 },
  edge: {
    minWidth: 1.5, maxWidth: 6.5, opacity: 0.9, labelFontSize: 10, labelBgPadding: 3, labelBgOpacity: 0.92,
    arrowWidthPx: 8, arrowHeightPx: 8, arrowRefX: 46, hitWidth: 14, haloPadding: 4, haloOpacity: 0.95
  },
  force: {
    linkBaseDistance: 18, linkScoreBonus: 10, linkStrength: 0.5,
    chargeStrength: -2500, chargeDistanceMin: 20, chargeDistanceMax: 5000,
    centerStrength: 0.08, collisionPadding: 12, collisionIterations: 4,
    driftX: 0.04, driftY: 0.04, tickMultiplier: 2
  },
  radial: { ringStep: 0.18, ringOffset: 40 },
  circular: { radiusFraction: 0.36 }
};

/** Decode a gufe GraphML conformer: a base-1-per-char .npy blob of float64s. */
function parseNpyCoords(confStr) {
  const buf = new ArrayBuffer(confStr.length);
  const bytes = new Uint8Array(buf);
  for (let i = 0; i < confStr.length; i++) bytes[i] = confStr.charCodeAt(i);
  const headerLen = new DataView(buf).getUint16(8, true);
  const header = String.fromCharCode.apply(null, new Uint8Array(buf, 10, headerLen));
  const nAtoms = parseInt(header.match(/'shape':[ ]*[(]([0-9]+),[ ]*([0-9]+)[)]/)[1], 10);
  const dv = new DataView(buf, 10 + headerLen);
  const coords = [];
  for (let i = 0; i < nAtoms; i++) {
    coords.push([dv.getFloat64(i * 24, true), dv.getFloat64(i * 24 + 8, true), dv.getFloat64(i * 24 + 16, true)]);
  }
  return coords;
}

/**
 * Rotate coordinates onto their principal axes, so the flattest projection
 * faces the viewer — what makes the node depictions legible.
 */
function alignFlat(coords) {
  const n = coords.length;
  const c = [0, 0, 0];
  coords.forEach(p => { for (let k = 0; k < 3; k++) c[k] += p[k]; });
  for (let k = 0; k < 3; k++) c[k] /= n;
  const pts = coords.map(p => [p[0] - c[0], p[1] - c[1], p[2] - c[2]]);

  const cov = [0, 0, 0, 0, 0, 0, 0, 0, 0];
  pts.forEach(p => {
    for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) cov[i * 3 + j] += p[i] * p[j];
  });
  const R = jacobiSym3(cov).vectors;   // columns = principal axes, descending
  return pts.map(p => [
    R[0] * p[0] + R[3] * p[1] + R[6] * p[2],
    R[1] * p[0] + R[4] * p[1] + R[7] * p[2],
    R[2] * p[0] + R[5] * p[1] + R[8] * p[2]
  ]);
}

/** Parse a GraphML document into `{ nodes, edges, molecules }`. */
function parseGraphML(xmlStr) {
  const doc = new DOMParser().parseFromString(xmlStr, 'text/xml');
  const parseErr = doc.getElementsByTagName('parsererror');
  if (parseErr.length > 0) throw new Error(parseErr[0].textContent);

  const keyNames = {};
  for (const k of doc.getElementsByTagName('key')) {
    keyNames[k.getAttribute('id')] = k.getAttribute('attr.name') || k.getAttribute('id');
  }

  const molecules = {};
  const nodes = [];
  for (const nodeEl of doc.getElementsByTagName('node')) {
    const id = nodeEl.getAttribute('id');
    try {
      const dataEl = Array.from(nodeEl.getElementsByTagName('data'))
        .find(d => d.getAttribute('key') === 'd0');
      if (!dataEl) throw new Error('no d0 data element');
      const json = JSON.parse(dataEl.textContent);
      const mol = {
        name: (json.molprops && json.molprops['ofe-name']) || id,
        symbols: json.atoms.map(a => ELEMENT_MAP[a[0]] || 'X'),
        bonds: json.bonds.map(b => [b[0], b[1], b[2]]),
        coords: alignFlat(parseNpyCoords(json.conformer[0]))
      };
      molecules[id] = mol;
      nodes.push({ id, name: mol.name });
    } catch (e) {
      if (typeof logStderr === 'function') logStderr('Failed to parse node ' + id + ': ' + errText(e));
    }
  }

  const edges = [];
  for (const edgeEl of doc.getElementsByTagName('edge')) {
    const source = edgeEl.getAttribute('source');
    const target = edgeEl.getAttribute('target');
    if (!molecules[source] || !molecules[target]) continue;

    let score = 0.5;
    let mapping = null;
    for (const d of edgeEl.getElementsByTagName('data')) {
      const keyName = (keyNames[d.getAttribute('key')] || '').toLowerCase();
      const text = d.textContent.trim();
      if (!text) continue;
      if (keyName.includes('score') || keyName.includes('lomap') || keyName === 'weight') {
        const v = parseFloat(text);
        if (!isNaN(v)) score = v;
      }
      if (mapping) continue;
      try { mapping = parseMappingValue(JSON.parse(text)); } catch (e) { /* not JSON */ }
    }
    edges.push({ source, target, score, mapping });
  }

  return { nodes, edges, molecules };
}

async function viewLigandNetwork(host, inputs) {
  const xml = await asText(inputs['network.graphml']);
  const split = el('div', 'flex:1;display:flex;flex-direction:row;min-height:0;overflow:hidden;');
  host.appendChild(split);

  const leftPane = el('div',
    'flex:1 1 50%;min-width:0;display:flex;flex-direction:column;background:' + T.netCanvasBg + ';');
  const rightPane = el('div',
    'flex:1 1 50%;min-width:0;display:flex;flex-direction:column;background:' + T.appBg + ';');
  split.appendChild(leftPane);
  split.appendChild(el('div', 'width:1px;flex-shrink:0;background:' + T.splitBorder + ';'));
  split.appendChild(rightPane);

  const mappingViewer = createMappingViewer(rightPane, 'plain');

  const svgWrap = el('div',
    'flex:1;position:relative;overflow:hidden;min-height:0;background:' + T.netCanvasBg + ';');
  leftPane.appendChild(svgWrap);

  const toolbar = el('div',
    'display:flex;align-items:center;gap:12px;flex-wrap:wrap;padding:10px 16px;flex-shrink:0;' +
    'background:' + T.toolbarBg + ';border-top:1px solid ' + T.toolbarBorder + ';');
  leftPane.appendChild(toolbar);

  const layoutSelect = el('select', SELECT_CSS);
  ['Force-directed', 'Circular', 'Radial Tree'].forEach(l => {
    const o = el('option', '', l);
    o.value = l;
    layoutSelect.appendChild(o);
  });
  toolbar.appendChild(el('label', 'font-size:12px;margin-left:auto;color:' + T.textMuted + ';', 'Layout:'));
  toolbar.appendChild(layoutSelect);

  const legend = el('div', 'display:flex;align-items:center;gap:6px;font-size:11px;color:' + T.textMuted + ';');
  legend.innerHTML =
    '<span>LOMAP score:</span><span style="display:flex;align-items:center;gap:4px;">' +
    '<span style="width:32px;height:4px;border-radius:2px;display:inline-block;' +
    'background:linear-gradient(to right,' + T.netEdgeRamp.join(',') + ');"></span><span>0 → 1</span></span>';
  toolbar.appendChild(legend);

  const tooltip = el('div',
    'position:absolute;pointer-events:none;opacity:0;max-width:220px;padding:10px 14px;border-radius:8px;' +
    'font-size:12px;z-index:10;transition:opacity .15s;box-shadow:0 4px 16px rgba(15,23,42,0.12);' +
    'background:' + T.tooltipBg + ';border:1px solid ' + T.tooltipBorder + ';color:' + T.textPrimary + ';');
  svgWrap.appendChild(tooltip);

  let graph;
  try {
    graph = parseGraphML(xml);
  } catch (e) {
    if (typeof logStderr === 'function') logStderr('GraphML parse failed: ' + errText(e));
    floatingWarning(svgWrap, 'GraphML parse error: ' + errText(e));
    mappingViewer.message('No network to show.');
    return { onResize: mappingViewer.resize, cleanup: mappingViewer.destroy };
  }

  // Depictions are RDKit-only, so start them without waiting on anything else.
  const rdkitReady = loadRDKit().catch(err => {
    if (typeof logStderr === 'function') logStderr('RDKit load failed: ' + errText(err));
    return null;
  });

  const d3 = await loadD3();
  const scoreColor = d3.scaleSequential().domain([0, 1])
    .interpolator(d3.interpolateRgbBasis(T.netEdgeRamp));

  let simulation = null;
  let svgEl = null;
  let activeEdgeIdx = -1;

  function selectEdge(idx) {
    const edge = graph.edges[idx];
    if (!edge) return;
    mappingViewer.show(graph.molecules[edge.source], graph.molecules[edge.target], edge.mapping);
  }

  function render(layout) {
    if (svgEl) svgEl.remove();
    if (simulation) { simulation.stop(); simulation = null; }

    const W = svgWrap.getBoundingClientRect().width || 800;
    const H = svgWrap.getBoundingClientRect().height || 600;
    const CX = W / 2, CY = H / 2;
    const NR = NETWORK_CONFIG.node.radius;
    const EC = NETWORK_CONFIG.edge;

    svgEl = d3.select(svgWrap).append('svg').attr('width', W).attr('height', H).style('display', 'block');
    const gRoot = svgEl.append('g');
    svgEl.call(d3.zoom().scaleExtent([0.15, 5]).on('zoom', e => gRoot.attr('transform', e.transform)));

    // Working copies (d3 mutates node objects with x/y/vx/vy).
    const nodes = graph.nodes.map(n => ({ ...n }));
    const nodeById = Object.fromEntries(nodes.map(n => [n.id, n]));
    const edges = graph.edges.map((e, i) => ({
      ...e, _idx: i, _color: scoreColor(e.score),
      source: nodeById[e.source], target: nodeById[e.target]
    }));

    // One arrow marker per distinct edge colour (SVG markers cannot inherit it).
    const defs = svgEl.append('defs');
    const markerIds = new Map();
    edges.forEach(e => {
      if (markerIds.has(e._color)) return;
      const id = 'arrow-' + e._color.replace(/[^a-zA-Z0-9]/g, '');
      markerIds.set(e._color, id);
      defs.append('marker')
        .attr('id', id).attr('viewBox', '0 -5 10 10')
        .attr('refX', EC.arrowRefX).attr('refY', 0)
        .attr('markerUnits', 'userSpaceOnUse')
        .attr('markerWidth', EC.arrowWidthPx).attr('markerHeight', EC.arrowHeightPx)
        .attr('orient', 'auto')
        .append('path').attr('d', 'M0,-5L10,0L0,5').attr('fill', e._color);
    });

    // ─── layouts ───
    if (layout === 'Circular') {
      const r = Math.min(W, H) * NETWORK_CONFIG.circular.radiusFraction;
      nodes.forEach((n, i) => {
        const a = (2 * Math.PI * i / nodes.length) - Math.PI / 2;
        n.x = n.fx = CX + r * Math.cos(a);
        n.y = n.fy = CY + r * Math.sin(a);
      });
    } else if (layout === 'Radial Tree') {
      const adj = Object.fromEntries(nodes.map(n => [n.id, []]));
      edges.forEach(e => { adj[e.source.id].push(e.target.id); adj[e.target.id].push(e.source.id); });
      const visited = new Set([nodes[0].id]);
      const levels = [[nodes[0].id]];
      for (;;) {
        const next = [];
        levels[levels.length - 1].forEach(id =>
          adj[id].forEach(nb => { if (!visited.has(nb)) { visited.add(nb); next.push(nb); } }));
        if (!next.length) break;
        levels.push(next);
      }
      const rStep = Math.min(W, H) * NETWORK_CONFIG.radial.ringStep;
      levels.forEach((level, li) => {
        const r = li === 0 ? 0 : li * rStep + NETWORK_CONFIG.radial.ringOffset;
        level.forEach((id, i) => {
          const a = (2 * Math.PI * i / level.length) - Math.PI / 2;
          const nd = nodeById[id];
          nd.x = nd.fx = CX + r * Math.cos(a);
          nd.y = nd.fy = CY + r * Math.sin(a);
        });
      });
      nodes.filter(n => !visited.has(n.id)).forEach((n, i) => {
        n.x = n.fx = 40 + i * 80;
        n.y = n.fy = 40;
      });
    }

    const edgeWidth = (d) => EC.minWidth + d.score * (EC.maxWidth - EC.minWidth);

    const linkHit = gRoot.append('g').selectAll('line').data(edges).join('line')
      .attr('stroke', 'transparent').attr('stroke-width', EC.hitWidth).style('cursor', 'pointer');

    const linkHalo = gRoot.append('g').selectAll('line').data(edges).join('line')
      .attr('stroke', T.netHaloColor)
      .attr('stroke-width', d => edgeWidth(d) + EC.haloPadding * 2)
      .attr('stroke-linecap', 'round').attr('opacity', 0).style('pointer-events', 'none');

    const link = gRoot.append('g').selectAll('line').data(edges).join('line')
      .attr('stroke', d => d._color).attr('stroke-width', edgeWidth)
      .attr('stroke-opacity', EC.opacity)
      .attr('marker-end', d => 'url(#' + markerIds.get(d._color) + ')')
      .style('pointer-events', 'none');

    // ─── edge score labels ───
    const linkLabelG = gRoot.append('g').selectAll('g').data(edges).join('g').attr('pointer-events', 'none');
    linkLabelG.append('rect')
      .attr('fill', T.netLabelBg).attr('opacity', EC.labelBgOpacity).attr('rx', 3).attr('ry', 3);
    linkLabelG.append('text')
      .attr('text-anchor', 'middle').attr('dominant-baseline', 'middle')
      .attr('font-size', EC.labelFontSize + 'px').attr('font-weight', '600')
      .attr('fill', T.netEdgeLabel)
      .text(d => d.score.toFixed(2))
      .each(function () {
        const bb = this.getBBox();
        const bg = this.previousSibling;
        if (!bg || !bg.setAttribute) return;
        bg.setAttribute('x', bb.x - EC.labelBgPadding);
        bg.setAttribute('y', bb.y - EC.labelBgPadding);
        bg.setAttribute('width', bb.width + 2 * EC.labelBgPadding);
        bg.setAttribute('height', bb.height + 2 * EC.labelBgPadding);
      });

    const moveTooltip = (ev) => {
      tooltip.style.left = (ev.offsetX + 14) + 'px';
      tooltip.style.top = (ev.offsetY - 10) + 'px';
    };
    const hideTooltip = () => { tooltip.style.opacity = '0'; };

    linkHit
      .on('mouseenter', (ev, d) => {
        tooltip.style.opacity = '1';
        tooltip.innerHTML =
          '<div style="font-weight:700;margin-bottom:4px;color:' + T.titleColor + ';">Transformation</div>' +
          '<div>' + esc(d.source.name) + ' → ' + esc(d.target.name) + '</div>' +
          '<div style="margin-top:6px;color:' + T.titleColor + ';">LOMAP score: <b>' + d.score.toFixed(3) + '</b></div>' +
          '<div style="margin-top:4px;font-size:10px;color:' + T.textMuted2 + ';">Click to select</div>';
      })
      .on('mousemove', moveTooltip)
      .on('mouseleave', hideTooltip)
      .on('click', (ev, d) => {
        ev.stopPropagation();
        hideTooltip();
        activeEdgeIdx = d._idx;
        refreshSelection();
        selectEdge(activeEdgeIdx);
      });

    // ─── nodes ───
    const nodeG = gRoot.append('g').selectAll('g').data(nodes).join('g')
      .style('cursor', 'grab')
      .call(d3.drag()
        .on('start', (ev, d) => { d.fx = d.x; d.fy = d.y; })
        .on('drag', (ev, d) => { d.fx = d.x = ev.x; d.fy = d.y = ev.y; ticked(); })
        .on('end', (ev, d) => { d.fx = ev.x; d.fy = ev.y; }))
      .on('mouseenter', (ev, d) => {
        tooltip.style.opacity = '1';
        tooltip.innerHTML =
          '<div style="font-weight:700;color:' + T.titleColor + ';">' + esc(d.name) + '</div>' +
          '<div style="font-size:10px;margin-top:3px;color:' + T.textMuted2 + ';">id: ' + esc(d.id) + '</div>';
      })
      .on('mousemove', moveTooltip)
      .on('mouseleave', hideTooltip)
      .on('click', ev => ev.stopPropagation());

    nodeG.append('circle')
      .attr('r', NR).attr('fill', T.netNodeFill)
      .attr('stroke', T.netNodeStroke).attr('stroke-width', NETWORK_CONFIG.node.strokeWidth);

    const depictionGroups = nodeG.append('g');
    const initials = nodeG.append('text')
      .attr('text-anchor', 'middle').attr('dominant-baseline', 'middle')
      .attr('font-size', NETWORK_CONFIG.node.initialsSize + 'px').attr('font-weight', '700')
      .attr('fill', T.netInitials).attr('pointer-events', 'none')
      .text(d => (d.name || '?').slice(0, 2).toUpperCase());

    nodeG.append('text')
      .attr('text-anchor', 'middle').attr('y', NR + 14)
      .attr('font-size', NETWORK_CONFIG.node.labelFontSize + 'px').attr('font-weight', '600')
      .attr('fill', T.netNodeLabel).attr('pointer-events', 'none')
      .text(d => d.name.length > NETWORK_CONFIG.node.labelMaxChars
        ? d.name.slice(0, NETWORK_CONFIG.node.labelMaxChars - 1) + '…'
        : d.name);

    // Inject the RDKit depictions into the node circles once RDKit is ready.
    rdkitReady.then(RDKit => {
      if (!RDKit) return;
      const RDKIT_SIZE = 200;
      const scale = ((NR - NETWORK_CONFIG.node.depictionPadding) * 2) / RDKIT_SIZE;
      const parser = new DOMParser();
      let injected = 0;

      depictionGroups.each(function (d) {
        const mol = graph.molecules[d.id];
        if (!mol) return;
        const svgStr = depictSVG(RDKit, buildMolBlock(mol), RDKIT_SIZE);
        if (!svgStr) return;
        const rdSvg = parser.parseFromString(svgStr, 'image/svg+xml').documentElement;
        if (!rdSvg || rdSvg.nodeName.toLowerCase() === 'parsererror') return;

        this.setAttribute('transform',
          'translate(' + (-scale * RDKIT_SIZE / 2) + ',' + (-scale * RDKIT_SIZE / 2) + ') scale(' + scale + ')');
        let appended = 0;
        for (const k of Array.from(rdSvg.childNodes)) {
          if (k.nodeType !== 1) continue;
          const tag = k.nodeName.toLowerCase();
          if (tag === 'defs' || tag === 'metadata' || tag === 'title') continue;
          // Drop RDKit's opaque white backing rect so the node fill shows through.
          if (tag === 'rect') {
            const fill = (k.getAttribute('fill') || '').toLowerCase();
            if (fill === '#ffffff' || fill === 'white' || fill === 'rgb(255,255,255)') continue;
          }
          this.appendChild(document.importNode(k, true));
          appended++;
        }
        if (appended > 0) injected++;
      });

      if (injected > 0) initials.style('display', 'none');
    });

    function ticked() {
      [link, linkHalo, linkHit].forEach(sel => sel
        .attr('x1', d => d.source.x).attr('y1', d => d.source.y)
        .attr('x2', d => d.target.x).attr('y2', d => d.target.y));
      linkLabelG.attr('transform', d =>
        'translate(' + (d.source.x + d.target.x) / 2 + ',' + ((d.source.y + d.target.y) / 2 - 6) + ')');
      nodeG.attr('transform', d => 'translate(' + (d.x ?? CX) + ',' + (d.y ?? CY) + ')');
    }

    function refreshSelection() {
      linkHalo.attr('opacity', d => d._idx === activeEdgeIdx ? EC.haloOpacity : 0);
    }

    if (layout === 'Force-directed') {
      const fc = NETWORK_CONFIG.force;
      simulation = d3.forceSimulation(nodes)
        .force('link', d3.forceLink(edges).id(d => d.id)
          .distance(d => fc.linkBaseDistance + (1 - d.score) * fc.linkScoreBonus)
          .strength(fc.linkStrength))
        .force('charge', d3.forceManyBody().strength(fc.chargeStrength)
          .distanceMin(fc.chargeDistanceMin).distanceMax(fc.chargeDistanceMax))
        .force('center', d3.forceCenter(CX, CY).strength(fc.centerStrength))
        .force('collision', d3.forceCollide(NR + fc.collisionPadding).strength(1).iterations(fc.collisionIterations))
        .force('x', d3.forceX(CX).strength(fc.driftX))
        .force('y', d3.forceY(CY).strength(fc.driftY))
        .stop();

      // Run the simulation to completion up front rather than animating: with a
      // few hundred nodes a per-frame DOM update is what melts the browser.
      const n = Math.ceil(Math.log(simulation.alphaMin()) / Math.log(1 - simulation.alphaDecay()));
      for (let i = 0; i < n * fc.tickMultiplier; i++) simulation.tick();
    }

    ticked();
    // Auto-select the first edge so the viewer pane is populated immediately.
    if (activeEdgeIdx < 0 && edges.length > 0) activeEdgeIdx = 0;
    refreshSelection();
    selectEdge(activeEdgeIdx);
  }

  layoutSelect.addEventListener('change', () => render(layoutSelect.value));
  render(layoutSelect.value);

  if (!graph.edges.length) mappingViewer.message('This network has no edges.');

  return {
    onResize() { render(layoutSelect.value); mappingViewer.resize(); },
    cleanup() {
      if (simulation) simulation.stop();
      mappingViewer.destroy();
    }
  };
}

// ============================================================================
// VIEW: AlchemicalNetwork — force graph of ChemicalSystems / Transformations
//
// The payload carries topology + composition only (no SDF/PDB), so nodes are
// glyphs coloured and shaped by their *component signature* — the sorted set of
// component types the ChemicalSystem contains. That is what distinguishes a
// complex leg (protein + ligand + solvent) from a solvent leg (ligand +
// solvent), which is the thing users are actually scanning for.
// ============================================================================

const SIG_COLORS = ['#0f766e', '#b45309', '#6d28d9', '#be123c', '#0369a1', '#4d7c0f', '#a21caf', '#78350f'];

const ALCHEMICAL_CONFIG = {
  node: { size: 420, radius: 18, strokeWidth: 1.5, labelFontSize: 10, labelMaxChars: 16, labelMaxNodes: 200 },
  edge: { width: 1.4, opacity: 0.75, haloWidth: 6, arrowRefX: 25, arrowWidthPx: 8, arrowHeightPx: 8, hitWidth: 12 },
  force: {
    linkDistance: 90, linkStrength: 0.35, chargeStrength: -420, chargeDistanceMax: 900,
    centerStrength: 0.08, collisionPadding: 6, collisionIterations: 2,
    driftX: 0.05, driftY: 0.05, maxTicks: 400
  },
  circular: { radiusFraction: 0.40 },
  dimOpacity: 0.12
};

const signatureOf = (node) => {
  const types = Array.from(new Set((node.components || []).map(c => shortType(c && c.type)))).sort();
  return types.length ? types.join(' + ') : '(no components)';
};

const nodeLabel = (node) => node.name || String(node.id || '').slice(0, 8);

const truncate = (s, n) => s.length > n ? s.slice(0, n - 1) + '…' : s;

async function viewAlchemicalNetwork(host, inputs) {
  const net = await asObject(inputs.alchemical_network);

  const split = el('div', 'flex:1;display:flex;flex-direction:row;min-height:0;overflow:hidden;');
  host.appendChild(split);

  const leftPane = el('div',
    'flex:1 1 auto;min-width:0;display:flex;flex-direction:column;background:' + T.netCanvasBg + ';');
  const panel = el('div',
    'flex:0 0 280px;width:280px;overflow:auto;padding:14px 16px;font-size:12px;box-sizing:border-box;' +
    'background:' + T.panelBg + ';color:' + T.textPrimary + ';');
  split.appendChild(leftPane);
  split.appendChild(el('div', 'width:1px;flex-shrink:0;background:' + T.splitBorder + ';'));
  split.appendChild(panel);

  const svgWrap = el('div',
    'flex:1;position:relative;overflow:hidden;min-height:0;background:' + T.netCanvasBg + ';');
  leftPane.appendChild(svgWrap);

  const toolbar = el('div',
    'display:flex;align-items:center;gap:12px;flex-wrap:wrap;padding:8px 14px;flex-shrink:0;' +
    'background:' + T.toolbarBg + ';border-top:1px solid ' + T.toolbarBorder + ';');
  leftPane.appendChild(toolbar);

  const titleEl = el('span',
    'font-weight:700;font-size:13px;letter-spacing:.02em;color:' + T.titleColor + ';', 'Alchemical Network');
  toolbar.appendChild(titleEl);

  const legendWrap = el('div',
    'display:flex;align-items:center;gap:10px;flex-wrap:wrap;font-size:11px;color:' + T.textMuted + ';');
  toolbar.appendChild(legendWrap);
  toolbar.appendChild(el('span', 'margin-left:auto;'));

  const filterInput = el('input', SELECT_CSS + 'width:140px;cursor:auto;');
  filterInput.type = 'search';
  filterInput.placeholder = 'Filter nodes…';
  toolbar.appendChild(filterInput);

  const layoutSelect = el('select', SELECT_CSS);
  ['Force-directed', 'Circular'].forEach(l => {
    const o = el('option', '', l);
    o.value = l;
    layoutSelect.appendChild(o);
  });
  toolbar.appendChild(el('label', 'font-size:12px;color:' + T.textMuted + ';', 'Layout:'));
  toolbar.appendChild(layoutSelect);

  const tooltip = el('div',
    'position:absolute;pointer-events:none;opacity:0;max-width:260px;padding:8px 12px;border-radius:8px;' +
    'font-size:11px;z-index:10;transition:opacity .12s;box-shadow:0 4px 16px rgba(15,23,42,0.18);' +
    'background:' + T.tooltipBg + ';border:1px solid ' + T.tooltipBorder + ';color:' + T.textPrimary + ';');
  svgWrap.appendChild(tooltip);

  const panelPlaceholder = (msg) => {
    panel.innerHTML = '<div style="color:' + T.textMuted2 + ';line-height:1.5;">' + esc(msg) + '</div>';
  };

  if (!net || !Array.isArray(net.nodes)) {
    panelPlaceholder('Waiting for an AlchemicalNetwork…');
    return {};
  }

  const nodesIn = net.nodes;
  const edgesIn = Array.isArray(net.edges) ? net.edges : [];
  const byId = Object.fromEntries(nodesIn.map(n => [n.id, n]));

  titleEl.textContent = (net.name || 'Alchemical Network') +
    '  ·  ' + nodesIn.length + ' systems, ' + edgesIn.length + ' transformations';

  const d3 = await loadD3();
  const SIG_SHAPES = [
    d3.symbolCircle, d3.symbolSquare, d3.symbolDiamond, d3.symbolTriangle,
    d3.symbolWye, d3.symbolStar, d3.symbolCross, d3.symbolCircle
  ];

  // Most common signature first, so the dominant leg gets the primary colour.
  const counts = new Map();
  nodesIn.forEach(n => {
    const sig = signatureOf(n);
    counts.set(sig, (counts.get(sig) || 0) + 1);
  });
  const sigIndex = new Map();
  Array.from(counts.entries())
    .sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0]))
    .forEach(([sig, count], i) => sigIndex.set(sig, {
      color: SIG_COLORS[i % SIG_COLORS.length],
      shape: SIG_SHAPES[i % SIG_SHAPES.length],
      count
    }));

  const colorForSig = (sig) => (sigIndex.get(sig) || {}).color || T.textMuted;
  const shapeForSig = (sig) => (sigIndex.get(sig) || {}).shape || d3.symbolCircle;

  sigIndex.forEach((entry, sig) => {
    const item = el('span', 'display:inline-flex;align-items:center;gap:5px;');
    const glyph = d3.create('svg').attr('width', 14).attr('height', 14);
    glyph.append('path')
      .attr('transform', 'translate(7,7)')
      .attr('d', d3.symbol().type(entry.shape).size(70)())
      .attr('fill', entry.color);
    item.appendChild(glyph.node());
    item.appendChild(el('span', '', sig + ' (' + entry.count + ')'));
    legendWrap.appendChild(item);
  });

  // ─── detail panel ───
  const panelHeading = (text) =>
    '<div style="font-weight:700;font-size:11px;letter-spacing:.08em;text-transform:uppercase;color:' +
    T.titleColor + ';margin-bottom:8px;">' + esc(text) + '</div>';

  const panelRow = (k, v) =>
    '<div style="margin-bottom:6px;"><div style="color:' + T.textMuted2 + ';font-size:10px;">' + esc(k) +
    '</div><div style="color:' + T.textPrimary + ';word-break:break-all;">' + esc(v) + '</div></div>';

  function showNodeDetail(node) {
    const comps = node.components || [];
    let html = panelHeading('Chemical System') +
      panelRow('name', node.name || '(unnamed)') +
      panelRow('key', node.id) +
      panelRow('signature', signatureOf(node)) +
      '<div style="margin-top:12px;font-weight:700;color:' + T.textMuted + ';">Components (' + comps.length + ')</div>';
    if (!comps.length) {
      html += '<div style="color:' + T.textMuted2 + ';margin-top:4px;">none</div>';
    } else {
      html += '<div style="margin-top:6px;">' + comps.map(c =>
        '<div style="padding:6px 8px;margin-bottom:4px;border-radius:0 4px 4px 0;' +
        'background:rgba(127,127,127,0.08);border-left:2px solid ' + colorForSig(signatureOf(node)) + ';">' +
        '<div style="font-weight:600;">' + esc(c.label) + '</div>' +
        '<div style="color:' + T.textMuted2 + ';font-size:10px;">' + esc(c.type) +
        (c.name ? ' — ' + esc(c.name) : '') + '</div></div>').join('') + '</div>';
    }
    panel.innerHTML = html;
  }

  function showEdgeDetail(edge) {
    const s = byId[edge.source.id ?? edge.source];
    const t = byId[edge.target.id ?? edge.target];
    panel.innerHTML = panelHeading('Transformation') +
      panelRow('name', edge.name || '(unnamed)') +
      panelRow('key', edge.id) +
      panelRow('protocol', edge.protocol || '(unknown)') +
      '<div style="margin-top:12px;font-weight:700;color:' + T.textMuted + ';">Endpoints</div>' +
      '<div style="margin-top:6px;">' +
      panelRow('state A (source)', s ? nodeLabel(s) : String(edge.source)) +
      panelRow('state B (target)', t ? nodeLabel(t) : String(edge.target)) + '</div>';
  }

  panelPlaceholder('Click a chemical system or a transformation to see its details.');

  // ─── graph ───
  let simulation = null;
  let svgSel = null;
  let selection = null;
  let applyFilter = () => {};

  function render(layout) {
    if (svgSel) { svgSel.remove(); svgSel = null; }
    if (simulation) { simulation.stop(); simulation = null; }

    const W = svgWrap.getBoundingClientRect().width || 800;
    const H = svgWrap.getBoundingClientRect().height || 600;
    const CX = W / 2, CY = H / 2;
    const NC = ALCHEMICAL_CONFIG.node, EC = ALCHEMICAL_CONFIG.edge;

    svgSel = d3.select(svgWrap).append('svg').attr('width', W).attr('height', H).style('display', 'block');
    const gRoot = svgSel.append('g');
    svgSel.call(d3.zoom().scaleExtent([0.1, 5]).on('zoom', e => gRoot.attr('transform', e.transform)));

    svgSel.append('defs').append('marker')
      .attr('id', 'an-arrow').attr('viewBox', '0 -5 10 10')
      .attr('refX', EC.arrowRefX).attr('refY', 0)
      .attr('markerUnits', 'userSpaceOnUse')
      .attr('markerWidth', EC.arrowWidthPx).attr('markerHeight', EC.arrowHeightPx)
      .attr('orient', 'auto')
      .append('path').attr('d', 'M0,-5L10,0L0,5').attr('fill', T.netEdgeLine);

    const nodes = nodesIn.map(n => ({ ...n, _sig: signatureOf(n), _label: nodeLabel(n) }));
    const nodeById = Object.fromEntries(nodes.map(n => [n.id, n]));
    const edges = edgesIn
      .filter(e => nodeById[e.source] && nodeById[e.target])
      .map(e => ({ ...e, source: nodeById[e.source], target: nodeById[e.target] }));

    if (layout === 'Circular') {
      const r = Math.min(W, H) * ALCHEMICAL_CONFIG.circular.radiusFraction;
      nodes.forEach((n, i) => {
        const a = (2 * Math.PI * i / nodes.length) - Math.PI / 2;
        n.x = n.fx = CX + r * Math.cos(a);
        n.y = n.fy = CY + r * Math.sin(a);
      });
    } else {
      nodes.forEach(n => { n.fx = null; n.fy = null; });
    }

    const linkHalo = gRoot.append('g').selectAll('line').data(edges).join('line')
      .attr('stroke', T.netHaloColor).attr('stroke-width', EC.haloWidth)
      .attr('stroke-linecap', 'round').attr('opacity', 0).style('pointer-events', 'none');

    const link = gRoot.append('g').selectAll('line').data(edges).join('line')
      .attr('stroke', T.netEdgeLine).attr('stroke-width', EC.width)
      .attr('stroke-opacity', EC.opacity).attr('marker-end', 'url(#an-arrow)')
      .style('pointer-events', 'none');

    const moveTooltip = (ev) => {
      tooltip.style.left = (ev.offsetX + 14) + 'px';
      tooltip.style.top = (ev.offsetY - 10) + 'px';
    };
    const hideTooltip = () => { tooltip.style.opacity = '0'; };

    const linkHit = gRoot.append('g').selectAll('line').data(edges).join('line')
      .attr('stroke', 'transparent').attr('stroke-width', EC.hitWidth).style('cursor', 'pointer')
      .on('mouseenter', (ev, d) => {
        tooltip.style.opacity = '1';
        tooltip.innerHTML =
          '<div style="font-weight:700;color:' + T.titleColor + ';">' + esc(d.name || 'Transformation') + '</div>' +
          '<div style="margin-top:3px;">' + esc(d.source._label) + ' → ' + esc(d.target._label) + '</div>' +
          '<div style="margin-top:3px;color:' + T.textMuted2 + ';">' + esc(d.protocol || '') + '</div>';
      })
      .on('mousemove', moveTooltip)
      .on('mouseleave', hideTooltip)
      .on('click', (ev, d) => {
        ev.stopPropagation();
        hideTooltip();
        selection = { kind: 'edge', id: d.id };
        showEdgeDetail(d);
        refreshSelection();
      });

    const nodeG = gRoot.append('g').selectAll('g').data(nodes).join('g')
      .style('cursor', 'grab')
      .call(d3.drag()
        .on('start', (ev, d) => { d.fx = d.x; d.fy = d.y; })
        .on('drag', (ev, d) => { d.fx = d.x = ev.x; d.fy = d.y = ev.y; ticked(); })
        .on('end', (ev, d) => { d.fx = ev.x; d.fy = ev.y; }))
      .on('mouseenter', (ev, d) => {
        tooltip.style.opacity = '1';
        const comps = d.components || [];
        tooltip.innerHTML =
          '<div style="font-weight:700;color:' + T.titleColor + ';">' + esc(d._label) + '</div>' +
          (comps.length
            ? '<div style="margin-top:4px;line-height:1.5;">' + comps.map(c =>
              esc(c.label) + ': ' + esc(c.type) + (c.name ? ' (' + esc(c.name) + ')' : '')).join('<br>') + '</div>'
            : '<div style="color:' + T.textMuted2 + ';margin-top:3px;">no components</div>');
      })
      .on('mousemove', moveTooltip)
      .on('mouseleave', hideTooltip)
      .on('click', (ev, d) => {
        ev.stopPropagation();
        selection = { kind: 'node', id: d.id };
        showNodeDetail(d);
        refreshSelection();
      });

    const nodeHalo = nodeG.append('path')
      .attr('d', d => d3.symbol().type(shapeForSig(d._sig)).size(NC.size * 2.1)())
      .attr('fill', T.netHaloColor).attr('opacity', 0).style('pointer-events', 'none');

    nodeG.append('path')
      .attr('d', d => d3.symbol().type(shapeForSig(d._sig)).size(NC.size)())
      .attr('fill', d => colorForSig(d._sig))
      .attr('stroke', T.netNodeStroke).attr('stroke-width', NC.strokeWidth);

    if (nodes.length <= NC.labelMaxNodes) {
      nodeG.append('text')
        .attr('text-anchor', 'middle').attr('y', NC.radius + 12)
        .attr('font-size', NC.labelFontSize + 'px').attr('font-weight', '600')
        .attr('fill', T.netNodeLabel).attr('pointer-events', 'none')
        .text(d => truncate(d._label, NC.labelMaxChars));
    }

    function ticked() {
      [link, linkHalo, linkHit].forEach(sel => sel
        .attr('x1', d => d.source.x).attr('y1', d => d.source.y)
        .attr('x2', d => d.target.x).attr('y2', d => d.target.y));
      nodeG.attr('transform', d => 'translate(' + (d.x ?? CX) + ',' + (d.y ?? CY) + ')');
    }

    function refreshSelection() {
      nodeHalo.attr('opacity', d => selection && selection.kind === 'node' && selection.id === d.id ? 0.9 : 0);
      linkHalo.attr('opacity', d => selection && selection.kind === 'edge' && selection.id === d.id ? 0.9 : 0);
    }

    applyFilter = () => {
      const q = filterInput.value.trim().toLowerCase();
      if (!q) {
        nodeG.attr('opacity', 1);
        link.attr('stroke-opacity', EC.opacity);
        return;
      }
      const matches = new Set(nodes
        .filter(n => (n._label + ' ' + n.id + ' ' + n._sig).toLowerCase().indexOf(q) !== -1)
        .map(n => n.id));
      nodeG.attr('opacity', d => matches.has(d.id) ? 1 : ALCHEMICAL_CONFIG.dimOpacity);
      link.attr('stroke-opacity', d => (matches.has(d.source.id) || matches.has(d.target.id))
        ? EC.opacity : ALCHEMICAL_CONFIG.dimOpacity);
    };

    if (layout !== 'Circular') {
      const fc = ALCHEMICAL_CONFIG.force;
      simulation = d3.forceSimulation(nodes)
        .force('link', d3.forceLink(edges).id(d => d.id).distance(fc.linkDistance).strength(fc.linkStrength))
        .force('charge', d3.forceManyBody().strength(fc.chargeStrength).distanceMax(fc.chargeDistanceMax))
        .force('center', d3.forceCenter(CX, CY).strength(fc.centerStrength))
        .force('collision', d3.forceCollide(NC.radius + fc.collisionPadding).iterations(fc.collisionIterations))
        .force('x', d3.forceX(CX).strength(fc.driftX))
        .force('y', d3.forceY(CY).strength(fc.driftY))
        .stop();

      // Run to completion up front rather than animating — see the ligand
      // network note: a per-frame DOM update is what melts the browser.
      const need = Math.ceil(Math.log(simulation.alphaMin()) / Math.log(1 - simulation.alphaDecay()));
      for (let i = 0; i < Math.min(need, fc.maxTicks); i++) simulation.tick();
    }

    ticked();
    refreshSelection();
    applyFilter();

    // Fit the drawn graph into view.
    if (nodes.length) {
      const pad = 60;
      const xs = nodes.map(n => n.x), ys = nodes.map(n => n.y);
      const minX = Math.min(...xs) - pad, maxX = Math.max(...xs) + pad;
      const minY = Math.min(...ys) - pad, maxY = Math.max(...ys) + pad;
      const k = Math.min(1, W / (maxX - minX), H / (maxY - minY));
      gRoot.attr('transform',
        'translate(' + (W / 2 - k * (minX + maxX) / 2) + ',' + (H / 2 - k * (minY + maxY) / 2) + ') scale(' + k + ')');
    }
  }

  layoutSelect.addEventListener('change', () => render(layoutSelect.value));
  filterInput.addEventListener('input', () => applyFilter());
  render(layoutSelect.value);

  return {
    onResize() { render(layoutSelect.value); },
    cleanup() { if (simulation) simulation.stop(); }
  };
}

// ============================================================================
// VIEW: LigandAtomMapping — the mapping viewer, standalone
// ============================================================================

async function viewLigandAtomMapping(host, inputs) {
  const sdfA = await asText(inputs['molA.sdf']);
  const sdfB = await asText(inputs['molB.sdf']);
  // gufe components may carry an empty name; fall back to the SDF title line
  // (handled in parseSDF) rather than clobbering it with a generic label.
  const nameA = (await asText(inputs.nameA)) || '';
  const nameB = (await asText(inputs.nameB)) || '';
  const annotations = await asObject(inputs.annotations);

  const header = headerStrip('Ligand atom mapping', '');
  host.appendChild(header);

  const body = el('div', 'flex:1;min-height:0;display:flex;flex-direction:column;');
  host.appendChild(body);
  const viewer = createMappingViewer(body, 'colored');

  if (!sdfA || !sdfB) {
    viewer.message('Payload is missing molA.sdf / molB.sdf.');
    return { onResize: viewer.resize, cleanup: viewer.destroy };
  }

  let molA, molB;
  try {
    molA = parseSDF(sdfA, nameA || 'molecule A');
    molB = parseSDF(sdfB, nameB || 'molecule B');
  } catch (e) {
    viewer.message('SDF parse error: ' + errText(e), true);
    return { onResize: viewer.resize, cleanup: viewer.destroy };
  }
  // Payload names win over whatever title line the SDF happens to carry.
  if (nameA) molA.name = nameA;
  if (nameB) molB.name = nameB;

  const mapping = viewer.show(molA, molB, inputs.mapping);

  header.titleEl.textContent = molA.name + ' → ' + molB.name;
  header.statsEl.appendChild(statChip('mapped', (mapping ? mapping.size : 0) + ' atoms'));
  header.statsEl.appendChild(statChip('unique ' + molA.name,
    molA.symbols.length - (mapping ? new Set(mapping.keys()).size : 0), T.chipUniqueA));
  header.statsEl.appendChild(statChip('unique ' + molB.name,
    molB.symbols.length - (mapping ? new Set(mapping.values()).size : 0), T.chipUniqueB));
  if (annotations && annotations.score != null) {
    const s = annotations.score;
    header.statsEl.appendChild(statChip('score', typeof s === 'number' ? s.toFixed(3) : String(s)));
  }

  if (!mapping) floatingWarning(body, 'No atom mapping in payload — showing plain 3D views.');

  return { onResize: viewer.resize, cleanup: viewer.destroy };
}

// ============================================================================
// DISPATCH
//
// Each view claims the payload by the key it needs; the first match wins, so
// the more specific pairs (molA.sdf + molB.sdf) come before the singles.
// ============================================================================

const VIEWS = [
  { key: 'network.graphml', render: viewLigandNetwork },
  { key: 'molA.sdf', render: viewLigandAtomMapping },
  { key: 'molecule.sdf', render: viewSmallMolecule },
  { key: 'protein.pdb', render: viewProtein },
  { key: 'solvent_component', render: viewSolvent },
  { key: 'chemical_system', render: viewChemicalSystem },
  { key: 'transformation', render: viewTransformation },
  { key: 'alchemical_network', render: viewAlchemicalNetwork }
];

const mount = (typeof root !== 'undefined' && root) ? root : document.body;
mount.style.background = T.appBg;
mount.style.fontFamily = "'Inter',system-ui,sans-serif";
mount.style.color = T.textPrimary;

let shell = null;    // the current view's container
let active = null;   // { key, handle } for the mounted view

/** Tear the mounted view down and hand back a fresh, empty shell. */
function resetShell() {
  if (active && active.handle && active.handle.cleanup) {
    try { active.handle.cleanup(); } catch (e) { console.warn('[gufe] cleanup failed:', e); }
  }
  active = null;
  mount.innerHTML = '';
  shell = el('div',
    'width:100%;height:100%;display:flex;flex-direction:column;overflow:hidden;background:' + T.appBg + ';');
  mount.appendChild(shell);
  return shell;
}

resetShell().appendChild(centredMessage('Waiting for data…'));

// A host can resize the frame without `onResize` ever firing, so watch the mount
// too. Debounced, because a graph view re-lays-out its whole simulation and a
// drag would otherwise fire this on every pixel.
let resizeTimer = null;
new ResizeObserver(() => {
  clearTimeout(resizeTimer);
  resizeTimer = setTimeout(() => onResize(), 150);
}).observe(mount);

// ============================================================================
// Module entry points
// ============================================================================

export async function onInputs(inputs) {
  if (!inputs) return;
  const view = VIEWS.find(v => inputs[v.key] != null);
  if (!view) {
    resetShell().appendChild(centredMessage(
      'No gufe visualization for these inputs (' + Object.keys(inputs).join(', ') + ').'));
    return;
  }

  const host = resetShell();
  try {
    active = { key: view.key, handle: await view.render(host, inputs) };
  } catch (e) {
    console.warn('[gufe] ' + view.key + ' failed:', e);
    if (typeof logStderr === 'function') logStderr('gufe viz ' + view.key + ' failed: ' + errText(e));
    host.appendChild(centredMessage('Failed to render ' + view.key + ': ' + errText(e), true));
  }
}

export function onResize() {
  if (active && active.handle && active.handle.onResize) active.handle.onResize();
}

export function cleanup() {
  resetShell();
}
