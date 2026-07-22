// ============================================================================
// Transformation — one alchemical edge, stateA → stateB.
//
// Top: a two-column component diff (stateA left, stateB right, aligned by
// label) marking each row unchanged / changed / added / removed.
// Bottom: if the edge carries atom mappings, a side-by-side 3Dmol viewer of
// molA and molB with core atoms grey and unique atoms red (A) / green (B).
// ============================================================================

var ThreeDmol = $3Dmol;

var mount = (typeof root !== 'undefined' && root) ? root : document.body;

// ============================================================================
// THEME — flip DARK_MODE to switch the entire app
// ============================================================================
var DARK_MODE = false;

var THEMES = {
  dark: {
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
    labelBg:         '#16213e',
    labelFg:         '#7ecfff',
    badgeBg:         '#0f3460',
    badgeFg:         '#7ecfff',
    cardBg:          '#16213e',
    cardBorder:      '#2a4a7f',

    // Diff status palette (CSS colors)
    diffUnchanged:   '#64748b',
    diffChanged:     '#d9a300',
    diffAdded:       '#2a9d4a',
    diffRemoved:     '#d62828',

    // Viewer panel
    viewerBg:        '0x1a1a2e',
    colorCore:       '0xaaaaaa',
    colorUniqueA:    '0xff4d4d',
    colorUniqueB:    '0x4dff88',
    chipCore:        '#aaaaaa',
    chipUniqueA:     '#ff4d4d',
    chipUniqueB:     '#4dff88',
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
    labelBg:         '#f0f4fa',
    labelFg:         '#1a4a8a',
    badgeBg:         '#e1e8f2',
    badgeFg:         '#1a4a8a',
    cardBg:          '#ffffff',
    cardBorder:      '#e2e8f0',

    diffUnchanged:   '#94a3b8',
    diffChanged:     '#b45309',
    diffAdded:       '#2a9d4a',
    diffRemoved:     '#d62828',

    viewerBg:        '0xffffff',
    colorCore:       '0x888888',
    colorUniqueA:    '0xd62828',
    colorUniqueB:    '0x2a9d4a',
    chipCore:        '#888888',
    chipUniqueA:     '#d62828',
    chipUniqueB:     '#2a9d4a',
    loadingFg:       '#888',
    errorFg:         '#c33'
  }
};

var T = DARK_MODE ? THEMES.dark : THEMES.light;

const NL = String.fromCharCode(10);
const DOLLAR = String.fromCharCode(36);

// ============================================================================
// INPUT NORMALIZATION
// ============================================================================

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
  if (v instanceof Blob || v instanceof ArrayBuffer || ArrayBuffer.isView(v)) {
    v = await asText(v);
  }
  if (typeof v === 'string') {
    try { return JSON.parse(v); } catch (e) { return null; }
  }
  return (typeof v === 'object') ? v : null;
}

// ============================================================================
// SDF PARSING (V2000) — elements, coordinates, bonds
// Internal molecule shape: { name, symbols, bonds: [[i,j,order]], coords }
// ============================================================================

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

// ============================================================================
// DOM SCAFFOLD
// ============================================================================

mount.innerHTML = '';
mount.style.background = T.appBg;
mount.style.fontFamily = "'Inter',system-ui,sans-serif";

const appWrap = document.createElement('div');
appWrap.style.cssText = 'width:100%;height:100%;display:flex;flex-direction:column;overflow:hidden;background:' + T.appBg + ';color:' + T.textPrimary + ';';
mount.appendChild(appWrap);

// ─── Header ───
const header = document.createElement('div');
header.style.cssText = 'display:flex;align-items:center;gap:12px;flex-wrap:wrap;padding:9px 14px;background:' + T.toolbarBg + ';border-bottom:1px solid ' + T.toolbarBorder + ';flex-shrink:0;';
appWrap.appendChild(header);

const headerTitle = document.createElement('span');
headerTitle.style.cssText = 'font-weight:700;font-size:15px;color:' + T.titleColor + ';';
header.appendChild(headerTitle);

const headerName = document.createElement('span');
headerName.style.cssText = 'font-size:12px;color:' + T.textMuted + ';';
header.appendChild(headerName);

const headerBadge = document.createElement('span');
headerBadge.style.cssText = 'margin-left:auto;font-size:11px;font-weight:600;padding:3px 9px;border-radius:10px;background:' + T.badgeBg + ';color:' + T.badgeFg + ';white-space:nowrap;';
header.appendChild(headerBadge);

// ─── Diff section ───
const diffSection = document.createElement('div');
diffSection.style.cssText = 'flex:0 0 auto;max-height:55%;overflow:auto;padding:12px 14px;border-bottom:1px solid ' + T.splitBorder + ';';
appWrap.appendChild(diffSection);

// ─── Mapping section ───
const mappingSection = document.createElement('div');
mappingSection.style.cssText = 'flex:1 1 auto;min-height:0;display:flex;flex-direction:column;position:relative;';
appWrap.appendChild(mappingSection);

const mappingBar = document.createElement('div');
mappingBar.style.cssText = 'display:flex;align-items:center;gap:12px;flex-wrap:wrap;padding:6px 14px;background:' + T.labelBg + ';flex-shrink:0;font-size:11px;color:' + T.textMuted + ';';
mappingSection.appendChild(mappingBar);

const viewerArea = document.createElement('div');
viewerArea.style.cssText = 'flex:1;display:flex;flex-direction:row;min-height:0;';
mappingSection.appendChild(viewerArea);

// ============================================================================
// DIFF RENDERING
// ============================================================================

const DIFF_STATUS = {
  unchanged: { label: 'unchanged', color: () => T.diffUnchanged },
  changed:   { label: 'changed',   color: () => T.diffChanged },
  added:     { label: 'added',     color: () => T.diffAdded },
  removed:   { label: 'removed',   color: () => T.diffRemoved }
};

// One-line summary of a component descriptor: what it actually is.
function componentDetail(c) {
  if (!c) return '';
  if (c.error) return 'error: ' + c.error;
  if (c.smiles && c.sdf) return c.smiles;
  if (c.pdb) return (c.pdb.split(NL).length) + ' PDB lines';
  if (c.smiles) {
    const bits = [c.smiles];
    if (c.ion_concentration) bits.push(c.ion_concentration);
    if (c.positive_ion && c.negative_ion) bits.push(c.positive_ion + '/' + c.negative_ion);
    return bits.join('  ·  ');
  }
  return '';
}

function diffStatus(a, b) {
  if (a && !b) return 'removed';
  if (b && !a) return 'added';
  if (a.type !== b.type) return 'changed';
  if ((a.name || '') !== (b.name || '')) return 'changed';
  return 'unchanged';
}

// A single component cell (one side of a diff row). `null` renders as an
// explicit empty slot so added/removed rows stay visually aligned.
function componentCell(c, status, side) {
  const cell = document.createElement('div');
  cell.style.cssText = 'flex:1 1 0;min-width:0;padding:7px 10px;border-radius:6px;border:1px solid ' + T.cardBorder + ';background:' + T.cardBg + ';';
  if (!c) {
    cell.style.borderStyle = 'dashed';
    cell.style.opacity = '0.5';
    cell.innerHTML = '<span style="font-size:11px;color:' + T.textMuted2 + ';">— absent in state' + side + ' —</span>';
    return cell;
  }
  const nameRow = document.createElement('div');
  nameRow.style.cssText = 'font-size:12px;font-weight:600;color:' + T.textPrimary + ';overflow:hidden;text-overflow:ellipsis;white-space:nowrap;';
  nameRow.textContent = c.name || '(unnamed)';
  nameRow.title = c.name || '';
  cell.appendChild(nameRow);

  const typeRow = document.createElement('div');
  typeRow.style.cssText = 'font-size:10px;color:' + T.textMuted2 + ';margin-top:2px;';
  typeRow.textContent = c.type || '';
  cell.appendChild(typeRow);

  const detail = componentDetail(c);
  if (detail) {
    const detailRow = document.createElement('div');
    detailRow.style.cssText = 'font-size:10px;font-family:ui-monospace,monospace;color:' + T.textMuted + ';margin-top:3px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;';
    detailRow.textContent = detail;
    detailRow.title = detail;
    cell.appendChild(detailRow);
  }
  return cell;
}

function renderLegend() {
  const legend = document.createElement('div');
  legend.style.cssText = 'display:flex;align-items:center;gap:14px;flex-wrap:wrap;font-size:11px;color:' + T.textMuted + ';margin-bottom:10px;';
  Object.keys(DIFF_STATUS).forEach(k => {
    const item = document.createElement('span');
    item.style.cssText = 'display:inline-flex;align-items:center;gap:5px;';
    item.innerHTML =
      '<span style="width:9px;height:9px;border-radius:2px;background:' + DIFF_STATUS[k].color() + ';display:inline-block;"></span>' +
      DIFF_STATUS[k].label;
    legend.appendChild(item);
  });
  return legend;
}

function renderDiff(stateA, stateB) {
  diffSection.innerHTML = '';

  const compsA = {}, compsB = {};
  (stateA.components || []).forEach(c => { compsA[c.label] = c; });
  (stateB.components || []).forEach(c => { compsB[c.label] = c; });
  const labels = Array.from(new Set(Object.keys(compsA).concat(Object.keys(compsB)))).sort();

  // Column headings
  const heads = document.createElement('div');
  heads.style.cssText = 'display:flex;gap:8px;align-items:center;margin-bottom:6px;font-size:11px;font-weight:700;color:' + T.textMuted + ';letter-spacing:.04em;text-transform:uppercase;';
  const spacer = document.createElement('div');
  spacer.style.cssText = 'flex:0 0 120px;';
  spacer.textContent = 'label';
  heads.appendChild(spacer);
  [stateA.name || 'state A', stateB.name || 'state B'].forEach(n => {
    const h = document.createElement('div');
    h.style.cssText = 'flex:1 1 0;min-width:0;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;';
    h.textContent = n;
    heads.appendChild(h);
  });
  diffSection.appendChild(heads);
  diffSection.appendChild(renderLegend());

  if (labels.length === 0) {
    const empty = document.createElement('div');
    empty.style.cssText = 'font-size:12px;color:' + T.textMuted2 + ';padding:8px 0;';
    empty.textContent = 'Neither chemical system declares any components.';
    diffSection.appendChild(empty);
    return;
  }

  labels.forEach(label => {
    const a = compsA[label] || null;
    const b = compsB[label] || null;
    const status = diffStatus(a, b);
    const color = DIFF_STATUS[status].color();

    const row = document.createElement('div');
    row.style.cssText = 'display:flex;gap:8px;align-items:stretch;margin-bottom:6px;';

    const labelCell = document.createElement('div');
    labelCell.style.cssText = 'flex:0 0 120px;display:flex;flex-direction:column;justify-content:center;gap:3px;border-left:3px solid ' + color + ';padding-left:8px;min-width:0;';
    labelCell.innerHTML =
      '<span style="font-size:12px;font-weight:600;color:' + T.textPrimary + ';overflow:hidden;text-overflow:ellipsis;white-space:nowrap;">' + label + '</span>' +
      '<span style="font-size:10px;font-weight:600;color:' + color + ';">' + DIFF_STATUS[status].label + '</span>';
    row.appendChild(labelCell);

    row.appendChild(componentCell(a, status, 'A'));
    row.appendChild(componentCell(b, status, 'B'));
    diffSection.appendChild(row);
  });
}

// ============================================================================
// MAPPING VIEWER — core atoms grey, unique atoms red (A) / green (B).
// Mirrors `ligand_network/code.js`'s `renderColored`.
// ============================================================================

let boxes = [];

function makeViewerBox(labelText) {
  const box = document.createElement('div');
  box.style.cssText = 'flex:1 1 0;display:flex;flex-direction:column;position:relative;min-width:0;min-height:0;';
  const label = document.createElement('div');
  label.style.cssText = 'padding:4px 10px;font-size:12px;font-weight:bold;color:' + T.labelFg + ';background:' + T.labelBg + ';overflow:hidden;text-overflow:ellipsis;white-space:nowrap;';
  label.textContent = labelText;
  box.appendChild(label);
  const container = document.createElement('div');
  container.style.cssText = 'flex:1;position:relative;';
  box.appendChild(container);
  viewerArea.appendChild(box);
  return { label, container, box, viewer: null };
}

function clearViewers() {
  boxes.forEach(b => { if (b.viewer) { try { b.viewer.clear(); } catch (e) {} } });
  viewerArea.innerHTML = '';
  boxes = [];
}

function showMappingMessage(text, isError) {
  clearViewers();
  const el = document.createElement('div');
  el.style.cssText = 'flex:1;display:flex;align-items:center;justify-content:center;text-align:center;padding:24px;font-size:13px;color:' +
    (isError ? T.errorFg : T.textMuted2) + ';';
  el.textContent = text;
  viewerArea.appendChild(el);
}

function statChip(label, value, dotColor) {
  const el = document.createElement('span');
  el.style.cssText = 'display:inline-flex;align-items:center;gap:5px;white-space:nowrap;';
  if (dotColor) {
    const dot = document.createElement('span');
    dot.style.cssText = 'width:8px;height:8px;border-radius:50%;background:' + dotColor + ';display:inline-block;';
    el.appendChild(dot);
  }
  const txt = document.createElement('span');
  txt.innerHTML = label + ' <b style="color:' + T.textPrimary + ';">' + value + '</b>';
  el.appendChild(txt);
  return el;
}

// Turn {"<int>": <int>} into a Map, flipping it if the indices only make sense
// the other way round (same guard `renderColored` uses).
function toMapping(mappingObj, nA, nB) {
  let m = new Map();
  if (mappingObj && typeof mappingObj === 'object') {
    Object.keys(mappingObj).forEach(k => {
      const ai = parseInt(k, 10);
      const bi = mappingObj[k];
      if (isFinite(ai) && typeof bi === 'number') m.set(ai, bi);
    });
  }
  if (m.size === 0) return null;

  let maxK = -1, maxV = -1;
  m.forEach((v, k) => { if (k > maxK) maxK = k; if (v > maxV) maxV = v; });
  if ((maxK >= nA || maxV >= nB) && maxK < nB && maxV < nA) {
    const flipped = new Map();
    m.forEach((v, k) => flipped.set(v, k));
    m = flipped;
  }
  return m;
}

function renderMapping(mappingPayload) {
  mappingBar.innerHTML = '';
  if (mappingSelect) mappingBar.appendChild(mappingSelect);

  let molA, molB;
  try {
    molA = parseSDF(mappingPayload['molA.sdf'], mappingPayload.nameA || 'molecule A');
    molB = parseSDF(mappingPayload['molB.sdf'], mappingPayload.nameB || 'molecule B');
  } catch (e) {
    showMappingMessage('SDF parse error: ' + e.message, true);
    return;
  }
  if (mappingPayload.nameA) molA.name = mappingPayload.nameA;
  if (mappingPayload.nameB) molB.name = mappingPayload.nameB;

  const mapping = toMapping(mappingPayload.mapping, molA.coords.length, molB.coords.length);
  const coreA = mapping ? new Set(mapping.keys())   : new Set();
  const coreB = mapping ? new Set(mapping.values()) : new Set();

  const stats = document.createElement('div');
  stats.style.cssText = 'display:flex;align-items:center;gap:12px;flex-wrap:wrap;margin-left:auto;';
  stats.appendChild(statChip('mapped', (mapping ? mapping.size : 0) + ' atoms', T.chipCore));
  stats.appendChild(statChip('unique A', molA.symbols.length - coreA.size, T.chipUniqueA));
  stats.appendChild(statChip('unique B', molB.symbols.length - coreB.size, T.chipUniqueB));
  mappingBar.appendChild(stats);

  clearViewers();
  const mols = [molA, molB];
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
    for (let ai = 0; ai < mols[i].symbols.length; ai++) {
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
}

// ============================================================================
// STATE + INIT
// ============================================================================

let currentMappings = [];
let mappingSelect = null;

function buildMappingSelect() {
  const sel = document.createElement('select');
  sel.style.cssText = 'background:' + T.selectBg + ';color:' + T.textPrimary + ';border:1px solid ' + T.selectBorder + ';border-radius:6px;padding:3px 8px;font-size:11px;cursor:pointer;font-family:inherit;';
  currentMappings.forEach((m, i) => {
    const o = document.createElement('option');
    o.value = String(i);
    o.textContent = 'mapping ' + (i + 1) + ': ' + (m.nameA || 'A') + ' → ' + (m.nameB || 'B');
    sel.appendChild(o);
  });
  sel.addEventListener('change', () => {
    const idx = parseInt(sel.value, 10);
    if (currentMappings[idx]) renderMapping(currentMappings[idx]);
  });
  return sel;
}

function init(tf) {
  const stateA = tf.stateA || { name: 'state A', components: [] };
  const stateB = tf.stateB || { name: 'state B', components: [] };

  headerTitle.textContent = (stateA.name || 'state A') + '  →  ' + (stateB.name || 'state B');
  headerName.textContent = tf.name ? tf.name : '(unnamed transformation)';
  headerBadge.textContent = tf.protocol || 'unknown protocol';

  renderDiff(stateA, stateB);

  currentMappings = Array.isArray(tf.mappings) ? tf.mappings : [];
  mappingSelect = null;
  mappingBar.innerHTML = '';

  if (currentMappings.length === 0) {
    showMappingMessage('No atom mapping on this transformation — nothing to superimpose.', false);
    mappingBar.textContent = 'atom mapping';
    return;
  }
  if (currentMappings.length > 1) mappingSelect = buildMappingSelect();
  renderMapping(currentMappings[0]);
}

function showFatal(msg) {
  diffSection.innerHTML = '';
  const warn = document.createElement('div');
  warn.textContent = '⚠ ' + msg;
  warn.style.cssText = 'background:#fee2e2;color:#991b1b;border:1px solid #fecaca;padding:6px 14px;border-radius:6px;font-size:12px;';
  diffSection.appendChild(warn);
  showMappingMessage('', false);
}

headerTitle.textContent = 'Transformation';
showMappingMessage('Waiting for transformation data…', false);

// ============================================================================
// RESIZE OBSERVER
// ============================================================================

const viewerRO = new ResizeObserver(() => {
  for (let i = 0; i < boxes.length; i++) {
    if (boxes[i].viewer) { boxes[i].viewer.resize(); boxes[i].viewer.render(); }
  }
});
viewerRO.observe(mappingSection);

// ============================================================================
// Module entry points
// ============================================================================

export async function onInputs(inputs) {
  if (!inputs) return;
  try {
    const tf = await asObject(inputs['transformation']);
    if (!tf || typeof tf !== 'object') {
      showMappingMessage('Payload is missing the `transformation` object.', false);
      return;
    }
    init(tf);
  } catch (e) {
    console.warn('[transformation] onInputs failed: ' + e.message);
    showFatal('Failed to render transformation: ' + e.message);
  }
}

export function onResize() {
  for (let i = 0; i < boxes.length; i++) {
    if (boxes[i].viewer) { boxes[i].viewer.resize(); boxes[i].viewer.render(); }
  }
}

export function cleanup() {
  try { viewerRO.disconnect(); } catch (e) {}
  clearViewers();
}
