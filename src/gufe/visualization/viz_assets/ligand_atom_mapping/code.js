// ============================================================================
// Ligand atom mapping — standalone 3D/2D viewer for a single
// `gufe.LigandAtomMapping`.
//
// This is the right-hand viewer pane of the `ligand_network` frame, lifted out
// to fill the whole page. Same five modes (plain / colored / lines / overlay /
// 2d), same Kabsch alignment, same view-sync loop. The only real differences:
// molecules arrive as SDF strings instead of GraphML nodes, and the mapping
// arrives directly as `inputs.mapping`.
// ============================================================================

var ThreeDmol = $3Dmol;

// ============================================================================
// THEME — flip DARK_MODE to switch the entire app
// ============================================================================
var DARK_MODE = false;

var THEMES = {
  dark: {
    appBg:           '#1a1a2e',
    toolbarBg:       '#16213e',
    toolbarBorder:   '#2a4a7f',
    titleColor:      '#7ecfff',
    textPrimary:     '#e6f3ff',
    textMuted:       '#9bb8d6',
    textMuted2:      '#7a96b8',

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
    chipUniqueA:     '#ff4d4d',
    chipUniqueB:     '#4dff88',
    loadingFg:       '#888',
    errorFg:         '#c33'
  },
  light: {
    appBg:           '#ffffff',
    toolbarBg:       '#f8fafc',
    toolbarBorder:   '#e2e8f0',
    titleColor:      '#0369a1',
    textPrimary:     '#1e293b',
    textMuted:       '#475569',
    textMuted2:      '#64748b',

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
    chipUniqueA:     '#d62828',
    chipUniqueB:     '#2a9d4a',
    loadingFg:       '#888',
    errorFg:         '#c33'
  }
};

var T = DARK_MODE ? THEMES.dark : THEMES.light;

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
function asObject(v) {
  if (v == null) return null;
  if (typeof v === 'string') {
    try { return JSON.parse(v); } catch (e) { return null; }
  }
  return (typeof v === 'object') ? v : null;
}

// ============================================================================
// MOLECULE HELPERS — SDF in, internal shape out, MOL block back out
//
// Internal molecule shape (same as the network frame's `parseMolecule`):
//   { name, symbols: ['C',...], bonds: [[i, j, order], ...], coords: [[x,y,z]] }
// ============================================================================

const NL = String.fromCharCode(10);
const DOLLAR = String.fromCharCode(36);

// Parse a V2000 SDF/MOL record into the internal molecule shape. Only the
// counts line + atom block + bond block are needed; property/data blocks are
// ignored. Atom indices in `bonds` are 0-based.
function parseSDF(sdfText, fallbackName) {
  if (!sdfText || !sdfText.trim()) throw new Error('empty SDF');
  const lines = sdfText.replace(/\r/g, '').split(NL);
  if (lines.length < 4) throw new Error('SDF too short');

  const countsLine = lines[3];
  if (countsLine.indexOf('V3000') !== -1) {
    throw new Error('V3000 molfiles are not supported');
  }
  const nAtoms = parseInt(countsLine.substring(0, 3), 10);
  const nBonds = parseInt(countsLine.substring(3, 6), 10);
  if (!isFinite(nAtoms) || nAtoms <= 0) throw new Error('bad counts line: ' + countsLine);

  const coords = [];
  const symbols = [];
  for (let i = 0; i < nAtoms; i++) {
    const ln = lines[4 + i];
    if (ln == null) throw new Error('truncated atom block');
    const x = parseFloat(ln.substring(0, 10));
    const y = parseFloat(ln.substring(10, 20));
    const z = parseFloat(ln.substring(20, 30));
    const sym = ln.substring(31, 34).trim() || 'X';
    coords.push([x || 0, y || 0, z || 0]);
    symbols.push(sym);
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

// MOL block (V2000) without the SDF terminator — rebuilt from the internal
// shape so the `lines` mode can feed 3Dmol translated coordinates.
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

// 2D depiction SVG (no highlights)
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
    console.warn('[mapping] depictMoleculeSVG threw for', mol.name, '-', e.message);
    return null;
  } finally {
    if (rdmol) { try { rdmol.delete(); } catch (e) {} }
  }
}

// ============================================================================
// DOM SCAFFOLD — header strip on top, viewer filling the rest
// ============================================================================

const root = (typeof globalThis.root !== 'undefined' && globalThis.root) || document.body;
root.innerHTML = '';
root.style.background = T.appBg;
root.style.fontFamily = "'Inter',system-ui,sans-serif";

const appWrap = document.createElement('div');
appWrap.style.cssText = 'width:100%;height:100%;display:flex;flex-direction:column;overflow:hidden;background:' + T.appBg + ';';
root.appendChild(appWrap);

// ─── Header strip ───
const header = document.createElement('div');
header.style.cssText = 'display:flex;align-items:center;gap:14px;flex-wrap:wrap;padding:8px 14px;background:' + T.toolbarBg + ';border-bottom:1px solid ' + T.toolbarBorder + ';flex-shrink:0;';
appWrap.appendChild(header);

const headerTitle = document.createElement('span');
headerTitle.style.cssText = 'font-weight:700;font-size:14px;color:' + T.titleColor + ';';
header.appendChild(headerTitle);

const headerStats = document.createElement('div');
headerStats.style.cssText = 'display:flex;align-items:center;gap:12px;font-size:11px;color:' + T.textMuted + ';margin-left:auto;flex-wrap:wrap;';
header.appendChild(headerStats);

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

function updateHeader() {
  headerTitle.textContent = '';
  headerStats.innerHTML = '';
  if (!currentMols || currentMols.length < 2) {
    headerTitle.textContent = 'Ligand atom mapping';
    return;
  }
  const molA = currentMols[0], molB = currentMols[1];
  headerTitle.textContent = molA.name + ' → ' + molB.name;

  const mapping = resolveMapping();
  const nMapped = mapping ? mapping.size : 0;
  const uniqueA = molA.symbols.length - (mapping ? new Set(mapping.keys()).size : 0);
  const uniqueB = molB.symbols.length - (mapping ? new Set(mapping.values()).size : 0);

  headerStats.appendChild(statChip('mapped', nMapped + ' atoms'));
  headerStats.appendChild(statChip('unique ' + molA.name, uniqueA, T.chipUniqueA));
  headerStats.appendChild(statChip('unique ' + molB.name, uniqueB, T.chipUniqueB));

  if (currentAnnotations && currentAnnotations.score != null) {
    const s = currentAnnotations.score;
    headerStats.appendChild(statChip('score', typeof s === 'number' ? s.toFixed(3) : String(s)));
  }
}

// ─── Viewer area ───
const viewerWrapper = document.createElement('div');
viewerWrapper.style.cssText = 'position:relative;flex:1;min-height:0;display:flex;flex-direction:column;';
appWrap.appendChild(viewerWrapper);

const viewerArea = document.createElement('div');
viewerArea.style.cssText = 'flex:1;display:flex;flex-direction:column;min-height:0;';
viewerWrapper.appendChild(viewerArea);

// ============================================================================
// VIEW MODES
// ============================================================================

const MODES = [
  { id: 'plain',   label: '3D',      title: 'Plain 3D view' },
  { id: 'colored', label: '3D-Map',  title: 'Color-coded by mapping' },
  { id: 'lines',   label: 'Pairs',   title: 'Dashed lines between mapped atoms' },
  { id: 'overlay', label: 'Overlay', title: 'Both molecules superimposed' },
  { id: '2d',      label: '2D',      title: '2D depictions with highlights' }
];

let currentMode = 'colored';
let currentMols = null;
let currentMapping = null;      // Map<aIdx, bIdx>
let currentAnnotations = null;

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

function resolveMapping() {
  return (currentMapping && currentMapping.size > 0) ? currentMapping : null;
}

// ============================================================================
// VIEWER BOX MANAGEMENT
// ============================================================================

let boxes = [];
let syncing = false;
let syncRunning = false;
let syncStopped = false;
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
    if (syncStopped) { syncRunning = false; return; }
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
    console.warn('[mapping] No mapping available for', mols[0].name, '->', mols[1].name);
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
    console.warn('[mapping] No mapping available; falling back to plain');
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

  const b = makeViewerBox(mol0.name + ' ↔ ' + mol1.name + '  (' + mapping.size + ' mapped pairs)');
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
    container.innerHTML = '<div style="color:' + T.loadingFg + ';font-size:12px;">Loading 2D depiction…</div>';
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

function renderMode() {
  if (!currentMols) return;
  switch (currentMode) {
    case 'colored': renderColored(currentMols); break;
    case 'lines':   renderLines(currentMols);   break;
    case 'overlay': renderOverlay(currentMols); break;
    case '2d':      render2D(currentMols);      break;
    default:        renderPlain(currentMols);
  }
}

// ============================================================================
// PLACEHOLDER / ERROR BANNERS
// ============================================================================

function showMessage(text, isError) {
  clearAllBoxes();
  const el = document.createElement('div');
  el.style.cssText = 'flex:1;display:flex;align-items:center;justify-content:center;text-align:center;padding:24px;font-size:13px;color:' +
    (isError ? T.errorFg : T.textMuted2) + ';';
  el.textContent = text;
  viewerArea.appendChild(el);
}

function showErrorBanner(msg) {
  const warn = document.createElement('div');
  warn.textContent = '⚠ ' + msg;
  warn.setAttribute('data-banner', '1');
  warn.style.cssText = 'position:absolute;top:10px;left:50%;transform:translateX(-50%);background:#fee2e2;color:#991b1b;border:1px solid #fecaca;padding:6px 14px;border-radius:6px;font-size:12px;z-index:20;max-width:90%;';
  viewerWrapper.appendChild(warn);
}

// ============================================================================
// RESIZE OBSERVER
// ============================================================================

const viewerRO = new ResizeObserver(() => {
  for (let i = 0; i < boxes.length; i++) {
    if (boxes[i].viewer) { boxes[i].viewer.resize(); boxes[i].viewer.render(); }
  }
});
viewerRO.observe(viewerWrapper);

// ============================================================================
// INIT
// ============================================================================

updateSwitcherHighlight();
updateHeader();
showMessage('Waiting for mapping data…', false);

function init(sdfA, sdfB, nameA, nameB, mappingObj, annotations) {
  // Clear any stale error banner from a previous payload.
  viewerWrapper.querySelectorAll('[data-banner]').forEach(n => n.remove());

  let molA, molB;
  try {
    molA = parseSDF(sdfA, nameA || 'molecule A');
    molB = parseSDF(sdfB, nameB || 'molecule B');
  } catch (e) {
    console.warn('[mapping] SDF parse failed: ' + e.message);
    currentMols = null;
    updateHeader();
    showMessage('SDF parse error: ' + e.message, true);
    return;
  }

  // Payload names win over whatever title line the SDF happens to carry.
  if (nameA) molA.name = nameA;
  if (nameB) molB.name = nameB;

  const m = new Map();
  if (mappingObj && typeof mappingObj === 'object') {
    Object.keys(mappingObj).forEach(k => {
      const ai = parseInt(k, 10);
      const bi = mappingObj[k];
      if (isFinite(ai) && typeof bi === 'number') m.set(ai, bi);
    });
  }

  currentMols = [molA, molB];
  currentMapping = m.size > 0 ? m : null;
  currentAnnotations = annotations || null;

  updateHeader();
  updateSwitcherHighlight();
  renderMode();

  if (!currentMapping) {
    showErrorBanner('No atom mapping in payload — showing plain 3D views.');
  }
}

// ============================================================================
// Module entry points
// ============================================================================

export async function onInputs(inputs) {
  if (!inputs) return;
  try {
    const sdfA = await asText(inputs['molA.sdf']);
    const sdfB = await asText(inputs['molB.sdf']);
    if (!sdfA || !sdfB) {
      showMessage('Payload is missing molA.sdf / molB.sdf.', false);
      return;
    }
    // gufe components may carry an empty name; fall back to the SDF title line
    // (handled in parseSDF) rather than clobbering it with a generic label.
    const nameA = (await asText(inputs['nameA'])) || '';
    const nameB = (await asText(inputs['nameB'])) || '';
    const mappingObj = asObject(inputs['mapping']);
    const annotations = asObject(inputs['annotations']);
    init(sdfA, sdfB, nameA, nameB, mappingObj, annotations);
  } catch (e) {
    console.warn('[mapping] onInputs failed: ' + e.message);
    showMessage('Failed to render mapping: ' + e.message, true);
  }
}

export function onResize() {
  for (let i = 0; i < boxes.length; i++) {
    if (boxes[i].viewer) { boxes[i].viewer.resize(); boxes[i].viewer.render(); }
  }
}

export function cleanup() {
  syncStopped = true;
  try { viewerRO.disconnect(); } catch (e) {}
  clearAllBoxes();
}
