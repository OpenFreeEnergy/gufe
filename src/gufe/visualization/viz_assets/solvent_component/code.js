// ============================================================================
// SolventComponent — a compact settings card.
//
// A gufe.SolventComponent has no coordinates; it is a *specification* of a
// solvent box (solvent SMILES, counter-ions, ion concentration, whether the
// system is neutralized). So we render:
//   • a title row with a small RDKit 2D depiction of the solvent molecule
//   • a row of labelled value chips
//   • a purely illustrative schematic of a solvent box with ion dots, the
//     number of dots scaled by the ion concentration
//
// Designed short-and-wide: the default frame height is only 320px.
// ============================================================================

// ============================================================================
// THEME — flip DARK_MODE to switch the entire app
// ============================================================================
var DARK_MODE = false;

var THEMES = {
  dark: {
    appBg:         '#1a1a2e',
    cardBg:        '#16213e',
    cardBorder:    '#2a4a7f',
    titleColor:    '#7ecfff',
    textPrimary:   '#e6f3ff',
    textMuted:     '#9bb8d6',
    textMuted2:    '#7a96b8',
    chipBg:        '#0f3460',
    chipBorder:    '#2a4a7f',
    boxFill:       '#0f2a4a',
    boxStroke:     '#2a4a7f',
    cationColor:   '#ff79c6',
    anionColor:    '#7ecfff',
    yesBg:         '#14532d',
    yesFg:         '#86efac',
    yesBorder:     '#166534',
    noBg:          '#4a1d1d',
    noFg:          '#fca5a5',
    noBorder:      '#7f1d1d',
    depictionBg:   '#ffffff',
    loadingFg:     '#888',
    errorFg:       '#c33'
  },
  light: {
    appBg:         '#ffffff',
    cardBg:        '#f8fafc',
    cardBorder:    '#e2e8f0',
    titleColor:    '#0369a1',
    textPrimary:   '#1e293b',
    textMuted:     '#475569',
    textMuted2:    '#64748b',
    chipBg:        '#ffffff',
    chipBorder:    '#cbd5e1',
    boxFill:       '#eff6ff',
    boxStroke:     '#bfdbfe',
    cationColor:   '#d62828',
    anionColor:    '#0f766e',
    yesBg:         '#dcfce7',
    yesFg:         '#166534',
    yesBorder:     '#bbf7d0',
    noBg:          '#fee2e2',
    noFg:          '#991b1b',
    noBorder:      '#fecaca',
    depictionBg:   '#ffffff',
    loadingFg:     '#888',
    errorFg:       '#c33'
  }
};

var T = DARK_MODE ? THEMES.dark : THEMES.light;

// ============================================================================
// CONFIG — geometry only; palette comes from T
// ============================================================================
const CONFIG = {
  depictionSize: 96,      // px, RDKit render size for the solvent molecule
  schematic: {
    width:        320,
    height:       120,
    padding:      14,
    dotRadius:    4,
    // molar → dot count. 0.15 M ≈ 9 ion pairs, capped so it stays legible.
    dotsPerMolar: 60,
    minDots:      0,
    maxDots:      26
  }
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

// Object-valued inputs may arrive as a JSON string.
async function asObject(v) {
  if (v == null) return null;
  if (typeof v === 'object' && !(v instanceof Blob) && !(v instanceof ArrayBuffer) && !ArrayBuffer.isView(v)) {
    return v;
  }
  const text = await asText(v);
  if (text == null) return null;
  return JSON.parse(text);
}

// RDKit loader (singleton) — copied verbatim from ligand_network/code.js
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

// 2D depiction SVG from a SMILES string.
function depictSmilesSVG(RDKit, smiles, size) {
  let rdmol = null;
  try {
    rdmol = RDKit.get_mol(smiles);
    if (!rdmol) return null;
    const svg = rdmol.get_svg(size, size);
    return svg || null;
  } catch (e) {
    console.warn('[solvent] depictSmilesSVG threw for', smiles, '-', e.message);
    return null;
  } finally {
    if (rdmol) { try { rdmol.delete(); } catch (e) {} }
  }
}

// "0.15 molar" / "0.15 mol/L" / 0.15 → 0.15. Returns null if unparseable.
function parseConcentration(value) {
  if (value == null) return null;
  if (typeof value === 'number') return isFinite(value) ? value : null;
  const m = String(value).match(/-?\d+(\.\d+)?([eE][-+]?\d+)?/);
  if (!m) return null;
  const n = parseFloat(m[0]);
  return isFinite(n) ? n : null;
}

const EM_DASH = '—';

function displayOrDash(v) {
  if (v == null) return EM_DASH;
  const s = String(v).trim();
  return s.length > 0 ? s : EM_DASH;
}

// ============================================================================
// DOM SCAFFOLD
// ============================================================================
var mount = (typeof root !== 'undefined' && root) ? root : document.body;

mount.innerHTML = '';
mount.style.background = T.appBg;
mount.style.fontFamily = "'Inter',system-ui,sans-serif";
mount.style.color = T.textPrimary;

const card = document.createElement('div');
card.style.cssText =
  'width:100%;height:100%;box-sizing:border-box;padding:14px 18px;' +
  'display:flex;flex-direction:column;gap:12px;overflow:auto;background:' + T.appBg + ';';
mount.appendChild(card);

// ─── Title row: "Solvent" + depiction + SMILES ───
const titleRow = document.createElement('div');
titleRow.style.cssText = 'display:flex;align-items:center;gap:14px;flex-shrink:0;';
card.appendChild(titleRow);

const depictionBox = document.createElement('div');
depictionBox.style.cssText =
  'width:64px;height:64px;flex-shrink:0;border:1px solid ' + T.cardBorder + ';border-radius:10px;' +
  'background:' + T.depictionBg + ';display:flex;align-items:center;justify-content:center;overflow:hidden;';
titleRow.appendChild(depictionBox);

const titleText = document.createElement('div');
titleText.style.cssText = 'display:flex;flex-direction:column;gap:2px;min-width:0;';
titleRow.appendChild(titleText);

const titleEl = document.createElement('div');
titleEl.textContent = 'Solvent';
titleEl.style.cssText =
  'font-weight:700;font-size:17px;color:' + T.titleColor + ';letter-spacing:.03em;';
titleText.appendChild(titleEl);

const smilesEl = document.createElement('div');
smilesEl.style.cssText =
  'font-size:12px;color:' + T.textMuted + ";font-family:ui-monospace,'SF Mono',Menlo,monospace;" +
  'white-space:nowrap;overflow:hidden;text-overflow:ellipsis;';
titleText.appendChild(smilesEl);

// ─── Chip row ───
const chipRow = document.createElement('div');
chipRow.style.cssText = 'display:flex;flex-wrap:wrap;gap:8px;flex-shrink:0;';
card.appendChild(chipRow);

function makeChip(label, valueNode) {
  const chip = document.createElement('div');
  chip.style.cssText =
    'display:flex;flex-direction:column;gap:2px;padding:6px 12px;border-radius:8px;' +
    'background:' + T.chipBg + ';border:1px solid ' + T.chipBorder + ';min-width:76px;';
  const l = document.createElement('div');
  l.textContent = label;
  l.style.cssText =
    'font-size:10px;text-transform:uppercase;letter-spacing:.06em;color:' + T.textMuted2 + ';';
  chip.appendChild(l);
  chip.appendChild(valueNode);
  chipRow.appendChild(chip);
  return chip;
}

function makeValueText(text) {
  const v = document.createElement('div');
  v.textContent = text;
  v.style.cssText = 'font-size:14px;font-weight:600;color:' + T.textPrimary + ';';
  return v;
}

function makePill(on) {
  const p = document.createElement('span');
  p.textContent = on ? 'yes' : 'no';
  const bg = on ? T.yesBg : T.noBg;
  const fg = on ? T.yesFg : T.noFg;
  const br = on ? T.yesBorder : T.noBorder;
  p.style.cssText =
    'display:inline-block;align-self:flex-start;padding:1px 10px;border-radius:999px;' +
    'font-size:12px;font-weight:700;background:' + bg + ';color:' + fg + ';border:1px solid ' + br + ';';
  return p;
}

// ─── Schematic ───
const schematicWrap = document.createElement('div');
schematicWrap.style.cssText =
  'flex:1;min-height:0;display:flex;flex-direction:column;gap:4px;';
card.appendChild(schematicWrap);

const schematicSvgWrap = document.createElement('div');
schematicSvgWrap.style.cssText = 'flex:1;min-height:0;';
schematicWrap.appendChild(schematicSvgWrap);

const schematicCaption = document.createElement('div');
schematicCaption.style.cssText = 'font-size:10px;color:' + T.textMuted2 + ';font-style:italic;';
schematicWrap.appendChild(schematicCaption);

// ─── Error banner (copied pattern from ligand_network's init()) ───
function showError(message) {
  const warn = document.createElement('div');
  warn.textContent = '⚠ ' + message;
  warn.style.cssText =
    'background:#fee2e2;color:#991b1b;border:1px solid #fecaca;padding:6px 14px;' +
    'border-radius:6px;font-size:12px;';
  card.insertBefore(warn, card.firstChild);
}

// ============================================================================
// SCHEMATIC — illustrative solvent box with a scattering of ion dots
// ============================================================================

// Deterministic PRNG so the scatter is stable across re-renders / resizes.
function makeRng(seed) {
  let s = seed >>> 0 || 1;
  return function() {
    s ^= s << 13; s >>>= 0;
    s ^= s >> 17;
    s ^= s << 5;  s >>>= 0;
    return s / 4294967296;
  };
}

function renderSchematic(solvent) {
  const S = CONFIG.schematic;
  const W = S.width, H = S.height;
  const pad = S.padding;

  const conc = parseConcentration(solvent.ion_concentration);
  const hasIons = solvent.positive_ion != null || solvent.negative_ion != null;

  let nPairs = 0;
  if (hasIons && conc != null && conc > 0) {
    nPairs = Math.round(conc * S.dotsPerMolar);
    nPairs = Math.max(S.minDots, Math.min(S.maxDots, nPairs));
  }

  const rng = makeRng(1337);
  const dots = [];
  for (let i = 0; i < nPairs; i++) {
    if (solvent.positive_ion != null) dots.push('+');
    if (solvent.negative_ion != null) dots.push('-');
  }

  const parts = [];
  parts.push(
    '<svg viewBox="0 0 ' + W + ' ' + H + '" preserveAspectRatio="xMidYMid meet" ' +
    'style="width:100%;height:100%;display:block;">'
  );
  // Solvent box
  parts.push(
    '<rect x="1" y="1" width="' + (W - 2) + '" height="' + (H - 2) + '" rx="12" ' +
    'fill="' + T.boxFill + '" stroke="' + T.boxStroke + '" stroke-width="1.5"/>'
  );

  // Ion dots, scattered inside the padded interior.
  const innerW = W - 2 * pad, innerH = H - 2 * pad;
  for (let i = 0; i < dots.length; i++) {
    const cx = pad + rng() * innerW;
    const cy = pad + rng() * innerH;
    const isCation = dots[i] === '+';
    const color = isCation ? T.cationColor : T.anionColor;
    parts.push(
      '<circle cx="' + cx.toFixed(1) + '" cy="' + cy.toFixed(1) + '" r="' + S.dotRadius + '" ' +
      'fill="' + color + '" fill-opacity="0.85"/>'
    );
    parts.push(
      '<text x="' + cx.toFixed(1) + '" y="' + (cy + 2.6).toFixed(1) + '" text-anchor="middle" ' +
      'font-size="7" font-weight="700" fill="#ffffff" pointer-events="none">' +
      (isCation ? '+' : '−') + '</text>'
    );
  }

  if (dots.length === 0) {
    parts.push(
      '<text x="' + (W / 2) + '" y="' + (H / 2 + 4) + '" text-anchor="middle" ' +
      'font-size="11" fill="' + T.textMuted2 + '">no ions</text>'
    );
  }

  parts.push('</svg>');
  schematicSvgWrap.innerHTML = parts.join('');

  const legend = [];
  if (solvent.positive_ion != null) legend.push(solvent.positive_ion + ' ● cation');
  if (solvent.negative_ion != null) legend.push(solvent.negative_ion + ' ● anion');
  schematicCaption.textContent =
    'Illustrative only — not a real solvent configuration. ' +
    'Ion dots are drawn in proportion to the requested concentration' +
    (legend.length ? ' (' + legend.join(', ') + ').' : '.');
}

// ============================================================================
// RENDER
// ============================================================================

let currentSolvent = null;

function renderPlaceholder(message) {
  smilesEl.textContent = '';
  depictionBox.innerHTML =
    '<div style="color:' + T.loadingFg + ';font-size:11px;">' + EM_DASH + '</div>';
  chipRow.innerHTML = '';
  schematicSvgWrap.innerHTML =
    '<div style="display:flex;align-items:center;justify-content:center;height:100%;' +
    'color:' + T.textMuted2 + ';font-size:12px;">' + message + '</div>';
  schematicCaption.textContent = '';
}

function render(solvent) {
  currentSolvent = solvent;

  // Title row
  const smiles = solvent.smiles != null ? String(solvent.smiles) : null;
  smilesEl.textContent = smiles ? 'SMILES  ' + smiles : 'SMILES  ' + EM_DASH;
  smilesEl.title = smiles || '';

  depictionBox.innerHTML =
    '<div style="color:' + T.loadingFg + ';font-size:10px;">…</div>';
  if (smiles) {
    loadRDKit().then(RDKit => {
      const svg = depictSmilesSVG(RDKit, smiles, CONFIG.depictionSize);
      if (!svg) {
        depictionBox.innerHTML =
          '<div style="color:' + T.errorFg + ';font-size:10px;">?</div>';
        return;
      }
      depictionBox.innerHTML = svg;
      const svgEl = depictionBox.querySelector('svg');
      if (svgEl) {
        svgEl.style.width = '100%';
        svgEl.style.height = '100%';
        svgEl.setAttribute('preserveAspectRatio', 'xMidYMid meet');
      }
    }).catch(err => {
      depictionBox.innerHTML =
        '<div style="color:' + T.errorFg + ';font-size:9px;text-align:center;">RDKit</div>';
      console.warn('[solvent] RDKit failed to load:', err.message);
    });
  } else {
    depictionBox.innerHTML =
      '<div style="color:' + T.textMuted2 + ';font-size:14px;">' + EM_DASH + '</div>';
  }

  // Chips
  chipRow.innerHTML = '';
  makeChip('Positive ion', makeValueText(displayOrDash(solvent.positive_ion)));
  makeChip('Negative ion', makeValueText(displayOrDash(solvent.negative_ion)));
  makeChip('Ion concentration', makeValueText(displayOrDash(solvent.ion_concentration)));
  makeChip('Neutralize', makePill(!!solvent.neutralize));

  // Schematic
  renderSchematic(solvent);
}

// ============================================================================
// Module entry points
// ============================================================================
export async function onInputs(inputs) {
  try {
    if (!inputs) { renderPlaceholder('Waiting for solvent data…'); return; }

    const payload = inputs['solvent_component'];
    if (payload == null) { renderPlaceholder('Waiting for solvent data…'); return; }

    const solvent = await asObject(payload);
    if (solvent == null || typeof solvent !== 'object') {
      renderPlaceholder('No solvent specification available.');
      return;
    }
    render(solvent);
  } catch (err) {
    console.warn('[solvent] onInputs failed:', err);
    renderPlaceholder('Could not read solvent data.');
    showError('Solvent payload error: ' + err.message);
  }
}

export function onResize() {
  // The card is pure CSS flex + a viewBox SVG, so it reflows on its own.
  // Re-running the schematic keeps the dot scatter crisp after big resizes.
  if (currentSolvent) renderSchematic(currentSolvent);
}

export function cleanup() {
  // No timers, simulations or viewers to tear down.
}

renderPlaceholder('Waiting for solvent data…');
