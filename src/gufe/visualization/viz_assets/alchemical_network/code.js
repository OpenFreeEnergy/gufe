// ============================================================================
// Alchemical network — force-directed graph of ChemicalSystems (nodes) joined
// by Transformations (edges).
//
// The payload carries topology + composition only (no SDF/PDB), so nodes are
// glyphs coloured/shaped by their *component signature* — the sorted set of
// component types the ChemicalSystem contains. That is what distinguishes a
// complex leg (protein + ligand + solvent) from a solvent leg (ligand +
// solvent), which is the thing users are actually scanning for.
// ============================================================================

import * as d3 from "https://cdn.jsdelivr.net/npm/d3@7/+esm";

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
    selectBg:      '#0f3460',
    selectBorder:  '#2a4a7f',
    tooltipBg:     '#16213e',
    tooltipBorder: '#2a4a7f',
    panelBg:       '#16213e',

    netNodeStroke: '#0b1020',
    netNodeLabel:  '#cfe6ff',
    netEdgeLine:   '#5f7ea8',
    netHaloColor:  '#ff79c6',
    netCanvasBg:   '#1a1a2e',

    errorFg:       '#f88'
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
    selectBg:      '#ffffff',
    selectBorder:  '#cbd5e1',
    tooltipBg:     '#ffffff',
    tooltipBorder: '#cbd5e1',
    panelBg:       '#f8fafc',

    netNodeStroke: '#ffffff',
    netNodeLabel:  '#334155',
    netEdgeLine:   '#94a3b8',
    netHaloColor:  '#fbcfe8',
    netCanvasBg:   '#ffffff',

    errorFg:       '#c33'
  }
};

var T = DARK_MODE ? THEMES.dark : THEMES.light;

// Categorical palette + glyph shapes, indexed by component signature.
const SIG_COLORS = [
  '#0f766e', '#b45309', '#6d28d9', '#be123c',
  '#0369a1', '#4d7c0f', '#a21caf', '#78350f'
];
const SIG_SHAPES = [
  d3.symbolCircle, d3.symbolSquare, d3.symbolDiamond, d3.symbolTriangle,
  d3.symbolWye, d3.symbolStar, d3.symbolCross, d3.symbolCircle
];

const CONFIG = {
  node: {
    size:           420,   // d3.symbol area
    radius:         18,    // approximate radius, used for arrow offset / collision
    strokeWidth:    1.5,
    labelFontSize:  10,
    labelMaxChars:  16,
    labelMaxNodes:  200,   // above this, labels are suppressed for legibility/perf
  },
  edge: {
    width:        1.4,
    opacity:      0.75,
    haloWidth:    6,
    arrowRefX:    25,      // viewBox units; * (markerWidth/viewBoxWidth) = user units
    arrowWidthPx: 8,
    arrowHeightPx: 8,
    hitWidth:     12,
  },
  force: {
    linkDistance:        90,
    linkStrength:        0.35,
    chargeStrength:      -420,
    chargeDistanceMax:   900,
    centerStrength:      0.08,
    collisionPadding:    6,
    collisionIterations: 2,
    driftX:              0.05,
    driftY:              0.05,
    maxTicks:            400,
  },
  circular: { radiusFraction: 0.40 },
  dimOpacity: 0.12,
};

// ============================================================================
// HELPERS
// ============================================================================

async function asText(v) {
  if (v == null) return null;
  if (v instanceof Blob) return await v.text();
  if (v instanceof ArrayBuffer) return new TextDecoder().decode(v);
  if (ArrayBuffer.isView(v)) return new TextDecoder().decode(v);
  return typeof v === 'string' ? v : String(v);
}

// Payload values may arrive as objects or as JSON strings/blobs.
async function asObject(v) {
  if (v == null) return null;
  if (typeof v === 'object' && !(v instanceof Blob) && !(v instanceof ArrayBuffer) && !ArrayBuffer.isView(v)) {
    return v;
  }
  const text = await asText(v);
  if (text == null) return null;
  return JSON.parse(text);
}

function esc(s) {
  return String(s == null ? '' : s)
    .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
}

// "SmallMoleculeComponent" -> "SmallMolecule"; keeps signatures readable.
function shortType(t) {
  const s = String(t || '');
  return s.endsWith('Component') && s.length > 9 ? s.slice(0, -9) : s;
}

function signatureOf(node) {
  const comps = Array.isArray(node.components) ? node.components : [];
  const types = Array.from(new Set(comps.map(c => shortType(c && c.type)))).sort();
  return types.length ? types.join(' + ') : '(no components)';
}

function nodeLabel(node) {
  if (node.name) return node.name;
  return String(node.id || '').slice(0, 8);
}

function truncate(s, n) {
  return s.length > n ? s.slice(0, n - 1) + '…' : s;
}

// ============================================================================
// DOM SCAFFOLD: graph (left) + detail panel (right)
// ============================================================================
var mount = (typeof root !== 'undefined' && root) ? root : document.body;
mount.innerHTML = '';
mount.style.background = T.appBg;
mount.style.fontFamily = "'Inter',system-ui,sans-serif";

const splitWrap = document.createElement('div');
splitWrap.style.cssText = 'width:100%;height:100%;display:flex;flex-direction:row;overflow:hidden;';
mount.appendChild(splitWrap);

const leftPane = document.createElement('div');
leftPane.style.cssText = 'flex:1 1 auto;min-width:0;height:100%;display:flex;flex-direction:column;background:' + T.netCanvasBg + ';';
splitWrap.appendChild(leftPane);

const divider = document.createElement('div');
divider.style.cssText = 'width:1px;background:' + T.splitBorder + ';flex-shrink:0;';
splitWrap.appendChild(divider);

const panel = document.createElement('div');
panel.style.cssText = 'flex:0 0 280px;width:280px;height:100%;overflow:auto;background:' + T.panelBg +
  ';color:' + T.textPrimary + ';padding:14px 16px;font-size:12px;box-sizing:border-box;';
splitWrap.appendChild(panel);

const svgWrap = document.createElement('div');
svgWrap.style.cssText = 'flex:1;position:relative;overflow:hidden;min-height:0;background:' + T.netCanvasBg + ';';
leftPane.appendChild(svgWrap);

const toolbar = document.createElement('div');
toolbar.style.cssText = 'display:flex;align-items:center;gap:12px;padding:8px 14px;background:' + T.toolbarBg +
  ';border-top:1px solid ' + T.toolbarBorder + ';flex-shrink:0;flex-wrap:wrap;';
leftPane.appendChild(toolbar);

const titleEl = document.createElement('span');
titleEl.style.cssText = 'font-weight:700;font-size:13px;color:' + T.titleColor + ';letter-spacing:.02em;';
toolbar.appendChild(titleEl);

const legendWrap = document.createElement('div');
legendWrap.style.cssText = 'display:flex;align-items:center;gap:10px;flex-wrap:wrap;font-size:11px;color:' + T.textMuted + ';';
toolbar.appendChild(legendWrap);

const spacer = document.createElement('span');
spacer.style.cssText = 'margin-left:auto;';
toolbar.appendChild(spacer);

const filterInput = document.createElement('input');
filterInput.type = 'search';
filterInput.placeholder = 'Filter nodes…';
filterInput.style.cssText = 'background:' + T.selectBg + ';color:' + T.textPrimary + ';border:1px solid ' +
  T.selectBorder + ';border-radius:6px;padding:4px 8px;font-size:12px;width:140px;font-family:inherit;';
toolbar.appendChild(filterInput);

const layoutLabel = document.createElement('label');
layoutLabel.textContent = 'Layout:';
layoutLabel.style.cssText = 'font-size:12px;color:' + T.textMuted + ';';
toolbar.appendChild(layoutLabel);

const layoutSelect = document.createElement('select');
layoutSelect.style.cssText = 'background:' + T.selectBg + ';color:' + T.textPrimary + ';border:1px solid ' +
  T.selectBorder + ';border-radius:6px;padding:4px 8px;font-size:12px;cursor:pointer;';
['Force-directed', 'Circular'].forEach(l => {
  const o = document.createElement('option');
  o.value = l; o.textContent = l;
  layoutSelect.appendChild(o);
});
toolbar.appendChild(layoutSelect);

const tooltip = document.createElement('div');
tooltip.style.cssText = 'position:absolute;pointer-events:none;opacity:0;background:' + T.tooltipBg + ';border:1px solid ' +
  T.tooltipBorder + ';border-radius:8px;padding:8px 12px;font-size:11px;color:' + T.textPrimary +
  ';max-width:260px;box-shadow:0 4px 16px rgba(15,23,42,0.18);transition:opacity .12s;z-index:10;';
svgWrap.appendChild(tooltip);

// ============================================================================
// STATE
// ============================================================================
let netData   = null;   // { name, nodes, edges }
let sigIndex  = null;   // signature -> {color, shape, count}
let simulation = null;
let svgSel    = null;
let selection = null;   // {kind:'node'|'edge', id}
let applyFilter = () => {};

// ============================================================================
// DETAIL PANEL
// ============================================================================
function panelPlaceholder(msg) {
  panel.innerHTML = '<div style="color:' + T.textMuted2 + ';line-height:1.5;">' + esc(msg) + '</div>';
}

function panelHeading(text) {
  return '<div style="font-weight:700;font-size:11px;letter-spacing:.08em;text-transform:uppercase;color:' +
    T.titleColor + ';margin-bottom:8px;">' + esc(text) + '</div>';
}

function panelRow(k, v) {
  return '<div style="margin-bottom:6px;">' +
    '<div style="color:' + T.textMuted2 + ';font-size:10px;">' + esc(k) + '</div>' +
    '<div style="color:' + T.textPrimary + ';word-break:break-all;">' + esc(v) + '</div>' +
    '</div>';
}

function showNodeDetail(node) {
  const comps = Array.isArray(node.components) ? node.components : [];
  let html = panelHeading('Chemical System');
  html += panelRow('name', node.name || '(unnamed)');
  html += panelRow('key', node.id);
  html += panelRow('signature', signatureOf(node));
  html += '<div style="margin-top:12px;font-weight:700;color:' + T.textMuted + ';">Components (' + comps.length + ')</div>';
  if (!comps.length) {
    html += '<div style="color:' + T.textMuted2 + ';margin-top:4px;">none</div>';
  } else {
    html += '<div style="margin-top:6px;">';
    comps.forEach(c => {
      html += '<div style="padding:6px 8px;margin-bottom:4px;border-left:2px solid ' +
        colorForSig(signatureOf(node)) + ';background:rgba(127,127,127,0.08);border-radius:0 4px 4px 0;">' +
        '<div style="font-weight:600;">' + esc(c.label) + '</div>' +
        '<div style="color:' + T.textMuted2 + ';font-size:10px;">' + esc(c.type) +
        (c.name ? ' — ' + esc(c.name) : '') + '</div>' +
        '</div>';
    });
    html += '</div>';
  }
  panel.innerHTML = html;
}

function showEdgeDetail(edge) {
  const byId = netData._byId;
  const s = byId[edge.source && edge.source.id ? edge.source.id : edge.source];
  const t = byId[edge.target && edge.target.id ? edge.target.id : edge.target];
  let html = panelHeading('Transformation');
  html += panelRow('name', edge.name || '(unnamed)');
  html += panelRow('key', edge.id);
  html += panelRow('protocol', edge.protocol || '(unknown)');
  html += '<div style="margin-top:12px;font-weight:700;color:' + T.textMuted + ';">Endpoints</div>';
  html += '<div style="margin-top:6px;">' +
    panelRow('state A (source)', s ? nodeLabel(s) : String(edge.source)) +
    panelRow('state B (target)', t ? nodeLabel(t) : String(edge.target)) +
    '</div>';
  panel.innerHTML = html;
}

// ============================================================================
// SIGNATURES + LEGEND
// ============================================================================
function buildSignatureIndex(nodes) {
  const counts = new Map();
  nodes.forEach(n => {
    const sig = signatureOf(n);
    counts.set(sig, (counts.get(sig) || 0) + 1);
  });
  // Most common signature first, so the dominant leg gets the primary colour.
  const sorted = Array.from(counts.entries()).sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0]));
  const index = new Map();
  sorted.forEach(([sig, count], i) => {
    index.set(sig, {
      color: SIG_COLORS[i % SIG_COLORS.length],
      shape: SIG_SHAPES[i % SIG_SHAPES.length],
      count
    });
  });
  return index;
}

function colorForSig(sig) {
  const e = sigIndex && sigIndex.get(sig);
  return e ? e.color : T.textMuted;
}

function shapeForSig(sig) {
  const e = sigIndex && sigIndex.get(sig);
  return e ? e.shape : d3.symbolCircle;
}

function renderLegend() {
  legendWrap.innerHTML = '';
  if (!sigIndex) return;
  sigIndex.forEach((entry, sig) => {
    const item = document.createElement('span');
    item.style.cssText = 'display:inline-flex;align-items:center;gap:5px;';
    const glyph = d3.create('svg').attr('width', 14).attr('height', 14);
    glyph.append('path')
      .attr('transform', 'translate(7,7)')
      .attr('d', d3.symbol().type(entry.shape).size(70)())
      .attr('fill', entry.color);
    item.appendChild(glyph.node());
    const txt = document.createElement('span');
    txt.textContent = sig + ' (' + entry.count + ')';
    item.appendChild(txt);
    legendWrap.appendChild(item);
  });
}

// ============================================================================
// GRAPH RENDER
// ============================================================================
function renderNetwork(data, layout) {
  if (svgSel) { svgSel.remove(); svgSel = null; }
  if (simulation) { simulation.stop(); simulation = null; }

  const W  = svgWrap.getBoundingClientRect().width  || 800;
  const H  = svgWrap.getBoundingClientRect().height || 600;
  const CX = W / 2, CY = H / 2;

  svgSel = d3.select(svgWrap).append('svg')
    .attr('width', W).attr('height', H).style('display', 'block');

  const gRoot = svgSel.append('g');
  svgSel.call(d3.zoom().scaleExtent([0.1, 5]).on('zoom', e => gRoot.attr('transform', e.transform)));

  svgSel.append('defs').append('marker')
    .attr('id', 'an-arrow')
    .attr('viewBox', '0 -5 10 10')
    .attr('refX', CONFIG.edge.arrowRefX)
    .attr('refY', 0)
    .attr('markerUnits', 'userSpaceOnUse')
    .attr('markerWidth',  CONFIG.edge.arrowWidthPx)
    .attr('markerHeight', CONFIG.edge.arrowHeightPx)
    .attr('orient', 'auto')
    .append('path').attr('d', 'M0,-5L10,0L0,5').attr('fill', T.netEdgeLine);

  // Working copies (d3 mutates node objects with x/y/vx/vy).
  const nodes = data.nodes.map(n => ({ ...n, _sig: signatureOf(n), _label: nodeLabel(n) }));
  const nodeById = Object.fromEntries(nodes.map(n => [n.id, n]));
  const edges = data.edges
    .filter(e => nodeById[e.source] && nodeById[e.target])
    .map(e => ({ ...e, source: nodeById[e.source], target: nodeById[e.target] }));

  if (layout === 'Circular') {
    const r = Math.min(W, H) * CONFIG.circular.radiusFraction;
    nodes.forEach((n, i) => {
      const a = (2 * Math.PI * i / nodes.length) - Math.PI / 2;
      n.x = CX + r * Math.cos(a); n.y = CY + r * Math.sin(a);
      n.fx = n.x; n.fy = n.y;
    });
  } else {
    nodes.forEach(n => { n.fx = null; n.fy = null; });
  }

  // ─── Edges ───
  const linkHalo = gRoot.append('g').selectAll('line').data(edges).join('line')
    .attr('stroke', T.netHaloColor)
    .attr('stroke-width', CONFIG.edge.haloWidth)
    .attr('stroke-linecap', 'round')
    .attr('opacity', 0)
    .style('pointer-events', 'none');

  const link = gRoot.append('g').selectAll('line').data(edges).join('line')
    .attr('stroke', T.netEdgeLine)
    .attr('stroke-width', CONFIG.edge.width)
    .attr('stroke-opacity', CONFIG.edge.opacity)
    .attr('marker-end', 'url(#an-arrow)')
    .style('pointer-events', 'none');

  const linkHit = gRoot.append('g').selectAll('line').data(edges).join('line')
    .attr('stroke', 'transparent')
    .attr('stroke-width', CONFIG.edge.hitWidth)
    .style('cursor', 'pointer')
    .on('mouseenter', (ev, d) => {
      tooltip.style.opacity = '1';
      tooltip.innerHTML =
        '<div style="font-weight:700;color:' + T.titleColor + ';">' + esc(d.name || 'Transformation') + '</div>' +
        '<div style="margin-top:3px;">' + esc(d.source._label) + ' → ' + esc(d.target._label) + '</div>' +
        '<div style="margin-top:3px;color:' + T.textMuted2 + ';">' + esc(d.protocol || '') + '</div>';
    })
    .on('mousemove', ev => {
      tooltip.style.left = (ev.offsetX + 14) + 'px';
      tooltip.style.top  = (ev.offsetY - 10) + 'px';
    })
    .on('mouseleave', () => { tooltip.style.opacity = '0'; })
    .on('click', (ev, d) => {
      ev.stopPropagation();
      tooltip.style.opacity = '0';
      selection = { kind: 'edge', id: d.id };
      showEdgeDetail(d);
      refreshSelection();
    });

  // ─── Nodes ───
  const showLabels = nodes.length <= CONFIG.node.labelMaxNodes;

  const nodeG = gRoot.append('g').selectAll('g').data(nodes).join('g')
    .style('cursor', 'grab')
    .call(d3.drag()
      .on('start', (ev, d) => { d.fx = d.x; d.fy = d.y; })
      .on('drag',  (ev, d) => { d.fx = ev.x; d.fy = ev.y; d.x = ev.x; d.y = ev.y; ticked(); })
      .on('end',   (ev, d) => { d.fx = ev.x; d.fy = ev.y; })
    )
    .on('mouseenter', (ev, d) => {
      tooltip.style.opacity = '1';
      const comps = Array.isArray(d.components) ? d.components : [];
      let html = '<div style="font-weight:700;color:' + T.titleColor + ';">' + esc(d._label) + '</div>';
      if (!comps.length) {
        html += '<div style="color:' + T.textMuted2 + ';margin-top:3px;">no components</div>';
      } else {
        html += '<div style="margin-top:4px;line-height:1.5;">' +
          comps.map(c =>
            esc(c.label) + ': ' + esc(c.type) + (c.name ? ' (' + esc(c.name) + ')' : '')
          ).join('<br>') + '</div>';
      }
      tooltip.innerHTML = html;
    })
    .on('mousemove', ev => {
      tooltip.style.left = (ev.offsetX + 14) + 'px';
      tooltip.style.top  = (ev.offsetY - 10) + 'px';
    })
    .on('mouseleave', () => { tooltip.style.opacity = '0'; })
    .on('click', (ev, d) => {
      ev.stopPropagation();
      selection = { kind: 'node', id: d.id };
      showNodeDetail(d);
      refreshSelection();
    });

  const nodeHalo = nodeG.append('path')
    .attr('d', d => d3.symbol().type(shapeForSig(d._sig)).size(CONFIG.node.size * 2.1)())
    .attr('fill', T.netHaloColor)
    .attr('opacity', 0)
    .style('pointer-events', 'none');

  nodeG.append('path')
    .attr('d', d => d3.symbol().type(shapeForSig(d._sig)).size(CONFIG.node.size)())
    .attr('fill', d => colorForSig(d._sig))
    .attr('stroke', T.netNodeStroke)
    .attr('stroke-width', CONFIG.node.strokeWidth);

  if (showLabels) {
    nodeG.append('text')
      .attr('text-anchor', 'middle')
      .attr('y', CONFIG.node.radius + 12)
      .attr('font-size', CONFIG.node.labelFontSize + 'px')
      .attr('font-weight', '600')
      .attr('fill', T.netNodeLabel)
      .attr('pointer-events', 'none')
      .text(d => truncate(d._label, CONFIG.node.labelMaxChars));
  }

  function ticked() {
    link
      .attr('x1', d => d.source.x).attr('y1', d => d.source.y)
      .attr('x2', d => d.target.x).attr('y2', d => d.target.y);
    linkHalo
      .attr('x1', d => d.source.x).attr('y1', d => d.source.y)
      .attr('x2', d => d.target.x).attr('y2', d => d.target.y);
    linkHit
      .attr('x1', d => d.source.x).attr('y1', d => d.source.y)
      .attr('x2', d => d.target.x).attr('y2', d => d.target.y);
    nodeG.attr('transform', d => 'translate(' + (d.x ?? CX) + ',' + (d.y ?? CY) + ')');
  }

  function refreshSelection() {
    nodeHalo.attr('opacity', d =>
      selection && selection.kind === 'node' && selection.id === d.id ? 0.9 : 0);
    linkHalo.attr('opacity', d =>
      selection && selection.kind === 'edge' && selection.id === d.id ? 0.9 : 0);
  }

  applyFilter = function() {
    const q = filterInput.value.trim().toLowerCase();
    if (!q) {
      nodeG.attr('opacity', 1);
      link.attr('stroke-opacity', CONFIG.edge.opacity);
      return;
    }
    const matches = new Set();
    nodes.forEach(n => {
      const hay = (n._label + ' ' + n.id + ' ' + n._sig).toLowerCase();
      if (hay.indexOf(q) !== -1) matches.add(n.id);
    });
    nodeG.attr('opacity', d => matches.has(d.id) ? 1 : CONFIG.dimOpacity);
    link.attr('stroke-opacity', d =>
      (matches.has(d.source.id) || matches.has(d.target.id))
        ? CONFIG.edge.opacity : CONFIG.dimOpacity);
  };

  // ─── Layout ───
  if (layout !== 'Circular') {
    const fc = CONFIG.force;
    simulation = d3.forceSimulation(nodes)
      .force('link', d3.forceLink(edges).id(d => d.id)
        .distance(fc.linkDistance).strength(fc.linkStrength))
      .force('charge', d3.forceManyBody()
        .strength(fc.chargeStrength).distanceMax(fc.chargeDistanceMax))
      .force('center', d3.forceCenter(CX, CY).strength(fc.centerStrength))
      .force('collision', d3.forceCollide(CONFIG.node.radius + fc.collisionPadding)
        .iterations(fc.collisionIterations))
      .force('x', d3.forceX(CX).strength(fc.driftX))
      .force('y', d3.forceY(CY).strength(fc.driftY))
      .stop();

    // Run the simulation to completion up front rather than animating: with a
    // few hundred nodes a per-frame DOM update is what melts the browser.
    const need = Math.ceil(Math.log(simulation.alphaMin()) / Math.log(1 - simulation.alphaDecay()));
    const n = Math.min(need, CONFIG.force.maxTicks);
    for (let i = 0; i < n; i++) simulation.tick();
  }

  ticked();
  refreshSelection();
  applyFilter();

  // Fit the drawn graph into view.
  const xs = nodes.map(n => n.x), ys = nodes.map(n => n.y);
  if (xs.length) {
    const pad = 60;
    const minX = Math.min(...xs) - pad, maxX = Math.max(...xs) + pad;
    const minY = Math.min(...ys) - pad, maxY = Math.max(...ys) + pad;
    const k = Math.min(1, Math.min(W / (maxX - minX), H / (maxY - minY)));
    const tx = W / 2 - k * (minX + maxX) / 2;
    const ty = H / 2 - k * (minY + maxY) / 2;
    gRoot.attr('transform', 'translate(' + tx + ',' + ty + ') scale(' + k + ')');
  }
}

// ============================================================================
// INIT
// ============================================================================
function showError(msg) {
  const warn = document.createElement('div');
  warn.textContent = '⚠ ' + msg;
  warn.style.cssText = 'position:absolute;top:20px;left:50%;transform:translateX(-50%);background:#fee2e2;color:#991b1b;' +
    'border:1px solid #fecaca;padding:6px 14px;border-radius:6px;font-size:12px;z-index:20;max-width:80%;';
  svgWrap.appendChild(warn);
}

function init(net) {
  if (!net || !Array.isArray(net.nodes)) {
    panelPlaceholder('Waiting for an AlchemicalNetwork…');
    titleEl.textContent = 'Alchemical Network';
    return;
  }
  netData = {
    name:  net.name || null,
    nodes: net.nodes,
    edges: Array.isArray(net.edges) ? net.edges : [],
  };
  netData._byId = Object.fromEntries(netData.nodes.map(n => [n.id, n]));

  sigIndex = buildSignatureIndex(netData.nodes);
  renderLegend();

  titleEl.textContent = (netData.name || 'Alchemical Network') +
    '  ·  ' + netData.nodes.length + ' systems, ' + netData.edges.length + ' transformations';

  selection = null;
  panelPlaceholder('Click a chemical system or a transformation to see its details.');
  renderNetwork(netData, layoutSelect.value);
}

layoutSelect.addEventListener('change', () => {
  if (netData) renderNetwork(netData, layoutSelect.value);
});
filterInput.addEventListener('input', () => applyFilter());

panelPlaceholder('Waiting for an AlchemicalNetwork…');
titleEl.textContent = 'Alchemical Network';

// ============================================================================
// Module entry points
// ============================================================================
export async function onInputs(inputs) {
  if (!inputs) return;
  try {
    const net = await asObject(inputs['alchemical_network']);
    if (!net) return;
    init(net);
  } catch (e) {
    if (typeof logStderr === 'function') logStderr('alchemical_network: ' + e.message);
    showError('Failed to read alchemical_network payload: ' + e.message);
  }
}

export function onResize() {
  if (netData) renderNetwork(netData, layoutSelect.value);
}

export function cleanup() {
  if (simulation) { simulation.stop(); simulation = null; }
}
