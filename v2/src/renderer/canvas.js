// Canvas renderer for the surface view.
//
// Phase 1 step 1+3: nucleotide concentrations as additive RGB +
// background temperature tint + translucent water overlay driven by waterLevel.
//
// We render the (W, H) field to an offscreen ImageData buffer (one pixel per
// grid cell), then drawImage upscales it to the visible canvas.

import { SPECIES_INDEX, SPECIES_COLOR } from '../simulation/constants.js';
import { getTemperature, getWaterLevel } from '../simulation/environment.js';

// Cached ImageData buffer (recreated when canvas size changes)
let cached = null;

export function renderWorld(canvas, world, cellPx) {
  const ctx = canvas.getContext('2d', { alpha: false });
  if (!ctx) return;

  const W = world.width;
  const H = world.height;
  const cw = W * cellPx;
  const ch = H * cellPx;

  if (canvas.width !== cw) canvas.width = cw;
  if (canvas.height !== ch) canvas.height = ch;

  // We render the W×H low-res image first, then upscale via drawImage.
  if (!cached || cached.width !== W || cached.height !== H) {
    const off = document.createElement('canvas');
    off.width = W;
    off.height = H;
    cached = {
      width: W,
      height: H,
      offCanvas: off,
      offCtx: off.getContext('2d', { alpha: false }),
      imgData: off.getContext('2d').createImageData(W, H),
    };
  }
  const { offCanvas, offCtx, imgData } = cached;

  const fA = world.fields[SPECIES_INDEX.A];
  const fU = world.fields[SPECIES_INDEX.U];
  const fG = world.fields[SPECIES_INDEX.G];
  const fC = world.fields[SPECIES_INDEX.C];

  const cA = SPECIES_COLOR.A;
  const cU = SPECIES_COLOR.U;
  const cG = SPECIES_COLOR.G;
  const cC = SPECIES_COLOR.C;

  // Phase 1 step 3: temperature tints background; water level adds blue veil.
  const temperature = getTemperature(world.tickCount);
  const waterLevel = getWaterLevel(world.tickCount);
  const bgR = 12 + Math.max(0, temperature) * 50;
  const bgG = 12 + Math.max(0, -temperature) * 5;
  const bgB = 18 + Math.max(0, -temperature) * 50;

  const data = imgData.data;
  const N = W * H;
  for (let k = 0; k < N; k++) {
    // Each species adds to RGB scaled by its concentration.
    // Concentration ~ [0, ~1.5]. Multiply by 1.5 to get reasonable brightness.
    const a = Math.min(1, fA[k] * 1.5);
    const u = Math.min(1, fU[k] * 1.5);
    const g = Math.min(1, fG[k] * 1.5);
    const c = Math.min(1, fC[k] * 1.5);

    let r = bgR + a * cA[0] + u * cU[0] + g * cG[0] + c * cC[0];
    let gg = bgG + a * cA[1] + u * cU[1] + g * cG[1] + c * cC[1];
    let b = bgB + a * cA[2] + u * cU[2] + g * cG[2] + c * cC[2];

    if (r > 255) r = 255;
    if (gg > 255) gg = 255;
    if (b > 255) b = 255;

    const o = k * 4;
    data[o] = r;
    data[o + 1] = gg;
    data[o + 2] = b;
    data[o + 3] = 255;
  }

  offCtx.putImageData(imgData, 0, 0);

  // Upscale (nearest-neighbor via image-rendering: pixelated on the canvas element)
  ctx.imageSmoothingEnabled = false;
  ctx.drawImage(offCanvas, 0, 0, cw, ch);

  // Translucent water overlay: stronger when waterLevel is positive
  const wet = Math.max(0, waterLevel);
  if (wet > 0) {
    ctx.fillStyle = `rgba(60, 130, 220, ${0.05 + wet * 0.18})`;
    ctx.fillRect(0, 0, cw, ch);
  }
  // Dry "exposed surface" tint when waterLevel is negative
  const dry = Math.max(0, -waterLevel);
  if (dry > 0) {
    ctx.fillStyle = `rgba(190, 130, 70, ${dry * 0.07})`;
    ctx.fillRect(0, 0, cw, ch);
  }

  // Phase 1 step 6: draw H-bond connections (dotted lines between paired strands).
  // Drawn FIRST so chain dots overlap them.
  drawHBonds(ctx, world, cellPx);

  // Phase 1 step 5: draw RNA chains as dots over the field.
  drawRnaChains(ctx, world, cellPx);
}

function drawHBonds(ctx, world, cellPx) {
  ctx.save();
  ctx.strokeStyle = 'rgba(150, 220, 255, 0.55)';
  ctx.lineWidth = 1;
  ctx.setLineDash([2, 2]);
  const drawn = new Set();
  for (const a of world.structures) {
    if (a.type !== 'rna' || a.hBondedTo == null) continue;
    if (drawn.has(a.id)) continue;
    let b = null;
    for (const s of world.structures) {
      if (s.id === a.hBondedTo) { b = s; break; }
    }
    if (!b) continue;
    drawn.add(a.id);
    drawn.add(b.id);
    const ax = a.position.x * cellPx + cellPx / 2;
    const ay = a.position.y * cellPx + cellPx / 2;
    const bx = b.position.x * cellPx + cellPx / 2;
    const by = b.position.y * cellPx + cellPx / 2;
    ctx.beginPath();
    ctx.moveTo(ax, ay);
    ctx.lineTo(bx, by);
    ctx.stroke();
  }
  ctx.restore();
}

function drawRnaChains(ctx, world, cellPx) {
  for (const st of world.structures) {
    if (st.type !== 'rna') continue;
    const L = st.sequence.length;
    const px = st.position.x * cellPx;
    const py = st.position.y * cellPx;

    // Color: shift from cyan (high AU) to magenta (high GC) based on composition
    let gc = 0;
    for (const b of st.sequence) if (b === 'G' || b === 'C') gc++;
    const gcRatio = gc / L;
    const r = Math.floor(120 + 120 * gcRatio);
    const g = Math.floor(220 - 100 * gcRatio);
    const b = Math.floor(220);

    // Size grows with chain length: 2-mer ≈ cellPx; longer ≈ cellPx * (1 + log)
    const size = Math.max(cellPx, Math.min(cellPx * 4, cellPx * (1 + Math.log2(L))));
    const off = (size - cellPx) / 2;

    // Halo for surface-bound (slightly brighter)
    if (st.surfaceBound) {
      ctx.fillStyle = `rgba(${r}, ${g}, ${b}, 0.95)`;
    } else {
      ctx.fillStyle = `rgba(${r}, ${g}, ${b}, 0.70)`;
    }
    ctx.fillRect(px - off, py - off, size, size);

    // For long chains, add a white core
    if (L >= 6) {
      ctx.fillStyle = 'rgba(255,255,255,0.7)';
      ctx.fillRect(px + cellPx * 0.25, py + cellPx * 0.25, cellPx * 0.5, cellPx * 0.5);
    }
  }
}
