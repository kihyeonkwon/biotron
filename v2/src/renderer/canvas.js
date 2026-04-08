// Canvas renderer for the surface view.
//
// Phase 1 step 1: just draws nucleotide concentrations as additive RGB.
// A is red, U is cyan, G is green, C is yellow. Higher concentration = brighter.
//
// We render straight to an ImageData buffer for speed (no per-pixel canvas calls).

import { SPECIES_INDEX, SPECIES_COLOR } from '../simulation/constants.js';

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

  const data = imgData.data;
  const N = W * H;
  for (let k = 0; k < N; k++) {
    // Each species adds to RGB scaled by its concentration.
    // Concentration ~ [0, ~0.3]. Multiply by 4 to get reasonable brightness.
    const a = Math.min(1, fA[k] * 4);
    const u = Math.min(1, fU[k] * 4);
    const g = Math.min(1, fG[k] * 4);
    const c = Math.min(1, fC[k] * 4);

    let r = a * cA[0] + u * cU[0] + g * cG[0] + c * cC[0];
    let gg = a * cA[1] + u * cU[1] + g * cG[1] + c * cC[1];
    let b = a * cA[2] + u * cU[2] + g * cG[2] + c * cC[2];

    // Soft saturate
    r = r > 255 ? 255 : r;
    gg = gg > 255 ? 255 : gg;
    b = b > 255 ? 255 : b;

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
}
