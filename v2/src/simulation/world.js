// World — grid + concentration field + tick loop.
//
// The concentration field is one Float32Array per species, of length W*H.
// Index = y * W + x. This is what makes the renderer + diffusion fast.
//
// Phase 1 step 1: only the field exists. Diffusion / cycles / chains
// are added in subsequent steps.

import { N_SPECIES, SPECIES, SPECIES_INDEX, INITIAL_CONCENTRATION } from './constants.js';
import { getTemperature, getWaterLevel } from './environment.js';

export class World {
  constructor({ width = 200, height = 200, seed = 1 } = {}) {
    this.width = width;
    this.height = height;
    this.size = width * height;
    this.tickCount = 0;

    // One Float32Array per species. fields[i] is the concentration grid for SPECIES[i].
    this.fields = new Array(N_SPECIES);
    for (let i = 0; i < N_SPECIES; i++) {
      this.fields[i] = new Float32Array(this.size);
    }

    // Structures (RNA chains, peptides, lipids, membranes) live here in later phases.
    this.structures = [];

    // Seeded RNG (mulberry32) for reproducible runs
    this._seed = seed >>> 0;

    this._initialFill();
  }

  _initialFill() {
    // Phase 1: a uniform random scatter of monomers across the surface so we
    // can see something on screen even before diffusion is wired in.
    for (let i = 0; i < N_SPECIES; i++) {
      const f = this.fields[i];
      for (let k = 0; k < this.size; k++) {
        f[k] = INITIAL_CONCENTRATION * this._rand();
      }
    }
  }

  // Mulberry32 — small fast deterministic PRNG
  _rand() {
    let t = (this._seed += 0x6D2B79F5);
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  }

  tick() {
    this.tickCount++;
    // Phase 1 step 1: tick is a no-op (concentrations exist but don't move).
    // Diffusion arrives in step 2; environment cycles in step 3.
  }

  // ─── helpers ───

  index(x, y) {
    return y * this.width + x;
  }

  getConcentration(speciesName, x, y) {
    const idx = SPECIES_INDEX[speciesName];
    return this.fields[idx][this.index(x, y)];
  }

  setConcentration(speciesName, x, y, value) {
    const idx = SPECIES_INDEX[speciesName];
    this.fields[idx][this.index(x, y)] = value;
  }

  getStats() {
    const temperature = getTemperature(this.tickCount);
    const waterLevel = getWaterLevel(this.tickCount);

    // Mean concentration per species (cheap reduction over the field)
    const concentrations = {};
    for (let i = 0; i < N_SPECIES; i++) {
      const f = this.fields[i];
      let sum = 0;
      for (let k = 0; k < this.size; k++) sum += f[k];
      concentrations[SPECIES[i]] = sum / this.size;
    }

    return {
      tick: this.tickCount,
      temperature,
      waterLevel,
      concentrations,
      structures: this.structures.length,
    };
  }
}
