// World — grid + concentration field + tick loop.
//
// The concentration field is one Float32Array per species, of length W*H.
// Index = y * W + x. This is what makes the renderer + diffusion fast.
//
// Phase 1 step 1: only the field exists. Diffusion / cycles / chains
// are added in subsequent steps.

import { N_SPECIES, SPECIES, SPECIES_INDEX, INITIAL_CONCENTRATION } from './constants.js';
import { getTemperature, getWaterLevel, getReactionRates } from './environment.js';
import {
  polymerizeRNA,
  processHydrogenBonds,
  templatedExtension,
  degradeAndSupply,
  activateRibozymes,
  attachAminoAcids,
  formPeptideBonds,
  nucleateLipids,
  assembleMembranes,
  updateMembranes,
  growMembranes,
  divideMembranes,
} from './reactions.js';
import { resetRNAIds } from './rna.js';
import { resetLipidIds } from './lipid.js';
import { createMilestoneTracker } from './metrics.js';

export class World {
  constructor({ width = 200, height = 200, seed = 1 } = {}) {
    resetRNAIds();
    resetLipidIds();
    this.width = width;
    this.height = height;
    this.size = width * height;
    this.tickCount = 0;

    // One Float32Array per species. fields[i] is the concentration grid for SPECIES[i].
    this.fields = new Array(N_SPECIES);
    // Double-buffer used by the diffusion step (avoids per-tick allocation).
    this._next = new Array(N_SPECIES);
    for (let i = 0; i < N_SPECIES; i++) {
      this.fields[i] = new Float32Array(this.size);
      this._next[i] = new Float32Array(this.size);
    }

    // Structures (RNA chains, peptides, lipids, membranes) live here in later phases.
    this.structures = [];

    // Milestone tracker (Phase 2 step 12)
    this.milestones = createMilestoneTracker();

    // Phase 4: per-cell membrane id (Int32, -1 = no membrane)
    this._cellMembrane = new Int32Array(this.size).fill(-1);

    // Seeded RNG (mulberry32) for reproducible runs
    this._seed = seed >>> 0;

    this._initialFill();
  }

  _initialFill() {
    // Phase 1: low uniform random base + a few high-concentration "drops" of
    // each nucleotide. The drops let you visually see the diffusion step
    // smear them out into smooth gradients over the first few hundred ticks.
    for (let i = 0; i < N_SPECIES; i++) {
      const f = this.fields[i];
      for (let k = 0; k < this.size; k++) {
        f[k] = 0.3 * INITIAL_CONCENTRATION * this._rand();
      }
    }

    // Drop one bright blob of each nucleotide (A, U, G, C) at four corners
    // of the grid so diffusion is immediately visible.
    const blobs = [
      ['A', this.width * 0.25, this.height * 0.25],
      ['U', this.width * 0.75, this.height * 0.25],
      ['G', this.width * 0.25, this.height * 0.75],
      ['C', this.width * 0.75, this.height * 0.75],
    ];
    for (const [name, cx, cy] of blobs) {
      this._addBlob(name, cx, cy, this.width * 0.05, 1.5);
    }

    // Phase 4: scatter some FA blobs so lipids have something to nucleate from.
    this._addBlob('FA', this.width * 0.5, this.height * 0.3, this.width * 0.06, 1.5);
    this._addBlob('FA', this.width * 0.3, this.height * 0.7, this.width * 0.06, 1.5);
    this._addBlob('FA', this.width * 0.7, this.height * 0.7, this.width * 0.06, 1.5);

    // Bonus seed: a few free amino acids so translation has something to chew on.
    for (const aa of ['Gly', 'Ala', 'Val', 'Asp', 'Glu']) {
      this._addBlob(aa, this.width * (0.2 + 0.6 * this._rand()),
                    this.height * (0.2 + 0.6 * this._rand()),
                    this.width * 0.04, 1.0);
    }
  }

  _addBlob(speciesName, cx, cy, radius, peak) {
    const idx = SPECIES_INDEX[speciesName];
    const f = this.fields[idx];
    const r2 = radius * radius;
    const ix = Math.floor(cx);
    const iy = Math.floor(cy);
    const r = Math.ceil(radius);
    for (let dy = -r; dy <= r; dy++) {
      const y = (iy + dy + this.height) % this.height;
      for (let dx = -r; dx <= r; dx++) {
        const x = (ix + dx + this.width) % this.width;
        const d2 = dx * dx + dy * dy;
        if (d2 > r2) continue;
        const fall = 1 - d2 / r2;
        f[y * this.width + x] += peak * fall * fall;
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
    const temperature = getTemperature(this.tickCount);
    const waterLevel = getWaterLevel(this.tickCount);
    const rates = getReactionRates(temperature, waterLevel);
    this._lastTemperature = temperature;
    this._lastWaterLevel = waterLevel;
    this._lastRates = rates;

    // Phase 1 step 2: diffusion (gated off when fully dry)
    if (waterLevel >= -0.8) {
      this._diffuse(rates.diffusion);
    }

    // Phase 1 step 4: evaporative concentration during dry phase.
    // When the surface is exposed (waterLevel < 0), water evaporates from cells
    // and the molecules left behind get more concentrated. This is the chemical
    // basis for surface polymerization (Deamer hot-pool model). Effect is small
    // per tick but accumulates and shows clearly during low tide.
    if (waterLevel < 0) {
      this._evaporativeConcentration(-waterLevel);
    }

    // Phase 1 step 5: non-templated RNA polymerization (R4)
    polymerizeRNA(this, rates);

    // Phase 1 step 6: hydrogen bonding between complementary strands (R3)
    if (this.structures.length > 1) {
      processHydrogenBonds(this, rates);
    }

    // Phase 2 step 8 + 10: templated extension (with mutation) on H-bonded pairs (R4 + R5)
    if (this.structures.length > 1) {
      templatedExtension(this, rates);
    }

    // Phase 3 step 14: ribozyme activation (R9). Refresh every 8 ticks
    // (sequence-changing events are infrequent compared to tick rate).
    if (this.tickCount % 8 === 0 && this.structures.length > 0) {
      activateRibozymes(this);
    }

    // Phase 3 step 13: codon → amino acid attachment (R7)
    if (this.structures.length > 0) {
      attachAminoAcids(this, rates);
    }

    // Phase 3 step 15: peptide bond formation (R8)
    if (this.structures.length > 0) {
      formPeptideBonds(this, rates);
    }

    // Phase 4 steps 17-19: lipid nucleation + membrane self-assembly
    nucleateLipids(this, rates);
    if (this.tickCount % 6 === 0) {
      assembleMembranes(this);
    }

    // Phase 4 step 21: vesicle growth + binary fission (cell reproduction)
    growMembranes(this);
    if (this.tickCount % 12 === 0) {
      divideMembranes(this);
    }

    // Phase 4 step 20+21: membrane state update (enclosed list, integrity, dissolution)
    updateMembranes(this, waterLevel);
    this._rebuildCellMembraneIndex();

    // Phase 2 step 11: degradation + ocean supply (R6)
    degradeAndSupply(this, rates);

    // Phase 1 step 4: structure adsorption / desorption
    if (this.structures.length > 0) {
      this._processAdsorption(rates);
    }

    // Phase 2 step 12: milestone observation (cheap, runs every tick)
    this.milestones.update(this);
  }

  _rebuildCellMembraneIndex() {
    // Mark each grid cell with its enclosing membrane id (-1 if none).
    // For overlapping membranes, the smaller one wins (more recent / inner).
    this._cellMembrane.fill(-1);
    const W = this.width;
    const H = this.height;
    let anyMembrane = false;
    for (const m of this.structures) {
      if (m.type !== 'membrane') continue;
      anyMembrane = true;
      const r = Math.ceil(m.radius);
      const r2 = m.radius * m.radius;
      for (let dy = -r; dy <= r; dy++) {
        for (let dx = -r; dx <= r; dx++) {
          if (dx * dx + dy * dy > r2) continue;
          const x = ((Math.round(m.center.x) + dx) % W + W) % W;
          const y = ((Math.round(m.center.y) + dy) % H + H) % H;
          const k = y * W + x;
          if (this._cellMembrane[k] === -1) this._cellMembrane[k] = m.id;
        }
      }
    }
    this._hasMembraneIndex = anyMembrane;
  }

  _evaporativeConcentration(dryFactor) {
    // Multiply each cell's concentration by a small factor (1 + 0.005 * dryFactor)
    // BUT cap the total per cell to avoid runaway. This models the fact that
    // water evaporates from exposed cells, leaving solutes behind.
    // dryFactor ∈ [0, 1] (1 = fully dry).
    const boost = 1 + 0.005 * dryFactor;
    for (let s = 0; s < N_SPECIES; s++) {
      const f = this.fields[s];
      for (let k = 0; k < this.size; k++) {
        const v = f[k] * boost;
        f[k] = v > 5.0 ? 5.0 : v;  // saturation cap
      }
    }
  }

  _processAdsorption(rates) {
    // Each structure may flip its surfaceBound state based on water level.
    // - Wet conditions favor desorption (structures float into solution).
    // - Dry conditions favor adsorption (structures stick to surface).
    // No-op until structures exist (Phase 1 step 5).
    for (const st of this.structures) {
      if (st.surfaceBound) {
        if (this._rand() < rates.desorption * 0.1) st.surfaceBound = false;
      } else {
        if (this._rand() < rates.adsorption * 0.1) st.surfaceBound = true;
      }
    }
  }

  _diffuse(D) {
    // Stability constraint: alpha must be ≤ 0.25 for stable explicit Laplacian
    const alpha = Math.min(0.24, D * 0.25);
    const W = this.width;
    const H = this.height;
    const cellMem = this._cellMembrane;
    const hasMembranes = this._hasMembraneIndex;
    // Phase 4: monomers crossing a membrane boundary are attenuated by this factor
    // (small molecules pass through but slowly — semi-permeable lipid bilayer)
    const SEAL = 0.30;

    for (let s = 0; s < N_SPECIES; s++) {
      const cur = this.fields[s];
      const nxt = this._next[s];

      for (let y = 0; y < H; y++) {
        const yUp = (y - 1 + H) % H;
        const yDn = (y + 1) % H;
        for (let x = 0; x < W; x++) {
          const xLf = (x - 1 + W) % W;
          const xRt = (x + 1) % W;
          const i = y * W + x;
          const c = cur[i];
          if (hasMembranes) {
            const myM = cellMem[i];
            const fU = (cellMem[yUp * W + x] === myM) ? 1.0 : SEAL;
            const fD = (cellMem[yDn * W + x] === myM) ? 1.0 : SEAL;
            const fL = (cellMem[y * W + xLf] === myM) ? 1.0 : SEAL;
            const fR = (cellMem[y * W + xRt] === myM) ? 1.0 : SEAL;
            const lap =
              fU * (cur[yUp * W + x] - c) +
              fD * (cur[yDn * W + x] - c) +
              fL * (cur[y * W + xLf] - c) +
              fR * (cur[y * W + xRt] - c);
            nxt[i] = c + alpha * lap;
          } else {
            const lap =
              cur[yUp * W + x] +
              cur[yDn * W + x] +
              cur[y * W + xLf] +
              cur[y * W + xRt] -
              4 * c;
            nxt[i] = c + alpha * lap;
          }
        }
      }
      this.fields[s] = nxt;
      this._next[s] = cur;
    }
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

    // Chain length histogram (Phase 1 step 5+)
    const rnaChains = this.structures.filter((s) => s.type === 'rna');
    const peptides = this.structures.filter((s) => s.type === 'peptide');
    const lipids = this.structures.filter((s) => s.type === 'lipid');
    const membranes = this.structures.filter((s) => s.type === 'membrane');
    const maxGen = membranes.reduce((m, v) => Math.max(m, v.generation), 0);
    const totalChildren = membranes.reduce((m, v) => m + v.childCount, 0);
    const lenHisto = {};
    for (const ch of rnaChains) {
      const L = ch.sequence.length;
      lenHisto[L] = (lenHisto[L] || 0) + 1;
    }
    const maxLen = rnaChains.reduce((m, ch) => Math.max(m, ch.sequence.length), 0);
    const hBondedCount = rnaChains.filter((c) => c.hBondedTo != null).length;
    const ribozymes = rnaChains.filter((c) => c.catalyticFunction != null);
    const ribByType = { peptidyl_transferase: 0, rna_replicase: 0, aminoacyl_transferase: 0 };
    for (const r of ribozymes) {
      for (const f of r.catalyticFunction) ribByType[f.type]++;
    }

    return {
      tick: this.tickCount,
      temperature,
      waterLevel,
      concentrations,
      structures: this.structures.length,
      rnaChains: rnaChains.length,
      maxRnaLen: maxLen,
      lenHisto,
      hBondedChains: hBondedCount,
      ribozymes: ribozymes.length,
      ribByType,
      peptides: peptides.length,
      lipids: lipids.length,
      membranes: membranes.length,
      maxVesicleGen: maxGen,
      vesicleDivisions: totalChildren,
      milestoneStatus: this.milestones.getStatus(),
    };
  }
}
