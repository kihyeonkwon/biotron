// Chemical reactions on the world.
//
// Each reaction reads/writes the world's concentration fields and structures.
// Reactions are batched per-tick by world.tick().

import { N_SPECIES, SPECIES_INDEX } from './constants.js';
import { makeRNA, NUCLEOTIDES, isComplementary } from './rna.js';

const A_IDX = SPECIES_INDEX.A;
const U_IDX = SPECIES_INDEX.U;
const G_IDX = SPECIES_INDEX.G;
const C_IDX = SPECIES_INDEX.C;
const NUC_IDX = [A_IDX, U_IDX, G_IDX, C_IDX];

const NUCLEOTIDE_THRESHOLD = 0.15;  // local concentration above which polymerization is feasible

/**
 * Phase 1 step 5: non-templated RNA polymerization (R4).
 *
 * Each tick we sample a fraction of the grid cells. In each sampled cell,
 * if local nucleotide concentration is high enough AND warm+dry conditions
 * are met, we either:
 *   (a) create a new 2-mer chain (sampling 2 nucleotides by abundance), or
 *   (b) extend an existing chain in that cell (consume one matching nucleotide).
 *
 * Surface-bound chains polymerize faster (clay catalysis effect).
 */
export function polymerizeRNA(world, rates) {
  const { width: W, height: H, fields, structures } = world;
  const formProb = rates.backboneFormSurface;
  if (formProb < 0.01) return 0;  // no polymerization in cold/wet conditions

  // Build a per-cell index of existing chains (cheap; chain count is small)
  const chainsByCell = new Map();
  for (const st of structures) {
    if (st.type !== 'rna') continue;
    const k = st.position.y * W + st.position.x;
    if (!chainsByCell.has(k)) chainsByCell.set(k, []);
    chainsByCell.get(k).push(st);
  }

  // Sample fraction of cells per tick (constant work per tick)
  const N_SAMPLES = 200;
  let formed = 0;
  for (let s = 0; s < N_SAMPLES; s++) {
    const x = Math.floor(world._rand() * W);
    const y = Math.floor(world._rand() * H);
    const k = y * W + x;

    // Local nucleotide concentration sum
    const cA = fields[A_IDX][k];
    const cU = fields[U_IDX][k];
    const cG = fields[G_IDX][k];
    const cC = fields[C_IDX][k];
    const totN = cA + cU + cG + cC;
    if (totN < NUCLEOTIDE_THRESHOLD) continue;

    // Existing chain(s) in this cell?
    const cellChains = chainsByCell.get(k);
    if (cellChains && cellChains.length > 0) {
      // Try to extend the first one
      const chain = cellChains[0];
      const ext = formProb * (chain.surfaceBound ? 1.0 : 0.4);
      if (world._rand() < ext) {
        const nuc = _pickNucleotide(world, k, fields);
        if (nuc) {
          chain.sequence.push(nuc);
          formed++;
        }
      }
    } else {
      // Spontaneous dimer formation
      if (world._rand() < formProb * 0.4) {
        const n1 = _pickNucleotide(world, k, fields);
        if (!n1) continue;
        const n2 = _pickNucleotide(world, k, fields);
        if (!n2) continue;
        const chain = makeRNA([n1, n2], x, y, true);
        structures.push(chain);
        formed++;
      }
    }
  }
  return formed;
}

// Sample one nucleotide from a cell weighted by its concentration; consume it.
function _pickNucleotide(world, k, fields) {
  const cA = fields[A_IDX][k];
  const cU = fields[U_IDX][k];
  const cG = fields[G_IDX][k];
  const cC = fields[C_IDX][k];
  const tot = cA + cU + cG + cC;
  if (tot <= 0) return null;
  let r = world._rand() * tot;
  if ((r -= cA) < 0) {
    fields[A_IDX][k] = Math.max(0, cA - 0.05);
    return 'A';
  }
  if ((r -= cU) < 0) {
    fields[U_IDX][k] = Math.max(0, cU - 0.05);
    return 'U';
  }
  if ((r -= cG) < 0) {
    fields[G_IDX][k] = Math.max(0, cG - 0.05);
    return 'G';
  }
  fields[C_IDX][k] = Math.max(0, cC - 0.05);
  return 'C';
}
