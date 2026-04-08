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

/**
 * Phase 1 step 6: Hydrogen bonding between complementary RNA strands (R3).
 *
 * For each cell with ≥2 unbonded chains, attempt antiparallel pairing.
 * Pair score = number of complementary positions when one strand is reversed.
 * Pairs form with probability rates.hBondForm * (1 + score/L).
 * Existing pairs break each tick with probability hBondBreakAU (weighted down
 * by GC content of the pairing — GC bonds are ~2x stronger).
 */
export function processHydrogenBonds(world, rates) {
  const W = world.width;
  const structures = world.structures;

  // ─── Step 1: break existing pairs ───
  // We walk pairs once. The lower-id chain in a pair owns the break decision
  // (so we don't double-roll).
  const pairBreakBaseAU = rates.hBondBreakAU;
  const pairBreakBaseGC = rates.hBondBreakGC;
  if (pairBreakBaseAU > 0.05 || pairBreakBaseGC > 0.05) {
    for (const a of structures) {
      if (a.type !== 'rna' || a.hBondedTo == null) continue;
      const b = _findById(structures, a.hBondedTo);
      if (!b || b.id < a.id) continue;
      // Pair score and GC fraction
      const { score, gcFrac } = _pairScore(a, b);
      // Effective break rate: weighted average of AU and GC break rates
      const breakRate = pairBreakBaseAU * (1 - gcFrac) + pairBreakBaseGC * gcFrac;
      // Stability inversely proportional to score
      const effective = breakRate / Math.max(1, score * 0.5);
      if (world._rand() < effective) {
        a.hBondedTo = null;
        b.hBondedTo = null;
      }
    }
  }

  // ─── Step 2: form new pairs ───
  if (rates.hBondForm < 0.05) return;
  // Index chains by cell
  const cellMap = new Map();
  for (const st of structures) {
    if (st.type !== 'rna') continue;
    if (st.hBondedTo != null) continue;
    const k = st.position.y * W + st.position.x;
    if (!cellMap.has(k)) cellMap.set(k, []);
    cellMap.get(k).push(st);
  }
  for (const list of cellMap.values()) {
    if (list.length < 2) continue;
    for (let i = 0; i < list.length - 1; i++) {
      const a = list[i];
      if (a.hBondedTo != null) continue;
      for (let j = i + 1; j < list.length; j++) {
        const b = list[j];
        if (b.hBondedTo != null) continue;
        if (a.sequence.length < 2 || b.sequence.length < 2) continue;
        const { score } = _pairScore(a, b);
        if (score === 0) continue;
        const minL = Math.min(a.sequence.length, b.sequence.length);
        const matchFrac = score / minL;
        // Probability scales with rates.hBondForm and match fraction
        if (matchFrac < 0.4) continue;  // need ≥40% complementary to nucleate
        if (world._rand() < rates.hBondForm * matchFrac) {
          a.hBondedTo = b.id;
          b.hBondedTo = a.id;
          break;
        }
      }
    }
  }
}

// Score the antiparallel pairing of two chains. Returns {score, gcFrac}.
// Aligns end-to-end (longer chain may "overhang" — only matching length scored).
function _pairScore(a, b) {
  const sa = a.sequence;
  const sb = b.sequence;
  const minL = Math.min(sa.length, sb.length);
  let score = 0;
  let gcMatches = 0;
  for (let i = 0; i < minL; i++) {
    const x = sa[i];
    const y = sb[sb.length - 1 - i];  // antiparallel: read b in reverse
    if (
      (x === 'A' && y === 'U') ||
      (x === 'U' && y === 'A')
    ) {
      score++;
    } else if (
      (x === 'G' && y === 'C') ||
      (x === 'C' && y === 'G')
    ) {
      score++;
      gcMatches++;
    }
  }
  return { score, gcFrac: score === 0 ? 0 : gcMatches / score };
}

function _findById(structures, id) {
  for (const s of structures) if (s.id === id) return s;
  return null;
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
