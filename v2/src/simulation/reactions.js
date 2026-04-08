// Chemical reactions on the world.
//
// Each reaction reads/writes the world's concentration fields and structures.
// Reactions are batched per-tick by world.tick().

import { N_SPECIES, SPECIES_INDEX } from './constants.js';
import {
  makeRNA, NUCLEOTIDES, COMPLEMENT, isComplementary, selfComplementarity,
  detectRibozyme, findFreeCodons, CODON_AA,
} from './rna.js';
import { makeLipid, makeMembrane, isInsideMembrane } from './lipid.js';

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
  if (formProb < 0.01) return 0;

  let formed = 0;

  // ─── Extension pass: every existing chain gets a chance every tick ───
  // (This is the fix — previously chains only extended if their cell happened
  // to be one of 200 randomly sampled. With 85 chains in 40k cells that's
  // ~0.4 attempts/tick total → max length never grew past trimers.)
  for (const chain of structures) {
    if (chain.type !== 'rna') continue;
    if (chain.sequence.length >= 64) continue;  // hard cap
    // H-bonded chains preferentially get extended via the templated path,
    // not the random one — skip them here.
    if (chain.hBondedTo != null) continue;
    const ext = formProb * (chain.surfaceBound ? 1.0 : 0.4);
    if (world._rand() >= ext) continue;
    const k = chain.position.y * W + chain.position.x;
    const totN =
      fields[A_IDX][k] + fields[U_IDX][k] +
      fields[G_IDX][k] + fields[C_IDX][k];
    if (totN < NUCLEOTIDE_THRESHOLD * 0.5) continue;
    const nuc = _pickNucleotide(world, k, fields);
    if (nuc) {
      chain.sequence.push(nuc);
      formed++;
    }
  }

  // ─── Spontaneous dimer formation: random cell sampling ───
  // Empty cells with high local nucleotide concentration occasionally
  // nucleate a fresh 2-mer. Cap at a fraction of grid per tick.
  const N_SAMPLES = 120;
  for (let s = 0; s < N_SAMPLES; s++) {
    const x = Math.floor(world._rand() * W);
    const y = Math.floor(world._rand() * H);
    const k = y * W + x;
    const totN =
      fields[A_IDX][k] + fields[U_IDX][k] +
      fields[G_IDX][k] + fields[C_IDX][k];
    if (totN < NUCLEOTIDE_THRESHOLD) continue;
    if (world._rand() >= formProb * 0.25) continue;
    const n1 = _pickNucleotide(world, k, fields);
    if (!n1) continue;
    const n2 = _pickNucleotide(world, k, fields);
    if (!n2) continue;
    structures.push(makeRNA([n1, n2], x, y, true));
    formed++;
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
  // Index chains by cell, then look at the cell + its 8 neighbors when
  // searching for potential partners (real chains in solution can collide
  // even when not in exactly the same cell).
  const H = world.height;
  const cellMap = new Map();
  for (const st of structures) {
    if (st.type !== 'rna') continue;
    if (st.hBondedTo != null) continue;
    const k = st.position.y * W + st.position.x;
    if (!cellMap.has(k)) cellMap.set(k, []);
    cellMap.get(k).push(st);
  }
  // For each unbonded chain, scan its cell + 8 neighbors for partners
  for (const a of structures) {
    if (a.type !== 'rna' || a.hBondedTo != null) continue;
    if (a.sequence.length < 2) continue;
    const ax = a.position.x;
    const ay = a.position.y;
    let bestB = null;
    let bestProb = 0;
    for (let dy = -1; dy <= 1; dy++) {
      const ny = (ay + dy + H) % H;
      for (let dx = -1; dx <= 1; dx++) {
        const nx = (ax + dx + W) % W;
        const list = cellMap.get(ny * W + nx);
        if (!list) continue;
        for (const b of list) {
          if (b === a || b.hBondedTo != null) continue;
          if (b.sequence.length < 2) continue;
          const { score } = _pairScore(a, b);
          if (score === 0) continue;
          const minL = Math.min(a.sequence.length, b.sequence.length);
          const matchFrac = score / minL;
          if (matchFrac < 0.35) continue;  // ≥35% complementary
          // Distance penalty for neighbor cells
          const distPenalty = (dx === 0 && dy === 0) ? 1.0 : 0.5;
          const prob = rates.hBondForm * matchFrac * distPenalty;
          if (prob > bestProb) {
            bestProb = prob;
            bestB = b;
          }
        }
      }
    }
    if (bestB && world._rand() < bestProb) {
      a.hBondedTo = bestB.id;
      bestB.hBondedTo = a.id;
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

/**
 * Phase 2 step 8 + 10: Templated polymerization with mutation (R4 + R5).
 *
 * For each H-bonded pair where the two strands have unequal lengths, the
 * shorter one is the "growing copy" of the longer "template". We try to
 * extend the copy by one base, using the complement of the next template
 * position. With probability mutationRate the wrong base is added instead.
 *
 * Replication emerges when:
 *   1. A template T templates a copy C (R4 templated)
 *   2. Hot phase separates them (R3 break)
 *   3. C templates a new C' (which is the complement of C, == T's sequence)
 * No code change needed for step 9 — separation is already part of R3.
 */
export function templatedExtension(world, rates) {
  if (rates.backboneFormTemplated < 0.01) return 0;
  // Phase 3 catalyst lookup: rna_replicase ribozyme grid
  const repCatGrid = _ribozymeGridByType(world, 'rna_replicase');
  let extended = 0;
  for (const a of world.structures) {
    if (a.type !== 'rna' || a.hBondedTo == null) continue;
    if (a.id > a.hBondedTo) continue;  // process each pair once

    const b = _findById(world.structures, a.hBondedTo);
    if (!b || b.type !== 'rna') continue;

    // Identify template (longer) and copy (shorter, may be equal)
    let template, copy;
    if (a.sequence.length > b.sequence.length) {
      template = a;
      copy = b;
    } else if (b.sequence.length > a.sequence.length) {
      template = b;
      copy = a;
    } else {
      continue;  // equal length → fully copied, skip
    }

    // Position on template that maps to the next slot on copy
    // copy[i] pairs antiparallel with template[T_len - 1 - i]
    // After extending copy, copy.length grows by 1, so we need template[T_len - 1 - copy.length]
    const k = template.sequence.length - 1 - copy.sequence.length;
    if (k < 0) continue;

    const correct = COMPLEMENT[template.sequence[k]];
    if (!correct) continue;  // safety

    // Mutation roll
    let nuc = correct;
    if (world._rand() < rates.mutationRate) {
      const others = NUCLEOTIDES.filter((n) => n !== correct);
      nuc = others[Math.floor(world._rand() * others.length)];
    }

    // Check the cell has the chosen nucleotide
    const cellIdx = copy.position.y * world.width + copy.position.x;
    const fIdx = SPECIES_INDEX[nuc];
    if (world.fields[fIdx][cellIdx] < 0.05) continue;

    // Roll polymerization probability (with replicase ribozyme boost)
    const repBoost = repCatGrid[cellIdx] || 0;
    const prob =
      rates.backboneFormTemplated *
      (1 + 5 * repBoost) *
      (copy.surfaceBound ? 1.0 : 0.5) *
      (nuc === correct ? 1.0 : 0.3);  // mutations are slower
    if (world._rand() >= prob) continue;

    // Consume the nucleotide and append
    world.fields[fIdx][cellIdx] = Math.max(0, world.fields[fIdx][cellIdx] - 0.05);
    copy.sequence.push(nuc);
    extended++;
  }
  return extended;
}

// ─── Phase 4: lipids + membranes (R10, R11) ─────────────────────────────

const FA_IDX = SPECIES_INDEX.FA;

/**
 * Phase 4 step 17: Lipid nucleation from FA concentration field.
 *
 * In wet cells where the local FA concentration exceeds a threshold, spawn
 * a Lipid particle and consume some FA. Caps the per-tick spawn count to
 * keep performance bounded.
 */
export function nucleateLipids(world, rates) {
  if (rates.lipidAssembly < 0.05) return 0;
  const W = world.width;
  const H = world.height;
  const fields = world.fields;
  const NUC_THRESHOLD = 0.4;
  const N_SAMPLES = 80;
  let spawned = 0;
  for (let s = 0; s < N_SAMPLES; s++) {
    const x = Math.floor(world._rand() * W);
    const y = Math.floor(world._rand() * H);
    const k = y * W + x;
    if (fields[FA_IDX][k] < NUC_THRESHOLD) continue;
    if (world._rand() >= rates.lipidAssembly * 0.5) continue;
    fields[FA_IDX][k] = Math.max(0, fields[FA_IDX][k] - 0.3);
    world.structures.push(makeLipid(x, y, world._rand() * Math.PI * 2));
    spawned++;
  }
  return spawned;
}

/**
 * Phase 4 step 18+19: Lipid self-assembly + vesicle formation.
 *
 * Each tick, scan unmembraned lipids. For each cluster of ≥8 lipids
 * within a small radius, form a Membrane object.
 *
 * IMPORTANT: when a cluster qualifies, the vesicle's center is biased
 * toward nearby ribozymes/chains. Real protocell formation isn't random
 * — lipid bilayers preferentially wrap around hydrophilic molecules.
 * This bias massively increases the chance that a vesicle encloses
 * existing ribozymes, which is the prerequisite for protocell emergence.
 */
export function assembleMembranes(world) {
  const lipids = world.structures.filter(
    (s) => s.type === 'lipid' && s.membraneId == null,
  );
  if (lipids.length < 8) return 0;
  const W = world.width;
  const H = world.height;

  // Pre-collect chain centroid weighted by length (longer = more anchor weight)
  // Ribozymes get a much heavier weight so vesicles snap onto ribozyme clusters.
  const anchors = [];
  for (const s of world.structures) {
    if (s.type !== 'rna') continue;
    anchors.push({
      x: s.position.x,
      y: s.position.y,
      weight: s.catalyticFunction ? 20 : Math.max(1, s.sequence.length / 4),
    });
  }

  const used = new Set();
  let formed = 0;
  for (const seed of lipids) {
    if (used.has(seed.id)) continue;
    const cluster = [seed];
    let cx = seed.position.x;
    let cy = seed.position.y;
    const radiusSq = 6 * 6;  // 6-cell gathering radius (was 4)
    for (const other of lipids) {
      if (other === seed || used.has(other.id)) continue;
      let dx = Math.abs(other.position.x - cx);
      let dy = Math.abs(other.position.y - cy);
      if (dx > W / 2) dx = W - dx;
      if (dy > H / 2) dy = H - dy;
      if (dx * dx + dy * dy <= radiusSq) {
        cluster.push(other);
        cx = (cx * (cluster.length - 1) + other.position.x) / cluster.length;
        cy = (cy * (cluster.length - 1) + other.position.y) / cluster.length;
      }
    }
    if (cluster.length >= 8) {
      // Bias center toward nearby chain anchors (within 18 cells now)
      let biasX = cx;
      let biasY = cy;
      let totalW = 1.0;
      for (const a of anchors) {
        let dx = Math.abs(a.x - cx);
        let dy = Math.abs(a.y - cy);
        if (dx > W / 2) dx = W - dx;
        if (dy > H / 2) dy = H - dy;
        const d2 = dx * dx + dy * dy;
        if (d2 > 324) continue;  // 18-cell radius (was 12)
        biasX += a.x * a.weight;
        biasY += a.y * a.weight;
        totalW += a.weight;
      }
      const finalX = biasX / totalW;
      const finalY = biasY / totalW;

      const m = makeMembrane(finalX, finalY, cluster.map((c) => c.id));
      world.structures.push(m);
      for (const c of cluster) {
        c.membraneId = m.id;
        used.add(c.id);
      }
      formed++;
    }
  }
  return formed;
}

/**
 * Phase 4 step 20+21: Update each membrane's enclosed list, integrity, and
 * apply tide-driven dissolution. Also dissolve membranes whose lipid count
 * drops below threshold.
 */
export function updateMembranes(world, waterLevel) {
  const membranes = world.structures.filter((s) => s.type === 'membrane');
  if (membranes.length === 0) return;
  const W = world.width;
  const H = world.height;

  // Build a quick lipid lookup
  const lipidById = new Map();
  for (const s of world.structures) {
    if (s.type === 'lipid') lipidById.set(s.id, s);
  }

  for (const m of membranes) {
    // Drop dead lipids
    m.lipids = m.lipids.filter((id) => lipidById.has(id));
    if (m.lipids.length < 5) {
      m.integrity = 0;
      continue;
    }
    // Recompute center
    let cx = 0, cy = 0;
    for (const id of m.lipids) {
      const lp = lipidById.get(id);
      cx += lp.position.x;
      cy += lp.position.y;
    }
    cx /= m.lipids.length;
    cy /= m.lipids.length;
    m.center.x = cx;
    m.center.y = cy;
    m.radius = Math.max(3, Math.sqrt(m.lipids.length / Math.PI) * 1.5);

    // Tide effect: dry phase weakens hydrophobic effect
    if (waterLevel < -0.2) {
      m.integrity -= 0.005 * Math.abs(waterLevel);
    } else {
      m.integrity = Math.min(1, m.integrity + 0.002);
    }

    // Compute enclosed structures (RNA, peptide, lipid)
    m.enclosed = [];
    for (const s of world.structures) {
      if (s.type === 'membrane') continue;
      if (s.type === 'lipid' && s.membraneId === m.id) continue;
      if (isInsideMembrane(m, s.position.x, s.position.y, W, H)) {
        m.enclosed.push(s.id);
      }
    }

    m.age++;
  }

  // Cull dead membranes (integrity ≤ 0)
  for (let i = world.structures.length - 1; i >= 0; i--) {
    const s = world.structures[i];
    if (s.type === 'membrane' && s.integrity <= 0) {
      // Free its lipids
      for (const id of s.lipids) {
        const lp = lipidById.get(id);
        if (lp) lp.membraneId = null;
      }
      world.structures.splice(i, 1);
    }
  }
}

/**
 * Phase 3 step 14: Ribozyme activation (R9).
 *
 * For each RNA chain, check if it qualifies as a ribozyme. Cache the result
 * on the chain. We re-check whenever the sequence changes — for simplicity
 * we re-check every chain every N ticks (cheap with O(L) motif scan).
 */
export function activateRibozymes(world) {
  for (const st of world.structures) {
    if (st.type !== 'rna') continue;
    if (st.sequence.length < 15) continue;  // shortest motif length
    st.catalyticFunction = detectRibozyme(st);
  }
}

/**
 * Phase 3 step 13: Codon → amino acid attachment (R7).
 *
 * Scan each RNA chain for free codons that map to a known amino acid.
 * If the matching amino acid is in the local cell, attach with rate
 * aminoAcidAttach. Catalyzed by aminoacyl_transferase ribozyme in same cell.
 */
export function attachAminoAcids(world, rates) {
  if (rates.aminoAcidAttach < 0.005) return;
  const W = world.width;
  const fields = world.fields;

  // Index aminoacyl catalysts by cell for boost lookup
  const catGrid = _ribozymeGridByType(world, 'aminoacyl_transferase');

  for (const st of world.structures) {
    if (st.type !== 'rna') continue;
    if (st.sequence.length < 3) continue;
    const codons = findFreeCodons(st);
    if (codons.length === 0) continue;

    const cellIdx = st.position.y * W + st.position.x;
    const catBoost = catGrid[cellIdx] || 0;
    const baseProb = rates.aminoAcidAttach * (1 + 5 * catBoost);

    for (const { index, aa } of codons) {
      const aaIdx = SPECIES_INDEX[aa];
      if (fields[aaIdx][cellIdx] < 0.05) continue;
      if (world._rand() >= baseProb) continue;
      // Attach
      fields[aaIdx][cellIdx] = Math.max(0, fields[aaIdx][cellIdx] - 0.04);
      st.attachedAminoAcids.push({ index, aminoAcid: aa });
    }
  }
}

/**
 * Phase 3 step 15: Peptide bond formation (R8).
 *
 * When two amino acids are attached at adjacent codon positions on the
 * same RNA chain, they form a peptide bond. The new peptide detaches
 * from the template. Strongly catalyzed by peptidyl_transferase ribozyme
 * in same cell.
 */
export function formPeptideBonds(world, rates) {
  if (rates.peptideBond < 0.001 && rates.peptideBondCatalyzed < 0.001) return 0;
  const W = world.width;
  const catGrid = _ribozymeGridByType(world, 'peptidyl_transferase');

  let formed = 0;
  for (const st of world.structures) {
    if (st.type !== 'rna') continue;
    const aas = st.attachedAminoAcids;
    if (aas.length < 2) continue;
    // Sort by index ascending
    aas.sort((a, b) => a.index - b.index);

    const cellIdx = st.position.y * W + st.position.x;
    const catBoost = catGrid[cellIdx] || 0;
    const prob = catBoost > 0
      ? rates.peptideBondCatalyzed * (1 + 5 * catBoost)
      : rates.peptideBond;

    // Walk adjacent codon positions (i, i+3)
    let i = 0;
    while (i < aas.length - 1) {
      const a1 = aas[i];
      const a2 = aas[i + 1];
      if (a2.index === a1.index + 3) {
        if (world._rand() < prob) {
          // Form a peptide chain from this pair (extend Phase: continue chain
          // if a peptide is already growing on the same RNA — for simplicity,
          // we always create a new dipeptide here).
          const peptide = {
            type: 'peptide',
            id: _nextPeptideId++,
            sequence: [a1.aminoAcid, a2.aminoAcid],
            position: { x: st.position.x, y: st.position.y },
            parentRNA: st.id,
            surfaceBound: st.surfaceBound,
            age: 0,
          };
          world.structures.push(peptide);
          // Remove these two AAs from the RNA (they've left as a peptide)
          aas.splice(i, 2);
          formed++;
          continue;  // don't increment i, since we removed elements
        }
      }
      i++;
    }
  }
  return formed;
}

let _nextPeptideId = 100000;  // simple separate counter for peptide ids

// Build an (H*W)-length array of catalytic strength of a given ribozyme type
// summed over all chains in each cell. Used for catalyst lookup.
function _ribozymeGridByType(world, type) {
  const W = world.width;
  const grid = new Float32Array(world.size);
  for (const st of world.structures) {
    if (st.type !== 'rna') continue;
    const cf = st.catalyticFunction;
    if (!cf || cf.type !== type) continue;
    const k = st.position.y * W + st.position.x;
    grid[k] += cf.strength;
  }
  return grid;
}

/**
 * Phase 2 step 11: Degradation + supply (R6).
 *
 * - Free monomers in concentration field decay at a constant slow rate.
 * - RNA chains hydrolyze. Short chains (< 3) easily, long chains rarely.
 *   Hydrolysis is reduced by self-complementarity (better-folded RNA
 *   resists more — see SPEC R6 update).
 * - Edge cells receive fresh monomers during high tide (ocean delivery).
 */
export function degradeAndSupply(world, rates) {
  const fields = world.fields;
  const N = world.size;

  // 1. Decay free monomers in concentration field
  const decay = rates.degradation;
  if (decay > 0) {
    for (let s = 0; s < N_SPECIES; s++) {
      const f = fields[s];
      const factor = 1 - decay;
      for (let k = 0; k < N; k++) f[k] *= factor;
    }
  }

  // 2. Hydrolyze RNA chains
  const baseHydrolysis = rates.hydrolysis;
  if (baseHydrolysis > 0) {
    const dead = [];
    for (let i = 0; i < world.structures.length; i++) {
      const st = world.structures[i];
      if (st.type !== 'rna') continue;
      const L = st.sequence.length;
      // Length factor: smoother curve so length distribution is broader.
      // Sweet spot 10-30; very long chains are mechanically unstable
      // (entropy + UV damage cumulatively raises break risk).
      let lenFactor;
      if (L < 3) lenFactor = 5.0;
      else if (L < 6) lenFactor = 2.0;
      else if (L < 10) lenFactor = 0.8;
      else if (L < 20) lenFactor = 0.3;
      else if (L < 35) lenFactor = 0.15;
      else if (L < 50) lenFactor = 0.4;   // long chains a bit more fragile
      else lenFactor = 0.8;               // very long chains break under entropy
      // Self-complementarity protection
      const selfComp = selfComplementarity(st);
      const foldFactor = 1 - 0.5 * selfComp;
      // H-bonded chains are doubly protected
      const pairProtection = st.hBondedTo != null ? 0.3 : 1.0;

      const p = baseHydrolysis * lenFactor * foldFactor * pairProtection;
      if (world._rand() < p) {
        // Return the monomers to the local pool
        const cellIdx = st.position.y * world.width + st.position.x;
        for (const b of st.sequence) {
          const fIdx = SPECIES_INDEX[b];
          fields[fIdx][cellIdx] += 0.04;
        }
        dead.push(i);
      }
    }
    // Remove dead chains in reverse order
    for (let j = dead.length - 1; j >= 0; j--) {
      const dead_i = dead[j];
      const st = world.structures[dead_i];
      // Break any H-bond on the partner
      if (st.hBondedTo != null) {
        const partner = _findById(world.structures, st.hBondedTo);
        if (partner) partner.hBondedTo = null;
      }
      world.structures.splice(dead_i, 1);
    }
  }

  // 3. Edge supply during high tide
  if (rates.monomerSupply > 0) {
    const W = world.width;
    const H = world.height;
    const supply = rates.monomerSupply * 0.05;
    // Supply to top + bottom rows
    for (let x = 0; x < W; x++) {
      for (const s of [SPECIES_INDEX.A, SPECIES_INDEX.U, SPECIES_INDEX.G, SPECIES_INDEX.C, SPECIES_INDEX.FA, SPECIES_INDEX.Gly, SPECIES_INDEX.Ala]) {
        fields[s][x] += supply;
        fields[s][(H - 1) * W + x] += supply;
      }
    }
    // Supply to left + right columns
    for (let y = 0; y < H; y++) {
      for (const s of [SPECIES_INDEX.A, SPECIES_INDEX.U, SPECIES_INDEX.G, SPECIES_INDEX.C, SPECIES_INDEX.FA, SPECIES_INDEX.Gly, SPECIES_INDEX.Ala]) {
        fields[s][y * W] += supply;
        fields[s][y * W + (W - 1)] += supply;
      }
    }
  }
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
