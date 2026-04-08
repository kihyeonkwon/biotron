// Milestone detection — observes the world state and records when each
// emergence milestone is first reached. Does NOT modify the world.

import { selfComplementarity } from './rna.js';

export const MILESTONES = [
  { id: 1, name: 'Base pairing',           description: 'A-U / G-C pairs form during cold phases' },
  { id: 2, name: 'Oligomers',              description: 'RNA chains length 3-10 exist' },
  { id: 3, name: 'Template copying',       description: 'A chain has >80% complementarity to another' },
  { id: 4, name: 'Self-replication',       description: 'A copied strand reaches its template length' },
  { id: 5, name: 'Ribozyme emergence',     description: 'An RNA chain gains catalytic function (Phase 3)' },
  { id: 6, name: 'Primitive translation',  description: 'Peptide chains form via codon → amino acid (Phase 3)' },
  { id: 7, name: 'Vesicle formation',      description: 'Lipids self-assemble into closed membranes (Phase 4)' },
  { id: 8, name: 'Compartmentalized RNA',  description: 'A vesicle contains a replicating RNA chain' },
  { id: 9, name: 'Protocell',              description: 'Vesicle with replication + translation + catalysis ★' },
];

export function createMilestoneTracker() {
  const reached = {}; // id -> tick when first reached
  let _completedCopies = 0;

  function update(world) {
    const tick = world.tickCount;
    const structures = world.structures;

    // M1 — base pairing (any H-bonded chain exists)
    if (!reached[1]) {
      for (const s of structures) {
        if (s.type === 'rna' && s.hBondedTo != null) {
          reached[1] = tick;
          break;
        }
      }
    }

    // M2 — oligomers (any chain length 3-10)
    if (!reached[2]) {
      for (const s of structures) {
        if (s.type === 'rna' && s.sequence.length >= 3 && s.sequence.length <= 10) {
          reached[2] = tick;
          break;
        }
      }
    }

    // M3 — template copying (an H-bonded pair where score / minL >= 0.8)
    if (!reached[3]) {
      for (const a of structures) {
        if (a.type !== 'rna' || a.hBondedTo == null) continue;
        if (a.id > a.hBondedTo) continue;
        let b = null;
        for (const s of structures) if (s.id === a.hBondedTo) { b = s; break; }
        if (!b) continue;
        const matchFrac = _antiparallelMatchFraction(a.sequence, b.sequence);
        if (matchFrac >= 0.8) {
          reached[3] = tick;
          break;
        }
      }
    }

    // M4 — self-replication (a copy reached its template length AND is no longer in progress)
    // Simple proxy: count the number of "completed" pairs where both chains have equal length
    // and score / L >= 0.7, AND have been bonded for at least one tick.
    // We track _completedCopies cumulatively; M4 trips when it crosses 5.
    if (!reached[4]) {
      let cur = 0;
      for (const a of structures) {
        if (a.type !== 'rna' || a.hBondedTo == null) continue;
        if (a.id > a.hBondedTo) continue;
        let b = null;
        for (const s of structures) if (s.id === a.hBondedTo) { b = s; break; }
        if (!b) continue;
        if (a.sequence.length !== b.sequence.length) continue;
        if (a.sequence.length < 4) continue;
        const matchFrac = _antiparallelMatchFraction(a.sequence, b.sequence);
        if (matchFrac >= 0.7) cur++;
      }
      _completedCopies = Math.max(_completedCopies, cur);
      if (_completedCopies >= 5) reached[4] = tick;
    }

    // M5 — ribozyme emergence (any RNA chain has catalyticFunction)
    if (!reached[5]) {
      for (const s of structures) {
        if (s.type === 'rna' && s.catalyticFunction != null) {
          reached[5] = tick;
          break;
        }
      }
    }

    // M6 — primitive translation (any peptide chain exists)
    if (!reached[6]) {
      for (const s of structures) {
        if (s.type === 'peptide') {
          reached[6] = tick;
          break;
        }
      }
    }

    // M7 — vesicle formation (any membrane exists)
    if (!reached[7]) {
      for (const s of structures) {
        if (s.type === 'membrane') {
          reached[7] = tick;
          break;
        }
      }
    }

    // M8 — compartmentalized RNA (a vesicle encloses an RNA chain)
    if (!reached[8]) {
      // Need a quick lookup of structure by id
      const idMap = new Map();
      for (const s of structures) idMap.set(s.id, s);
      for (const m of structures) {
        if (m.type !== 'membrane') continue;
        for (const eid of m.enclosed) {
          const en = idMap.get(eid);
          if (en && en.type === 'rna') {
            reached[8] = tick;
            break;
          }
        }
        if (reached[8]) break;
      }
    }

    // M9 — PROTOCELL: vesicle containing all 3 ribozyme types
    if (!reached[9]) {
      const idMap = new Map();
      for (const s of structures) idMap.set(s.id, s);
      for (const m of structures) {
        if (m.type !== 'membrane') continue;
        const types = new Set();
        for (const eid of m.enclosed) {
          const en = idMap.get(eid);
          if (en && en.type === 'rna' && en.catalyticFunction) {
            types.add(en.catalyticFunction.type);
          }
        }
        if (
          types.has('rna_replicase') &&
          types.has('peptidyl_transferase') &&
          types.has('aminoacyl_transferase')
        ) {
          reached[9] = tick;
          break;
        }
      }
    }
  }

  function getReached() {
    return { ...reached };
  }

  function getStatus() {
    return MILESTONES.map((m) => ({
      ...m,
      tickReached: reached[m.id] ?? null,
    }));
  }

  return { update, getReached, getStatus };
}

function _antiparallelMatchFraction(seqA, seqB) {
  const minL = Math.min(seqA.length, seqB.length);
  if (minL === 0) return 0;
  let score = 0;
  for (let i = 0; i < minL; i++) {
    const x = seqA[i];
    const y = seqB[seqB.length - 1 - i];
    if (
      (x === 'A' && y === 'U') || (x === 'U' && y === 'A') ||
      (x === 'G' && y === 'C') || (x === 'C' && y === 'G')
    ) score++;
  }
  return score / minL;
}
