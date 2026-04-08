// RNA chain — sequence of A/U/G/C bases.
//
// Chains are tracked individually (in world.structures[]). They have a
// position on the grid, may be surface-bound, may be H-bonded to a
// complementary partner, and may carry a catalytic function (Phase 3).

let _nextId = 1;

export const NUCLEOTIDES = ['A', 'U', 'G', 'C'];

export const COMPLEMENT = {
  A: 'U',
  U: 'A',
  G: 'C',
  C: 'G',
};

export function isComplementary(a, b) {
  return COMPLEMENT[a] === b;
}

export function makeRNA(sequence, x, y, surfaceBound = true) {
  return {
    type: 'rna',
    id: _nextId++,
    sequence: sequence.slice(),
    position: { x, y },
    surfaceBound,
    hBondedTo: null,           // chain id of paired strand, if any
    catalyticFunction: null,   // {type, strength} | null  (Phase 3)
    attachedAminoAcids: [],    // [{index, aminoAcid}]      (Phase 3)
    age: 0,
  };
}

export function rnaLength(chain) {
  return chain.sequence.length;
}

export function gcRatio(chain) {
  let n = 0;
  for (const b of chain.sequence) if (b === 'G' || b === 'C') n++;
  return n / chain.sequence.length;
}

// Self-complementarity score: count of intramolecular complementary pairs
// at distance >= 3, normalized by length. Used for ribozyme strength
// (Phase 3) and hydrolysis resistance (Phase 2). O(L²) per call.
export function selfComplementarity(chain) {
  const seq = chain.sequence;
  const L = seq.length;
  if (L < 6) return 0;
  let score = 0;
  for (let i = 0; i < L; i++) {
    for (let j = i + 3; j < L; j++) {
      if (isComplementary(seq[i], seq[j])) score++;
    }
  }
  return score / L;
}

export function resetRNAIds() {
  _nextId = 1;
}
