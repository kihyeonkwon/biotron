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

// Codon → amino acid (Yarus stereochemical hypothesis, 5 simplest amino acids)
// These triplets have measurable physical-chemistry affinity for these AAs
// independent of any biological machinery.
export const CODON_AA = {
  GCC: 'Ala',
  GGC: 'Gly',
  GUC: 'Val',
  GAC: 'Asp',
  GAG: 'Glu',
};

// Ribozyme motifs (R9 — encoded chemistry rule)
// All motifs are 6 nucleotides (~1/4096 per window) and use the same
// length threshold so all 3 types appear at comparable rates.
// Length threshold 12 instead of 15-20 to broaden the ribozyme population
// (was 2.7% of chains, target ~10%).
export const RIBOZYME_MOTIFS = {
  peptidyl_transferase:  { motif: 'GGCGCC', minLength: 12 },
  rna_replicase:         { motif: 'CCCUUU', minLength: 12 },
  aminoacyl_transferase: { motif: 'GGGAAA', minLength: 12 },
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

// Check whether an RNA chain qualifies as a ribozyme. Returns
// {type, strength} or null. The combined strength formula is
//   strength = 0.5 * GC_ratio + 0.5 * selfComplementarity
// per the SPEC R9 update.
export function detectRibozyme(chain) {
  const seqStr = chain.sequence.join('');
  const L = chain.sequence.length;
  for (const [type, def] of Object.entries(RIBOZYME_MOTIFS)) {
    if (L < def.minLength) continue;
    if (seqStr.indexOf(def.motif) === -1) continue;
    const strength = 0.5 * gcRatio(chain) + 0.5 * selfComplementarity(chain);
    return { type, strength };
  }
  return null;
}

// Find all codon positions in a chain that are NOT yet attached to an AA.
// Returns array of {index, aa} for codons that map to a known amino acid.
export function findFreeCodons(chain) {
  const seq = chain.sequence;
  const L = seq.length;
  const taken = new Set(chain.attachedAminoAcids.map((x) => x.index));
  const out = [];
  for (let i = 0; i + 3 <= L; i += 3) {
    if (taken.has(i)) continue;
    const codon = seq[i] + seq[i + 1] + seq[i + 2];
    const aa = CODON_AA[codon];
    if (aa) out.push({ index: i, aa });
  }
  return out;
}
