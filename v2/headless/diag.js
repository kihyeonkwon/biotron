// Diagnostic: run + dump per-vesicle content breakdown.
//
// Tells us EXACTLY why M9 isn't tripping by showing how ribozymes are
// distributed across vesicles.

import { World } from '../src/simulation/world.js';

const TICKS = 30000;
const SEED = 1;

const w = new World({ width: 200, height: 200, seed: SEED });
for (let t = 0; t < TICKS; t++) w.tick();

const stats = w.getStats();
const membranes = w.structures.filter((s) => s.type === 'membrane');
const idMap = new Map();
for (const s of w.structures) idMap.set(s.id, s);

console.log(`After ${TICKS} ticks (seed ${SEED}):`);
console.log(`  total chains: ${stats.rnaChains}`);
console.log(`  total ribozymes: ${stats.ribozymes}`);
console.log(`  by type: rep=${stats.ribByType.rna_replicase} aa=${stats.ribByType.aminoacyl_transferase} pT=${stats.ribByType.peptidyl_transferase}`);
console.log(`  total vesicles: ${stats.membranes}`);
console.log();

// Per-vesicle breakdown
let vesiclesWithAnyRibozyme = 0;
let vesiclesWith2Types = 0;
let vesiclesWith3Types = 0;
const histogram = { 0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, '6+': 0 };
const typeMixHistogram = { 0: 0, 1: 0, 2: 0, 3: 0 };

for (const m of membranes) {
  let nRibz = 0;
  const types = new Set();
  for (const eid of m.enclosed) {
    const en = idMap.get(eid);
    if (!en || en.type !== 'rna' || !en.catalyticFunction) continue;
    nRibz++;
    for (const f of en.catalyticFunction) types.add(f.type);
  }
  if (nRibz > 0) vesiclesWithAnyRibozyme++;
  if (types.size >= 2) vesiclesWith2Types++;
  if (types.size >= 3) vesiclesWith3Types++;
  const bin = nRibz >= 6 ? '6+' : nRibz;
  histogram[bin]++;
  typeMixHistogram[types.size]++;
}

console.log('VESICLE → ENCLOSED RIBOZYME COUNT:');
for (const k of [0, 1, 2, 3, 4, 5, '6+']) {
  console.log(`  ${k} ribozymes: ${histogram[k]} vesicles`);
}
console.log();
console.log('VESICLE → ENCLOSED RIBOZYME TYPE COUNT:');
for (const k of [0, 1, 2, 3]) {
  console.log(`  ${k} distinct types: ${typeMixHistogram[k]} vesicles`);
}
console.log();
console.log(`Vesicles with ANY ribozyme:    ${vesiclesWithAnyRibozyme} / ${membranes.length}`);
console.log(`Vesicles with 2 distinct types: ${vesiclesWith2Types}`);
console.log(`Vesicles with 3 distinct types: ${vesiclesWith3Types}  ← M9 needs ≥1`);

// Sample: dump a few vesicles that have the most ribozymes
const sortedVesicles = membranes
  .map((m) => {
    let n = 0;
    const types = new Set();
    const seqs = [];
    for (const eid of m.enclosed) {
      const en = idMap.get(eid);
      if (!en || en.type !== 'rna' || !en.catalyticFunction) continue;
      n++;
      for (const f of en.catalyticFunction) types.add(f.type);
      seqs.push({ len: en.sequence.length, types: en.catalyticFunction.map((f) => f.type) });
    }
    return { id: m.id, n, types, seqs, total: m.enclosed.length, radius: m.radius };
  })
  .filter((x) => x.n > 0)
  .sort((a, b) => b.n - a.n)
  .slice(0, 5);

console.log();
console.log('TOP 5 VESICLES BY RIBOZYME COUNT:');
for (const v of sortedVesicles) {
  console.log(`  vesicle #${v.id}  r=${v.radius.toFixed(1)}  enclosed=${v.total}  ribozymes=${v.n}  types={${[...v.types].join(',')}}`);
  for (const s of v.seqs.slice(0, 8)) {
    console.log(`    chain L=${s.len} types=[${s.types.join(',')}]`);
  }
}
