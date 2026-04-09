// Run several seeds in sequence and report milestone outcomes.
// Used for sanity checking parameter changes across RNG variation.

import { World } from '../src/simulation/world.js';

const SEEDS = [1, 2, 3, 7, 42];
const TICKS = 30000;

console.log(`Sweep: ${SEEDS.length} seeds × ${TICKS} ticks each`);
console.log('─'.repeat(80));
console.log(
  ['seed', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9',
   'chains', 'maxL', 'ribz', 'rep', 'aa', 'pT', 'vesi', 'gen', 'div'].join('\t'),
);

const aggregated = { passed: [0, 0, 0, 0, 0, 0, 0, 0, 0] };

for (const seed of SEEDS) {
  const w = new World({ width: 200, height: 200, seed });
  for (let t = 0; t < TICKS; t++) w.tick();
  const stats = w.getStats();
  const m = stats.milestoneStatus;
  const reach = (id) => {
    const x = m.find((mm) => mm.id === id);
    return x?.tickReached != null ? `t${x.tickReached}` : '—';
  };
  for (let i = 0; i < 9; i++) {
    if (m[i].tickReached != null) aggregated.passed[i]++;
  }
  console.log(
    [
      `s${seed}`,
      reach(1), reach(2), reach(3), reach(4), reach(5),
      reach(6), reach(7), reach(8), reach(9),
      stats.rnaChains, stats.maxRnaLen, stats.ribozymes,
      stats.ribByType?.rna_replicase ?? 0,
      stats.ribByType?.aminoacyl_transferase ?? 0,
      stats.ribByType?.peptidyl_transferase ?? 0,
      stats.membranes, stats.maxVesicleGen, stats.vesicleDivisions,
    ].join('\t'),
  );
}

console.log('─'.repeat(80));
console.log('Pass rates across seeds:');
for (let i = 0; i < 9; i++) {
  const rate = aggregated.passed[i] / SEEDS.length;
  const bar = '█'.repeat(Math.round(rate * 10));
  console.log(`  M${i + 1}  ${aggregated.passed[i]}/${SEEDS.length}  ${bar}`);
}
