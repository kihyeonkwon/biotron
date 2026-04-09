// Headless BIOTRON v2 runner.
//
// Imports the same simulation code the React app uses, runs N ticks,
// and prints periodic stats + a final summary. Lets Claude (or anyone)
// iterate on parameter tuning without firing up a browser.
//
// Usage:
//   node headless/run.js                   # default 50000 ticks
//   node headless/run.js --ticks 20000     # custom tick count
//   node headless/run.js --ticks 30000 --seed 7
//   node headless/run.js --quiet           # only final summary
//   node headless/run.js --json            # machine-readable JSON output
//   node headless/run.js --interval 500    # sample every 500 ticks (default 1000)

import { World } from '../src/simulation/world.js';

// ─── Argument parsing ───
const args = process.argv.slice(2);
function getArg(name, def) {
  const i = args.indexOf(`--${name}`);
  if (i === -1) return def;
  return args[i + 1];
}
function hasFlag(name) {
  return args.includes(`--${name}`);
}

const TICKS = parseInt(getArg('ticks', '50000'));
const SEED = parseInt(getArg('seed', '1'));
const INTERVAL = parseInt(getArg('interval', '1000'));
const QUIET = hasFlag('quiet');
const JSON_OUT = hasFlag('json');
const W = parseInt(getArg('width', '200'));
const H = parseInt(getArg('height', '200'));

// ─── Run ───
const t0 = Date.now();
const world = new World({ width: W, height: H, seed: SEED });

const samples = [];
const HEADER = [
  'tick', 'chains', 'maxL', 'hbond',
  'ribz', 'rep', 'aa', 'pT',
  'peps', 'lipid', 'vesi', 'gen', 'div',
  'temp', 'water',
];

function rowFromStats(stats) {
  return {
    tick: stats.tick,
    chains: stats.rnaChains ?? 0,
    maxL: stats.maxRnaLen ?? 0,
    hbond: stats.hBondedChains ?? 0,
    ribz: stats.ribozymes ?? 0,
    rep: stats.ribByType?.rna_replicase ?? 0,
    aa: stats.ribByType?.aminoacyl_transferase ?? 0,
    pT: stats.ribByType?.peptidyl_transferase ?? 0,
    peps: stats.peptides ?? 0,
    lipid: stats.lipids ?? 0,
    vesi: stats.membranes ?? 0,
    gen: stats.maxVesicleGen ?? 0,
    div: stats.vesicleDivisions ?? 0,
    temp: stats.temperature?.toFixed(2) ?? '0',
    water: stats.waterLevel?.toFixed(2) ?? '0',
  };
}

function printHeader() {
  if (JSON_OUT) return;
  const widths = { tick: 6, chains: 6, maxL: 5, hbond: 6,
                   ribz: 5, rep: 4, aa: 4, pT: 4,
                   peps: 5, lipid: 6, vesi: 5, gen: 4, div: 4,
                   temp: 6, water: 6 };
  const cols = HEADER.map((k) => k.padStart(widths[k]));
  console.log(cols.join(' '));
  console.log(HEADER.map((k) => '─'.repeat(widths[k])).join(' '));
}

function printRow(row) {
  if (JSON_OUT) return;
  const widths = { tick: 6, chains: 6, maxL: 5, hbond: 6,
                   ribz: 5, rep: 4, aa: 4, pT: 4,
                   peps: 5, lipid: 6, vesi: 5, gen: 4, div: 4,
                   temp: 6, water: 6 };
  const cells = HEADER.map((k) => String(row[k]).padStart(widths[k]));
  console.log(cells.join(' '));
}

if (!QUIET) printHeader();

for (let t = 0; t < TICKS; t++) {
  world.tick();
  if ((t + 1) % INTERVAL === 0 || t === TICKS - 1) {
    const stats = world.getStats();
    const row = rowFromStats(stats);
    samples.push(row);
    if (!QUIET) printRow(row);
  }
}

const dt = (Date.now() - t0) / 1000;

// ─── Final summary ───
const finalStats = world.getStats();
const finalMilestones = finalStats.milestoneStatus;

if (JSON_OUT) {
  console.log(JSON.stringify({
    config: { ticks: TICKS, seed: SEED, width: W, height: H, interval: INTERVAL },
    runtimeSeconds: dt,
    samples,
    final: rowFromStats(finalStats),
    milestones: finalMilestones,
  }, null, 2));
} else {
  console.log();
  console.log('─'.repeat(60));
  console.log(`Done in ${dt.toFixed(1)}s  (${(TICKS / dt).toFixed(0)} ticks/s)`);
  console.log();
  console.log('MILESTONES:');
  for (const m of finalMilestones) {
    const mark = m.tickReached != null ? '✓' : '·';
    const t = m.tickReached != null ? `t${m.tickReached.toLocaleString()}` : '—';
    console.log(`  ${mark} M${m.id} ${m.name.padEnd(24)} ${t}`);
  }
  console.log();
  const passed = finalMilestones.filter((m) => m.tickReached != null).length;
  console.log(`SCORE: ${passed} / 9 milestones`);
}
