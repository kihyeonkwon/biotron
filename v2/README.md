# BIOTRON v2 — Abiogenesis Simulator

A simulation of life's origin from first principles. AUGC nucleotides, amino acids, and lipids interact via physicochemical rules on a prehistoric Earth surface, subject to day/night and tidal cycles. **No biology is hardcoded.** Replication, catalysis, translation, and compartmentalization must emerge from simple chemical rules.

> **Philosophical premise**: Life is a recursive algorithm that emerged from simple chemical rules at planetary scale. The same way token prediction generates intelligence from a simple principle, complementary base pairing and thermodynamic cycling generated life from a simple principle.

## Quick start

```bash
cd v2
npm install
npm run dev
```

Open http://localhost:5173 in your browser. Click **▶** to start. Drag the speed slider for time compression.

## What you should see

The 200×200 grid is a tidal pool. Background tint shifts with temperature; a translucent blue overlay rises and falls with the tide.

| Tick range | Expected emergence |
|---|---|
| 0–100 | Diffusion: nucleotide blobs at corners spread out into smooth gradients |
| 100–500 | Warm + dry phases trigger spontaneous RNA polymerization (small colored squares appear) |
| 500–2000 | H-bonds form between complementary chains during cold (dotted cyan lines) |
| 2000–10000 | Templated extension on bonded pairs → replication |
| 10000+ | Long chains reach ribozyme length → golden glow appears |
| 10000+ | Lipids nucleate from FA blobs → aggregate → form vesicles (orange circles) |
| Eventually | A vesicle encloses replicating + translating ribozymes → **PROTOCELL** |

The sidebar tracks all 9 emergence milestones with the tick they were first reached.

## Architecture

```
v2/
├── SPEC.md                      # full design spec
├── package.json
├── vite.config.js
├── index.html
├── src/
│   ├── App.jsx                  # React UI: canvas + sidebar + controls
│   ├── main.jsx
│   ├── styles.css
│   ├── simulation/
│   │   ├── world.js             # World class + tick loop + Float32Array fields
│   │   ├── environment.js       # Temperature + waterLevel cycles + reaction rates
│   │   ├── constants.js         # 12 species (AUGC + 5 amino acids + P/Mg/FA)
│   │   ├── rna.js               # RNA chain factory + helpers (ribozyme detection,
│   │   │                        # codon scanning, self-complementarity score)
│   │   ├── lipid.js             # Lipid + Membrane factories
│   │   ├── reactions.js         # All R1-R11 chemical rules
│   │   └── metrics.js           # Milestone tracker
│   └── renderer/
│       └── canvas.js            # ImageData-based renderer (concentration field
│                                # heatmap + chain/peptide/lipid/membrane overlays)
```

## The 5 build phases

1. **Foundation** (steps 1-7): grid, diffusion, environment cycles, surface adsorption, RNA polymerization, H-bonding, basic stats
2. **Replication** (steps 8-12): templated polymerization, mutation, degradation + supply, milestone detection 1-4
3. **Translation** (steps 13-16): codon → amino acid attachment, ribozyme activation, peptide bond formation, milestones 5-6
4. **Compartmentalization** (steps 17-22): lipid nucleation, self-assembly, vesicles, semi-permeable membranes, milestones 7-9 (protocell goal)
5. **Polish** (steps 23-27): performance, presets, README

Each step is committed separately so you can `git checkout` any earlier stage and see what it looked like.

## The encoded chemistry rules (the only "magic")

Two rules represent physical-chemistry laws and are explicitly encoded. Everything else emerges.

1. **Codon → amino acid affinity** (R7) — based on Yarus stereochemical hypothesis (2009):
   ```
   GCC → Ala     GGC → Gly     GUC → Val
   GAC → Asp     GAG → Glu
   ```
   These specific RNA triplets have measurable physical-chemistry affinity for these amino acids, independent of any biological machinery.

2. **Ribozyme motifs** (R9) — sequence determines structure determines function:
   ```
   length ≥ 15 + contains "GGCGCC"      → peptidyl_transferase
   length ≥ 15 + contains "CCCUUU"      → rna_replicase
   length ≥ 20 + contains "GGGAAACCC"   → aminoacyl_transferase
   ```
   Ribozyme strength is `0.5 × GC_ratio + 0.5 × self-complementarity` — a chain is a better catalyst when it has more stable bonds AND more potential for hairpin folding.

Everything else — base pairing, polymerization, replication, peptide formation, vesicle self-assembly, compartmentalized selection — has to *emerge* from these plus the basic R1-R6/R8/R10/R11 dynamics.

## Controls

- **▶ / ⏸** — play / pause
- **step** — advance 1 tick (when paused)
- **↻** — reset to fresh world
- **speed slider** — 1× to 32× ticks per animation frame (real-time → time-lapse)

## Scientific references

- Miller-Urey (1953) — Prebiotic amino acid synthesis
- Szostak lab — RNA replication and protocell formation
- Deamer — Lipid vesicles and tidal cycling in abiogenesis
- Yarus (2009) — Stereochemical basis for codon-amino acid affinity
- Lincoln & Joyce (2009) — Self-replicating RNA systems
- Ferris (2006) — Montmorillonite clay catalysis of RNA polymerization
