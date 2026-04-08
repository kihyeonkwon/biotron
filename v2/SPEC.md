# Biotron v2: Abiogenesis Simulator
## From Chemistry to the First Cell

### IMPORTANT: Fresh start. All new code goes in `v2/` directory. Do not modify existing files.

---

## Vision

A simulation of life's origin from first principles. AUGC nucleotides, amino acids, and lipids interact via physicochemical rules on a prehistoric Earth surface, subject to day/night temperature cycles and tidal wet/dry cycles. No biology is hardcoded. Complex behaviors — replication, catalysis, translation, membrane formation, and compartmentalization — must **emerge** from simple chemical rules.

Two rules are explicitly encoded because they represent physical chemistry laws (like gravity):
1. **Codon-amino acid affinity** — specific RNA triplets have chemical affinity for specific amino acids
2. **Ribozyme motifs** — RNA sequences above a length threshold with specific motifs gain catalytic function

Everything else must emerge.

### Philosophical Premise
"Life is a recursive algorithm that emerged from simple chemical rules at planetary scale. The same way token prediction generates intelligence from a simple principle, complementary base pairing and thermodynamic cycling generated life from a simple principle."

---

## World Model

### Space: Surface + Solution

The world represents a **tidal pool on prehistoric Earth** — a mineral surface partially covered by water that rises and falls.

```
Cross-section view:

   Water surface (rises/falls with tide)
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  ← waterLevel
   |  SOLUTION ZONE                  |  ← 3D free diffusion
   |  (dissolved monomers float)     |
   |                                 |
   ==================================  ← SURFACE (mineral/clay)
      adsorbed molecules here            ← 2D constrained movement
      chains grow here
      catalysis happens here
```

Implementation:
- **Surface**: 2D grid (200x200 default, toroidal)
- **Solution**: Represented as a concentration field above the surface
- Particles can **adsorb** (solution → surface) and **desorb** (surface → solution)
- Surface = where polymerization and catalysis primarily happen
- Solution = where mixing, meeting, and transport happen

### Rendering
- Top-down view of the surface (primary view)
- Water level shown as translucent blue overlay that rises/falls
- When water covers a cell: slightly blue-tinted, high mobility
- When water recedes: exposed surface, low mobility, concentration increases

---

## Particle System: Hybrid Model

### Background: Concentration Fields (not individually tracked)

Each grid cell holds continuous concentration values for dissolved monomers:

```javascript
cell.concentrations = {
  // Nucleotides (RNA world — U not T)
  A: 0.0,   // Adenine
  U: 0.0,   // Uracil (NOT Thymine — RNA world predates DNA)
  G: 0.0,   // Guanine
  C: 0.0,   // Cytosine

  // Amino acids (5 simplest, found in meteorites)
  Gly: 0.0, // Glycine
  Ala: 0.0, // Alanine
  Val: 0.0, // Valine
  Asp: 0.0, // Aspartate
  Glu: 0.0, // Glutamate

  // Other essential molecules
  P: 0.0,   // Phosphate (energy + backbone material)
  Mg: 0.0,  // Magnesium ions (catalyst)
  FA: 0.0,  // Fatty acids (membrane material)
}
```

Concentration dynamics:
- Diffuse according to Fick's law (gradient-driven)
- Diffusion rate depends on water level (no water = no diffusion)
- Consumed when incorporated into structures
- Returned when structures degrade
- Supplied at edges (simulating ocean input during high tide)

### Structures: Individually Tracked Entities

```javascript
// RNA Chain
{
  type: "rna",
  id: unique,
  sequence: ["A", "U", "G", "C", "G", "C", "A", ...],
  position: {x, y},
  surfaceBound: true/false,
  hBondedTo: chainId | null,     // complementary strand
  catalyticFunction: null | {type, strength},
  attachedAminoAcids: [{index: 3, aminoAcid: "Ala"}, ...],
  age: 0,
}

// Peptide Chain
{
  type: "peptide",
  id: unique,
  sequence: ["Gly", "Ala", "Val", ...],
  position: {x, y},
  parentRNA: chainId,           // which RNA templated this
  surfaceBound: true/false,
  age: 0,
}

// Lipid
{
  type: "lipid",
  id: unique,
  position: {x, y},
  headAngle: 0-2π,              // direction head faces
  membraneId: null | membraneId, // part of which membrane
  age: 0,
}

// Membrane (vesicle)
{
  type: "membrane",
  id: unique,
  lipids: [lipidId, ...],       // constituent lipids
  center: {x, y},
  radius: float,
  enclosed: [chainId, ...],     // what's inside
  integrity: 0-1,               // structural health
}
```

Multiple structures can occupy the same grid cell. The grid cell is a spatial index, not a physical container.

---

## Environment Cycles

Temperature and water level are **independent continuous variables**, each driven by multiple overlapping sinusoidal cycles.

### Prehistoric Earth Parameters
- Day length: ~8 hours (Earth rotated faster 4 billion years ago)
- Tidal period: ~4 hours (Moon was closer, stronger tides)
- Year: ~365 modern days, but ~1095 "prehistoric days" at 8hr/day
- Tidal range: much larger than today (Moon 2-3x closer)

### Time Scale
```
1 tick = 1 prehistoric hour
1 day = 8 ticks
1 tidal cycle = 4 ticks (2 tides per day)
1 year = 8760 ticks (1095 days × 8 hours)
1 glacial cycle = 876000 ticks (100 years, compressed)
```

### Temperature Function
```javascript
function getTemperature(tick) {
  const dayLength = 8;
  const yearLength = 8760;
  const glacialLength = 876000;

  // Diurnal: day/night cycle
  const diurnal = Math.sin(2 * Math.PI * tick / dayLength);

  // Seasonal: summer/winter
  const seasonal = Math.sin(2 * Math.PI * tick / yearLength);

  // Glacial: ice age / interglacial (compressed)
  const glacial = Math.sin(2 * Math.PI * tick / glacialLength);

  // Noise: volcanic activity, random weather
  const noise = 0.05 * (Math.random() - 0.5);

  // Weighted sum → range [-1, +1]
  // -1 = extreme cold (glacial winter night)
  // +1 = extreme heat (interglacial summer noon)
  return 0.30 * diurnal
       + 0.35 * seasonal
       + 0.25 * glacial
       + noise;
}
```

### Water Level Function
```javascript
function getWaterLevel(tick) {
  const tidalLength = 4;       // ~4 hours per tidal cycle
  const yearLength = 8760;
  const glacialLength = 876000;

  // Tidal: strongest and fastest cycle
  const tidal = Math.sin(2 * Math.PI * tick / tidalLength);

  // Seasonal rain: offset from temperature seasonal
  const rain = Math.sin(2 * Math.PI * tick / yearLength + Math.PI / 3);

  // Glacial: ice ages lock up water → lower sea level
  const glacial = -Math.sin(2 * Math.PI * tick / glacialLength);

  const noise = 0.05 * (Math.random() - 0.5);

  // range [-1, +1]
  // -1 = fully exposed (extreme low tide + dry season + ice age)
  // +1 = fully submerged (high tide + rainy season + interglacial)
  return 0.55 * tidal
       + 0.25 * rain
       + 0.15 * glacial
       + noise;
}
```

### Environment → Reaction Rates

All reaction probabilities are continuous functions of temperature and waterLevel:

```javascript
function getReactionRates(temperature, waterLevel) {
  // Helpers: clamp to [0,1]
  const warm = Math.max(0, temperature);      // 0 when cold, up to 1
  const cold = Math.max(0, -temperature);     // 0 when warm, up to 1
  const wet = Math.max(0, waterLevel);        // 0 when dry, up to 1
  const dry = Math.max(0, -waterLevel);       // 0 when wet, up to 1

  return {
    // Movement
    diffusion: 0.1 + 0.7 * wet + 0.2 * warm,
    // High when wet+warm, low when dry+cold

    // Surface interaction
    adsorption: 0.2 + 0.6 * dry,
    // Particles stick to surface when dry
    desorption: 0.1 + 0.7 * wet + 0.2 * warm,
    // Particles release when wet+warm

    // Hydrogen bonding (complementary pairing)
    hBondForm: 0.2 + 0.7 * cold,
    // Forms in cold (annealing)
    hBondBreakAU: 0.1 + 0.8 * warm,
    // A-U bonds break easily when warm
    hBondBreakGC: 0.05 + 0.5 * warm,
    // G-C bonds are stronger, need more heat

    // Backbone bond (polymerization / chain growth)
    backboneForm: 0.4 * warm * dry,
    // Dehydration condensation: needs heat AND dryness
    backboneFormTemplated: 0.6 * warm * dry,
    // Higher when on a template strand
    backboneFormSurface: 0.5 * warm * dry,
    // Higher on mineral surface (clay catalysis)

    // Hydrolysis (chain breaking)
    hydrolysis: 0.002 + 0.005 * wet * warm,
    // Slow, but faster in hot water

    // Amino acid attachment to RNA
    aminoAcidAttach: 0.15 * cold * wet,
    // Needs solution (wet) and stability (cold)

    // Peptide bond formation
    peptideBond: 0.02 * warm * dry,
    // Dehydration condensation, like backbone
    peptideBondCatalyzed: 0.25 * warm * dry,
    // Much higher near ribozyme

    // Lipid self-assembly
    lipidAssembly: 0.3 * wet,
    // Hydrophobic effect requires water

    // Mutation during replication
    mutationRate: 0.02 * warm,
    // Higher temperature = more errors

    // Degradation of free monomers
    degradation: 0.001,
    // Constant low rate

    // Supply of new monomers (ocean input during high tide)
    monomerSupply: 0.3 * Math.max(0, wet - 0.5),
    // Only when significantly wet (high tide brings material)
  };
}
```

---

## Reaction Rules

### R1: Diffusion
```
Free monomers (concentration field):
  - Diffuse along concentration gradient
  - Rate proportional to diffusion coefficient
  - No diffusion when waterLevel < -0.8 (completely dry)

Structures (chains, peptides, lipids):
  - Move to random adjacent cell
  - Probability = diffusion rate / sqrt(structure_length)
  - Longer structures move slower
  - Surface-bound structures don't move (until desorbed)
```

### R2: Surface Adsorption / Desorption
```
Adsorption (solution → surface):
  - Any structure in a cell with waterLevel dropping below threshold
  - Probability = adsorption rate
  - Once adsorbed: movement stops, reaction rates change

Desorption (surface → solution):
  - When waterLevel rises (tide comes in)
  - Probability = desorption rate
  - Longer chains are harder to desorb (prob / sqrt(length))
```

### R3: Hydrogen Bonding (Complementary Pairing)
```
Formation:
  - A-U and G-C only
  - Two RNA strands in same cell, antiparallel alignment
  - Match complementary bases along the length
  - Probability per base pair = hBondForm rate
  - Partial pairing OK (not all bases need to match)

Breaking:
  - A-U bonds: probability = hBondBreakAU rate
  - G-C bonds: probability = hBondBreakGC rate (lower)
  - If >50% of bonds in a duplex break → strands separate
```

### R4: Backbone Bond (RNA Polymerization)
```
Non-templated (random polymerization):
  - Two nucleotides in same cell
  - At least one is surface-bound OR local concentration > threshold
  - Consumes 1 nucleotide from concentration field
  - Appends to existing chain (or creates new 2-mer)
  - Probability = backboneFormSurface (if surface) or backboneForm * 0.1

Templated (copying):
  - Free nucleotide in same cell as a template strand
  - Nucleotide is complementary to the next uncopied base on template
  - Consumes 1 nucleotide from concentration field
  - Appends to growing complementary strand
  - Probability = backboneFormTemplated
  - On wrong base: probability drops to backboneFormTemplated * mutationRate
```

### R5: Mutation
```
During templated polymerization (R4):
  - Instead of correct complementary base, wrong base is added
  - Probability = mutationRate
  - Random selection among the 3 wrong bases
```

### R6: Degradation and Supply
```
Degradation:
  - Free monomers (in concentration field): reduce by degradation rate
  - Unbound short chains (length < 3): degrade with probability 0.01
  - Bound/long chains: protected (probability 0.0001)

Supply (ocean input):
  - At grid edges, when waterLevel > 0.5 (high tide):
  - Add random nucleotides, amino acids, phosphate, fatty acids
  - Rate = monomerSupply
  - Simulates tidal delivery of raw materials
```

### R7: Codon-Amino Acid Affinity [ENCODED RULE]
```
This is a physical chemistry law, not emergent behavior.

When an RNA chain has a free (non-backbone-bonded) trinucleotide sequence
AND the matching amino acid is available in the concentration field:

Codon → Amino Acid mapping:
  GCC → Ala (Alanine)
  GGC → Gly (Glycine)
  GUC → Val (Valine)
  GAC → Asp (Aspartate)
  GAG → Glu (Glutamate)

  (Based on stereochemical hypothesis — Yarus et al. 2009
   These specific triplets have measurable chemical affinity
   for these amino acids, independent of any biological machinery)

Process:
  1. Scan RNA chain for codon triplets
  2. If matching amino acid exists in cell concentration
  3. Probability = aminoAcidAttach rate
  4. Amino acid attaches to RNA at that codon position
  5. Consumes amino acid from concentration field
```

### R8: Peptide Bond Formation
```
When two amino acids are attached to the SAME RNA chain
at adjacent codon positions:

Without catalyst:
  - Probability = peptideBond rate (very low)

With ribozyme catalyst nearby (same cell):
  - Probability = peptideBondCatalyzed rate (much higher)

Process:
  1. Peptide bond forms between adjacent amino acids
  2. New peptide chain is created (or extended)
  3. Peptide detaches from RNA template
  4. This is primitive translation
```

### R9: Ribozyme Activation [ENCODED RULE]
```
This is a physical chemistry law — sequence determines structure determines function.

When an RNA chain reaches sufficient length AND contains specific motifs:

Ribozyme types:
  Length >= 15 AND contains "GGCGCC":
    → peptidyl_transferase
    → Catalyzes peptide bond formation (R8) in same cell
    → strength = GC_ratio of full sequence

  Length >= 15 AND contains "CCCUUU":
    → rna_replicase
    → Catalyzes backbone bond formation (R4) in same cell
    → strength = GC_ratio of full sequence

  Length >= 20 AND contains "GGGAAACCC":
    → aminoacyl_transferase
    → Catalyzes amino acid attachment (R7) in same cell
    → strength = GC_ratio of full sequence

Catalytic effect:
  - Multiply the base reaction probability by (1 + 5 * strength)
  - Range: 1x (no catalyst) to ~6x (strong ribozyme)

  Higher GC content = more stable structure = better catalyst
  This is chemically justified: GC-rich RNA folds more stably
```

### R10: Lipid Self-Assembly
```
Fatty acids (FA) in concentration field exhibit hydrophobic effect:

Step 1: Lipid nucleation
  - When FA concentration > threshold in a wet cell
  - Spontaneously create individual Lipid particles
  - Head angle = toward nearest water (perpendicular to surface)

Step 2: Lipid aggregation
  - Lipids in adjacent cells attract each other (tail-to-tail)
  - Align: tails face inward, heads face outward (toward water)
  - When >= 8 lipids form a connected ring → Membrane (vesicle)

Step 3: Vesicle dynamics
  - Vesicle has center, radius, integrity
  - Interior = enclosed space
  - Molecules inside cannot diffuse out (unless integrity drops)
  - Vesicle grows by incorporating more lipids
  - Vesicle splits if radius > threshold (mechanical instability)

Tide interaction:
  - Wet phase: vesicles stable (hydrophobic effect strong)
  - Dry phase: vesicles destabilize (no water = no hydrophobic driving force)
  - This creates natural "breathing" — vesicles form when wet,
    release contents when dry, reform when wet again
  - Replicating molecules can spread to new vesicles this way
```

### R11: Compartmentalization Effects
```
When a vesicle (membrane) encloses structures:

  - Internal concentrations are isolated from external
  - Catalytic effects are concentrated (ribozyme only helps what's inside)
  - Favorable mutations benefit only the compartment they're in
  - → Natural selection at the protocell level

  - If a vesicle contains:
    - Self-replicating RNA (via rna_replicase ribozyme)
    - Translation machinery (aminoacyl_transferase + peptidyl_transferase)
    - Growing peptides
    - Lipid production (some peptides could catalyze FA synthesis — stretch)
  - → This is a protocell
  - → The simulation's ultimate goal
```

---

## Emergence Milestones (observe, don't code)

```
Milestone 1: Base pairing
  A-U and G-C pairs form during cold phases
  Metric: count of hydrogen-bonded pairs

Milestone 2: Oligomers
  Short RNA chains (3-10 nucleotides) form
  Metric: chain length distribution

Milestone 3: Template copying
  A strand serves as template → complementary strand grows
  Metric: count of chains with >80% complementarity to another chain

Milestone 4: Self-replication
  Strand A templates strand B, strand B templates strand A
  Metric: detect complementary pair cycles

Milestone 5: Ribozyme emergence
  An RNA chain gains catalytic function
  Metric: count of active ribozymes by type

Milestone 6: Primitive translation
  Amino acids attach to RNA via codon matching → peptide forms
  Metric: count of peptides, length distribution

Milestone 7: Vesicle formation
  Lipids self-assemble into closed membranes
  Metric: count of vesicles, average radius

Milestone 8: Compartmentalized replication
  A vesicle encloses replicating RNA
  Metric: vesicles containing RNA chains

Milestone 9: Protocell ★
  A vesicle with replication + translation + catalysis
  Metric: vesicles containing all three ribozyme types
  THIS IS THE GOAL
```

---

## Technical Implementation

### Recommended: Pure JavaScript + React

Single codebase, immediate visual feedback, no WebSocket complexity.
Can port simulation core to Web Workers for performance.

### File Structure
```
v2/
├── index.html
├── package.json
├── src/
│   ├── App.jsx            # Main React component
│   ├── simulation/
│   │   ├── world.js       # World state, grid, tick loop
│   │   ├── environment.js # Temperature, waterLevel calculations
│   │   ├── concentrations.js # Concentration field diffusion
│   │   ├── rna.js         # RNA chain operations
│   │   ├── peptide.js     # Peptide chain operations
│   │   ├── lipid.js       # Lipid particles and membranes
│   │   ├── reactions.js   # All reaction rules (R1-R11)
│   │   ├── ribozyme.js    # Ribozyme detection and catalysis
│   │   └── metrics.js     # Milestone tracking, statistics
│   ├── renderer/
│   │   ├── canvas.js      # Main canvas rendering
│   │   ├── colors.js      # Color schemes for particles/phases
│   │   └── overlays.js    # Water level overlay, stats display
│   └── ui/
│       ├── controls.js    # Play/pause, speed, parameter sliders
│       ├── stats.js       # Real-time statistics panel
│       └── milestones.js  # Milestone achievement display
├── styles.css
└── README.md
```

### Visualization

```
Main Canvas (top-down surface view):
  Background: shifts with temperature (blue-cold to red-warm)
  Water overlay: translucent blue, height = waterLevel

  Nucleotides in concentration field:
    Subtle colored dots proportional to concentration
    A=red, U=cyan, G=green, C=yellow

  RNA chains: colored line segments following sequence
    Backbone bonds: thick solid lines
    Hydrogen bonds: thin dotted lines to paired strand
    Ribozymes: glowing outline

  Amino acids attached to RNA: small squares on the chain

  Peptides: different colored line segments
    Gly=white, Ala=pink, Val=orange, Asp=red, Glu=purple

  Lipids: small arrows (head=circle, tail=line)
    Membrane: connected lipids forming arc/circle
    Vesicle interior: slightly highlighted area

  Protocell: golden glow around vesicle containing all components

Stats Panel:
  Environment: temperature bar, water level bar, current "time of day"
  Population counts: RNA chains, peptides, lipids, membranes
  Chain lengths: histogram
  Ribozymes: count by type
  Milestones: checklist with achievement timestamps
  Cycle visualization: circular diagram showing current position in day/tide cycles

Controls:
  Play / Pause / Step (1 tick)
  Speed: ticks per second (1-1000)
  Time compression: skip to next interesting event
  Camera: zoom, pan
  Spawn rates (edge sliders for each monomer type)
  Grid size
  Reset

  Advanced:
    Cycle parameters (day length, tidal period, etc.)
    Reaction rate multipliers
    Toggle individual rules on/off (for debugging)
```

---

## Build Order

### Phase 1: Foundation (MVP)
```
1. Grid + concentration field + rendering
   See monomers on a 2D surface with colored dots

2. Diffusion
   Concentrations spread along gradients

3. Environment cycles
   Temperature and waterLevel oscillating
   Background color shifting, water overlay rising/falling

4. Surface adsorption/desorption
   Monomers stick when dry, release when wet
   Visual: dots get "fixed" to surface vs floating

5. RNA polymerization (R4 non-templated)
   Short chains start forming on surface during warm+dry
   Visual: colored line segments appear

6. Hydrogen bonding (R3)
   Complementary strands pair during cold
   Visual: dotted lines between paired strands

7. Stats panel
   Basic counts and chain length histogram
```

### Phase 2: Replication
```
8. Templated polymerization (R4 templated)
   Complementary strand grows along template

9. Warm/cold strand separation
   Duplexes separate during warm phase
   Each strand becomes template
   = REPLICATION (emergence, not coded)

10. Mutation (R5)
    Variation in sequences

11. Degradation + Supply (R6)
    Population dynamics, selection pressure

12. Milestone detection (1-4)
```

### Phase 3: Translation
```
13. Codon-amino acid attachment (R7)
    Amino acids bind to RNA at codon positions

14. Ribozyme activation (R9)
    Long RNA with motifs gain catalytic function
    Visual: glowing outline on ribozymes

15. Peptide bond formation (R8)
    Amino acids link into peptide chains
    = PRIMITIVE TRANSLATION (emergence from rules)

16. Milestone detection (5-6)
```

### Phase 4: Compartmentalization
```
17. Lipid particle creation from concentration field

18. Lipid self-assembly (R10)
    Hydrophobic clustering, alignment

19. Vesicle formation
    Lipids form closed membranes
    Visual: circular structures with interior

20. Compartmentalization effects (R11)
    Interior isolation, concentrated catalysis

21. Vesicle dynamics
    Growth, splitting, tidal breathing

22. Milestone detection (7-9)
    PROTOCELL = victory condition
```

### Phase 5: Polish
```
23. Performance optimization (Web Workers)
24. Export/save interesting runs
25. Time-lapse recording
26. Parameter presets ("Earth-like", "Mars-like", "Europa-like")
27. README with scientific references
```

---

## Key Reminders for Claude Code

1. **Fresh start in v2/.** Do not touch existing files.
2. **Build in order.** Each phase depends on the previous one working.
3. **Visual feedback at every step.** If you can't see it, it's not working.
4. **Commit after each numbered step** with descriptive messages.
5. **Test emergence, not correctness.** There's no "right answer." Watch what happens.
6. **Performance matters.** 200x200 grid with hundreds of structures must run smoothly.
7. **The simulation should be mesmerizing to watch.** If it's boring, something is wrong.
8. **Use typed arrays for concentration fields** (Float32Array for speed).
9. **Spatial hashing for structure queries** (which structures are in this cell?).
10. **The fun is in the emergence.** Don't skip to Phase 4. Let each phase surprise you.

---

## Scientific References (for context, not implementation)

- Miller-Urey experiment (1953): Prebiotic amino acid synthesis
- Szostak lab: RNA replication and protocell formation
- Deamer: Lipid vesicles and tidal cycling in abiogenesis
- Yarus (2009): Stereochemical basis for codon-amino acid affinity
- Lincoln & Joyce (2009): Self-replicating RNA systems
- Ferris (2006): Montmorillonite clay catalysis of RNA polymerization
- Oparin-Haldane hypothesis: Primordial soup
- RNA World hypothesis: RNA preceded DNA and proteins
