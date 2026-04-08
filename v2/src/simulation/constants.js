// BIOTRON v2 — shared constants

// Concentration field species (in this canonical order, used as Float32Array indices)
export const SPECIES = [
  'A',   // Adenine
  'U',   // Uracil (RNA world — not Thymine)
  'G',   // Guanine
  'C',   // Cytosine
  'Gly', // Glycine
  'Ala', // Alanine
  'Val', // Valine
  'Asp', // Aspartate
  'Glu', // Glutamate
  'P',   // Phosphate
  'Mg',  // Magnesium
  'FA',  // Fatty acid
];

export const SPECIES_INDEX = Object.fromEntries(SPECIES.map((s, i) => [s, i]));
export const N_SPECIES = SPECIES.length;

// Display colors per species (used by the renderer)
export const SPECIES_COLOR = {
  A:   [255, 90, 90],     // red
  U:   [120, 220, 230],   // cyan
  G:   [110, 220, 120],   // green
  C:   [255, 215, 90],    // yellow
  Gly: [240, 240, 240],   // white
  Ala: [255, 180, 200],   // pink
  Val: [255, 160, 70],    // orange
  Asp: [220, 80, 80],     // dark red
  Glu: [180, 100, 220],   // purple
  P:   [200, 180, 100],   // tan
  Mg:  [120, 200, 200],   // teal
  FA:  [200, 160, 100],   // sandy
};

// Time scale (Phase 1: just oscillation defaults — full cycles in environment.js)
export const TICKS_PER_DAY  = 8;
export const TICKS_PER_TIDE = 4;

// Default initial concentrations (uniform random fill at world creation)
export const INITIAL_CONCENTRATION = 0.10;
