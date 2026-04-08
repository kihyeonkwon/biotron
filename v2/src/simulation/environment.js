// Environment cycles — temperature + water level
//
// Both range over [-1, +1]. Driven by overlapping sinusoidal cycles
// (diurnal / tidal / seasonal / glacial). Phase 1 uses a simplified
// stack; full multi-scale cycles arrive in Phase 1 step 3.

const TWO_PI = Math.PI * 2;

// Cycle lengths (in ticks)
export const CYCLES = {
  day: 8,
  tide: 4,
  year: 8760,
  glacial: 876_000,
};

export function getTemperature(tick) {
  const diurnal  = Math.sin(TWO_PI * tick / CYCLES.day);
  const seasonal = Math.sin(TWO_PI * tick / CYCLES.year);
  const glacial  = Math.sin(TWO_PI * tick / CYCLES.glacial);
  const noise    = 0.05 * (Math.random() - 0.5);
  return clamp(0.30 * diurnal + 0.35 * seasonal + 0.25 * glacial + noise, -1, 1);
}

export function getWaterLevel(tick) {
  const tidal   = Math.sin(TWO_PI * tick / CYCLES.tide);
  const rain    = Math.sin(TWO_PI * tick / CYCLES.year + Math.PI / 3);
  const glacial = -Math.sin(TWO_PI * tick / CYCLES.glacial);
  const noise   = 0.05 * (Math.random() - 0.5);
  return clamp(0.55 * tidal + 0.25 * rain + 0.15 * glacial + noise, -1, 1);
}

export function getReactionRates(temperature, waterLevel) {
  const warm = Math.max(0, temperature);
  const cold = Math.max(0, -temperature);
  const wet  = Math.max(0, waterLevel);
  const dry  = Math.max(0, -waterLevel);

  return {
    diffusion:               0.10 + 0.70 * wet + 0.20 * warm,
    adsorption:              0.20 + 0.60 * dry,
    desorption:              0.10 + 0.70 * wet + 0.20 * warm,
    hBondForm:               0.20 + 0.70 * cold,
    hBondBreakAU:            0.10 + 0.80 * warm,
    hBondBreakGC:            0.05 + 0.50 * warm,
    backboneForm:            0.40 * warm * dry,
    backboneFormTemplated:   0.60 * warm * dry,
    backboneFormSurface:     0.50 * warm * dry,
    hydrolysis:              0.002 + 0.005 * wet * warm,
    aminoAcidAttach:         0.15 * cold * wet,
    peptideBond:             0.02 * warm * dry,
    peptideBondCatalyzed:    0.25 * warm * dry,
    lipidAssembly:           0.30 * wet,
    mutationRate:            0.02 * warm,
    degradation:             0.001,
    monomerSupply:           0.30 * Math.max(0, wet - 0.5),
  };
}

function clamp(x, lo, hi) {
  return x < lo ? lo : x > hi ? hi : x;
}
