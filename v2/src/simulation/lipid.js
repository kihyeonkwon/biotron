// Lipid particles + Membrane (vesicle) structures.
//
// A Lipid is an individual particle on the surface. It has a head angle
// (orientation toward water), and may belong to a membrane.
//
// A Membrane is a closed cluster of lipids forming a vesicle. We model it
// as a center + radius (circular) for simplicity. A structure is "inside"
// the membrane if its position is within radius. Free monomer diffusion
// across the membrane boundary is attenuated (semi-permeable, R11).

let _nextLipidId = 200000;
let _nextMembraneId = 300000;

export function makeLipid(x, y, headAngle = 0) {
  return {
    type: 'lipid',
    id: _nextLipidId++,
    position: { x, y },
    headAngle,
    membraneId: null,
    age: 0,
  };
}

export function makeMembrane(centerX, centerY, lipidIds, generation = 0, parentId = null) {
  // Larger radius so each vesicle can capture multiple ribozymes.
  // 8 lipids → ~8 cell radius (~200 cells interior),
  // 40 lipids → ~18 cell radius (~1000 cells interior).
  const radius = Math.max(6, Math.sqrt(lipidIds.length / Math.PI) * 5);
  return {
    type: 'membrane',
    id: _nextMembraneId++,
    lipids: lipidIds.slice(),
    center: { x: centerX, y: centerY },
    radius,
    enclosed: [],          // populated each tick from world.structures
    integrity: 1.0,
    age: 0,
    generation,            // vesicle lineage depth
    parentId,              // parent vesicle id (null = spontaneous)
    childCount: 0,         // how many times this vesicle has divided
  };
}

export function isInsideMembrane(membrane, x, y, width, height) {
  // Toroidal-aware distance
  let dx = Math.abs(x - membrane.center.x);
  let dy = Math.abs(y - membrane.center.y);
  if (dx > width / 2) dx = width - dx;
  if (dy > height / 2) dy = height - dy;
  return dx * dx + dy * dy <= membrane.radius * membrane.radius;
}

export function resetLipidIds() {
  _nextLipidId = 200000;
  _nextMembraneId = 300000;
}
