"""BIOTRON world: 2D grid of cells with diffusing monomer pool and strands.

Step pipeline (in order):
  1. Diffuse free monomers across the grid (random walk)
  2. (optional) background monomer inflow
  3. Spontaneous dimer formation in occupied cells
  4. For each strand:
       a. Form H-bonds at empty positions (cold-favored, stacking-boosted)
       b. Polymerize daughter cov bonds between adjacent anchored monomers
       c. Boltzmann break of H-bonds
       d. Boltzmann break of daughter cov bonds
       e. Boltzmann break of parent cov bonds
       f. Lift-off: cov-connected daughter components with no anchors detach
       g. Fragmentation: split parent if any parent cov broke
"""
from dataclasses import dataclass
import numpy as np

from . import physics
from .strand import Strand, reset_ids


@dataclass
class WorldConfig:
    grid_h: int = 64
    grid_w: int = 64

    # Initial monomer pool
    initial_mono_per_cell: int = 12   # per type → 48 monomers per cell total
    inflow_per_cell: int = 0          # background source per type per step
    max_per_cell: int = 60            # per type cap (only enforced if inflow > 0)

    # Diffusion
    diffuse_p: float = 0.10           # prob per monomer per step to move

    # Spontaneous dimer formation
    spont_rate: float = 0.0008        # attempts per cell per step

    # H-bond formation
    hbond_base_rate: float = 0.6      # base prob per attempt
    stacking_bonus: float = 3.0       # multiplier per anchored neighbor

    # Daughter polymerization (covalent bond formation between adjacent daughter monomers)
    poly_base_rate: float = 0.04

    # End extension: free monomer joins a strand end (template-free covalent growth)
    extend_rate: float = 0.02

    # Bond energies (None = use physics module defaults; set to override per-experiment)
    # Regime α (default): weak H-bond, constant warm works.
    # Regime β: strong H-bond (~0.80), cycling required.
    e_hb: float           = physics.E_HB
    e_stack: float        = physics.E_STACK
    e_cov_parent: float   = physics.E_COV_PARENT
    e_cov_daughter: float = physics.E_COV_DAUGHTER

    # Birth detection
    min_birth_length: int = 3         # min daughter length to count as birth

    # Stage 1: mutations
    # Per-position substitution probability applied to daughter monomers at liftoff.
    # Default 0 → byte-identical regression with v1 behavior.
    # Recommended: 0.01 (gives ~10% of L=10 daughters carrying ≥1 mutation,
    # safely below the error catastrophe threshold ~1/L = 0.1).
    p_mut: float = 0.0

    # Stage 2: catalysis
    # Per-cell catalytic boost on daughter polymerization.
    # Effective rate = poly_base_rate * (1 + k_poly_boost * cell_polymerize_strength)
    # Default 0 → no catalytic effect (regression).
    k_poly_boost: float = 0.0

    # Stage 3: lipid synthesis (catalysis-driven monomer → lipid conversion)
    # k_lipid: per (cell × lipid-catalyst-strength) Bernoulli rate of converting
    # one free monomer in that cell into one lipid unit.
    # Default 0 → no lipid synthesis (regression with Stage 2).
    k_lipid: float = 0.0
    # Lipid diffuses much slower than free monomers (lets blobs accumulate).
    lipid_diffuse_p: float = 0.02
    # Chemotactic bias: fraction of moves that go toward the highest-lipid neighbor
    # instead of uniform-random. 0 = pure random walk; 1 = pure attraction.
    lipid_chemotaxis: float = 0.3
    # Lipid bond energy for Boltzmann decay. Higher = more stable.
    e_lipid: float = 1.0
    # Per-step base rate of lipid → free monomer conversion attempt
    # (multiplied by exp(-e_lipid/kT)). Keeps lipid budget bounded.
    lipid_decay_attempt_rate: float = 1.0

    # Stage 4: membrane self-assembly
    # A cell becomes "lipid-filled" when its lipid count crosses HIGH; it stays
    # filled until lipid drops below LOW (hysteresis prevents flicker).
    # is_membrane[x,y] = lipid_filled[x,y] AND any 4-neighbor is NOT lipid_filled.
    # This "surface-only" rule automatically gives ring topology around lipid blobs.
    # Default thresholds disabled (HIGH=999) → no membranes (regression).
    membrane_threshold_high: int = 999
    membrane_threshold_low: int = 999

    # Stage 5: membrane as diffusion barrier
    # Fraction of monomer/lipid moves blocked when destination is a membrane cell.
    # 0 = membranes don't impede diffusion (regression with Stage 4).
    # 0.95 = realistic seal with 5% leak (avoids frozen states, biologically reasonable).
    membrane_seal: float = 0.0

    # Limits
    max_strand_length: int = 24
    max_strands: int = 4000


class World:
    def __init__(self, cfg: WorldConfig = None, seed: int = None):
        self.cfg = cfg or WorldConfig()
        if seed is not None:
            np.random.seed(seed)
        reset_ids()

        H, W = self.cfg.grid_h, self.cfg.grid_w
        # Free monomer pool: (H, W, 4 types)
        self.free = np.full((H, W, physics.N_TYPES),
                            self.cfg.initial_mono_per_cell, dtype=np.int32)
        self.strands: list[Strand] = []
        self.step_count = 0
        self.env_E = 0.5

        # Stage 2: per-cell catalytic strength grids, refreshed once per step.
        self._cat_grid_poly = np.zeros((H, W), dtype=np.int32)
        self._cat_grid_lipid = np.zeros((H, W), dtype=np.int32)

        # Stage 3: lipid pool per cell.
        self.lipid = np.zeros((H, W), dtype=np.int32)

        # Stage 4: membrane state
        # _lipid_filled tracks "filled" with hysteresis (HIGH crosses to True,
        # LOW crosses to False). is_membrane is the surface of filled regions.
        self._lipid_filled = np.zeros((H, W), dtype=bool)
        self.is_membrane = np.zeros((H, W), dtype=bool)

        self.total_births = 0
        self.total_deaths = 0
        self.history = {
            "step": [], "env_E": [], "n_strands": [],
            "total_free": [], "births": [], "deaths": [],
            "max_len": [], "avg_len": [], "n_anchored": [],
        }

    # ─────────── Helpers ───────────

    def total_free(self) -> int:
        return int(self.free.sum())

    def n_strands(self) -> int:
        return len(self.strands)

    def _take(self, x, y, t) -> bool:
        if self.free[x, y, t] > 0:
            self.free[x, y, t] -= 1
            return True
        return False

    def _put(self, x, y, t):
        self.free[x, y, t] += 1

    # ─────────── Step ───────────

    def step(self, env_E: float = None):
        if env_E is not None:
            self.env_E = float(env_E)
        E = self.env_E

        # 1. Diffusion
        self._diffuse()

        # 2. Background inflow (skip if 0)
        if self.cfg.inflow_per_cell > 0:
            self.free += self.cfg.inflow_per_cell
            np.minimum(self.free, self.cfg.max_per_cell, out=self.free)

        # 3. Spontaneous dimer formation
        self._spontaneous_dimers()

        # 3.5 Stage 2: rebuild per-cell catalyst grids from current strand population
        self._rebuild_catalyst_grids()

        # 3.6 Stage 3: lipid subsystem (synthesis + diffusion + decay)
        # Order: synthesize first (catalyst-driven), then diffuse, then decay.
        # Skipping when k_lipid == 0 keeps regression with Stage 2.
        if self.cfg.k_lipid > 0:
            self._synthesize_lipids()
        if self.lipid.any():
            self._diffuse_lipids()
            self._decay_lipids(E)

        # 3.7 Stage 4: update membrane state from current lipid distribution
        if self.cfg.membrane_threshold_high < 999 or self._lipid_filled.any():
            self._update_membranes()

        # 4. Per-strand step
        births = 0
        new_strands: list[Strand] = []
        dead_strands: list[Strand] = []

        for s in list(self.strands):  # iterate snapshot
            self._strand_end_extend(s, E)
            self._strand_form_hbonds(s, E)
            self._strand_polymerize_daughter(s, E)
            self._strand_break_hbonds(s, E)
            self._strand_break_daughter_cov(s, E)
            broken = self._strand_break_parent_cov(s, E)
            lifted, b = self._strand_liftoff(s)
            new_strands.extend(lifted)
            births += b
            s.repl_count += b
            if broken:
                fragments = self._strand_fragment(s, broken)
                if fragments:
                    new_strands.extend(fragments)
                dead_strands.append(s)
            s.age += 1

        # Apply deaths
        if dead_strands:
            dead_ids = {id(d) for d in dead_strands}
            self.strands = [s for s in self.strands if id(s) not in dead_ids]
            self.total_deaths += len(dead_strands)

        # Apply births (capped)
        for ns in new_strands:
            if len(self.strands) < self.cfg.max_strands:
                self.strands.append(ns)

        self.total_births += births
        self.step_count += 1
        self._record(births, len(dead_strands))

    # ─────────── Step components ───────────

    # ─────────── Stage 3: lipid subsystem ───────────

    def _synthesize_lipids(self):
        """Catalyst-driven free monomer → lipid conversion.

        For each cell with lipid catalyst strength > 0, attempt
        `int(k_lipid * cat_strength)` Bernoulli conversions of one
        random free monomer (weighted by abundance) into one lipid unit.
        Consumes the monomer (negative feedback on local pool).
        """
        cfg = self.cfg
        H, W = cfg.grid_h, cfg.grid_w
        # Iterate only over cells with non-zero lipid catalysts
        ys, xs = np.where(self._cat_grid_lipid > 0)  # numpy returns row, col
        # Note: my (x,y) convention treats x as row and y as col throughout the file
        for x, y in zip(ys, xs):
            cs = int(self._cat_grid_lipid[x, y])
            n_attempts = max(1, int(cfg.k_lipid * cs))
            for _ in range(n_attempts):
                cell = self.free[x, y]
                tot = int(cell.sum())
                if tot < 1:
                    break
                probs = cell.astype(float) / tot
                t = int(np.random.choice(physics.N_TYPES, p=probs))
                if not self._take(x, y, t):
                    break
                self.lipid[x, y] += 1

    def _diffuse_lipids(self):
        """Slow random walk with optional chemotactic bias toward lipid-rich neighbors."""
        cfg = self.cfg
        p = cfg.lipid_diffuse_p
        if p <= 0:
            return
        H, W = cfg.grid_h, cfg.grid_w

        # Step 1: how many lipid units in each cell try to move
        moving = np.random.binomial(self.lipid, p)
        if not moving.any():
            return

        # Step 2: split into 4 directions. Without chemotaxis: uniform 1/4.
        # With chemotaxis: a fraction `lipid_chemotaxis` go preferentially toward
        # the highest-lipid neighbor; the rest are uniform.
        chemo = float(cfg.lipid_chemotaxis)
        if chemo > 0:
            # Compute neighbor-lipid for each direction (use current self.lipid as bias)
            # up neighbor of cell (i,j) is (i-1,j); we use np.roll-style with padding 0
            up_n   = np.zeros_like(self.lipid); up_n[1:, :]   = self.lipid[:-1, :]
            dn_n   = np.zeros_like(self.lipid); dn_n[:-1, :]  = self.lipid[1:, :]
            lf_n   = np.zeros_like(self.lipid); lf_n[:, 1:]   = self.lipid[:, :-1]
            rt_n   = np.zeros_like(self.lipid); rt_n[:, :-1]  = self.lipid[:, 1:]
            # For each cell, find the direction with the highest neighbor lipid
            # (ties broken by argmax, which is deterministic).
            stack = np.stack([up_n, dn_n, lf_n, rt_n], axis=-1)  # (H,W,4)
            best_dir = np.argmax(stack, axis=-1)  # (H,W) ∈ {0..3}

        # For simplicity: split moving into chemo_count + random_count, then
        # the chemo_count all go in best_dir, the random_count split 1/4 each.
        if chemo > 0:
            chemo_n = np.random.binomial(moving, chemo)
            random_n = moving - chemo_n
        else:
            chemo_n = np.zeros_like(moving)
            random_n = moving

        # Random part: split 4 ways via successive binomial
        n_up = np.random.binomial(random_n, 0.25)
        rest = random_n - n_up
        n_dn = np.random.binomial(rest, 1.0 / 3.0)
        rest -= n_dn
        n_lf = np.random.binomial(rest, 0.5)
        n_rt = rest - n_lf

        # Add chemotactic part to the corresponding direction
        if chemo > 0:
            n_up += chemo_n * (best_dir == 0)
            n_dn += chemo_n * (best_dir == 1)
            n_lf += chemo_n * (best_dir == 2)
            n_rt += chemo_n * (best_dir == 3)

        # Apply moves (subtract from origin, add to destination, edge cells keep their own)
        self.lipid -= moving  # everyone tried to move
        L = self.lipid
        L[:-1, :] += n_up[1:, :]
        L[0, :]   += n_up[0, :]
        L[1:, :]  += n_dn[:-1, :]
        L[-1, :]  += n_dn[-1, :]
        L[:, :-1] += n_lf[:, 1:]
        L[:, 0]   += n_lf[:, 0]
        L[:, 1:]  += n_rt[:, :-1]
        L[:, -1]  += n_rt[:, -1]

    def _update_membranes(self):
        """Stage 4: refresh _lipid_filled (with hysteresis) then is_membrane (surface-only).

        Hysteresis avoids flickering when lipid count fluctuates around a single
        threshold:
          - lipid >= HIGH      → filled (regardless of previous state)
          - lipid < LOW        → not filled
          - LOW <= lipid < HIGH → keep previous state

        is_membrane is the boundary of any filled region: a filled cell that has
        at least one non-filled 4-neighbor (or sits on the grid edge — those count
        as having an "outside" non-filled neighbor).
        """
        cfg = self.cfg
        HIGH = cfg.membrane_threshold_high
        LOW = cfg.membrane_threshold_low

        prev_filled = self._lipid_filled
        new_filled = np.where(
            self.lipid >= HIGH, True,
            np.where(self.lipid < LOW, False, prev_filled)
        )
        self._lipid_filled = new_filled

        # Surface-only rule: filled AND has any non-filled 4-neighbor
        # Treat off-grid as non-filled (so filled cells on the boundary are membrane)
        H, W = new_filled.shape
        # Build "any non-filled neighbor" mask
        has_non_filled_n = np.zeros_like(new_filled)
        # Up neighbor (i-1): off-grid for row 0 → counts as non-filled
        has_non_filled_n[1:, :] |= ~new_filled[:-1, :]
        has_non_filled_n[0, :] = True
        # Down neighbor
        has_non_filled_n[:-1, :] |= ~new_filled[1:, :]
        has_non_filled_n[-1, :] = True
        # Left
        has_non_filled_n[:, 1:] |= ~new_filled[:, :-1]
        has_non_filled_n[:, 0] = True
        # Right
        has_non_filled_n[:, :-1] |= ~new_filled[:, 1:]
        has_non_filled_n[:, -1] = True

        self.is_membrane = new_filled & has_non_filled_n

    def _decay_lipids(self, E):
        """Boltzmann decay: each lipid unit has prob exp(-e_lipid/kT) of dissolving
        per attempt. On dissolution, returns one random monomer to the local cell."""
        cfg = self.cfg
        if cfg.lipid_decay_attempt_rate <= 0:
            return
        p_break = physics.break_prob(cfg.e_lipid, E) * cfg.lipid_decay_attempt_rate
        if p_break <= 0:
            return
        p_break = min(0.95, p_break)
        # Vectorized: sample binomial(lipid_count, p_break) per cell to count decays
        decays = np.random.binomial(self.lipid, p_break)
        if not decays.any():
            return
        self.lipid -= decays
        # Return decayed lipids as random monomer types
        # For simplicity: return as a uniform random distribution across the 4 types
        # (could make it depend on local pool composition, but uniform is fine for v1)
        H, W = cfg.grid_h, cfg.grid_w
        # Use multinomial-style: each decayed lipid becomes a random type
        # Vectorize via 4 binomial splits
        n_t0 = np.random.binomial(decays, 0.25)
        rest = decays - n_t0
        n_t1 = np.random.binomial(rest, 1.0 / 3.0)
        rest -= n_t1
        n_t2 = np.random.binomial(rest, 0.5)
        n_t3 = rest - n_t2
        self.free[..., 0] += n_t0
        self.free[..., 1] += n_t1
        self.free[..., 2] += n_t2
        self.free[..., 3] += n_t3

    def _rebuild_catalyst_grids(self):
        """Stage 2 + Stage 3: aggregate per-cell catalytic strengths from all strands.
        Each grid is skipped when its corresponding rate is 0 (regression-friendly)."""
        cfg = self.cfg
        need_poly = cfg.k_poly_boost > 0
        need_lipid = cfg.k_lipid > 0
        if not need_poly and not need_lipid:
            return
        if need_poly:
            self._cat_grid_poly.fill(0)
        if need_lipid:
            self._cat_grid_lipid.fill(0)
        for s in self.strands:
            if need_poly and s._cat_poly:
                self._cat_grid_poly[s.x, s.y] += s._cat_poly
            if need_lipid and s._cat_lipid:
                self._cat_grid_lipid[s.x, s.y] += s._cat_lipid

    def _diffuse(self):
        """Random walk: each free monomer moves with prob diffuse_p, equally to 4 neighbors.
        Boundary: monomers that would leave the grid stay in place.

        Stage 5: when a destination cell is a membrane, the move is blocked with
        probability `membrane_seal`. Blocked monomers stay in their origin cell.
        """
        p = self.cfg.diffuse_p
        if p <= 0:
            return
        seal = self.cfg.membrane_seal
        gate = seal > 0 and self.is_membrane.any()

        moving = np.random.binomial(self.free, p)
        self.free -= moving
        # Split into 4 directions via successive binomial (each direction 1/4)
        n_up = np.random.binomial(moving, 0.25)
        rest = moving - n_up
        n_dn = np.random.binomial(rest, 1.0 / 3.0)
        rest -= n_dn
        n_lf = np.random.binomial(rest, 0.5)
        n_rt = rest - n_lf

        if gate:
            # For each direction, find moves whose destination is a membrane cell
            # and apply membrane_seal as a Bernoulli rejection per monomer.
            # We approximate by computing per-cell pass fractions.
            # Origin cell (i,j) sends "up" to (i-1,j). The destination is membrane if
            # is_membrane[i-1, j]. Blocked moves stay at (i,j).
            mem = self.is_membrane
            # For each direction, build a (H,W,K)-broadcastable mask of "destination
            # is membrane" indexed by ORIGIN. We need to align dimensions.
            # up: dest of origin (i,j) = (i-1,j); dest_membrane[i,j] = mem[i-1,j] for i>=1
            up_dest_mem = np.zeros_like(mem)
            up_dest_mem[1:, :] = mem[:-1, :]
            dn_dest_mem = np.zeros_like(mem)
            dn_dest_mem[:-1, :] = mem[1:, :]
            lf_dest_mem = np.zeros_like(mem)
            lf_dest_mem[:, 1:] = mem[:, :-1]
            rt_dest_mem = np.zeros_like(mem)
            rt_dest_mem[:, :-1] = mem[:, 1:]
            # Apply seal: blocked moves return to origin (subtract from outgoing,
            # add back to free at origin).
            # We use binomial(n, seal) to count blocked moves per direction.
            for n_arr, dest_mem in (
                (n_up, up_dest_mem),
                (n_dn, dn_dest_mem),
                (n_lf, lf_dest_mem),
                (n_rt, rt_dest_mem),
            ):
                # Blocked count per cell per type
                # Broadcast dest_mem (H,W) over the K-axis
                dest_mem3 = dest_mem[:, :, None]
                blocked = np.where(dest_mem3, np.random.binomial(n_arr, seal), 0)
                n_arr -= blocked
                self.free += blocked  # return to origin

        # Apply shifts. "up" = row index decreases.
        # Cells that can't move (edge) keep their monomers.
        f = self.free
        f[:-1, :, :] += n_up[1:, :, :]
        f[0, :, :]   += n_up[0, :, :]   # row 0 keeps its own "up"
        f[1:, :, :]  += n_dn[:-1, :, :]
        f[-1, :, :]  += n_dn[-1, :, :]
        f[:, :-1, :] += n_lf[:, 1:, :]
        f[:, 0, :]   += n_lf[:, 0, :]
        f[:, 1:, :]  += n_rt[:, :-1, :]
        f[:, -1, :]  += n_rt[:, -1, :]

    def _spontaneous_dimers(self):
        """Pick random cells; if they have >= 2 monomers, form a covalent dimer.
        Sequence is random (weighted by local abundance) — no template involvement."""
        cfg = self.cfg
        if len(self.strands) >= cfg.max_strands:
            return
        H, W = cfg.grid_h, cfg.grid_w
        expected = H * W * cfg.spont_rate
        n_attempts = int(np.random.poisson(expected))
        n_attempts = min(n_attempts, cfg.max_strands - len(self.strands))
        for _ in range(n_attempts):
            x = np.random.randint(H)
            y = np.random.randint(W)
            cell = self.free[x, y]
            tot = cell.sum()
            if tot < 2:
                continue
            probs = cell / tot
            t1 = int(np.random.choice(physics.N_TYPES, p=probs))
            if not self._take(x, y, t1):
                continue
            cell2 = self.free[x, y]
            tot2 = cell2.sum()
            if tot2 < 1:
                self._put(x, y, t1)
                continue
            probs2 = cell2 / tot2
            t2 = int(np.random.choice(physics.N_TYPES, p=probs2))
            if not self._take(x, y, t2):
                self._put(x, y, t1)
                continue
            self.strands.append(Strand([t1, t2], x, y, gen=0))

    def _strand_end_extend(self, s: Strand, E):
        """Template-free covalent growth: a free monomer joins one end of the strand.

        This is just the same physics as spontaneous dimer formation, applied to
        a strand end and a same-cell free monomer. Slow, cold-favored. Does not
        provide any template information — only seeds longer parents."""
        cfg = self.cfg
        cold = max(0.0, 1.0 - E)
        if cold < 0.05:
            return
        if s.length >= cfg.max_strand_length:
            return
        if np.random.random() >= cfg.extend_rate * cold:
            return
        cell = self.free[s.x, s.y]
        tot = cell.sum()
        if tot < 1:
            return
        probs = cell / tot
        t = int(np.random.choice(physics.N_TYPES, p=probs))
        if not self._take(s.x, s.y, t):
            return
        # Append at left or right end (50/50)
        if np.random.random() < 0.5:
            s.mono = np.concatenate([[t], s.mono]).astype(np.int8)
            s.cov = np.concatenate([[True], s.cov]).astype(bool)
            s.d_present = np.concatenate([[False], s.d_present]).astype(bool)
            s.d_mono = np.concatenate([[0], s.d_mono]).astype(np.int8)
            s.d_anchor = np.concatenate([[False], s.d_anchor]).astype(bool)
            s.d_cov = np.concatenate([[False], s.d_cov]).astype(bool)
        else:
            s.mono = np.concatenate([s.mono, [t]]).astype(np.int8)
            s.cov = np.concatenate([s.cov, [True]]).astype(bool)
            s.d_present = np.concatenate([s.d_present, [False]]).astype(bool)
            s.d_mono = np.concatenate([s.d_mono, [0]]).astype(np.int8)
            s.d_anchor = np.concatenate([s.d_anchor, [False]]).astype(bool)
            s.d_cov = np.concatenate([s.d_cov, [False]]).astype(bool)
        # Stage 2: sequence changed → recompute catalytic cache
        s.recompute_catalytic_cache()

    def _strand_form_hbonds(self, s: Strand, E):
        """Try to grab a complementary monomer from the local cell at each empty position."""
        cfg = self.cfg
        cold = max(0.0, 1.0 - E)
        if cold < 0.05:
            return
        L = s.length
        for i in range(L):
            if s.d_present[i]:
                continue
            need = int(physics.COMPLEMENT[int(s.mono[i])])
            local = int(self.free[s.x, s.y, need])
            if local <= 0:
                continue
            stack_n = 0
            if i > 0 and s.d_anchor[i - 1]:
                stack_n += 1
            if i < L - 1 and s.d_anchor[i + 1]:
                stack_n += 1
            stack_factor = (cfg.stacking_bonus ** stack_n) if stack_n > 0 else 1.0
            # Concentration factor (saturating)
            conc = min(1.0, local / 8.0)
            prob = cfg.hbond_base_rate * cold * conc * stack_factor
            if prob > 0.95:
                prob = 0.95
            if np.random.random() < prob and self._take(s.x, s.y, need):
                s.d_present[i] = True
                s.d_mono[i] = need
                s.d_anchor[i] = True

    def _strand_polymerize_daughter(self, s: Strand, E):
        """Form covalent bonds between adjacent daughter monomers.
        Stage 2: rate is boosted by the per-cell polymerase catalyst strength."""
        cfg = self.cfg
        cold = max(0.0, 1.0 - E)
        if cold < 0.05:
            return
        # Stage 2 catalyst boost (1 + k * cell_strength). For k=0 this is a no-op.
        cat_boost = 1.0
        if cfg.k_poly_boost > 0:
            cat_boost = 1.0 + cfg.k_poly_boost * float(self._cat_grid_poly[s.x, s.y])
        L = s.length
        for i in range(L - 1):
            if s.d_cov[i]:
                continue
            if not (s.d_present[i] and s.d_present[i + 1]):
                continue
            anchor_factor = 1.5 if (s.d_anchor[i] and s.d_anchor[i + 1]) else 1.0
            prob = cfg.poly_base_rate * cold * anchor_factor * cat_boost
            if prob > 0.95:
                prob = 0.95
            if np.random.random() < prob:
                s.d_cov[i] = True

    def _strand_break_hbonds(self, s: Strand, E):
        """Boltzmann break of H-bonds. If the daughter monomer is not held by daughter cov,
        it returns to the local pool."""
        L = s.length
        e_hb = self.cfg.e_hb
        e_stack = self.cfg.e_stack
        for i in range(L):
            if not s.d_anchor[i]:
                continue
            stack = 0
            if i > 0 and s.d_anchor[i - 1]:
                stack += 1
            if i < L - 1 and s.d_anchor[i + 1]:
                stack += 1
            e = e_hb + stack * e_stack
            if np.random.random() < physics.break_prob(e, E):
                s.d_anchor[i] = False
                left_cov  = i > 0     and s.d_cov[i - 1] and s.d_present[i - 1]
                right_cov = i < L - 1 and s.d_cov[i]     and s.d_present[i + 1]
                if not left_cov and not right_cov:
                    self._put(s.x, s.y, int(s.d_mono[i]))
                    s.d_present[i] = False
                    s.d_mono[i] = 0

    def _strand_break_daughter_cov(self, s: Strand, E):
        """Boltzmann break of daughter cov bonds (weakly reversible).
        If a now-cov-orphaned monomer also has no anchor, return it to pool."""
        L = s.length
        bp = physics.break_prob(self.cfg.e_cov_daughter, E)
        for i in range(max(0, L - 1)):
            if not s.d_cov[i]:
                continue
            if np.random.random() < bp:
                s.d_cov[i] = False
                # Check left side (position i)
                if s.d_present[i] and not s.d_anchor[i]:
                    has_left = i > 0 and s.d_cov[i - 1] and s.d_present[i - 1]
                    if not has_left:
                        self._put(s.x, s.y, int(s.d_mono[i]))
                        s.d_present[i] = False
                        s.d_mono[i] = 0
                # Check right side (position i+1)
                j = i + 1
                if s.d_present[j] and not s.d_anchor[j]:
                    has_right = j < L - 1 and s.d_cov[j] and s.d_present[j + 1]
                    if not has_right:
                        self._put(s.x, s.y, int(s.d_mono[j]))
                        s.d_present[j] = False
                        s.d_mono[j] = 0

    def _strand_break_parent_cov(self, s: Strand, E):
        """Boltzmann break of parent cov bonds. Returns list of broken indices."""
        broken = []
        L = s.length
        bp = physics.break_prob(self.cfg.e_cov_parent, E)
        for i in range(max(0, L - 1)):
            if not s.cov[i]:
                continue
            if np.random.random() < bp:
                s.cov[i] = False
                broken.append(i)
        return broken

    def _strand_liftoff(self, s: Strand):
        """Detach maximal cov-connected daughter components that have no remaining anchors."""
        L = s.length
        if L == 0:
            return [], 0
        new_strands = []
        births = 0
        visited = np.zeros(L, dtype=bool)
        for i in range(L):
            if not s.d_present[i] or visited[i]:
                continue
            comp = [i]
            visited[i] = True
            # Walk left
            j = i
            while j > 0 and s.d_cov[j - 1] and s.d_present[j - 1] and not visited[j - 1]:
                j -= 1
                comp.append(j)
                visited[j] = True
            # Walk right
            j = i
            while j < L - 1 and s.d_cov[j] and s.d_present[j + 1] and not visited[j + 1]:
                j += 1
                comp.append(j)
                visited[j] = True
            comp.sort()
            # Anchored anywhere?
            if any(bool(s.d_anchor[k]) for k in comp):
                continue
            # Lift off
            child_mono = [int(s.d_mono[k]) for k in comp]
            for k in comp:
                s.d_present[k] = False
                s.d_mono[k] = 0
            for k in range(len(comp) - 1):
                if comp[k + 1] - comp[k] == 1:
                    s.d_cov[comp[k]] = False
            if len(child_mono) >= self.cfg.min_birth_length:
                # Stage 1: substitution mutations at birth.
                # For each position, with prob p_mut replace with a uniform-random
                # *non-original* monomer type (3 choices).
                n_mut = 0
                p_mut = self.cfg.p_mut
                if p_mut > 0:
                    for k in range(len(child_mono)):
                        if np.random.random() < p_mut:
                            orig = child_mono[k]
                            # pick a random other type (0..3 excluding orig)
                            new_t = np.random.randint(physics.N_TYPES - 1)
                            if new_t >= orig:
                                new_t += 1
                            child_mono[k] = int(new_t)
                            n_mut += 1
                child = Strand(
                    child_mono, s.x, s.y,
                    gen=s.gen + 1, parent_id=s.id,
                    mutation_count=n_mut,
                    total_mutations_inherited=s.total_mutations_inherited + n_mut,
                )
                child.cov[:] = True  # the lifted segment is fully cov-linked by construction
                new_strands.append(child)
                births += 1
            else:
                for m in child_mono:
                    self._put(s.x, s.y, m)
        return new_strands, births

    def _strand_fragment(self, s: Strand, broken_indices):
        """Split parent strand at broken cov positions. Returns list of new fragment strands.
        The original strand should be removed by the caller."""
        L = s.length
        cuts = sorted(set(broken_indices))
        # Build half-open segments [start, end)
        segments = []
        start = 0
        for c in cuts:
            end = c + 1  # mono[start..c] inclusive → half-open [start, c+1)
            if end > start:
                segments.append((start, end))
            start = c + 1
        if start < L:
            segments.append((start, L))

        fragments = []
        for (a, b) in segments:
            seg_len = b - a
            if seg_len < 1:
                continue
            if seg_len == 1:
                # Single monomer: return parent monomer + any daughter at this position
                self._put(s.x, s.y, int(s.mono[a]))
                if s.d_present[a]:
                    self._put(s.x, s.y, int(s.d_mono[a]))
                continue
            new = Strand(s.mono[a:b].tolist(), s.x, s.y, gen=s.gen, parent_id=s.parent_id)
            new.cov[:] = s.cov[a:b - 1]
            new.d_present[:] = s.d_present[a:b]
            new.d_mono[:] = s.d_mono[a:b]
            new.d_anchor[:] = s.d_anchor[a:b]
            new.d_cov[:] = s.d_cov[a:b - 1]
            fragments.append(new)
        return fragments

    # ─────────── History ───────────

    def _record(self, births, deaths):
        h = self.history
        n = len(self.strands)
        if n > 0:
            lens = [s.length for s in self.strands]
            avg_len = float(np.mean(lens))
            max_len = int(max(lens))
            n_anchored = sum(1 for s in self.strands if s.d_anchor.any())
        else:
            avg_len = 0.0
            max_len = 0
            n_anchored = 0
        h["step"].append(self.step_count)
        h["env_E"].append(self.env_E)
        h["n_strands"].append(n)
        h["total_free"].append(self.total_free())
        h["births"].append(births)
        h["deaths"].append(deaths)
        h["max_len"].append(max_len)
        h["avg_len"].append(avg_len)
        h["n_anchored"].append(n_anchored)
