"""Strand: a polymer of monomers, with optional H-bonded daughter chain."""
import numpy as np
from . import physics
from . import catalysis


_next_id = 0
def _new_id():
    global _next_id
    _next_id += 1
    return _next_id

def reset_ids():
    global _next_id
    _next_id = 0


class Strand:
    """A linear polymer at grid cell (x, y).

    Parent strand:
      mono       (L,)   int8  — parent monomer types (0..3)
      cov        (L-1,) bool  — parent backbone covalent bonds

    Daughter (H-bonded chain on the parent template):
      d_present  (L,)   bool  — daughter monomer present at this position
      d_mono     (L,)   int8  — daughter monomer type (valid where d_present)
      d_anchor   (L,)   bool  — currently H-bonded to parent (subset of d_present)
      d_cov      (L-1,) bool  — daughter backbone cov bond between d[i] and d[i+1]

    Invariant: a daughter monomer remains in d_present iff it is anchored OR
    cov-linked to a present neighbor. Once both supports vanish, it goes back
    to the local monomer pool.

    Lift-off: a maximal cov-linked daughter component with no remaining
    anchors detaches as a new free strand.
    """

    __slots__ = (
        "id", "mono", "cov",
        "d_present", "d_mono", "d_anchor", "d_cov",
        "x", "y", "gen", "parent_id", "age", "repl_count",
        "mutation_count", "total_mutations_inherited",
        # Stage 2: cached catalytic strengths (recomputed only when mono changes)
        "_cat_lipid", "_cat_poly", "_cat_extend",
    )

    def __init__(self, mono, x, y, gen=0, parent_id=None,
                 mutation_count=0, total_mutations_inherited=0):
        self.id = _new_id()
        self.mono = np.asarray(mono, dtype=np.int8)
        L = len(self.mono)
        self.cov = np.ones(max(0, L - 1), dtype=bool)
        self.d_present = np.zeros(L, dtype=bool)
        self.d_mono    = np.zeros(L, dtype=np.int8)
        self.d_anchor  = np.zeros(L, dtype=bool)
        self.d_cov     = np.zeros(max(0, L - 1), dtype=bool)
        self.x = int(x)
        self.y = int(y)
        self.gen = int(gen)
        self.parent_id = parent_id
        self.age = 0
        self.repl_count = 0
        # Mutations introduced at this strand's own birth (during liftoff)
        self.mutation_count = int(mutation_count)
        # Cumulative mutations along ancestral chain (this + all parents)
        self.total_mutations_inherited = int(total_mutations_inherited)
        # Catalytic strengths cached from sequence (Stage 2)
        self._cat_lipid = 0
        self._cat_poly = 0
        self._cat_extend = 0
        self.recompute_catalytic_cache()

    def recompute_catalytic_cache(self):
        """Recompute cached catalytic strengths from current self.mono.
        Call this any time self.mono is mutated in place (extend, fragment)."""
        m = self.mono
        self._cat_lipid  = catalysis.catalytic_strength(m, catalysis.DEFAULT_MOTIFS["lipid_synth"])
        self._cat_poly   = catalysis.catalytic_strength(m, catalysis.DEFAULT_MOTIFS["polymerize"])
        self._cat_extend = catalysis.catalytic_strength(m, catalysis.DEFAULT_MOTIFS["extend"])

    @property
    def length(self):
        return len(self.mono)

    def sequence_str(self):
        return "".join(physics.LABELS[int(m)] for m in self.mono)

    def __repr__(self):
        return f"<Strand #{self.id} L={self.length} g={self.gen} seq={self.sequence_str()}>"
