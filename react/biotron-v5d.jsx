import { useState, useEffect, useRef, useCallback } from "react";

// ============================================================
// BIOTRON v5d: Fluctuating Environment
// "Birth sensitive, survival resilient"
// Parameters drift over time like geological cycles
// Life emerges when sweet spot is randomly hit
//
// Changes from v5:
// - Decay removed as separate rule → strand survives if
//   total bond energy > environment thermal energy
// - Parameters aggressively tuned for time compression
// - Millions of years → hundreds of steps
// ============================================================

const SITES = 2;

function randomMono() {
  return Array.from({ length: SITES }, () => Math.random() < 0.5 ? 1 : 0);
}
function compMono(m) { return m.map(s => 1 - s); }
function monoKey(m) { return m.join(""); }
function monoMatch(a, b) {
  for (let i = 0; i < a.length; i++) if (a[i] + b[i] !== 1) return false;
  return true;
}

const COLORS = { "00": "#FF6B6B", "11": "#4ECDC4", "01": "#FFEAA7", "10": "#DDA0DD" };
const LABELS = { "00": "α", "01": "β", "10": "β̄", "11": "ᾱ" };

let _id = 0;

function makeStrand(monomers) {
  return {
    id: ++_id,
    mono: monomers.map(m => m.slice()),
    hbonds: new Array(monomers.length).fill(null),
    age: 0,
    gen: 0,
    parentId: null,
    replCount: 0,
  };
}

// Bond energy of a strand: sum of covalent bonds + H-bonds
// Covalent bond = 1.0 energy unit per bond (strong)
// H-bond = 0.3 energy unit per bond (weak)
// Total determines survival against thermal energy
function strandBondEnergy(strand) {
  // Covalent bonds = (length - 1) backbone bonds
  const covalent = (strand.mono.length - 1) * 1.0;
  // H-bonds
  let hbondCount = 0;
  let consecutiveBonus = 0;
  let run = 0;
  for (let i = 0; i < strand.hbonds.length; i++) {
    if (strand.hbonds[i] !== null) {
      hbondCount++;
      run++;
      if (run >= 2) consecutiveBonus += 0.15; // stacking energy
    } else {
      run = 0;
    }
  }
  const hbondE = hbondCount * 0.3 + consecutiveBonus;
  return covalent + hbondE;
}

// Thermal energy that a strand must withstand
// Shorter strands are more vulnerable (surface area ratio)
function thermalChallenge(envEnergy, strandLen) {
  // Smaller molecules get hit harder by thermal fluctuations
  return envEnergy * (3.0 / (1 + strandLen * 0.3));
}

class World {
  constructor(cfg) {
    this.cfg = cfg;
    this.strands = [];
    this.freeMonos = {};
    this.step = 0;
    this.phase = 0;

    // Each fluctuating param has its OWN independent oscillators
    // So they move independently — not all up/down together
    this.fluctuatingParams = {
      hRate: { base: cfg.hRate, range: 0.4, oscs: [
        { freq: 0.0008, amp: 0.5, ph: Math.random() * 6.28 },
        { freq: 0.005, amp: 0.3, ph: Math.random() * 6.28 },
        { freq: 0.02, amp: 0.2, ph: Math.random() * 6.28 },
      ]},
      stackB: { base: cfg.stackB, range: 2.5, oscs: [
        { freq: 0.0006, amp: 0.4, ph: Math.random() * 6.28 },
        { freq: 0.004, amp: 0.35, ph: Math.random() * 6.28 },
        { freq: 0.015, amp: 0.25, ph: Math.random() * 6.28 },
      ]},
      detachR: { base: cfg.detachR, range: 0.03, oscs: [
        { freq: 0.001, amp: 0.45, ph: Math.random() * 6.28 },
        { freq: 0.007, amp: 0.3, ph: Math.random() * 6.28 },
        { freq: 0.025, amp: 0.25, ph: Math.random() * 6.28 },
      ]},
      spontRate: { base: cfg.spontRate, range: 0.06, oscs: [
        { freq: 0.0005, amp: 0.5, ph: Math.random() * 6.28 },
        { freq: 0.003, amp: 0.3, ph: Math.random() * 6.28 },
        { freq: 0.018, amp: 0.2, ph: Math.random() * 6.28 },
      ]},
      survivalFactor: { base: cfg.survivalFactor, range: 0.35, oscs: [
        { freq: 0.0009, amp: 0.4, ph: Math.random() * 6.28 },
        { freq: 0.006, amp: 0.35, ph: Math.random() * 6.28 },
        { freq: 0.022, amp: 0.25, ph: Math.random() * 6.28 },
      ]},
      sepR: { base: cfg.sepR, range: 0.6, oscs: [
        { freq: 0.0007, amp: 0.45, ph: Math.random() * 6.28 },
        { freq: 0.005, amp: 0.3, ph: Math.random() * 6.28 },
        { freq: 0.02, amp: 0.25, ph: Math.random() * 6.28 },
      ]},
    };

    this.popHist = [];
    this.monoHist = [];
    this.energyHist = [];
    this.birthHist = [];
    this.maxLenHist = [];
    this.avgLenHist = [];
    this.hbondedHist = [];
    this.deathHist = [];
    this.envFitnessHist = []; // how close to "sweet spot"
    this.events = [];
    this.totalBirths = 0;
    this.sustained = false;
    this.sustainedStep = -1;
    this.consRepl = 0;
    this.effectiveParams = {};

    this.init();
  }

  init() {
    ["00", "01", "10", "11"].forEach(k => {
      this.freeMonos[k] = Math.floor(this.cfg.monoPool / 4);
    });
    this.updateEffective();
    this.snap();
  }

  updateEffective() {
    this.effectiveParams = {};
    for (const [key, spec] of Object.entries(this.fluctuatingParams)) {
      // Each param has its own oscillator sum → independent movement
      let mod = 0;
      for (const osc of spec.oscs) {
        mod += Math.sin(this.step * osc.freq + osc.ph) * osc.amp;
      }
      mod = Math.max(-1, Math.min(1, mod));
      const val = spec.base + mod * spec.range;
      this.effectiveParams[key] = Math.max(0.001, val);
    }
  }

  // Environment "fitness" for life: how favorable are current conditions
  getEnvFitness() {
    const ep = this.effectiveParams;
    // Sweet spot: high hRate, high stackB, low detachR, moderate survival
    const score = (ep.hRate / 0.8) * (ep.stackB / 5) * (1 - ep.detachR / 0.05) * (ep.spontRate / 0.1);
    return Math.max(0, Math.min(1, score));
  }

  energy() { return (Math.sin(this.phase) + 1) / 2; }
  totalFree() { return Object.values(this.freeMonos).reduce((a, b) => a + b, 0); }

  takeMono(type) {
    const k = monoKey(type);
    if ((this.freeMonos[k] || 0) > 0) { this.freeMonos[k]--; return true; }
    return false;
  }
  returnMono(type) {
    this.freeMonos[monoKey(type)] = (this.freeMonos[monoKey(type)] || 0) + 1;
  }

  runStep() {
    this.step++;
    this.phase += this.cfg.cycleSpeed;
    this.updateEffective();
    const E = this.energy();
    const ep = this.effectiveParams;
    let births = 0;
    let deaths = 0;

    // Monomer inflow
    ["00", "01", "10", "11"].forEach(k => {
      this.freeMonos[k] = Math.min(
        (this.freeMonos[k] || 0) + this.cfg.monoInflow,
        this.cfg.maxMono
      );
    });

    // === Spontaneous polymerization: free monomers → dimers ===
    if (this.strands.length < this.cfg.maxStrands) {
      for (let i = 0; i < this.cfg.spontAttempts; i++) {
        if (Math.random() < ep.spontRate) {
          const m1 = randomMono();
          const m2 = randomMono();
          if (this.takeMono(m1) && this.takeMono(m2)) {
            this.strands.push(makeStrand([m1, m2]));
          }
        }
      }
    }

    // === End extension: strands grow by capturing monomers ===
    for (const s of this.strands) {
      if (s.mono.length >= this.cfg.maxLen) continue;
      if (Math.random() < this.cfg.extendRate && E < 0.5) {
        const m = randomMono();
        if (this.takeMono(m)) {
          if (Math.random() < 0.5) { s.mono.push(m); s.hbonds.push(null); }
          else { s.mono.unshift(m); s.hbonds.unshift(null); }
        }
      }
    }

    // === COLD PHASE: H-bonding ===
    if (E < this.cfg.annealT) {
      for (const s of this.strands) {
        for (let i = 0; i < s.mono.length; i++) {
          if (s.hbonds[i] !== null) continue;

          const left = i > 0 && s.hbonds[i - 1] !== null;
          const right = i < s.mono.length - 1 && s.hbonds[i + 1] !== null;
          const stack = (left ? ep.stackB : 1) * (right ? ep.stackB : 1);

          const needed = compMono(s.mono[i]);
          const conc = (this.freeMonos[monoKey(needed)] || 0) / this.cfg.maxMono;
          const prob = ep.hRate * conc * stack * (1 - E);

          if (Math.random() < prob && this.takeMono(needed)) {
            s.hbonds[i] = needed;
          }
        }
      }

      // Spontaneous H-bond detach
      for (const s of this.strands) {
        for (let i = 0; i < s.hbonds.length; i++) {
          if (s.hbonds[i] === null) continue;
          const left = i > 0 && s.hbonds[i - 1] !== null;
          const right = i < s.hbonds.length - 1 && s.hbonds[i + 1] !== null;
          const hold = (left ? 0.4 : 0) + (right ? 0.4 : 0);
          if (Math.random() < ep.detachR * (1 - hold)) {
            this.returnMono(s.hbonds[i]);
            s.hbonds[i] = null;
          }
        }
      }
    }

    // === HOT PHASE: Strand separation from template ===
    if (E > this.cfg.denatT) {
      const newStrands = [];
      for (const s of this.strands) {
        const runs = [];
        let rs = -1;
        for (let i = 0; i <= s.mono.length; i++) {
          const b = i < s.mono.length && s.hbonds[i] !== null;
          if (b && rs === -1) rs = i;
          if (!b && rs !== -1) { runs.push([rs, i]); rs = -1; }
        }

        for (const [start, end] of runs) {
          const len = end - start;
          const strength = len * this.cfg.bondStr;
          const sep = (E - this.cfg.denatT) * ep.sepR / (1 + strength);

          if (len >= 2 && Math.random() < sep) {
            const childMono = [];
            for (let i = start; i < end; i++) {
              childMono.push(s.hbonds[i]);
              s.hbonds[i] = null;
            }
            if (this.strands.length + newStrands.length < this.cfg.maxStrands) {
              const child = makeStrand(childMono);
              child.gen = s.gen + 1;
              child.parentId = s.id;
              newStrands.push(child);
              s.replCount++;
              births++;
              this.totalBirths++;

              if (this.totalBirths <= 30) {
                this.events.push({
                  step: this.step,
                  msg: `🧬 Birth: #${s.id}(len${s.mono.length},gen${s.gen}) → #${child.id}(len${len},gen${child.gen})`,
                });
              }
            }
          } else if (E > 0.85) {
            for (let i = start; i < end; i++) {
              if (s.hbonds[i]) { this.returnMono(s.hbonds[i]); s.hbonds[i] = null; }
            }
          }
        }
      }
      this.strands.push(...newStrands);
    }

    // === BOLTZMANN SURVIVAL (replaces decay) ===
    // P(break) = exp(-bondEnergy / thermalEnergy)
    // Like real thermodynamics: even stable molecules can break
    // by thermal fluctuation, but it's exponentially unlikely
    const before = this.strands.length;
    this.strands = this.strands.filter(s => {
      const bondE = strandBondEnergy(s);
      const thermal = thermalChallenge(E, s.mono.length);

      // Boltzmann: probability of breaking ∝ exp(-ΔE/kT)
      // ΔE = bondE, kT = thermal * survivalFactor
      const kT = Math.max(0.01, thermal * ep.survivalFactor);
      const breakProb = Math.exp(-bondE / kT) * 0.3; // 0.3 = attempt frequency per step

      if (Math.random() < breakProb) {
        for (const m of s.mono) this.returnMono(m);
        for (const m of s.hbonds) if (m) this.returnMono(m);
        return false;
      }
      s.age++;
      return true;
    });
    deaths = before - this.strands.length;

    // Detection
    this.birthHist.push(births);
    this.deathHist.push(deaths);
    if (births > 0) this.consRepl++;
    else this.consRepl = Math.max(0, this.consRepl - 1);
    if (this.consRepl >= 6 && !this.sustained) {
      this.sustained = true;
      this.sustainedStep = this.step;
      this.events.push({ step: this.step, msg: "🧬🧬🧬 SUSTAINED EMERGENT REPLICATION!" });
    }
    this.snap();
  }

  snap() {
    this.popHist.push(this.strands.length);
    this.monoHist.push(this.totalFree());
    this.energyHist.push(this.energy());
    this.envFitnessHist.push(this.getEnvFitness());
    const lens = this.strands.map(s => s.mono.length);
    this.avgLenHist.push(lens.length > 0 ? lens.reduce((a, b) => a + b, 0) / lens.length : 0);
    this.maxLenHist.push(lens.length > 0 ? Math.max(...lens) : 0);
    this.hbondedHist.push(this.strands.filter(s => s.hbonds.some(h => h !== null)).length);
  }

  getReplicators() {
    return [...this.strands].filter(s => s.replCount > 0).sort((a, b) => b.replCount - a.replCount).slice(0, 8);
  }
  getHbonded() {
    return this.strands.filter(s => s.hbonds.some(h => h !== null)).slice(0, 12);
  }
  getLongest() {
    return [...this.strands].sort((a, b) => b.mono.length - a.mono.length).slice(0, 6);
  }
}

// === VIS ===
function Spark({ data, color, w = 95, h = 20, label }) {
  if (!data || data.length < 2) return null;
  const d = data.slice(-250);
  const max = Math.max(...d, 0.01);
  const pts = d.map((v, i) => `${(i / (d.length - 1)) * w},${h - (v / max) * (h - 3) - 1.5}`).join(" ");
  const last = d[d.length - 1];
  return (
    <div>
      {label && <div style={{ fontSize: 7, color: "#333" }}>{label} <span style={{ color: "#555" }}>{typeof last === 'number' ? (last < 10 ? last.toFixed(1) : Math.round(last)) : ''}</span></div>}
      <svg width={w} height={h} style={{ display: "block" }}>
        <polyline points={pts} fill="none" stroke={color} strokeWidth="1.5" strokeLinecap="round" />
      </svg>
    </div>
  );
}

function SV({ strand }) {
  const len = strand.mono.length;
  const cw = Math.min(13, 120 / len);
  const hasH = strand.hbonds.some(h => h !== null);
  const h = hasH ? 19 : 8;
  return (
    <svg width={cw * len} height={h} style={{ display: "block" }}>
      {strand.mono.map((m, i) => (
        <g key={i}>
          <rect x={i * cw} y={0} width={cw - 0.8} height={7} rx={1.5} fill={COLORS[monoKey(m)]} opacity={0.85} />
          {strand.hbonds[i] !== null && (
            <>
              <line x1={i * cw + cw / 2} y1={7} x2={i * cw + cw / 2} y2={11}
                stroke={monoMatch(m, strand.hbonds[i]) ? "#4ECDC466" : "#FF6B6B66"} strokeWidth={0.6} />
              <rect x={i * cw} y={11} width={cw - 0.8} height={7} rx={1.5}
                fill={COLORS[monoKey(strand.hbonds[i])]} opacity={0.85}
                stroke={monoMatch(m, strand.hbonds[i]) ? "#4ECDC4" : "#FF6B6B"} strokeWidth={0.4} />
            </>
          )}
          {strand.hbonds[i] === null && hasH && (
            <rect x={i * cw} y={11} width={cw - 0.8} height={7} rx={1.5}
              fill="none" stroke="#151520" strokeWidth={0.3} strokeDasharray="1.5,1.5" />
          )}
        </g>
      ))}
    </svg>
  );
}

function EBar({ v, w = 200 }) {
  const c = v < 0.35 ? "#4ECDC4" : v > 0.65 ? "#FF6B6B" : "#FFEAA7";
  const label = v < 0.35 ? "COLD·bind" : v > 0.65 ? "HOT·split" : "transition";
  return (
    <div style={{ display: "flex", alignItems: "center", gap: 5 }}>
      <div style={{ position: "relative", width: w, height: 8, background: "#111", borderRadius: 4 }}>
        <div style={{ position: "absolute", left: 0, top: 0, width: `${v * 100}%`, height: "100%",
          background: "linear-gradient(90deg, #4ECDC4, #FFEAA7, #FF6B6B)", borderRadius: 4, opacity: 0.6 }} />
        <div style={{ position: "absolute", left: `${v * 100}%`, top: -2, width: 7, height: 12,
          background: c, borderRadius: 3, transform: "translateX(-3px)", border: "1px solid #000" }} />
      </div>
      <span style={{ fontSize: 8, color: c, fontWeight: 600 }}>{label}</span>
    </div>
  );
}

// TIME-COMPRESSED DEFAULTS
const DEF = {
  monoPool: 800,
  maxMono: 400,
  monoInflow: 5,
  maxStrands: 800,
  maxLen: 16,
  cycleSpeed: 0.025,   // slower cycle → more time per phase
  annealT: 0.45,
  denatT: 0.55,
  hRate: 0.5,           // 2x faster H-bonding
  stackB: 4,            // strong stacking
  detachR: 0.015,       // very low detach
  bondStr: 0.12,
  sepR: 1.2,            // aggressive separation
  spontRate: 0.08,      // 4x faster spontaneous poly
  spontAttempts: 8,
  extendRate: 0.03,
  survivalFactor: 0.7,  // how much bond energy helps survival
};

export default function BiotronV5d() {
  const [cfg, setCfg] = useState(DEF);
  const [world, setWorld] = useState(null);
  const [running, setRunning] = useState(false);
  const [speed, setSpeed] = useState(40);
  const ref = useRef(null);
  const iRef = useRef(null);

  const init = useCallback(() => {
    _id = 0;
    const w = new World(cfg);
    ref.current = w;
    setWorld({ ...w });
    setRunning(false);
    if (iRef.current) clearInterval(iRef.current);
  }, [cfg]);

  useEffect(() => { init(); }, [init]);
  useEffect(() => {
    if (running && ref.current) {
      iRef.current = setInterval(() => {
        for (let i = 0; i < 5; i++) ref.current.runStep();
        setWorld({ ...ref.current, strands: [...ref.current.strands] });
      }, speed);
    }
    return () => { if (iRef.current) clearInterval(iRef.current); };
  }, [running, speed]);

  const E = world ? ref.current?.energy() || 0 : 0;
  const hbonded = world ? ref.current?.getHbonded() : [];
  const replicators = world ? ref.current?.getReplicators() : [];
  const longest = world ? ref.current?.getLongest() : [];

  return (
    <div style={{ minHeight: "100vh", background: "#030306", color: "#bbb",
      fontFamily: "'JetBrains Mono', 'SF Mono', monospace", padding: 14, boxSizing: "border-box" }}>

      <div style={{ marginBottom: 8 }}>
        <h1 style={{ fontSize: 18, fontWeight: 800, color: "#4ECDC4", margin: 0, letterSpacing: "0.1em" }}>BIOTRON v5d</h1>
        <p style={{ fontSize: 7, color: "#1a1a1a", margin: "2px 0" }}>
          Fluctuating environment · Birth sensitive, survival resilient · Boltzmann thermodynamics
        </p>
      </div>

      {world?.sustained && (
        <div style={{ background: "#4ECDC408", border: "2px solid #4ECDC4", borderRadius: 8, padding: "10px 14px", marginBottom: 8 }}>
          <div style={{ fontSize: 14, fontWeight: 700, color: "#4ECDC4" }}>🧬 EMERGENT REPLICATION — Step {world.sustainedStep}</div>
        </div>
      )}

      <EBar v={E} />

      {/* Environment fitness indicator */}
      {world?.effectiveParams && (
        <div style={{ display: "flex", alignItems: "center", gap: 8, margin: "4px 0 8px", fontSize: 8 }}>
          <span style={{ color: "#333" }}>Env fitness:</span>
          <div style={{ width: 100, height: 6, background: "#111", borderRadius: 3, overflow: "hidden" }}>
            <div style={{ width: `${(ref.current?.getEnvFitness() || 0) * 100}%`, height: "100%",
              background: `hsl(${(ref.current?.getEnvFitness() || 0) * 120}, 70%, 45%)`, borderRadius: 3 }} />
          </div>
          {Object.entries(world.effectiveParams || {}).map(([k, v]) => (
            <span key={k} style={{ color: "#333" }}>{k}:<span style={{ color: "#666" }}>{v.toFixed(2)}</span></span>
          ))}
        </div>
      )}

      <div style={{ display: "grid", gridTemplateColumns: "repeat(9, 1fr)", gap: 4, margin: "4px 0" }}>
        <Spark data={world?.popHist} color="#4ECDC4" label="Strands" />
        <Spark data={world?.monoHist} color="#F0B27A" label="Monomers" />
        <Spark data={world?.energyHist} color="#FF6B6B" label="Energy" />
        <Spark data={world?.envFitnessHist} color="#82E0AA" label="Env Fitness" />
        <Spark data={world?.birthHist} color="#FFEAA7" label="Births" />
        <Spark data={world?.deathHist} color="#FF6B6B" label="Deaths" />
        <Spark data={world?.avgLenHist} color="#DDA0DD" label="Avg len" />
        <Spark data={world?.maxLenHist} color="#85C1E9" label="Max len" />
        <Spark data={world?.hbondedHist} color="#96CEB4" label="H-bonded" />
      </div>

      <div style={{ display: "grid", gridTemplateColumns: "repeat(auto-fit, minmax(85px, 1fr))", gap: 3, marginBottom: 8 }}>
        {[
          { l: "Cycle Spd", k: "cycleSpeed", min: 0.005, max: 0.1, s: 0.005 },
          { l: "H-bond Rate", k: "hRate", min: 0.1, max: 0.8, s: 0.05 },
          { l: "Stack Bonus", k: "stackB", min: 1, max: 6, s: 0.5 },
          { l: "Detach", k: "detachR", min: 0.005, max: 0.08, s: 0.005 },
          { l: "Sep Rate", k: "sepR", min: 0.3, max: 2, s: 0.1 },
          { l: "Spont Poly", k: "spontRate", min: 0.01, max: 0.15, s: 0.01 },
          { l: "Extend", k: "extendRate", min: 0, max: 0.08, s: 0.005 },
          { l: "Survival", k: "survivalFactor", min: 0.3, max: 1.2, s: 0.05 },
          { l: "Mono Inflow", k: "monoInflow", min: 1, max: 12, s: 1 },
          { l: "Bond Str", k: "bondStr", min: 0.05, max: 0.3, s: 0.01 },
        ].map(({ l, k, min, max, s }) => (
          <div key={k} style={{ background: "#080810", borderRadius: 3, padding: "2px 5px", border: "1px solid #0c0c14" }}>
            <div style={{ fontSize: 6, color: "#222" }}>{l}</div>
            <div style={{ fontSize: 10, fontWeight: 700, color: "#4ECDC4" }}>{cfg[k] < 1 ? cfg[k].toFixed(3) : cfg[k]}</div>
            <input type="range" min={min} max={max} step={s} value={cfg[k]}
              onChange={e => setCfg(p => ({ ...p, [k]: parseFloat(e.target.value) }))} style={{ width: "100%", accentColor: "#4ECDC4" }} />
          </div>
        ))}
      </div>

      <div style={{ display: "flex", gap: 5, marginBottom: 8, alignItems: "center" }}>
        <button onClick={init} style={{ background: "#080810", color: "#4ECDC4", border: "1px solid #4ECDC4", borderRadius: 3, padding: "3px 8px", cursor: "pointer", fontSize: 9, fontFamily: "inherit" }}>↻</button>
        <button onClick={() => setRunning(!running)} style={{ background: running ? "#FF6B6B10" : "#4ECDC410", color: running ? "#FF6B6B" : "#4ECDC4", border: `1px solid ${running ? "#FF6B6B" : "#4ECDC4"}`, borderRadius: 3, padding: "3px 8px", cursor: "pointer", fontSize: 9, fontFamily: "inherit" }}>{running ? "⏸" : "▶"}</button>
        <span style={{ fontSize: 7, color: "#1a1a1a" }}>Speed</span>
        <input type="range" min={10} max={200} value={200 - speed} onChange={e => setSpeed(200 - parseInt(e.target.value))} style={{ width: 35, accentColor: "#4ECDC4" }} />
        <span style={{ fontSize: 8, color: "#1a1a1a", marginLeft: "auto" }}>
          Step <b style={{ color: "#888" }}>{world?.step || 0}</b> · Births <b style={{ color: "#FFEAA7" }}>{world?.totalBirths || 0}</b>
        </span>
      </div>

      {/* H-bonded = template assembly happening */}
      {hbonded.length > 0 && (
        <div style={{ marginBottom: 8 }}>
          <div style={{ fontSize: 8, fontWeight: 600, color: "#96CEB4", marginBottom: 3 }}>🔗 Template Assembly ({hbonded.length})</div>
          <div style={{ display: "flex", gap: 4, flexWrap: "wrap" }}>
            {hbonded.map(s => {
              const f = s.hbonds.filter(h => h !== null).length;
              return (
                <div key={s.id} style={{ background: "#080810", border: `1px solid ${f === s.mono.length ? "#4ECDC4" : "#0e0e16"}`, borderRadius: 4, padding: "3px 5px" }}>
                  <div style={{ fontSize: 6, color: "#444" }}>#{s.id} len:{s.mono.length} gen:{s.gen} {f}/{s.mono.length}{f === s.mono.length ? " ✓" : ""}</div>
                  <SV strand={s} />
                </div>
              );
            })}
          </div>
        </div>
      )}

      {/* Replicators */}
      {replicators.length > 0 && (
        <div style={{ marginBottom: 8 }}>
          <div style={{ fontSize: 8, fontWeight: 600, color: "#FF6B6B", marginBottom: 3 }}>🏆 Replicators</div>
          <div style={{ display: "flex", gap: 4, flexWrap: "wrap" }}>
            {replicators.map(s => (
              <div key={s.id} style={{ background: "#080810", border: "1px solid #0e0e16", borderRadius: 4, padding: "3px 5px" }}>
                <div style={{ fontSize: 6, color: "#555" }}>#{s.id} <b style={{ color: "#FF6B6B" }}>×{s.replCount}</b> len:{s.mono.length} gen:{s.gen}</div>
                <SV strand={s} />
              </div>
            ))}
          </div>
        </div>
      )}

      {/* Longest */}
      {longest.length > 0 && longest[0].mono.length > 2 && (
        <div style={{ marginBottom: 8 }}>
          <div style={{ fontSize: 8, fontWeight: 600, color: "#85C1E9", marginBottom: 3 }}>📏 Longest</div>
          <div style={{ display: "flex", gap: 4, flexWrap: "wrap" }}>
            {longest.filter(s => s.mono.length > 2).map(s => (
              <div key={s.id} style={{ background: "#080810", border: "1px solid #0e0e16", borderRadius: 4, padding: "3px 5px" }}>
                <div style={{ fontSize: 6, color: "#555" }}>#{s.id} len:<b style={{ color: "#85C1E9" }}>{s.mono.length}</b> gen:{s.gen}</div>
                <SV strand={s} />
              </div>
            ))}
          </div>
        </div>
      )}

      {world?.events?.length > 0 && (
        <div style={{ background: "#080810", border: "1px solid #0c0c14", borderRadius: 3, padding: "4px 7px", marginBottom: 8, maxHeight: 45, overflowY: "auto" }}>
          {world.events.slice(-6).reverse().map((e, i) => (
            <div key={i} style={{ fontSize: 7, color: "#333" }}><span style={{ color: "#4ECDC4" }}>[{e.step}]</span> {e.msg}</div>
          ))}
        </div>
      )}

      <div style={{ background: "#080810", border: "1px solid #0c0c14", borderRadius: 3, padding: 7, fontSize: 7, color: "#1a1a1a", lineHeight: 1.6 }}>
        <b style={{ color: "#FFEAA7" }}>v5d: Fluctuating Environment</b><br />
        • 환경 파라미터가 4개 주기(기후/화산/조석/기상)로 fluctuate<br />
        • Env Fitness 높을 때 = sweet spot → birth 가능. 낮을 때 = 혹독한 환경<br />
        • Birth sensitive: 특정 조건에서만 탄생. Survival resilient: 한번 태어나면 버팀<br />
        • Boltzmann: P(break) = exp(-bondE/kT) — 물리적 열역학
      </div>

      <style>{`
        input[type="range"]{height:2px;-webkit-appearance:none;appearance:none;background:#151515;border-radius:1px;outline:none}
        input[type="range"]::-webkit-slider-thumb{-webkit-appearance:none;width:8px;height:8px;border-radius:50%;background:#4ECDC4;cursor:pointer}
        ::-webkit-scrollbar{width:2px}::-webkit-scrollbar-thumb{background:#111;border-radius:1px}
      `}</style>
    </div>
  );
}
