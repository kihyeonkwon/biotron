import { useEffect, useRef, useState } from 'react';
import { World } from './simulation/world.js';
import { renderWorld } from './renderer/canvas.js';

function ChainHistogram({ lenHisto }) {
  if (!lenHisto) return null;
  // Bin: 2, 3, 4, 5, 6-9, 10-19, 20+
  const bins = [
    { label: '2', test: (L) => L === 2 },
    { label: '3', test: (L) => L === 3 },
    { label: '4', test: (L) => L === 4 },
    { label: '5', test: (L) => L === 5 },
    { label: '6-9', test: (L) => L >= 6 && L <= 9 },
    { label: '10-19', test: (L) => L >= 10 && L <= 19 },
    { label: '20+', test: (L) => L >= 20 },
  ];
  const counts = bins.map((b) => {
    let n = 0;
    for (const L in lenHisto) {
      if (b.test(parseInt(L))) n += lenHisto[L];
    }
    return n;
  });
  const max = Math.max(1, ...counts);
  return (
    <div style={{ marginTop: 6 }}>
      <div style={{ fontSize: 8, color: '#444', marginBottom: 3 }}>length histogram</div>
      <svg width="100%" height="40" viewBox="0 0 280 40" preserveAspectRatio="none">
        {counts.map((c, i) => {
          const w = 280 / counts.length;
          const h = (c / max) * 32;
          return (
            <g key={i}>
              <rect
                x={i * w + 2}
                y={32 - h}
                width={w - 4}
                height={h}
                fill="#4ECDC4"
                opacity={c > 0 ? 0.8 : 0.2}
              />
              <text
                x={i * w + w / 2}
                y={39}
                fill="#666"
                fontSize="6"
                textAnchor="middle"
              >
                {bins[i].label}
              </text>
            </g>
          );
        })}
      </svg>
    </div>
  );
}

const GRID_W = 200;
const GRID_H = 200;
const CELL_PX = 3;          // canvas cell size
const TICK_INTERVAL_MS = 16; // ~60 fps

export default function App() {
  const canvasRef = useRef(null);
  const worldRef = useRef(null);
  const rafRef = useRef(null);
  const [running, setRunning] = useState(false);
  const [tick, setTick] = useState(0);
  const [stats, setStats] = useState({});

  // Initialize world once
  useEffect(() => {
    worldRef.current = new World({ width: GRID_W, height: GRID_H });
    setStats(worldRef.current.getStats());
    drawNow();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  // Tick loop
  useEffect(() => {
    if (!running) return;
    let last = performance.now();
    const loop = (now) => {
      if (now - last >= TICK_INTERVAL_MS) {
        worldRef.current.tick();
        last = now;
        if (worldRef.current.tickCount % 4 === 0) {
          drawNow();
          setTick(worldRef.current.tickCount);
          setStats(worldRef.current.getStats());
        }
      }
      rafRef.current = requestAnimationFrame(loop);
    };
    rafRef.current = requestAnimationFrame(loop);
    return () => cancelAnimationFrame(rafRef.current);
  }, [running]);

  function drawNow() {
    const canvas = canvasRef.current;
    if (!canvas || !worldRef.current) return;
    renderWorld(canvas, worldRef.current, CELL_PX);
  }

  function reset() {
    setRunning(false);
    worldRef.current = new World({ width: GRID_W, height: GRID_H });
    setTick(0);
    setStats(worldRef.current.getStats());
    drawNow();
  }

  function step() {
    worldRef.current.tick();
    setTick(worldRef.current.tickCount);
    setStats(worldRef.current.getStats());
    drawNow();
  }

  return (
    <div className="app">
      <div className="canvas-pane">
        <canvas
          ref={canvasRef}
          className="world"
          width={GRID_W * CELL_PX}
          height={GRID_H * CELL_PX}
        />
      </div>

      <aside className="sidebar">
        <h1>
          BIOTRON v2
          <span className="sub">abiogenesis simulator</span>
        </h1>

        <div className="section">
          <div className="section-title">control</div>
          <div className="controls">
            <button className="btn" onClick={reset}>↻</button>
            <button
              className={`btn ${running ? 'active' : ''}`}
              onClick={() => setRunning((r) => !r)}
            >
              {running ? '⏸' : '▶'}
            </button>
            <button className="btn" onClick={step} disabled={running}>
              step
            </button>
          </div>
          <div className="row">
            <span className="label">tick</span>
            <span className="value">{tick.toLocaleString()}</span>
          </div>
        </div>

        <div className="section">
          <div className="section-title">environment</div>
          <div className="row">
            <span className="label">temperature</span>
            <span className="value">{(stats.temperature ?? 0).toFixed(2)}</span>
          </div>
          <div className="bar">
            <div
              className="bar-fill"
              style={{ width: `${(((stats.temperature ?? 0) + 1) / 2) * 100}%` }}
            />
          </div>
          <div className="row">
            <span className="label">water level</span>
            <span className="value">{(stats.waterLevel ?? 0).toFixed(2)}</span>
          </div>
          <div className="bar">
            <div
              className="bar-fill"
              style={{
                width: `${(((stats.waterLevel ?? 0) + 1) / 2) * 100}%`,
                background: 'linear-gradient(90deg, #c98c5e, #5a8fb8, #2a4a8a)',
              }}
            />
          </div>
        </div>

        <div className="section">
          <div className="section-title">concentrations (mean)</div>
          {['A', 'U', 'G', 'C'].map((k) => (
            <div className="row" key={k}>
              <span className="label">{k}</span>
              <span className="value">
                {(stats.concentrations?.[k] ?? 0).toFixed(3)}
              </span>
            </div>
          ))}
        </div>

        <div className="section">
          <div className="section-title">structures</div>
          <div className="row">
            <span className="label">RNA chains</span>
            <span className="value">{stats.rnaChains ?? 0}</span>
          </div>
          <div className="row">
            <span className="label">max length</span>
            <span className="value">{stats.maxRnaLen ?? 0}</span>
          </div>
          <div className="row">
            <span className="label">H-bonded</span>
            <span className="value">{stats.hBondedChains ?? 0}</span>
          </div>
          <ChainHistogram lenHisto={stats.lenHisto} />
        </div>

        <div className="section">
          <div className="section-title">milestones</div>
          {(stats.milestoneStatus ?? []).map((m) => (
            <div className="row" key={m.id} title={m.description}>
              <span className="label">
                {m.tickReached != null ? '✓' : '·'} M{m.id} {m.name}
              </span>
              <span className="value">
                {m.tickReached != null ? `t${m.tickReached.toLocaleString()}` : '—'}
              </span>
            </div>
          ))}
        </div>

        <div className="section">
          <div className="section-title">phase</div>
          <div className="row">
            <span className="label">build status</span>
            <span className="value">phase 2 / step 12</span>
          </div>
          <div className="row">
            <span className="label">grid</span>
            <span className="value">{GRID_W}×{GRID_H}</span>
          </div>
        </div>
      </aside>
    </div>
  );
}
