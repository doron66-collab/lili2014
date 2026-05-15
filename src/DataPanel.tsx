import { useState, useEffect, useRef } from 'react';

// ─── Trials data ───────────────────────────────────────────────────────────────
const TRIALS = [
  {
    mutation: 'TP53', variant: 'Y220C', mColor: '#3399ff',
    drugs: [
      { drug: 'Rezatapopt (PC14586)', sponsor: 'PMV Pharma', progress: 50, phase: 'Phase II', trial: 'PYNNACLE · NCT04585750', note: 'Y220C cryptic pocket binder · NDA target Q1 2027' },
      { drug: 'Eprenetapopt (APR-246)', sponsor: 'Aprea Therapeutics', progress: 50, phase: 'Phase II', trial: 'NCT02999893', note: 'Pan-mutation covalent p53 reactivator' },
      { drug: 'Siremadlin (HDM201)', sponsor: 'Novartis', progress: 37, phase: 'Phase I/II', trial: 'NCT02143635', note: 'MDM2 inhibitor · TP53-WT tumors only' },
    ],
  },
  {
    mutation: 'KEAP1', variant: 'LOF', mColor: '#33ffaa',
    drugs: [
      { drug: 'VVD-065', sponsor: 'Vividion Therapeutics', progress: 25, phase: 'Phase I', trial: 'NCT05954312', note: 'First-in-class covalent NRF2 inhibitor' },
      { drug: 'Sapanisertib (TAK-228)', sponsor: 'Takeda', progress: 55, phase: 'Phase II', trial: 'Completed 2023', note: 'ORR 25% in NFE2L2-mutant LUSC' },
    ],
  },
  {
    mutation: 'CDKN2A', variant: 'p16 loss', mColor: '#ffaa33',
    drugs: [
      { drug: 'Palbociclib + Olaparib', sponsor: 'Pfizer / Academic', progress: 37, phase: 'Phase I/II', trial: 'Multiple NCTs', note: 'CDK4/6 + PARP1 degradation synergy (2024)' },
      { drug: 'Abemaciclib + Paclitaxel', sponsor: 'Eli Lilly', progress: 37, phase: 'Phase IB/II', trial: 'NCT03310879', note: 'DCR 66.7% · ORR 7.4%' },
    ],
  },
  {
    mutation: 'STK11', variant: 'LKB1 loss', mColor: '#ff3366',
    drugs: [
      { drug: 'Ceralasertib + Durvalumab', sponsor: 'AstraZeneca', progress: 75, phase: 'Phase III', trial: 'LATIFY · NCT05450692', note: 'ATR+PD-L1 · STK11/KEAP1-enriched biomarker' },
      { drug: 'Bemcentinib (BGB324)', sponsor: 'Blueprint Medicines', progress: 50, phase: 'Phase II', trial: 'AXL-STK11 cohort', note: 'AXL inhibitor · FDA Fast Track designation' },
      { drug: 'Telaglenastat (CB-839)', sponsor: 'Calithera Biosciences', progress: 50, phase: 'Phase II', trial: 'BeGIN · NCT03872427', note: 'Glutaminase inhibitor · LKB1/KEAP1/NRF2 cohort' },
    ],
  },
];

// ─── Co-mutation network ───────────────────────────────────────────────────────
const NET_NODES = [
  { id: 'TP53',   x: 170, y: 42,  r: 22, color: '#3399ff', freq: '~45%' },
  { id: 'KEAP1',  x: 308, y: 128, r: 18, color: '#33ffaa', freq: '~12%' },
  { id: 'CDKN2A', x: 255, y: 248, r: 20, color: '#ffaa33', freq: '~35%' },
  { id: 'STK11',  x: 82,  y: 248, r: 20, color: '#ff3366', freq: '~22%' },
  { id: 'KRAS',   x: 32,  y: 128, r: 16, color: '#aaaaaa', freq: '~28%' },
];
const NET_EDGES = [
  { a: 'STK11',   b: 'KRAS',   pct: 54, lx: 48,  ly: 185 },
  { a: 'STK11',   b: 'TP53',   pct: 44, lx: 115, ly: 138 },
  { a: 'STK11',   b: 'CDKN2A', pct: 37, lx: 168, ly: 255 },
  { a: 'STK11',   b: 'KEAP1',  pct: 27, lx: 200, ly: 198 },
  { a: 'KEAP1',   b: 'KRAS',   pct: 35, lx: 168, ly: 118 },
  { a: 'TP53',    b: 'KRAS',   pct: 25, lx: 94,  ly: 75  },
  { a: 'TP53',    b: 'KEAP1',  pct: 20, lx: 244, ly: 75  },
  { a: 'CDKN2A',  b: 'KRAS',   pct: 20, lx: 138, ly: 195 },
  { a: 'CDKN2A',  b: 'TP53',   pct: 15, lx: 218, ly: 145 },
];

// ─── Global trial sites ────────────────────────────────────────────────────────
const SITES = [
  { name: 'AstraZeneca',         loc: 'Cambridge, UK',       lon:   0.1, lat: 52.2, color: '#3399ff', targets: 'STK11 · KEAP1', drug: 'Ceralasertib + Durvalumab' },
  { name: 'IBM / Cleveland Clinic', loc: 'Cleveland, USA',   lon: -81.7, lat: 41.5, color: '#cc44ff', targets: 'sqDRIFT · QC',  drug: 'Protein Simulation' },
  { name: 'PMV Pharma',          loc: 'Cranbury, NJ',        lon: -74.5, lat: 40.3, color: '#3399ff', targets: 'TP53',           drug: 'Rezatapopt' },
  { name: 'Vividion Therapeutics', loc: 'San Diego, USA',    lon:-117.2, lat: 32.7, color: '#33ffaa', targets: 'KEAP1',          drug: 'VVD-065' },
  { name: 'Daiichi Sankyo',      loc: 'Tokyo, Japan',        lon: 139.7, lat: 35.7, color: '#3399ff', targets: 'TP53',           drug: 'Milademetan' },
  { name: 'Novartis',            loc: 'Basel, Switzerland',  lon:   7.6, lat: 47.6, color: '#3399ff', targets: 'TP53',           drug: 'Siremadlin' },
  { name: 'Calithera Biosciences', loc: 'San Francisco, USA',lon:-122.4, lat: 37.8, color: '#ff3366', targets: 'STK11 · KEAP1', drug: 'Telaglenastat' },
  { name: 'Google Quantum AI',   loc: 'Mountain View, USA',  lon:-122.1, lat: 37.4, color: '#cc44ff', targets: 'Quantum · KRAS', drug: 'Willow Chip' },
  { name: 'QuEra Computing',     loc: 'Boston, USA',         lon: -71.1, lat: 42.4, color: '#cc44ff', targets: 'Quantum · QC',  drug: 'Neutral Atom QC' },
  { name: 'Pfizer',              loc: 'New York, USA',       lon: -74.0, lat: 40.7, color: '#ffaa33', targets: 'CDKN2A',        drug: 'Palbociclib + Olaparib' },
  { name: 'Blueprint Medicines', loc: 'Cambridge, MA',       lon: -71.0, lat: 42.3, color: '#ff3366', targets: 'STK11',          drug: 'Bemcentinib' },
  { name: 'Wellcome Leap / Q4Bio', loc: 'Washington DC',     lon: -77.0, lat: 38.9, color: '#cc44ff', targets: 'Quantum · Bio', drug: 'Q4Bio Challenge' },
  { name: 'Eli Lilly',           loc: 'Indianapolis, USA',   lon: -86.1, lat: 39.8, color: '#ffaa33', targets: 'CDKN2A',        drug: 'Abemaciclib' },
];

// Continent lon/lat polygons for equirectangular SVG (viewBox 480×240)
const CONTINENTS: [number, number][][] = [
  [[-168,71],[-52,71],[-52,50],[-65,44],[-66,25],[-90,16],[-105,16],[-118,22],[-125,38],[-168,60]], // North America
  [[-82,12],[-60,15],[-35,5],[-35,-6],[-48,-28],[-52,-34],[-65,-55],[-75,-50]],                     // South America
  [[-10,36],[36,36],[36,50],[28,70],[10,72],[-10,60]],                                               // Europe
  [[-18,38],[52,38],[52,12],[42,-12],[32,-28],[18,-35],[-18,-35]],                                   // Africa
  [[26,10],[145,10],[145,75],[60,75],[26,70]],                                                       // Asia
  [[114,-44],[154,-44],[154,-10],[130,-10],[114,-22]],                                               // Australia
  [[-72,76],[-18,76],[-18,83],[-72,83]],                                                            // Greenland
];
const MAP_W = 480, MAP_H = 240;
const sx = (lon: number) => (lon + 180) * (MAP_W / 360);
const sy = (lat: number) => (90 - lat) * (MAP_H / 180);

// ─── Sub-components ────────────────────────────────────────────────────────────
const PHASES = ['PRE', 'I', 'II', 'III', 'FDA'];

function ProgressBar({ progress, color }: { progress: number; color: string }) {
  const [w, setW] = useState(0);
  useEffect(() => { const t = setTimeout(() => setW(progress), 120); return () => clearTimeout(t); }, [progress]);
  return (
    <div style={{ margin: '6px 0 10px' }}>
      <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: 5 }}>
        {PHASES.map((p, i) => (
          <span key={p} style={{ fontSize: 7.5, letterSpacing: 1, fontFamily: 'Courier New', color: i * 25 <= progress ? color : 'rgba(160,185,230,.65)' }}>{p}</span>
        ))}
      </div>
      <div style={{ position: 'relative', height: 6, background: 'rgba(40,65,120,.35)', borderRadius: 3 }}>
        <div style={{ position: 'absolute', inset: 0, height: '100%', width: `${w}%`, background: `linear-gradient(90deg,${color}66,${color})`, borderRadius: 3, transition: 'width 1.3s cubic-bezier(.4,0,.2,1)', boxShadow: `0 0 8px ${color}77` }} />
        {PHASES.map((_, i) => (
          <div key={i} style={{ position: 'absolute', top: -1, left: `${i * 25}%`, width: 8, height: 8, borderRadius: '50%', background: i * 25 <= w ? color : 'rgba(40,65,120,.6)', border: `1.5px solid ${i * 25 <= w ? color : 'rgba(60,90,150,.4)'}`, transform: 'translateX(-50%)', boxShadow: i * 25 <= w ? `0 0 6px ${color}` : 'none', transition: `background 0.3s ${i * 0.15}s, box-shadow 0.3s ${i * 0.15}s` }} />
        ))}
      </div>
    </div>
  );
}

function TrialsTab() {
  return (
    <div>
      {TRIALS.map(({ mutation, variant, mColor, drugs }) => (
        <div key={mutation} style={{ marginBottom: 18, padding: '12px 14px', background: 'rgba(8,16,45,.55)', borderRadius: 8, border: `1px solid ${mColor}22` }}>
          <div style={{ display: 'flex', alignItems: 'baseline', gap: 8, marginBottom: 10 }}>
            <span style={{ color: mColor, fontWeight: 'bold', fontSize: 15, letterSpacing: 2 }}>{mutation}</span>
            <span style={{ color: 'rgba(200,220,255,.85)', fontSize: 9, letterSpacing: 1 }}>· {variant}</span>
          </div>
          {drugs.map(({ drug, sponsor, progress, phase, trial, note }) => (
            <div key={drug} style={{ marginBottom: 12, paddingBottom: 12, borderBottom: '1px solid rgba(40,70,130,.25)' }}>
              <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'baseline', marginBottom: 1 }}>
                <span style={{ color: 'rgba(200,218,255,.88)', fontSize: 10.5, fontWeight: 'bold' }}>{drug}</span>
                <span style={{ color: `${mColor}cc`, fontSize: 8.5, letterSpacing: 1, background: `${mColor}18`, padding: '1px 7px', borderRadius: 10, border: `1px solid ${mColor}33` }}>{phase}</span>
              </div>
              <div style={{ color: 'rgba(190,210,255,.9)', fontSize: 8.5, marginBottom: 5 }}>{sponsor} · {trial}</div>
              <ProgressBar progress={progress} color={mColor} />
              <div style={{ color: `${mColor}`, fontSize: 9 }}>{note}</div>
            </div>
          ))}
        </div>
      ))}
    </div>
  );
}

function NetworkTab() {
  const [animated, setAnimated] = useState(false);
  useEffect(() => { const t = setTimeout(() => setAnimated(true), 150); return () => clearTimeout(t); }, []);

  function getNode(id: string) { return NET_NODES.find(n => n.id === id)!; }

  return (
    <div>
      <div style={{ color: 'rgba(185,210,255,.9)', fontSize: 8.5, letterSpacing: 2, marginBottom: 10, textAlign: 'center' }}>
        CO-MUTATION FREQUENCY IN NSCLC ADENOCARCINOMA
      </div>
      <svg viewBox="0 0 340 295" width="100%" style={{ display: 'block', overflow: 'visible' }}>
        <defs>
          {NET_NODES.map(n => (
            <radialGradient key={n.id} id={`rg-${n.id}`} cx="50%" cy="50%" r="50%">
              <stop offset="0%" stopColor={n.color} stopOpacity="0.95" />
              <stop offset="100%" stopColor={n.color} stopOpacity="0.15" />
            </radialGradient>
          ))}
          <filter id="glow">
            <feGaussianBlur stdDeviation="3" result="blur" />
            <feMerge><feMergeNode in="blur" /><feMergeNode in="SourceGraphic" /></feMerge>
          </filter>
        </defs>

        {/* Edges */}
        {NET_EDGES.map(({ a, b, pct, lx, ly }) => {
          const na = getNode(a), nb = getNode(b);
          const len = Math.hypot(nb.x - na.x, nb.y - na.y);
          const sw = pct / 13;
          return (
            <g key={`${a}-${b}`}>
              <line x1={na.x} y1={na.y} x2={nb.x} y2={nb.y}
                stroke={na.color} strokeWidth={sw} strokeOpacity={0.18}
                strokeDasharray={len} strokeDashoffset={animated ? 0 : len}
                style={{ transition: 'stroke-dashoffset 1.1s ease' }} />
              <line x1={na.x} y1={na.y} x2={nb.x} y2={nb.y}
                stroke={na.color} strokeWidth={sw * 0.4} strokeOpacity={0.55}
                strokeDasharray={len} strokeDashoffset={animated ? 0 : len}
                style={{ transition: 'stroke-dashoffset 1.1s ease 0.1s' }} />
              <text x={lx} y={ly} fill={na.color} fontSize={7.5} opacity={animated ? 0.75 : 0}
                textAnchor="middle" fontFamily="Courier New"
                style={{ transition: 'opacity 0.8s ease 0.9s' }}>
                {pct}%
              </text>
            </g>
          );
        })}

        {/* Nodes */}
        {NET_NODES.map(n => (
          <g key={n.id} filter="url(#glow)">
            <circle cx={n.x} cy={n.y} r={n.r + 10} fill={n.color} opacity={0.05} />
            <circle cx={n.x} cy={n.y} r={n.r + 5}  fill={n.color} opacity={0.1}  />
            <circle cx={n.x} cy={n.y} r={n.r}       fill={`url(#rg-${n.id})`} stroke={n.color} strokeWidth={1.5} strokeOpacity={0.85} />
            <text x={n.x} y={n.y - n.r - 5} fill={n.color} fontSize={9} fontWeight="bold"
              textAnchor="middle" fontFamily="Courier New" letterSpacing="1">{n.id}</text>
            <text x={n.x} y={n.y + 4} fill="rgba(230,242,255,.95)" fontSize={7}
              textAnchor="middle" fontFamily="Courier New">{n.freq}</text>
          </g>
        ))}
      </svg>

      <div style={{ display: 'flex', flexWrap: 'wrap', gap: 10, justifyContent: 'center', marginTop: 6 }}>
        {NET_NODES.map(n => (
          <div key={n.id} style={{ display: 'flex', alignItems: 'center', gap: 5, fontSize: 8.5, fontFamily: 'Courier New', color: n.color }}>
            <div style={{ width: 7, height: 7, borderRadius: '50%', background: n.color, boxShadow: `0 0 5px ${n.color}` }} />
            {n.id} <span style={{ color: 'rgba(180,205,255,.85)' }}>{n.freq}</span>
          </div>
        ))}
      </div>

      <div style={{ marginTop: 14, padding: '10px 12px', background: 'rgba(8,16,45,.55)', borderRadius: 7, border: '1px solid rgba(255,50,100,.2)' }}>
        <div style={{ color: '#ff3366', fontSize: 8.5, letterSpacing: 1.5, marginBottom: 5 }}>● KEY INSIGHT</div>
        <div style={{ color: 'rgba(225,238,255,.95)', fontSize: 9.5, lineHeight: 1.65 }}>
          STK11 is the highest-connectivity node — co-mutating with KRAS (54%), TP53 (44%), CDKN2A (37%), and KEAP1 (27%). The STK11/KEAP1 combination drives ATR dependency exploited by Ceralasertib in the Phase III LATIFY trial.
        </div>
      </div>
    </div>
  );
}

function MapTab() {
  const [hovered, setHovered] = useState<typeof SITES[0] | null>(null);
  const tickRef = useRef(0);
  const [tick, setTick] = useState(0);

  useEffect(() => {
    const iv = setInterval(() => { tickRef.current += 1; setTick(tickRef.current); }, 120);
    return () => clearInterval(iv);
  }, []);

  return (
    <div>
      <div style={{ color: 'rgba(185,210,255,.9)', fontSize: 8.5, letterSpacing: 2, marginBottom: 8, textAlign: 'center' }}>
        GLOBAL RESEARCH · CLINICAL TRIALS + QUANTUM COMPUTING
      </div>

      <div style={{ position: 'relative', borderRadius: 8, overflow: 'hidden', border: '1px solid rgba(50,90,190,.25)' }}>
        <svg viewBox={`0 0 ${MAP_W} ${MAP_H}`} width="100%" style={{ display: 'block', background: 'rgba(4,10,32,.85)' }}>
          {/* Grid */}
          {[-60,-30,0,30,60].map(lat => (
            <line key={`lat${lat}`} x1={0} y1={sy(lat)} x2={MAP_W} y2={sy(lat)} stroke="rgba(40,70,140,.18)" strokeWidth={0.5} />
          ))}
          {[-120,-60,0,60,120].map(lon => (
            <line key={`lon${lon}`} x1={sx(lon)} y1={0} x2={sx(lon)} y2={MAP_H} stroke="rgba(40,70,140,.18)" strokeWidth={0.5} />
          ))}
          {/* Equator */}
          <line x1={0} y1={sy(0)} x2={MAP_W} y2={sy(0)} stroke="rgba(60,100,200,.28)" strokeWidth={0.8} strokeDasharray="4,4" />

          {/* Continents */}
          {CONTINENTS.map((coords, ci) => (
            <polygon key={ci}
              points={coords.map(([lo, la]) => `${sx(lo).toFixed(1)},${sy(la).toFixed(1)}`).join(' ')}
              fill="rgba(25,48,110,.6)" stroke="rgba(60,100,200,.38)" strokeWidth={0.7} />
          ))}

          {/* Trial sites */}
          {SITES.map((site, i) => {
            const x = sx(site.lon), y = sy(site.lat);
            const phase = ((tick * 0.08 + i * 0.55) % 1);
            const pr = 3 + phase * 12;
            const po = (1 - phase) * 0.65;
            return (
              <g key={i} style={{ cursor: 'pointer' }}
                onMouseEnter={() => setHovered(site)}
                onMouseLeave={() => setHovered(null)}>
                {/* Pulse ring */}
                <circle cx={x} cy={y} r={pr} fill="none" stroke={site.color} strokeWidth={0.8} opacity={po} />
                {/* Core */}
                <circle cx={x} cy={y} r={3.2} fill={site.color} opacity={0.88} />
                <circle cx={x} cy={y} r={3.2} fill="none" stroke="rgba(255,255,255,.4)" strokeWidth={0.6} />
              </g>
            );
          })}
        </svg>

        {/* Hover tooltip */}
        {hovered && (
          <div style={{ position: 'absolute', bottom: 8, left: '50%', transform: 'translateX(-50%)', background: 'rgba(4,10,32,.97)', border: `1px solid ${hovered.color}55`, borderRadius: 7, padding: '7px 13px', fontFamily: 'Courier New', fontSize: 9, whiteSpace: 'nowrap', pointerEvents: 'none', zIndex: 10 }}>
            <div style={{ color: hovered.color, fontWeight: 'bold', letterSpacing: 1.5, marginBottom: 2 }}>{hovered.name}</div>
            <div style={{ color: 'rgba(200,220,255,.92)', marginBottom: 1 }}>{hovered.loc}</div>
            <div style={{ color: 'rgba(210,228,255,.92)' }}>
              <span style={{ color: hovered.color }}>Target: </span>{hovered.targets}
              <span style={{ color: 'rgba(180,210,255,.85)' }}> · </span>
              {hovered.drug}
            </div>
          </div>
        )}
      </div>

      {/* Legend */}
      <div style={{ display: 'flex', flexWrap: 'wrap', gap: 10, justifyContent: 'center', marginTop: 10 }}>
        {[
          { color: '#3399ff', label: 'TP53' },
          { color: '#33ffaa', label: 'KEAP1' },
          { color: '#ffaa33', label: 'CDKN2A' },
          { color: '#ff3366', label: 'STK11' },
          { color: '#cc44ff', label: 'Quantum' },
        ].map(({ color, label }) => (
          <div key={label} style={{ display: 'flex', alignItems: 'center', gap: 5, fontSize: 8.5, fontFamily: 'Courier New', color }}>
            <div style={{ width: 7, height: 7, borderRadius: '50%', background: color, boxShadow: `0 0 5px ${color}` }} />
            {label}
          </div>
        ))}
      </div>

      <div style={{ textAlign: 'center', marginTop: 8, fontSize: 8, color: 'rgba(170,200,255,.8)', fontFamily: 'Courier New', letterSpacing: 1.5 }}>
        {SITES.length} ACTIVE SITES · 10 COUNTRIES · 4 CONTINENTS
      </div>

      {/* Site list */}
      <div style={{ marginTop: 14, display: 'flex', flexDirection: 'column', gap: 6 }}>
        {SITES.map((site, i) => (
          <div key={i} onMouseEnter={() => setHovered(site)} onMouseLeave={() => setHovered(null)}
            style={{ display: 'flex', gap: 10, alignItems: 'flex-start', padding: '7px 10px', borderRadius: 6, background: hovered?.name === site.name ? 'rgba(20,40,100,.6)' : 'rgba(8,16,45,.4)', border: `1px solid ${hovered?.name === site.name ? site.color + '44' : 'rgba(40,70,130,.2)'}`, cursor: 'default', transition: 'all 0.2s' }}>
            <div style={{ width: 7, height: 7, borderRadius: '50%', background: site.color, boxShadow: `0 0 6px ${site.color}`, marginTop: 2, flexShrink: 0 }} />
            <div>
              <div style={{ color: site.color, fontSize: 9.5, fontFamily: 'Courier New', letterSpacing: 0.5 }}>{site.name}</div>
              <div style={{ color: 'rgba(185,210,255,.88)', fontSize: 8, fontFamily: 'Courier New' }}>{site.loc} · {site.drug}</div>
            </div>
          </div>
        ))}
      </div>
    </div>
  );
}

// ─── Main DataPanel ────────────────────────────────────────────────────────────
const TABS = [
  { id: 'trials'  as const, icon: '⬡', label: 'TRIALS'  },
  { id: 'network' as const, icon: '◈', label: 'NETWORK' },
  { id: 'map'     as const, icon: '◉', label: 'MAP'     },
];

export default function DataPanel() {
  const [open, setOpen] = useState(false);
  const [tab, setTab] = useState<'trials' | 'network' | 'map'>('trials');

  return (
    <>
      {/* Slide toggle tab */}
      <div onClick={() => setOpen(o => !o)} style={{
        position: 'fixed', right: open ? 382 : 0, top: '50%', transform: 'translateY(-50%)',
        background: 'rgba(4,10,32,.95)', border: '1px solid rgba(70,130,255,.32)',
        borderRight: open ? '1px solid rgba(70,130,255,.32)' : 'none',
        borderRadius: '8px 0 0 8px', padding: '16px 7px', cursor: 'pointer',
        writingMode: 'vertical-rl', fontFamily: 'Courier New', fontSize: 9,
        letterSpacing: 3, color: 'rgba(100,155,255,.8)', zIndex: 1001,
        transition: 'right 0.36s cubic-bezier(.4,0,.2,1)', userSelect: 'none',
        boxShadow: '-4px 0 18px rgba(0,0,0,.5)',
      }}>
        {open ? '◀ CLOSE' : '▶ DATA'}
      </div>

      {/* Panel */}
      <div style={{
        position: 'fixed', right: 0, top: 0, height: '100vh', width: 382,
        background: 'rgba(2,6,20,.97)', borderLeft: '1px solid rgba(60,110,230,.22)',
        zIndex: 1000, display: 'flex', flexDirection: 'column',
        transform: open ? 'translateX(0)' : 'translateX(100%)',
        transition: 'transform 0.36s cubic-bezier(.4,0,.2,1)',
        fontFamily: 'Courier New',
      }}>

        {/* Header */}
        <div style={{ padding: '14px 16px 0', borderBottom: '1px solid rgba(60,110,230,.18)', flexShrink: 0 }}>
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: 10 }}>
            <div style={{ fontSize: 10.5, letterSpacing: 3, color: 'rgba(180,210,255,.95)' }}>NSCLC · DATA PANEL</div>
            <div style={{ fontSize: 8, color: 'rgba(160,195,255,.8)', letterSpacing: 1 }}>IBM sqDRIFT 2026</div>
          </div>
          {/* Tab buttons */}
          <div style={{ display: 'flex', gap: 2 }}>
            {TABS.map(t => (
              <button key={t.id} onClick={() => setTab(t.id)} style={{
                flex: 1, padding: '8px 4px', fontSize: 9, letterSpacing: 1.5, fontFamily: 'Courier New', cursor: 'pointer',
                background: tab === t.id ? 'rgba(40,90,210,.28)' : 'transparent',
                border: `1px solid ${tab === t.id ? 'rgba(70,140,255,.5)' : 'rgba(60,100,190,.18)'}`,
                borderBottom: 'none', borderRadius: '6px 6px 0 0',
                color: tab === t.id ? 'rgba(220,235,255,1)' : 'rgba(160,190,255,.75)',
                transition: 'all 0.2s',
              }}>
                {t.icon} {t.label}
              </button>
            ))}
          </div>
        </div>

        {/* Scrollable content */}
        <div style={{ flex: 1, overflowY: 'auto', overflowX: 'hidden', padding: '14px 14px 20px' }}>
          {tab === 'trials'  && <TrialsTab />}
          {tab === 'network' && <NetworkTab />}
          {tab === 'map'     && <MapTab />}
        </div>

        {/* Footer */}
        <div style={{ padding: '8px 14px', borderTop: '1px solid rgba(60,100,190,.25)', fontSize: 7.5, letterSpacing: 1.5, color: 'rgba(150,185,255,.7)', textAlign: 'center', flexShrink: 0 }}>
          IBM sqDRIFT · NSCLC QUANTUM SIMULATION · CGU DISSERTATION 2026
        </div>
      </div>
    </>
  );
}
