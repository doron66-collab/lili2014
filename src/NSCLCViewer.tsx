import { useEffect, useRef, useState } from 'react';
import TP53LoopsViewer from './TP53LoopsViewer';
import PDBMolViewer from './PDBMolViewer';
import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';
import { EffectComposer } from 'three/examples/jsm/postprocessing/EffectComposer.js';
import { RenderPass } from 'three/examples/jsm/postprocessing/RenderPass.js';
import { UnrealBloomPass } from 'three/examples/jsm/postprocessing/UnrealBloomPass.js';

// ── Seeded RNG ────────────────────────────────────────────────────────────────
function seededRNG(seed: number) {
  let s = (seed | 0) >>> 0 || 1;
  return () => { s = (Math.imul(1664525, s) + 1013904223) | 0; return (s >>> 0) / 4294967296; };
}
function genPocket(seed: number, cx: number, cy: number, cz: number, sp: number) {
  const rng = seededRNG(seed + 8888);
  return Array.from({ length: 18 }, () => ({
    x: cx + (rng() - 0.5) * sp, y: cy + (rng() - 0.5) * sp, z: cz + (rng() - 0.5) * sp, r: 0.07 + rng() * 0.09,
  }));
}

// ── Alpha-helix ribbon path ───────────────────────────────────────────────────
function makeHelixPath(seed: number, turns: number, radius: number, rise: number, offsetX: number, offsetY: number, offsetZ: number) {
  const rng = seededRNG(seed + 1000);
  const pts: THREE.Vector3[] = [];
  const totalPoints = turns * 20;
  const tiltX = (rng() - 0.5) * 0.6, tiltZ = (rng() - 0.5) * 0.6;
  for (let i = 0; i <= totalPoints; i++) {
    const t = i / totalPoints, angle = t * turns * Math.PI * 2;
    const x = Math.cos(angle) * radius + offsetX;
    const y = (t - 0.5) * rise * turns + offsetY;
    const z = Math.sin(angle) * radius + offsetZ;
    pts.push(new THREE.Vector3(x + tiltX * y, y, z + tiltZ * y));
  }
  return pts;
}

// ── Beta-sheet path ───────────────────────────────────────────────────────────
function makeSheetPath(seed: number, length: number, ox: number, oy: number, oz: number) {
  const rng = seededRNG(seed + 2000);
  const pts: THREE.Vector3[] = [];
  const angle = rng() * Math.PI;
  for (let i = 0; i <= 12; i++) {
    const t = (i / 12 - 0.5) * length;
    pts.push(new THREE.Vector3(ox + Math.cos(angle) * t, oy + (rng() - 0.5) * 0.3, oz + Math.sin(angle) * t));
  }
  return pts;
}

// ── Mutation data ─────────────────────────────────────────────────────────────
const MUT = [
  {
    id: 'TP53', variant: 'Y220C', drug: 'Rezatapopt', sub: 'PC14586',
    phase: 'Phase II — PYNNACLE', color: 0x3399ff, pCol: 0xff6633, dCol: 0xffbb00,
    mech: 'Binds Y220C cryptic pocket → restores wild-type p53 conformation → reactivates tumor suppression',
    pocket: genPocket(101, 0.5, 0.65, 0.95, 0.85),
    atoms: [[0,0,0,.22,0xffbb00],[.52,.28,.1,.17,0xff9900],[-.52,.28,-.1,.17,0xffdd33],[.15,-.50,.25,.18,0xffcc00],[-.15,-.50,-.2,.16,0xff8800],[.70,-.15,-.30,.13,0xffaa00],[-.65,.10,.35,.13,0xffdd00],[0,.60,-.15,.15,0xffcc44]],
    helices: [[101,3.2,1.0,1.4,-1.4,0.2,0.3],[202,2.5,0.8,1.2,1.2,-0.5,-0.4],[303,2.0,0.7,1.0,0.1,1.1,-1.0],[404,1.8,0.65,0.9,-0.5,-1.2,0.7]],
    sheets:  [[111,2.2,-0.8,0.4,0.6],[222,1.8,0.9,-0.6,-0.5],[333,2.0,-0.2,-0.9,0.8]],
  },
  {
    id: 'KEAP1', variant: 'LOF', drug: 'VVD-065', sub: 'Covalent NRF2 inhibitor',
    phase: 'Phase I — NCT05954312', color: 0x33ffaa, pCol: 0xaaff33, dCol: 0xff33aa,
    mech: 'Covalent NRF2 inhibition in KEAP1-deficient cells → restores chemosensitivity',
    pocket: genPocket(202, -0.3, 0.72, 0.82, 0.92),
    atoms: [[0,0,0,.21,0xff33aa],[.48,.22,.15,.16,0xff55cc],[-.48,.22,-.15,.16,0xff11aa],[0,.55,.22,.18,0xff77dd],[.32,-.40,.28,.15,0xee33bb],[-.32,-.40,-.18,.15,0xff44cc],[.62,.08,-.32,.12,0xff22bb],[-.62,.08,.22,.12,0xff66dd]],
    helices: [[505,3.0,0.9,1.3,1.3,0.3,0.2],[606,2.2,0.75,1.1,-1.1,-0.4,0.5],[707,2.8,0.85,1.2,0.2,1.0,-0.9],[808,1.9,0.7,1.0,-0.4,-1.0,0.6]],
    sheets:  [[515,2.1,0.7,0.5,-0.7],[626,1.9,-0.9,-0.5,0.5],[737,2.3,0.1,-1.0,-0.4]],
  },
  {
    id: 'CDKN2A', variant: 'p16 loss', drug: 'Palbociclib + Olaparib', sub: 'CDK4/6 + PARP synergy',
    phase: 'Phase I/II — PARP1 degradation 2024', color: 0xffaa33, pCol: 0xffee22, dCol: 0x33aaff,
    mech: 'CDK4/6 inhibition restores G1 checkpoint lost by p16 deletion + PARP1 synthetic lethality',
    pocket: genPocket(303, 0.6, -0.42, 0.82, 0.88),
    atoms: [[0,0,0,.22,0x33aaff],[.50,.20,0,.17,0x1188ff],[-.42,.32,.18,.17,0x55ccff],[.08,-.48,.32,.16,0x33aaff],[-.28,-.32,-.32,.15,0x1166ff],[.58,-.22,-.22,.13,0x77ddff],[-.55,.12,.38,.13,0x44bbff],[0,.58,-.18,.16,0x22aaff]],
    helices: [[909,2.8,0.85,1.25,1.1,-0.3,0.4],[110,2.4,0.8,1.15,-1.0,0.4,-0.3],[211,2.1,0.7,1.05,0.3,-1.1,0.8],[312,1.7,0.6,0.95,-0.6,1.0,-0.6]],
    sheets:  [[919,2.0,-0.7,0.3,0.7],[120,1.7,0.8,-0.7,-0.5],[231,2.2,-0.1,-0.8,0.6]],
  },
  {
    id: 'STK11', variant: 'LKB1 loss', drug: 'Ceralasertib', sub: '+ Durvalumab (ATR inhibitor)',
    phase: 'Phase III — LATIFY NCT05450692', color: 0xff3366, pCol: 0xff6688, dCol: 0x99ff33,
    mech: 'ATR inhibition targets replication stress created by LKB1 loss (sqDRIFT-modellable synthetic lethality)',
    pocket: genPocket(404, -0.5, -0.32, 0.88, 0.92),
    atoms: [[0,0,0,.22,0x99ff33],[.50,.25,.12,.16,0x77ee11],[-.50,.25,-.12,.16,0xbbff55],[.18,-.50,.22,.17,0x99ff33],[-.18,-.50,-.22,.17,0x77ee22],[.68,-.12,-.28,.13,0xaaff44],[-.68,-.12,.28,.13,0x88ff22],[.28,.52,.30,.15,0xccffaa]],
    helices: [[413,3.1,0.95,1.35,-1.2,0.4,0.2],[514,2.6,0.82,1.18,1.0,-0.5,-0.4],[615,2.2,0.72,1.08,0.2,1.2,-0.7],[716,1.8,0.65,0.98,-0.5,-1.1,0.6]],
    sheets:  [[423,2.1,0.6,0.4,-0.8],[534,1.8,-0.8,-0.6,0.4],[645,2.0,0.1,-0.9,-0.5]],
  },
  {
    // TP53 C275F — structural β-sandwich core mutant (p.Cys275Phe)
    // Detected: Emek Medical Center MI25-0349 · June 2025
    // sqDRIFT active space: ~44 active electrons · ~44 orbitals · ~88 qubits
    id: 'TP53', variant: 'C275F', drug: 'Eprenetapopt', sub: 'APR-246 · Pan-mutant p53 reactivator',
    phase: 'No C275F-specific trial · Research stage', color: 0xffdd22, pCol: 0xff7700, dCol: 0x44ffee,
    mech: 'p.Cys275Phe replaces core β-sandwich Cys with bulky Phe → hydrophobic core collapse → global misfolding. No mutation-specific drug exists. Quantum simulation (sqDRIFT) of the π-electron pocket is the frontier approach. Active space: ~44e · ~88 qubits',
    pocket: genPocket(505, 0.12, 0.22, 0.30, 0.70),
    atoms: [[0,0,0,.21,0x44ffee],[.50,.24,.08,.16,0x22ddcc],[-.44,.28,-.12,.16,0x66ffee],[.12,-.46,.28,.17,0x44ffee],[-.30,-.28,-.28,.15,0x22ccdd],[.62,-.16,-.24,.13,0x55eecc],[-.60,.12,.32,.13,0x44ffdd],[.24,.50,.24,.15,0x66eedd],[.38,.10,-.46,.13,0x33eedd]],
    helices: [[551,2.6,0.80,1.15, 1.1, 0.2, 0.2],[662,2.0,0.68,1.05,-0.8,-0.4, 0.5],[773,1.6,0.60,0.90, 0.1,-0.9, 0.7]],
    sheets:  [[551,2.5, 0.7, 0.4, 0.5],[662,2.3,-0.8, 0.4,-0.4],[773,2.1, 0.2,-0.7, 0.6],[884,1.9,-0.4, 0.6,-0.5],[995,2.0, 0.6,-0.3,-0.7]],
  },
] as const;

// PDB crystallographic data per mutation (parallel to MUT array)
const PDB_MAP = [
  { pdb: '2VUK', chain: 'A', highlightRes: [220] },        // TP53 Y220C
  { pdb: '2FLU', chain: 'X', highlightRes: [] as number[] }, // KEAP1 LOF
  { pdb: '2A5E', chain: 'A', highlightRes: [] as number[] }, // CDKN2A p16 INK4a
  { pdb: '2QK7', chain: 'A', highlightRes: [] as number[] }, // STK11 LKB1
  { pdb: '2OCJ', chain: 'A', highlightRes: [275] },         // TP53 C275F
] as const;

// ── Patient-report bridge ─────────────────────────────────────────────────────
// The main platform (Assignment10_Prototype.html, served at the same origin)
// writes the active NGS/TEMPUS variant set to localStorage under this key when a
// report is loaded. This viewer reads it on mount and renders THOSE mutations
// instead of the built-in demo set above. If the key is absent (viewer opened
// standalone) we fall back to all five demo models, preserving prior behaviour.
const BRIDGE_KEY = 'solange_3d_variants';

// simulation_id → index into MUT / PDB_MAP (a bespoke, hand-built 3D model)
const MODEL_INDEX: Record<string, number> = {
  TP53_Y220C: 0,
  KEAP1_LOF:  1,
  CDKN2A_P16: 2,
  STK11_LKB1: 3,
  TP53_C275F: 4,
};
// Genes with exactly ONE bespoke model — allow a gene-level match when the
// simulation_id doesn't line up exactly (e.g. an uploaded KEAP1 frameshift or a
// CDKN2A deletion still maps to the single KEAP1 / CDKN2A model). TP53 is
// deliberately excluded: it has two models, so an unmatched TP53 variant falls
// through to the real-PDB path rather than guessing.
const GENE_MODEL_INDEX: Record<string, number> = {
  KEAP1: 1, CDKN2A: 2, STK11: 3,
};

interface BridgeVariant {
  gene: string;
  mutation: string;
  hgvs?: string;
  simulation_id: string;
  tier?: string;
  badge_type?: string;
  allele_frequency?: number | null;
  active_electrons?: number | null;
  full_qubits?: number | null;
  color?: string;     // CSS hex string, e.g. '#E8A020'
  source?: string;    // 'ngs' | 'research' | …
}

type ViewEntry =
  | { kind: 'full'; modelIndex: number; meta: BridgeVariant }
  | { kind: 'pdb';  meta: BridgeVariant };

// Build the ordered list the viewer should show, from the patient report if
// present, otherwise the built-in demo set.
function buildViewList(): { entries: ViewEntry[]; reportId: string | null } {
  try {
    const raw = typeof localStorage !== 'undefined' ? localStorage.getItem(BRIDGE_KEY) : null;
    if (raw) {
      const data = JSON.parse(raw);
      const vars: BridgeVariant[] = Array.isArray(data?.variants) ? data.variants : [];
      if (vars.length) {
        const entries: ViewEntry[] = vars.map(v => {
          const simId = (v.simulation_id || '').toUpperCase();
          const gene  = (v.gene || '').toUpperCase();
          let mi = MODEL_INDEX[simId];
          if (mi === undefined && GENE_MODEL_INDEX[gene] !== undefined) mi = GENE_MODEL_INDEX[gene];
          return mi !== undefined
            ? { kind: 'full', modelIndex: mi, meta: v }
            : { kind: 'pdb', meta: v };
        });
        return { entries, reportId: data.report_id || null };
      }
    }
  } catch {
    // Corrupt/blocked storage — fall through to the demo set.
  }
  // Fallback: the original five demo models.
  const simIds = Object.keys(MODEL_INDEX); // insertion order matches MUT indices
  const entries: ViewEntry[] = MUT.map((m, i) => ({
    kind: 'full',
    modelIndex: i,
    meta: { gene: m.id, mutation: m.variant, simulation_id: simIds[i], source: 'demo' },
  }));
  return { entries, reportId: null };
}

// Resolve a real structure for an unknown gene via the same backend endpoint the
// platform's "Structure Sources" card uses (UniProt → RCSB PDB → AlphaFold).
async function resolveStructure(gene: string): Promise<{ pdb: string; chain: string; url?: string } | null> {
  const base = (typeof window !== 'undefined' && (window as any).QCAIHPC_API_BASE)
    || 'https://qcaihpc-simulation-api.onrender.com';
  try {
    const r = await fetch(`${base}/api/pdb/lookup/gene/${encodeURIComponent(gene)}`).then(res => res.json());
    if (r?.pdb_ids?.length) return { pdb: r.pdb_ids[0], chain: 'A' };
    if (r?.alphafold_model_url) return { pdb: `AF-${r.uniprot_id || gene}`, chain: 'A', url: r.alphafold_model_url };
  } catch {
    // Backend asleep / offline — caller surfaces a message.
  }
  return null;
}

function cssHexToInt(hex?: string): number {
  if (!hex) return 0x06b6d4;
  const n = parseInt(hex.replace('#', ''), 16);
  return Number.isFinite(n) ? n : 0x06b6d4;
}

// ── Electron-density particle system ─────────────────────────────────────────
function makeElectronCloud(cx: number, cy: number, cz: number, color: number, n = 180) {
  const rng = seededRNG(cx * 1000 + cy * 100 + cz * 10 + color);
  const positions = new Float32Array(n * 3);
  const phases    = new Float32Array(n);
  for (let i = 0; i < n; i++) {
    const r = 0.55 + rng() * 0.65;
    const th = rng() * Math.PI * 2, ph = Math.acos(2 * rng() - 1);
    positions[i*3]   = cx + r * Math.sin(ph) * Math.cos(th);
    positions[i*3+1] = cy + r * Math.sin(ph) * Math.sin(th);
    positions[i*3+2] = cz + r * Math.cos(ph);
    phases[i] = rng() * Math.PI * 2;
  }
  const geo = new THREE.BufferGeometry();
  geo.setAttribute('position', new THREE.BufferAttribute(positions, 3));
  const mat = new THREE.PointsMaterial({ color: new THREE.Color(color), size: 0.055, transparent: true, opacity: 0, blending: THREE.AdditiveBlending, depthWrite: false });
  const pts = new THREE.Points(geo, mat);
  (pts as any).phases = phases;
  (pts as any).basePositions = positions.slice();
  return pts;
}

// ── sqDRIFT wavefunction rings ────────────────────────────────────────────────
function makeWavefunctionRings(color: number) {
  const grp = new THREE.Group();
  const ringCount = 5;
  for (let r = 0; r < ringCount; r++) {
    const pts: THREE.Vector3[] = [];
    const segments = 80;
    for (let i = 0; i <= segments; i++) {
      const angle = (i / segments) * Math.PI * 2;
      pts.push(new THREE.Vector3(Math.cos(angle), 0, Math.sin(angle)));
    }
    const geo = new THREE.BufferGeometry().setFromPoints(pts);
    const mat = new THREE.LineBasicMaterial({ color: new THREE.Color(color), transparent: true, opacity: 0, blending: THREE.AdditiveBlending });
    const ring = new THREE.Line(geo, mat);
    ring.scale.setScalar(0.35 + r * 0.22);
    ring.rotation.x = (r * Math.PI) / ringCount + 0.3;
    ring.rotation.z = (r * Math.PI * 0.7) / ringCount;
    grp.add(ring);
  }
  grp.visible = false;
  return grp;
}

interface PdbMutInfo {
  id: string; variant: string; pdb: string; chain: string;
  highlightRes: number[]; color: number; drug: string; phase: string;
  url?: string;   // optional explicit structure URL (e.g. AlphaFold model)
}

export default function NSCLCViewer() {
  const mountRef = useRef<HTMLDivElement>(null);
  const [showLoops, setShowLoops] = useState(false);
  const [pdbMut, setPdbMut] = useState<PdbMutInfo | null>(null);
  // Read the patient report (or demo fallback) once, at mount.
  const [view] = useState(buildViewList);

  useEffect(() => {
    const el = mountRef.current;
    if (!el) return;

    // ── Resolve the view list into scene-drivable entries ───────────────────
    // 'full' entries (genes we have a hand-built 3D model for) drive the Three.js
    // scene below. 'pdb' entries (everything else in the report) are reachable
    // from the nav and open the real crystal/AlphaFold structure on demand.
    const SCENE = view.entries
      .filter((e): e is Extract<ViewEntry, { kind: 'full' }> => e.kind === 'full')
      .map(e => ({ model: MUT[e.modelIndex], pdb: PDB_MAP[e.modelIndex], meta: e.meta }));
    const W = el.clientWidth, H = el.clientHeight;

    // ── Renderer ──────────────────────────────────────────────────────────────
    const renderer = new THREE.WebGLRenderer({ antialias: true });
    renderer.setSize(W, H);
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    renderer.toneMapping = THREE.ACESFilmicToneMapping;
    renderer.toneMappingExposure = 1.2;
    el.appendChild(renderer.domElement);

    const scene = new THREE.Scene();
    scene.background = new THREE.Color(0x020610);
    scene.fog = new THREE.FogExp2(0x020610, 0.032);

    const camera = new THREE.PerspectiveCamera(52, W / H, 0.1, 200);
    camera.position.set(0, 0.5, 10.5);

    const composer = new EffectComposer(renderer);
    composer.addPass(new RenderPass(scene, camera));
    const bloom = new UnrealBloomPass(new THREE.Vector2(W, H), 1.1, 0.55, 0.45);
    composer.addPass(bloom);

    const controls = new OrbitControls(camera, renderer.domElement);
    controls.enableDamping = true; controls.dampingFactor = 0.07;
    controls.minDistance = 4; controls.maxDistance = 20;

    scene.add(new THREE.AmbientLight(0x0a1428, 5));
    const L1 = new THREE.PointLight(0x4488ff, 90, 35);
    const L2 = new THREE.PointLight(0xff4488, 60, 28);
    const L3 = new THREE.PointLight(0x44ffaa, 50, 25);
    L1.position.set(6,5,5); L2.position.set(-5,-3,3); L3.position.set(0,6,-5);
    [L1,L2,L3].forEach(l => scene.add(l));

    // ── Helpers ───────────────────────────────────────────────────────────────
    const hx = (h: number) => new THREE.Color(h);

    function mkBond(p1: THREE.Vector3, p2: THREE.Vector3, col: number) {
      const d = new THREE.Vector3().subVectors(p2, p1), len = d.length();
      const m = new THREE.Mesh(
        new THREE.CylinderGeometry(0.033, 0.033, len, 8),
        new THREE.MeshStandardMaterial({ color: hx(col), emissive: hx(col), emissiveIntensity: 0.35, transparent: true, opacity: 0, roughness: 0.35 })
      );
      m.position.copy(p1).add(p2).multiplyScalar(0.5);
      m.quaternion.setFromUnitVectors(new THREE.Vector3(0,1,0), d.normalize());
      return m;
    }

    // ── Build scene objects ───────────────────────────────────────────────────
    const roots: THREE.Group[] = [];
    const pocketGrps: THREE.Group[] = [];
    const drugGrps: THREE.Group[] = [];
    const pocketCenters: THREE.Vector3[] = [];
    const electronClouds: THREE.Points[] = [];
    const wavefunctionGrps: THREE.Group[] = [];
    const heatmapMats: THREE.MeshStandardMaterial[][] = [];

    for (let i = 0; i < SCENE.length; i++) {
      const m = SCENE[i].model;
      const root = new THREE.Group();
      root.visible = (i === 0);
      scene.add(root);

      // ── Alpha helices ─────────────────────────────────────────────────────
      const helixMat = new THREE.MeshStandardMaterial({ color: hx(m.color), emissive: hx(m.color), emissiveIntensity: 0.12, roughness: 0.35, metalness: 0.25 });
      (m.helices as number[][]).forEach(([seed, turns, radius, rise, ox, oy, oz]) => {
        const pts = makeHelixPath(seed, turns, radius, rise, ox, oy, oz);
        const curve = new THREE.CatmullRomCurve3(pts);
        const tubeGeo = new THREE.TubeGeometry(curve, pts.length * 3, 0.075, 8, false);
        root.add(new THREE.Mesh(tubeGeo, helixMat.clone()));
      });

      // ── Beta sheets ───────────────────────────────────────────────────────
      const sheetMat = new THREE.MeshStandardMaterial({ color: hx(m.color), emissive: hx(m.color), emissiveIntensity: 0.08, roughness: 0.5, metalness: 0.15, side: THREE.DoubleSide });
      (m.sheets as number[][]).forEach(([seed, length, ox, oy, oz]) => {
        const pts = makeSheetPath(seed, length, ox, oy, oz);
        const curve = new THREE.CatmullRomCurve3(pts);
        // Flat ribbon via TubeGeometry with flat profile
        const tubeGeo = new THREE.TubeGeometry(curve, 20, 0.18, 4, false);
        const mesh = new THREE.Mesh(tubeGeo, sheetMat.clone());
        mesh.scale.y = 0.22; // flatten into ribbon
        root.add(mesh);
      });

      // ── Binding pocket (heatmap residues) ─────────────────────────────────
      const pkGrp = new THREE.Group();
      const pktMats: THREE.MeshStandardMaterial[] = [];
      (m.pocket as {x:number,y:number,z:number,r:number}[]).forEach(s => {
        const mat = new THREE.MeshStandardMaterial({ color: hx(m.pCol), emissive: hx(m.pCol), emissiveIntensity: 0.6, roughness: 0.15, metalness: 0.4 });
        const ms = new THREE.Mesh(new THREE.SphereGeometry(s.r, 12, 12), mat);
        ms.position.set(s.x, s.y, s.z);
        pkGrp.add(ms); pktMats.push(mat);
      });
      root.add(pkGrp); pocketGrps.push(pkGrp); heatmapMats.push(pktMats);

      // ── Pocket center ─────────────────────────────────────────────────────
      const pcArr = m.pocket as {x:number,y:number,z:number,r:number}[];
      const pc = pcArr.reduce((a, s) => ({ x:a.x+s.x, y:a.y+s.y, z:a.z+s.z }), {x:0,y:0,z:0});
      pc.x /= pcArr.length; pc.y /= pcArr.length; pc.z /= pcArr.length;
      pocketCenters.push(new THREE.Vector3(pc.x, pc.y, pc.z));

      // ── Electron density cloud ────────────────────────────────────────────
      const cloud = makeElectronCloud(pc.x, pc.y, pc.z, m.pCol);
      root.add(cloud); electronClouds.push(cloud);

      // ── Drug molecule (ball-and-stick) ────────────────────────────────────
      const dgGrp = new THREE.Group();
      (m.atoms as number[][]).forEach(([px,py,pz,r,c], j) => {
        const ms = new THREE.Mesh(new THREE.SphereGeometry(r, 16, 16),
          new THREE.MeshStandardMaterial({ color:hx(c), emissive:hx(c), emissiveIntensity:0.6, roughness:0.15, metalness:0.45, transparent:true, opacity:0 }));
        ms.position.set(px,py,pz); dgGrp.add(ms);
        if (j < m.atoms.length - 1) {
          const [nx,ny,nz] = m.atoms[j+1] as number[];
          dgGrp.add(mkBond(new THREE.Vector3(px,py,pz), new THREE.Vector3(nx,ny,nz), c));
        }
      });

      // ── sqDRIFT wavefunction rings ────────────────────────────────────────
      const wfGrp = makeWavefunctionRings(m.dCol);
      dgGrp.add(wfGrp); wavefunctionGrps.push(wfGrp);

      // ── C275F: Phe275 aromatic π-ring (benzene 6e/6-orbital system) ──────
      if (m.variant === 'C275F') {
        const piGrp = new THREE.Group();
        // Hexagonal ring (benzene)
        for (let ring = 0; ring < 2; ring++) {
          const hexPts: THREE.Vector3[] = [];
          const hexR = ring === 0 ? 0.42 : 0.22;
          for (let k = 0; k <= 6; k++) {
            const a = (k / 6) * Math.PI * 2;
            hexPts.push(new THREE.Vector3(Math.cos(a) * hexR, 0, Math.sin(a) * hexR));
          }
          const hexGeo = new THREE.BufferGeometry().setFromPoints(hexPts);
          const hexMat = new THREE.LineBasicMaterial({ color: ring === 0 ? 0xffdd22 : 0xff9900, transparent: true, opacity: ring === 0 ? 0.85 : 0.5, blending: THREE.AdditiveBlending });
          piGrp.add(new THREE.Line(hexGeo, hexMat));
        }
        // π electron lobes (above/below the ring plane)
        for (const yOff of [-0.18, 0.18]) {
          const lobePts: THREE.Vector3[] = [];
          for (let k = 0; k <= 60; k++) {
            const a = (k / 60) * Math.PI * 2;
            lobePts.push(new THREE.Vector3(Math.cos(a) * 0.30, yOff, Math.sin(a) * 0.30));
          }
          const lobeGeo = new THREE.BufferGeometry().setFromPoints(lobePts);
          const lobeMat = new THREE.LineBasicMaterial({ color: 0xffcc00, transparent: true, opacity: 0.3, blending: THREE.AdditiveBlending });
          piGrp.add(new THREE.Line(lobeGeo, lobeMat));
        }
        piGrp.position.set(pc.x, pc.y, pc.z);
        piGrp.rotation.x = 0.4; piGrp.rotation.z = 0.3;
        piGrp.visible = false;
        root.add(piGrp);
        (root as any).piRing = piGrp;
      }

      const sx = pc.x+5.5, sy = pc.y+3.5, sz = pc.z+2.0;
      dgGrp.position.set(sx, sy, sz);
      (dgGrp as any).startPos = new THREE.Vector3(sx, sy, sz);
      (dgGrp as any).targetPos = pocketCenters[i].clone();
      root.add(dgGrp); drugGrps.push(dgGrp); roots.push(root);
    }

    // ── Stars ─────────────────────────────────────────────────────────────────
    const geo = new THREE.BufferGeometry();
    const N = 600, pos = new Float32Array(N*3);
    for (let i = 0; i < N; i++) { pos[i*3]=(Math.random()-.5)*50; pos[i*3+1]=(Math.random()-.5)*50; pos[i*3+2]=(Math.random()-.5)*50; }
    geo.setAttribute('position', new THREE.BufferAttribute(pos, 3));
    scene.add(new THREE.Points(geo, new THREE.PointsMaterial({ color:0x1a2d4a, size:0.045, transparent:true, opacity:0.6 })));

    // ── State ─────────────────────────────────────────────────────────────────
    let cur = 0, animState = 0, animT = 0, T = 0;

    // ── UI ────────────────────────────────────────────────────────────────────
    const ui = document.createElement('div');
    ui.style.cssText = 'position:absolute;top:0;left:0;width:100%;height:100%;pointer-events:none;';
    el.style.position = 'relative'; el.appendChild(ui);

    const infoBox = document.createElement('div');
    infoBox.style.cssText = 'position:absolute;top:14px;left:14px;background:rgba(2,6,18,.92);border:1px solid rgba(70,150,255,.3);border-radius:12px;padding:18px 22px;max-width:275px;font-family:Courier New,monospace;box-shadow:0 0 32px rgba(40,100,255,.2);backdrop-filter:blur(12px);';
    ui.appendChild(infoBox);

    const rightLabel = document.createElement('div');
    rightLabel.style.cssText = 'position:absolute;top:14px;right:14px;font-family:Courier New,monospace;font-size:10px;letter-spacing:2.5px;color:rgba(160,200,255,.85);text-align:right;';
    ui.appendChild(rightLabel);

    const statusBox = document.createElement('div');
    statusBox.style.cssText = 'position:absolute;bottom:16px;right:14px;background:rgba(2,6,18,.85);border:1px solid rgba(70,240,160,.28);border-radius:6px;padding:7px 13px;font-family:Courier New,monospace;font-size:10px;letter-spacing:1px;';
    ui.appendChild(statusBox);

    // Legend
    const legend = document.createElement('div');
    legend.style.cssText = 'position:absolute;bottom:60px;right:14px;background:rgba(2,6,18,.82);border:1px solid rgba(70,140,255,.2);border-radius:8px;padding:10px 14px;font-family:Courier New,monospace;font-size:9px;letter-spacing:1px;color:rgba(150,180,255,.6);line-height:1.9;';
    legend.innerHTML = `
      <div style="color:rgba(200,220,255,.85);margin-bottom:4px;letter-spacing:2px;">LEGEND</div>
      <div>▬▬ <span style="opacity:.8">α-HELIX / β-SHEET</span></div>
      <div>● <span style="color:rgba(255,150,80,.8)">BINDING POCKET</span></div>
      <div>✦ <span style="color:rgba(180,255,180,.8)">ELECTRON DENSITY</span></div>
      <div>◯ <span style="color:rgba(200,220,255,.8)">sqDRIFT WAVEFUNCTION</span></div>`;
    ui.appendChild(legend);

    // ── Back to Platform — only when launched from a parent window ──────────
    if (typeof window !== 'undefined' && window.opener && !window.opener.closed) {
      const backBtn = document.createElement('button');
      backBtn.innerHTML = '← Back to Platform';
      backBtn.style.cssText = 'position:absolute;top:14px;left:50%;transform:translateX(-50%);background:rgba(8,18,50,.88);border:1px solid rgba(70,140,255,.4);color:rgba(150,200,255,.95);border-radius:18px;padding:7px 18px;cursor:pointer;font-size:11px;font-family:Courier New,monospace;letter-spacing:1px;pointer-events:all;transition:background .2s;z-index:20;';
      backBtn.onmouseover = () => backBtn.style.background = 'rgba(25,65,155,.9)';
      backBtn.onmouseout  = () => backBtn.style.background = 'rgba(8,18,50,.88)';
      backBtn.onclick = () => window.close();
      ui.appendChild(backBtn);
    }

    const nav = document.createElement('div');
    nav.style.cssText = 'position:absolute;bottom:16px;left:50%;transform:translateX(-50%);display:flex;gap:9px;align-items:center;pointer-events:all;';
    ui.appendChild(nav);

    function mkBtn(txt: string, fn: () => void) {
      const b = document.createElement('button');
      b.textContent = txt;
      b.style.cssText = 'background:rgba(8,18,50,.88);border:1px solid rgba(70,140,255,.4);color:rgba(150,200,255,.9);border-radius:18px;padding:6px 15px;cursor:pointer;font-size:11px;font-family:Courier New,monospace;letter-spacing:1px;pointer-events:all;transition:background .2s;';
      b.onmouseover = () => b.style.background = 'rgba(25,65,155,.9)';
      b.onmouseout  = () => b.style.background = 'rgba(8,18,50,.88)';
      b.onclick = fn; return b;
    }

    const dockBtn = mkBtn('⬡ DOCK MOLECULE', triggerDock);
    const prevBtn = mkBtn('◀', () => { if (SCENE.length) switchTo((cur - 1 + SCENE.length) % SCENE.length); });
    const nextBtn = mkBtn('▶', () => { if (SCENE.length) switchTo((cur + 1) % SCENE.length); });

    // One marker per detected variant. Solid ring = built-in 3D model (click →
    // animated scene). Dashed ring = real-PDB / AlphaFold fallback (click →
    // load the crystallographic structure for that gene).
    let sceneCounter = 0;
    const dots = view.entries.map(entry => {
      const d = document.createElement('button');
      if (entry.kind === 'full') {
        const sceneIndex = sceneCounter++;
        const c = '#' + MUT[entry.modelIndex].color.toString(16).padStart(6, '0');
        d.style.cssText = `width:9px;height:9px;border-radius:50%;border:2px solid ${c};background:transparent;cursor:pointer;pointer-events:all;padding:0;transition:all .22s;`;
        d.title = `${entry.meta.gene} ${entry.meta.mutation || ''}`.trim();
        d.onclick = () => switchTo(sceneIndex);
        return { el: d, sceneIndex, color: c };
      }
      const c = entry.meta.color || '#7fc7ff';
      d.style.cssText = `width:9px;height:9px;border-radius:50%;border:2px dashed ${c};background:transparent;cursor:pointer;pointer-events:all;padding:0;transition:all .22s;`;
      d.title = `${entry.meta.gene} ${entry.meta.mutation || ''} — real PDB structure`.trim();
      d.onclick = () => openPdbFallback(entry.meta);
      return { el: d, sceneIndex: -1, color: c };
    });
    // 6th dot — Loops 3D structural viewer
    const loopsBtn = document.createElement('button');
    loopsBtn.textContent = '🧬';
    loopsBtn.title = 'TP53 L1·L2·L3 Loops 3D';
    loopsBtn.style.cssText = 'background:rgba(0,180,180,.12);border:1.5px solid rgba(0,220,220,.5);color:rgba(0,240,220,.9);border-radius:18px;padding:5px 11px;cursor:pointer;font-size:13px;pointer-events:all;transition:background .2s;letter-spacing:1px;';
    loopsBtn.onmouseover = () => loopsBtn.style.background = 'rgba(0,180,180,.28)';
    loopsBtn.onmouseout  = () => loopsBtn.style.background = 'rgba(0,180,180,.12)';
    loopsBtn.onclick = () => setShowLoops(true);

    // PDB crystallographic viewer button
    const pdbBtn = document.createElement('button');
    pdbBtn.textContent = '🔬';
    pdbBtn.title = 'View PDB Crystal Structure';
    pdbBtn.style.cssText = 'background:rgba(80,160,255,.12);border:1.5px solid rgba(100,180,255,.5);color:rgba(140,210,255,.9);border-radius:18px;padding:5px 11px;cursor:pointer;font-size:13px;pointer-events:all;transition:background .2s;letter-spacing:1px;';
    pdbBtn.onmouseover = () => { pdbBtn.style.background = 'rgba(80,160,255,.28)'; };
    pdbBtn.onmouseout  = () => { pdbBtn.style.background = 'rgba(80,160,255,.12)'; };
    pdbBtn.onclick = () => {
      if (!SCENE.length) return;
      const m = SCENE[cur].model;
      const p = SCENE[cur].pdb;
      setPdbMut({
        id: SCENE[cur].meta.gene || m.id, variant: SCENE[cur].meta.mutation || m.variant,
        pdb: p.pdb, chain: p.chain, highlightRes: [...p.highlightRes],
        color: m.color, drug: m.drug, phase: m.phase,
      });
    };
    nav.append(prevBtn, ...dots.map(d => d.el), nextBtn, dockBtn, loopsBtn, pdbBtn);

    // Open the real crystallographic / AlphaFold structure for a variant that
    // has no built-in 3D model (resolved live via the platform backend).
    async function openPdbFallback(meta: BridgeVariant) {
      const gene = (meta.gene || '').toUpperCase();
      statusBox.textContent = `◌ Resolving ${gene} structure…`;
      statusBox.style.color = 'rgba(255,200,80,.9)';
      const s = await resolveStructure(gene);
      if (!s) {
        statusBox.textContent = `✗ No structure found for ${gene}`;
        statusBox.style.color = 'rgba(255,120,120,.9)';
        return;
      }
      statusBox.textContent = `🔬 ${gene} · ${s.pdb}`;
      statusBox.style.color = 'rgba(110,200,255,.9)';
      setPdbMut({
        id: gene,
        variant: meta.mutation || '',
        pdb: s.pdb, chain: s.chain, highlightRes: [],
        color: cssHexToInt(meta.color),
        drug: meta.source === 'research' ? 'Research literature target' : 'NGS-detected variant',
        phase: meta.tier ? `Tier ${meta.tier}` : '—',
        url: s.url,
      });
    }

    function refreshUI() {
      if (!SCENE.length) {
        infoBox.innerHTML = `
          <div style="color:#7fc7ff;font-size:15px;font-weight:bold;letter-spacing:2px;margin-bottom:8px;">No stylized model</div>
          <div style="color:rgba(200,220,255,.85);font-size:10px;line-height:1.7;">None of the detected variants has a built-in animated model. Use the dashed markers below to open the real PDB / AlphaFold structure for each gene.</div>`;
        dots.forEach(d => { d.el.style.transform = 'scale(1)'; if (d.sceneIndex >= 0) d.el.style.background = 'transparent'; });
        rightLabel.innerHTML = `${view.entries.length} VARIANT(S)<br><span style="opacity:.4;font-size:9px;">IBM sqDRIFT · NSCLC</span>`;
        return;
      }
      const m = SCENE[cur].model;
      const meta = SCENE[cur].meta;
      const cc = '#' + m.color.toString(16).padStart(6,'0');
      const dc = '#' + m.dCol.toString(16).padStart(6,'0');
      const pc = '#' + m.pCol.toString(16).padStart(6,'0');
      const isC275F = m.variant === 'C275F';
      const bqpBox = isC275F ? `
        <div style="margin-top:12px;padding:9px 11px;background:rgba(255,221,34,.06);border:1px solid rgba(255,221,34,.25);border-radius:7px;">
          <div style="color:#ffdd22;font-size:8.5px;letter-spacing:2px;margin-bottom:5px;">◈ BQP QUANTUM TARGET</div>
          <div style="color:rgba(255,240,180,.92);font-size:9px;line-height:1.7;">
            Phe275 π-system: 6e / 6 orbitals<br>
            Full pocket active space: ~44e · ~44 orbitals<br>
            Jordan-Wigner encoding: <span style="color:#ffdd22;font-weight:bold;">~88 qubits</span><br>
            Status: <span style="color:#44ffee;">within IBM Heron r3 range</span>
          </div>
        </div>
        <div style="margin-top:8px;padding:5px 8px;background:rgba(255,120,0,.06);border-left:2px solid rgba(255,119,0,.5);border-radius:0 4px 4px 0;">
          <div style="color:rgba(255,180,80,.9);font-size:8px;letter-spacing:1px;">STRUCTURAL MUTANT · No approved drug · No C275F-specific trial · Quantum simulation frontier</div>
        </div>` : '';
      infoBox.innerHTML = `
        <div style="color:${cc};font-size:20px;font-weight:bold;letter-spacing:3px;margin-bottom:2px;">${meta.gene || m.id}</div>
        <div style="color:rgba(200,220,255,.85);font-size:9.5px;letter-spacing:2px;margin-bottom:14px;">MUTATION · ${meta.mutation || m.variant}</div>
        <div style="color:${dc};font-size:13px;font-weight:bold;margin-bottom:3px;">⬡ ${m.drug}</div>
        <div style="color:rgba(190,215,255,.9);font-size:9.5px;margin-bottom:3px;">${m.sub}</div>
        <div style="color:rgba(180,210,255,.9);font-size:9.5px;margin-bottom:14px;">${m.phase}</div>
        <div style="border-top:1px solid rgba(70,150,255,.15);padding-top:12px;">
          <div style="color:${pc};font-size:9px;letter-spacing:2px;margin-bottom:5px;">● BINDING MECHANISM</div>
          <div style="color:rgba(220,235,255,.95);font-size:10.5px;line-height:1.7;">${m.mech}</div>
        </div>
        ${bqpBox}`;
      dots.forEach(d => {
        const active = d.sceneIndex === cur;
        if (d.sceneIndex >= 0) d.el.style.background = active ? d.color : 'transparent';
        d.el.style.transform = active ? 'scale(1.5)' : 'scale(1)';
      });
      const labels = ['● POCKET ACTIVE', '● LIGAND APPROACHING →', '✓ DOCKED  |  sqDRIFT ACTIVE'];
      statusBox.textContent = labels[Math.min(animState, 2)];
      statusBox.style.color = animState >= 2 ? 'rgba(80,255,160,.9)'
                            : animState === 1 ? 'rgba(255,200,80,.9)'
                            : 'rgba(110,200,255,.75)';
      rightLabel.innerHTML = `MUTATION ${cur+1} / ${SCENE.length}<br><span style="opacity:.4;font-size:9px;">IBM sqDRIFT · NSCLC</span>`;
    }

    function resetDrug(idx: number) {
      const dg = drugGrps[idx], pc = pocketCenters[idx];
      dg.position.set(pc.x+5.5, pc.y+3.5, pc.z+2.0);
      (dg as any).startPos = dg.position.clone();
      dg.rotation.set(0,0,0);
      // Hide π-ring if present
      const piRing = (roots[idx] as any).piRing;
      if (piRing) piRing.visible = false;
      // hide drug atoms/bonds
      dg.children.forEach((c: any) => {
        if (c.isGroup) { // wavefunction rings
          c.visible = false;
          c.children.forEach((l: any) => { if (l.material) l.material.opacity = 0; });
        } else if (c.material) c.material.opacity = 0;
      });
      // reset heatmap
      heatmapMats[idx].forEach(mat => {
        mat.color.set(new THREE.Color(SCENE[idx].model.pCol));
        mat.emissive.set(new THREE.Color(SCENE[idx].model.pCol));
        mat.emissiveIntensity = 0.6;
      });
      // hide electron cloud
      (electronClouds[idx].material as THREE.PointsMaterial).opacity = 0;
    }

    function switchTo(idx: number) {
      if (idx === cur || !roots[idx]) return;
      if (roots[cur]) roots[cur].visible = false;
      cur = idx; animState = 0; animT = 0;
      dockBtn.textContent = '⬡ DOCK MOLECULE';
      resetDrug(idx); roots[idx].visible = true; refreshUI();
    }

    function triggerDock() {
      if (!SCENE.length) return;
      if (animState >= 2) {
        animState = 0; animT = 0;
        dockBtn.textContent = '⬡ DOCK MOLECULE';
        resetDrug(cur);
      } else { animState = 1; }
      refreshUI();
    }

    const eioq = (t: number) => t < .5 ? 8*t*t*t*t : 1 - Math.pow(-2*t+2, 4)/2;

    // ── Heatmap colors ────────────────────────────────────────────────────────
    const HEAT = [0x0000ff, 0x0088ff, 0x00ffaa, 0xaaff00, 0xffff00, 0xff8800, 0xff2200];
    function heatColor(t: number) {
      const idx = Math.min(Math.floor(t * (HEAT.length - 1)), HEAT.length - 2);
      const f = t * (HEAT.length - 1) - idx;
      return new THREE.Color(HEAT[idx]).lerp(new THREE.Color(HEAT[idx + 1]), f);
    }

    let rafId: number;
    function animate() {
      rafId = requestAnimationFrame(animate);
      T += 0.011;

      const root = roots[cur];
      if (!root) {
        // No built-in 3D scene for this report (all variants use PDB fallback).
        // Keep the camera/lights alive; structures open in the PDB overlay.
        L1.position.set(Math.sin(T * .22) * 8, Math.cos(T * .18) * 7, 6);
        L2.position.set(Math.cos(T * .19) * 7, -4, Math.sin(T * .27) * 6);
        controls.update();
        composer.render();
        return;
      }
      const pkGrp = pocketGrps[cur];
      const dgGrp = drugGrps[cur];
      const cloud = electronClouds[cur];
      const wfGrp = wavefunctionGrps[cur];
      const hmMats = heatmapMats[cur];

      // Slow protein rotation
      root.rotation.y = T * 0.12;
      root.rotation.x = Math.sin(T * 0.08) * 0.05;

      // ── Pocket pulse ──────────────────────────────────────────────────────
      const pulse = 0.5 + 0.5 * Math.sin(T * 2.8);
      if (animState < 2) {
        pkGrp.children.forEach((c: any) => { if (c.material) c.material.emissiveIntensity = 0.3 + 0.45 * pulse; });
      }

      // ── Electron density cloud animation ──────────────────────────────────
      const cloudMat = cloud.material as THREE.PointsMaterial;
      if (animState >= 1) {
        const targetOpacity = animState >= 2 ? 0.55 : Math.min((animT) * 1.5, 0.45);
        cloudMat.opacity += (targetOpacity - cloudMat.opacity) * 0.06;
        // Animate individual particle positions
        const posAttr = cloud.geometry.getAttribute('position') as THREE.BufferAttribute;
        const base = (cloud as any).basePositions as Float32Array;
        const phases = (cloud as any).phases as Float32Array;
        for (let j = 0; j < phases.length; j++) {
          const jitter = Math.sin(T * 2.2 + phases[j]) * 0.04;
          posAttr.setXYZ(j, base[j*3] + jitter, base[j*3+1] + jitter * 0.7, base[j*3+2] + jitter * 0.9);
        }
        posAttr.needsUpdate = true;
      } else {
        cloudMat.opacity *= 0.93;
      }

      // ── Drug approach animation ───────────────────────────────────────────
      if (animState >= 1) {
        animT = Math.min(animT + 0.005, 1.0);
        const et = eioq(animT);
        dgGrp.position.lerpVectors((dgGrp as any).startPos, (dgGrp as any).targetPos, et);
        if (animState === 1) {
          dgGrp.rotation.y = T * 1.6 * (1 - animT);
          dgGrp.rotation.x = T * 1.1 * (1 - animT);
        }
        const opac = Math.min(animT * 2.5, 1.0);
        dgGrp.children.forEach((c: any) => { if (!c.isGroup && c.material) c.material.opacity = opac; });

        // ── Heatmap activation as drug approaches ─────────────────────────
        hmMats.forEach((mat, j) => {
          const phase = j / hmMats.length;
          const t = Math.max(0, Math.min((animT - phase * 0.3) * 2.5, 1.0));
          if (t > 0) {
            const hc = heatColor(t);
            mat.color.set(hc); mat.emissive.set(hc);
            mat.emissiveIntensity = 0.5 + t * 1.2;
          }
        });

        // ── Docked ────────────────────────────────────────────────────────
        if (animT >= 1.0 && animState === 1) {
          animState = 2;
          dockBtn.textContent = '↺ RESET';
          bloom.strength = 2.5;
          setTimeout(() => { bloom.strength = 1.1; }, 700);
          // Show wavefunction
          wfGrp.visible = true;
          wfGrp.children.forEach((l: any) => { l.material.opacity = 0.5; });
          // Show Phe275 π-ring for C275F
          const piRing = (roots[cur] as any).piRing;
          if (piRing) piRing.visible = true;
          refreshUI();
        }
      }

      // ── Docked idle animations ────────────────────────────────────────────
      // ── C275F π-ring animation ────────────────────────────────────────────
      const piRing = (roots[cur] as any).piRing;
      if (piRing && animState >= 2) {
        piRing.rotation.y = T * 0.7;
        piRing.children.forEach((l: any, idx: number) => {
          if (l.material) l.material.opacity = idx < 2
            ? (idx === 0 ? 0.7 + 0.3 * Math.sin(T * 3.2) : 0.35 + 0.2 * Math.sin(T * 2.8 + 1))
            : 0.18 + 0.15 * Math.sin(T * 4 + idx);
        });
      } else if (piRing) {
        piRing.visible = false;
      }

      if (animState >= 2) {
        dgGrp.position.y += Math.sin(T * 1.8) * 0.0008;
        dgGrp.rotation.set(0, 0, 0);

        // Wavefunction rings spin + pulse
        wfGrp.visible = true;
        wfGrp.rotation.y = T * 0.9;
        wfGrp.rotation.x = Math.sin(T * 0.5) * 0.3;
        wfGrp.children.forEach((ring: any, idx: number) => {
          ring.material.opacity = 0.25 + 0.25 * Math.sin(T * 2.5 + idx * 1.1);
          ring.scale.setScalar((0.35 + idx * 0.22) * (1 + 0.06 * Math.sin(T * 3 + idx)));
        });

        // Pocket heatmap breathing
        hmMats.forEach((mat, j) => {
          mat.emissiveIntensity = 1.2 + 0.4 * Math.sin(T * 2.4 + j * 0.4);
        });
      }

      // Moving lights
      L1.position.set(Math.sin(T * .22) * 8, Math.cos(T * .18) * 7, 6);
      L2.position.set(Math.cos(T * .19) * 7, -4, Math.sin(T * .27) * 6);

      controls.update();
      composer.render();
    }

    refreshUI(); animate();

    const onResize = () => {
      const w = el.clientWidth, h = el.clientHeight;
      camera.aspect = w / h; camera.updateProjectionMatrix();
      renderer.setSize(w, h); composer.setSize(w, h);
    };
    window.addEventListener('resize', onResize);

    return () => {
      cancelAnimationFrame(rafId);
      window.removeEventListener('resize', onResize);
      controls.dispose(); renderer.dispose();
      if (renderer.domElement.parentElement) renderer.domElement.parentElement.removeChild(renderer.domElement);
      if (ui.parentElement) ui.parentElement.removeChild(ui);
    };
  }, []);

  return (
    <div style={{ width: '100%', height: '100vh', position: 'relative' }}>
      <div ref={mountRef} style={{ width: '100%', height: '100%', background: '#020610' }} />
      {showLoops && (
        <div style={{ position: 'absolute', inset: 0, zIndex: 50 }}>
          <TP53LoopsViewer onBack={() => setShowLoops(false)} />
        </div>
      )}
      {pdbMut && (
        <PDBMolViewer mutation={pdbMut} onBack={() => setPdbMut(null)} />
      )}
    </div>
  );
}
