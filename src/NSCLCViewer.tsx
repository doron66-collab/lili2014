import { useEffect, useState } from 'react';
import PDBMolViewer from './PDBMolViewer';

// ── Mutation info — text only, no fabricated geometry ────────────────────────
// Previously this file also built a hand-crafted, entirely fictional 3D model
// per mutation (procedurally generated "pocket" spheres, invented atom
// coordinates for a fake drug molecule, a scripted docking animation) with no
// basis in any real structural or computational data. Removed 2026-09-09: it
// duplicated the real PDB crystal-structure viewer below for every single one
// of these five targets (all five already have a real PDB entry in PDB_MAP),
// and a docking animation implying a real drug binds at a site is actively
// misleading for several of these where the note itself says no
// mutation-specific therapy exists. What's kept is exactly the reproducibility
// info: drug name, trial phase, and mechanism text — attached to the real
// structure instead of a fabricated one.
const MUT_INFO = [
  {
    id: 'TP53', variant: 'Y220C', drug: 'Rezatapopt', sub: 'PC14586',
    phase: 'Phase II — PYNNACLE', color: 0x3399ff,
    mech: 'Binds Y220C cryptic pocket → restores wild-type p53 conformation → reactivates tumor suppression',
  },
  {
    id: 'KEAP1', variant: 'LOF', drug: 'VVD-065', sub: 'Covalent NRF2 inhibitor',
    phase: 'Phase I — NCT05954312', color: 0x33ffaa,
    mech: 'Covalent NRF2 inhibition in KEAP1-deficient cells → restores chemosensitivity',
  },
  {
    id: 'CDKN2A', variant: 'p16 loss', drug: 'Palbociclib + Olaparib', sub: 'CDK4/6 + PARP synergy',
    phase: 'Phase I/II — PARP1 degradation 2024', color: 0xffaa33,
    mech: 'CDK4/6 inhibition restores G1 checkpoint lost by p16 deletion + PARP1 synthetic lethality',
  },
  {
    id: 'STK11', variant: 'LKB1 loss', drug: 'Ceralasertib', sub: '+ Durvalumab (ATR inhibitor)',
    phase: 'Phase III — LATIFY NCT05450692', color: 0xff3366,
    mech: 'ATR inhibition targets replication stress created by LKB1 loss (sqDRIFT-modellable synthetic lethality)',
  },
  {
    // TP53 C275F — structural β-sandwich core mutant (p.Cys275Phe)
    // Detected: Emek Medical Center MI25-0349 · June 2025
    // sqDRIFT active space: 18–48 active electrons (AVAS-tested, criterion-dependent) · 24–56 qubits
    id: 'TP53', variant: 'C275F', drug: 'Eprenetapopt', sub: 'APR-246 · Pan-mutant p53 reactivator',
    phase: 'No C275F-specific trial · Research stage', color: 0xffdd22,
    mech: 'p.Cys275Phe replaces core β-sandwich Cys with bulky Phe → hydrophobic core collapse → global misfolding. No mutation-specific drug exists. Quantum simulation (sqDRIFT) of the π-electron pocket is the frontier approach. Active space: 18–48e · 24–56 qubits',
  },
] as const;

// PDB crystallographic data per mutation (parallel to MUT_INFO)
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

// simulation_id → index into MUT_INFO / PDB_MAP
const MODEL_INDEX: Record<string, number> = {
  TP53_Y220C: 0,
  KEAP1_LOF:  1,
  CDKN2A_P16: 2,
  STK11_LKB1: 3,
  TP53_C275F: 4,
};
// Genes with exactly ONE known entry — allow a gene-level match when the
// simulation_id doesn't line up exactly (e.g. an uploaded KEAP1 frameshift or a
// CDKN2A deletion still maps to the single KEAP1 / CDKN2A entry). TP53 is
// deliberately excluded: it has two entries, so an unmatched TP53 variant falls
// through to the real-PDB resolution path rather than guessing.
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
  | { kind: 'known'; modelIndex: number; meta: BridgeVariant }
  | { kind: 'pdb';   meta: BridgeVariant };

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
            ? { kind: 'known', modelIndex: mi, meta: v }
            : { kind: 'pdb', meta: v };
        });
        return { entries, reportId: data.report_id || null };
      }
    }
  } catch {
    // Corrupt/blocked storage — fall through to the demo set.
  }
  // Fallback: the original five demo entries.
  const simIds = Object.keys(MODEL_INDEX); // insertion order matches MUT_INFO indices
  const entries: ViewEntry[] = MUT_INFO.map((m, i) => ({
    kind: 'known',
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
  // Render's free/starter tier can take 30-60s+ to wake a sleeping backend from
  // cold, and this fetch previously had NO client-side timeout at all - if the
  // backend never answers (asleep, unreachable, DNS/CORS issue), the caller's
  // "Resolving <gene> structure…" status just sat there forever with no
  // fallback message, indistinguishable from the viewer being stuck. Found
  // live 2026-09-09 (ARID2, an AlphaFold-only entry with no RCSB PDB, so it
  // always takes this fallback path). 25s covers a real cold-start without
  // making a genuinely offline backend look hung for longer than necessary.
  const controller = new AbortController();
  const timeout = setTimeout(() => controller.abort(), 25000);
  try {
    const r = await fetch(`${base}/api/pdb/lookup/gene/${encodeURIComponent(gene)}`, { signal: controller.signal })
      .then(res => res.json());
    if (r?.pdb_ids?.length) return { pdb: r.pdb_ids[0], chain: 'A' };
    if (r?.alphafold_model_url) return { pdb: `AF-${r.uniprot_id || gene}`, chain: 'A', url: r.alphafold_model_url };
  } catch {
    // Backend asleep / offline / timed out — caller surfaces a message.
  } finally {
    clearTimeout(timeout);
  }
  return null;
}

function cssHexToInt(hex?: string): number {
  if (!hex) return 0x06b6d4;
  const n = parseInt(hex.replace('#', ''), 16);
  return Number.isFinite(n) ? n : 0x06b6d4;
}

interface PdbMutInfo {
  id: string; variant: string; pdb: string; chain: string;
  highlightRes: number[]; color: number; drug: string; phase: string;
  sub?: string; mech?: string;
  url?: string;   // optional explicit structure URL (e.g. AlphaFold model)
}

export default function NSCLCViewer() {
  const [view] = useState(buildViewList);
  const [cur, setCur] = useState(0);
  const [resolved, setResolved] = useState<PdbMutInfo | null>(null);
  const [status, setStatus] = useState<string | null>(null);

  const entry = view.entries[cur];

  useEffect(() => {
    if (!entry) { setResolved(null); return; }
    let cancelled = false;

    if (entry.kind === 'known') {
      const m = MUT_INFO[entry.modelIndex];
      const p = PDB_MAP[entry.modelIndex];
      setStatus(null);
      setResolved({
        id: entry.meta.gene || m.id, variant: entry.meta.mutation || m.variant,
        pdb: p.pdb, chain: p.chain, highlightRes: [...p.highlightRes],
        color: m.color, drug: m.drug, phase: m.phase, sub: m.sub, mech: m.mech,
      });
      return;
    }

    // 'pdb' entry — no known model, resolve a real structure live.
    const gene = (entry.meta.gene || '').toUpperCase();
    setResolved(null);
    setStatus(`◌ Resolving ${gene} structure…`);
    resolveStructure(gene).then(s => {
      if (cancelled) return;
      if (!s) {
        setStatus(`✗ No structure found for ${gene}`);
        return;
      }
      setStatus(null);
      setResolved({
        id: gene, variant: entry.meta.mutation || '',
        pdb: s.pdb, chain: s.chain, highlightRes: [],
        color: cssHexToInt(entry.meta.color),
        drug: entry.meta.source === 'research' ? 'Research literature target' : 'NGS-detected variant',
        phase: entry.meta.tier ? `Tier ${entry.meta.tier}` : '—',
        url: s.url,
      });
    });

    return () => { cancelled = true; };
  }, [cur, entry]);

  function closeOrBack() {
    if (typeof window !== 'undefined' && window.opener && !window.opener.closed) window.close();
  }

  if (!view.entries.length) {
    return (
      <div style={{
        width: '100%', height: '100vh', display: 'flex', alignItems: 'center', justifyContent: 'center',
        background: '#020610', color: 'rgba(200,220,255,.85)', fontFamily: 'Courier New, monospace',
      }}>
        No variants to display.
      </div>
    );
  }

  if (!resolved) {
    return (
      <div style={{
        width: '100%', height: '100vh', display: 'flex', alignItems: 'center', justifyContent: 'center',
        background: '#020610', color: 'rgba(200,220,255,.85)', fontFamily: 'Courier New, monospace', fontSize: 13,
      }}>
        {status || 'Loading…'}
      </div>
    );
  }

  return (
    <PDBMolViewer
      key={`${cur}-${resolved.pdb}`}
      mutation={resolved}
      onBack={closeOrBack}
      onPrev={view.entries.length > 1 ? () => setCur((cur - 1 + view.entries.length) % view.entries.length) : undefined}
      onNext={view.entries.length > 1 ? () => setCur((cur + 1) % view.entries.length) : undefined}
      navPosition={view.entries.length > 1 ? `${cur + 1} / ${view.entries.length}` : undefined}
    />
  );
}
