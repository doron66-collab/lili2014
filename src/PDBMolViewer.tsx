import { useEffect, useRef, useState } from 'react';
import * as NGL from 'ngl';

interface MutInfo {
  id: string;
  variant: string;
  pdb: string;
  chain: string;
  highlightRes?: number[];
  color: number;
  drug: string;
  phase: string;
  sub?: string;    // drug detail line, e.g. "APR-246 · Pan-mutant p53 reactivator"
  mech?: string;   // mechanism-of-action text
  url?: string;   // explicit structure URL (e.g. AlphaFold); defaults to RCSB
}

interface Props {
  mutation: MutInfo;
  onBack: () => void;
  onPrev?: () => void;
  onNext?: () => void;
  navPosition?: string;   // e.g. "2 / 5" — shown next to the prev/next controls
}

// ── Pocket detection (fpocket, backend/routes/pocket.py) ─────────────────────
// The SECOND, independent axis of target assessment alongside SOLANGE's
// existing electronic-structure/DMRG classification: does a 3D cavity even
// exist here for a drug-like ligand, regardless of Class A/B/C. Triggered by
// the user's own choice of structure — not a fixed list — per Doron's own
// framing: "משתמש קצה הגיע לתלת מימד ... צריך להוסיף שם כפתור מצא כיסים".
interface PocketResult {
  pocket_id: number;
  score?: number;
  druggability_score?: number;
  volume?: number;
  x?: number; y?: number; z?: number; radius?: number;
  chains: string[];
  single_chain: boolean;
  includes_target_residue?: boolean;
}
interface PocketResponse {
  pdb: string;
  near_residue?: string | null;
  n_pockets_total: number;
  n_single_chain_druggable: number;
  best_single_chain_pocket: PocketResult | null;
  all_pockets: PocketResult[];
}

function apiBase(): string {
  return (typeof window !== 'undefined' && (window as any).QCAIHPC_API_BASE)
    || 'https://qcaihpc-simulation-api.onrender.com';
}

// Map hex color → NGL color string
function hexToNGLColor(hex: number): string {
  return '#' + hex.toString(16).padStart(6, '0');
}

interface PdbMeta {
  title: string;
  method?: string;
  resolution?: number;
}

// Real, sourced structure data from RCSB's own data API — the one piece of
// substantive information available for ANY real PDB entry, not just the
// five curated demo mutations. Curated mech/sub/drug text only exists for
// those five (it's chemist-authored, per-mutation content); a gene pulled
// through the generic search box has none of that, so without this the
// manual-search path showed a bare rotating structure with no data at all
// (reported live 2026-09-09). Deliberately factual and RCSB-sourced only —
// no per-mutation mechanism is invented for an arbitrary structure.
async function fetchPdbMeta(pdbId: string): Promise<PdbMeta | null> {
  const controller = new AbortController();
  const timeout = setTimeout(() => controller.abort(), 10000);
  try {
    const res = await fetch(`https://data.rcsb.org/rest/v1/core/entry/${pdbId.toUpperCase()}`, {
      signal: controller.signal,
    });
    if (!res.ok) return null;
    const r = await res.json();
    const title = r?.struct?.title;
    if (!title) return null;
    return {
      title,
      method: r?.exptl?.[0]?.method,
      resolution: r?.rcsb_entry_info?.resolution_combined?.[0],
    };
  } catch {
    return null;
  } finally {
    clearTimeout(timeout);
  }
}

// Whole-protein display style — user-selectable (Doron: "כדורים, קווים,
// קפסולות טיקטק" / balls, wires, tic-tac capsules — and explicitly for the
// WHOLE protein, not just the mutation residue, after the first version
// only let this control the mutation site). English-only labels (repeated
// feedback: no Hebrew in UI controls). Per-style NGL representation params
// that otherwise look visually inconsistent at one shared default (e.g.
// spacefill needs a much smaller radiusScale than licorice to read as
// atoms rather than a solid blob over an entire chain).
const STRUCT_STYLES = [
  { key: 'cartoon',   label: 'RIBBON',   repr: 'cartoon',   params: { colorScheme: 'residueindex', smoothSheet: true, opacity: 0.92 } },
  { key: 'licorice',  label: 'WIRE',     repr: 'licorice',  params: { colorScheme: 'element', opacity: 0.9, radiusScale: 0.35 } },
  { key: 'spacefill', label: 'BALLS',    repr: 'spacefill', params: { colorScheme: 'element', opacity: 0.9, radiusScale: 0.5 } },
  { key: 'hyperball', label: 'CAPSULES', repr: 'hyperball', params: { colorScheme: 'element', opacity: 0.9, radiusScale: 0.4, shrink: 0.3 } },
] as const;
type StructStyleKey = typeof STRUCT_STYLES[number]['key'];

export default function PDBMolViewer({ mutation, onBack, onPrev, onNext, navPosition }: Props) {
  const mountRef = useRef<HTMLDivElement>(null);
  const stageRef = useRef<any>(null);
  const structComponentRef = useRef<any>(null);
  const baseReprRef = useRef<any>(null);
  const otherChainsReprRef = useRef<any>(null);
  const mutSiteReprRef = useRef<any>(null);
  const pocketShapeRef = useRef<any>(null);
  const [structReady, setStructReady] = useState(0);
  const [structStyle, setStructStyle] = useState<StructStyleKey>('cartoon');
  const [spinning, setSpinning] = useState(true);
  const [pdbMeta, setPdbMeta] = useState<PdbMeta | null>(null);
  const [chainInfo, setChainInfo] = useState<{ chainname: string; description: string }[]>([]);
  const [hoverInfo, setHoverInfo] = useState<{ chain: string; resno: number; resname: string; x: number; y: number } | null>(null);
  const [pocketStatus, setPocketStatus] = useState<string | null>(null);
  const [pocketResult, setPocketResult] = useState<PocketResponse | null>(null);
  const [selectedPocketId, setSelectedPocketId] = useState<number | null>(null);
  const [showAllPockets, setShowAllPockets] = useState(false);
  const isAlphaFold = mutation.pdb.startsWith('AF-');

  // Reset pocket results whenever the user navigates to a different structure
  // — a stale pocket overlay from the previous mutation rendered on the new
  // one would misattribute a cavity to the wrong target.
  useEffect(() => {
    setPocketResult(null);
    setPocketStatus(null);
    setSelectedPocketId(null);
    setShowAllPockets(false);
  }, [mutation.pdb]);

  // fpocket reports every geometric cavity on the WHOLE protein surface —
  // for TP53 2OCJ that was 59 of them. Rendering all of them at once (what
  // the first version did) reads as meaningless clutter, not information:
  // reported live by Doron ("סתם שרטוט תלת מימדי שמסתובב - לא רלוונטי").
  // Default to the handful that actually matter: the best druggable
  // candidate, anything else single-chain that clears the bar, and —
  // separately — anything (even multi-chain or low-score) that happens to
  // sit AT the mutation residue, since that's the one question a user
  // looking at a specific mutation actually has. "Show all" is still one
  // click away for anyone who wants the raw landscape.
  const allPockets = pocketResult?.all_pockets || [];
  const withCoords = allPockets.filter(p => p.x != null && p.y != null && p.z != null);
  const atMutationSite = withCoords.filter(p => p.includes_target_residue);
  const topCandidates = withCoords
    .filter(p => p.single_chain)
    .sort((a, b) => (b.score ?? -999) - (a.score ?? -999))
    .slice(0, 8);
  const shownPockets = showAllPockets
    ? withCoords
    : Array.from(new Map([...atMutationSite, ...topCandidates].map(p => [p.pocket_id, p])).values());
  const bestId = pocketResult?.best_single_chain_pocket?.pocket_id;

  function goToMutation() {
    const component = structComponentRef.current;
    if (!component || !mutation.highlightRes || mutation.highlightRes.length === 0) return;
    const sele = mutation.highlightRes.map(r => `${r}:${mutation.chain}`).join(' or ');
    component.autoView(sele, 700);
  }

  function focusPocket(p: PocketResult) {
    setSelectedPocketId(p.pocket_id);
    const stage = stageRef.current;
    if (stage && p.x != null && p.y != null && p.z != null) {
      stage.animationControls.move([p.x, p.y, p.z], 700);
      stage.animationControls.zoom(Math.max((p.radius || 5) * 3.2, 18), 700);
    }
  }

  async function findPockets() {
    if (isAlphaFold) {
      setPocketStatus('✗ fpocket needs a real RCSB crystal structure — not available for AlphaFold models');
      return;
    }
    setPocketStatus('◌ Running fpocket…');
    setPocketResult(null);
    const nearResidue = mutation.highlightRes && mutation.highlightRes.length > 0
      ? `${mutation.chain}:${mutation.highlightRes[0]}`
      : undefined;
    const controller = new AbortController();
    const timeout = setTimeout(() => controller.abort(), 130000); // fpocket itself is capped at 120s server-side
    try {
      const params = new URLSearchParams({ pdb: mutation.pdb, min_druggability: '0.5' });
      if (nearResidue) params.set('near_residue', nearResidue);
      const res = await fetch(`${apiBase()}/api/pocket/detect?${params}`, { signal: controller.signal });
      if (!res.ok) {
        const body = await res.json().catch(() => null);
        setPocketStatus(`✗ ${body?.detail || `pocket detection failed (HTTP ${res.status})`}`);
        return;
      }
      const data: PocketResponse = await res.json();
      setPocketResult(data);
      setPocketStatus(null);
    } catch (e: any) {
      setPocketStatus(e?.name === 'AbortError'
        ? '✗ pocket detection timed out'
        : '✗ could not reach pocket-detection backend (may be cold-starting — try again in ~30s)');
    } finally {
      clearTimeout(timeout);
    }
  }

  // Render pocket spheres as a separate NGL shape layer, independent of the
  // structure component above — added/removed on top of whatever structure
  // is already loaded, in the same PDB coordinate frame so they align
  // automatically. Colour encodes what the backend already decided (see
  // pocket.py / run_pocket_detect.py's own docstring on why the combined
  // score, not druggability alone, is the trustworthy signal): gold = the
  // best single-chain druggable candidate, cyan = other single-chain
  // candidates, dim red = multi-chain (crystal-contact artifact, not a real
  // binding site). A white halo marks any pocket that actually overlaps the
  // mutation residue — the one thing a user looking at a specific mutation
  // most wants to know at a glance, separate from which pocket ranks best
  // overall. The selected pocket (from the list panel) gets a brighter
  // pulsing-scale halo instead, as the click-to-focus feedback.
  useEffect(() => {
    const stage = stageRef.current;
    if (pocketShapeRef.current) {
      for (const c of pocketShapeRef.current) stage?.removeComponent(c);
      pocketShapeRef.current = null;
    }
    if (!stage || !shownPockets.length) return;

    // disableImpostor: raycasted screen-space impostor spheres (NGL's
    // default, faster to render) go dark/black once the camera crosses
    // into or very near the sphere's own volume — exactly what happens
    // when "focus" zooms the camera in on a small pocket, or the user
    // scrolls in close manually (reported live 2026-10-07: "הצבע הזהוב
    // שנעלם בפוזיציה מסויימת כאשר מגיעים אליה"). Real mesh geometry
    // doesn't have that failure mode, at the (here, negligible — a
    // handful of low-poly spheres) cost of being somewhat more expensive
    // to render than impostors.
    const shape = new (NGL as any).Shape('pockets', { disableImpostor: true });
    const haloShape = new (NGL as any).Shape('pocket-halos', { disableImpostor: true });
    for (const p of shownPockets) {
      if (p.x == null || p.y == null || p.z == null) continue;
      const color: [number, number, number] = !p.single_chain
        ? [0.6, 0.18, 0.18]
        : p.pocket_id === bestId
        ? [1.0, 0.82, 0.1]
        : [0.85, 0.25, 0.78]; // magenta — stays distinct from the cartoon's own rainbow (residueindex) coloring
      const radius = Math.max(p.radius || 4, 2.5);
      shape.addSphere([p.x, p.y, p.z], color, radius);
      if (p.pocket_id === selectedPocketId) {
        haloShape.addSphere([p.x, p.y, p.z], [1, 1, 1], radius * 1.6);
      } else if (p.includes_target_residue) {
        haloShape.addSphere([p.x, p.y, p.z], [1, 1, 1], radius * 1.35);
      }
    }
    // side:'double' + diffuseInterior:true: disableImpostor alone wasn't
    // enough (confirmed live — still went dark up close). Root cause is
    // standard mesh lighting: a sphere's front faces get culled once the
    // camera passes inside it, and the back faces it would show instead
    // have no normal-facing light hitting them, rendering solid black.
    // diffuseInterior makes NGL light the interior ignoring normal facing;
    // side:'double' ensures the interior faces are even drawn in the first
    // place. Verified both are real, accepted repr params (not guessed)
    // in a real headless-Chromium run before using them.
    const shapeComp = stage.addComponentFromObject(shape);
    shapeComp.addRepresentation('buffer', { opacity: 0.55, side: 'double', diffuseInterior: true });
    const haloComp = stage.addComponentFromObject(haloShape);
    haloComp.addRepresentation('buffer', { opacity: 0.16, side: 'double', diffuseInterior: true });
    pocketShapeRef.current = [shapeComp, haloComp];

    return () => {
      stage?.removeComponent(shapeComp);
      stage?.removeComponent(haloComp);
      pocketShapeRef.current = null;
    };
  }, [pocketResult, showAllPockets, selectedPocketId]);

  useEffect(() => {
    setPdbMeta(null);
    setChainInfo([]);
    if (isAlphaFold) return; // predicted model, not an RCSB experimental entry
    let cancelled = false;
    fetchPdbMeta(mutation.pdb).then(m => { if (!cancelled) setPdbMeta(m); });
    return () => { cancelled = true; };
  }, [mutation.pdb, isAlphaFold]);

  useEffect(() => {
    const el = mountRef.current;
    if (!el) return;
    // Defensive: NGL.Stage appends its own <canvas> into `el` and this effect
    // re-runs (new Stage) every time mutation.pdb changes, since PDBMolViewer
    // is now the ONLY view and stays mounted across navigation instead of
    // being conditionally unmounted. If the previous Stage's own dispose()
    // ever leaves its canvas behind for any reason, canvases would silently
    // accumulate here on every "next" click - each holding its own live
    // WebGL context, of which a browser allows only a small fixed number
    // before it starts failing/evicting contexts (surfacing as exactly what
    // was reported live 2026-09-09: navigation appearing stuck on the first
    // mutation, and the whole machine feeling sluggish after a few clicks).
    // Clearing the container explicitly before creating the new Stage makes
    // that impossible regardless of what dispose() itself guarantees.
    el.innerHTML = '';
    let cancelled = false;

    const stage = new NGL.Stage(el, {
      backgroundColor: '#020d1f',
      quality: 'high',
      antialias: true,
      impostor: true,
      tooltip: false, // NGL's default raw hover tooltip ("sphere: 113 (pockets)")
      // is meaningless to a non-technical viewer — the pocket list panel is
      // the real UI for inspecting a candidate (reported live 2026-10-07).
    });
    stageRef.current = stage;

    const pdbUrl = mutation.url || `https://files.rcsb.org/download/${mutation.pdb}.pdb`;

    stage.loadFile(pdbUrl, { ext: 'pdb', defaultRepresentation: false }).then((component: any) => {
      // `stage.dispose()` in this effect's cleanup does NOT reliably run
      // before a slow network load resolves — clicking prev/next quickly
      // (or a slow RCSB response racing a fast click) let this callback add
      // representations to, and call autoView on, a component whose stage
      // had already been torn down. That's consistent with what was reported
      // live 2026-09-09: the viewer intermittently not showing a mutation at
      // all, and going back not re-showing one that displayed fine moments
      // earlier — a disposed/half-torn-down stage rendering nothing rather
      // than throwing a visible error.
      if (cancelled) return;
      const ch = mutation.chain;
      structComponentRef.current = component;

      // Per-chain identity, straight from the loaded file's own entity
      // records (COMPND/entity description — no extra network call) —
      // answers "what are these chains" for any multi-chain complex, after
      // Doron asked this directly for 2WTK's LKB1/STRADalpha/MO25alpha trio.
      const seen = new Set<string>();
      const chains: { chainname: string; description: string }[] = [];
      component.structure.eachChain((cp: any) => {
        if (seen.has(cp.chainname)) return;
        seen.add(cp.chainname);
        chains.push({ chainname: cp.chainname, description: cp.entity?.description || '(no entity description in file)' });
      });
      chains.sort((a, b) => a.chainname.localeCompare(b.chainname));
      setChainInfo(chains);

      // The whole-chain base representation (ribbon / wire / balls /
      // capsules) is added by the dedicated effect below, keyed on
      // structStyle — not hardcoded here — so switching style swaps one
      // representation instead of reloading the whole PDB file.
      setStructReady(v => v + 1);

      if (mutation.highlightRes && mutation.highlightRes.length > 0) {
        // Pocket neighbourhood licorice — chain-qualified
        const pocketSele = mutation.highlightRes
          .flatMap(r => Array.from({ length: 11 }, (_, i) => r - 5 + i))
          .filter(r => r > 0)
          .map(r => `${r}:${ch}`)
          .join(' or ');
        component.addRepresentation('licorice', {
          sele: pocketSele,
          colorScheme: 'element',
          opacity: 0.85,
          radiusScale: 0.6,
        });
      }

      // Zinc ion — no chain filter needed (heteroatom)
      component.addRepresentation('spacefill', {
        sele: 'ZN',
        color: '#aaffdd',
        opacity: 1.0,
        radiusScale: 1.2,
      });

      component.autoView(800);
    }).catch((err: any) => {
      if (!cancelled) console.error('NGL load error:', err);
    });

    stage.setSpin([0, 1, 0], 0.006);

    // Custom hover tooltip (chain + residue) — NGL's own built-in tooltip is
    // disabled above because its raw form ("sphere: 113 (pockets)") was
    // meaningless for a pocket sphere; this replaces it with something
    // useful for the STRUCTURE itself: which chain/residue is under the
    // cursor. Answers "how do I know where the mutation is" by letting
    // Doron hover around and read off chain+residue directly, instead of
    // needing a pre-set highlight for every possible residue of interest.
    const onHover = (pickingProxy: any) => {
      if (pickingProxy && pickingProxy.atom) {
        const a = pickingProxy.atom;
        const pos = pickingProxy.mouse?.position;
        setHoverInfo({
          chain: a.chainname, resno: a.resno, resname: a.resname,
          x: pos?.x ?? 0, y: pos?.y ?? 0,
        });
      } else {
        setHoverInfo(null);
      }
    };
    stage.signals.hovered.add(onHover);

    const onResize = () => stage.handleResize();
    window.addEventListener('resize', onResize);

    return () => {
      cancelled = true;
      window.removeEventListener('resize', onResize);
      stage.signals.hovered.remove(onHover);
      setHoverInfo(null);
      stageRef.current = null;
      structComponentRef.current = null;
      baseReprRef.current = null;
      otherChainsReprRef.current = null;
      mutSiteReprRef.current = null;
      // stage.dispose() removes NGL's own bookkeeping and the canvas element,
      // but (confirmed against NGL's own dispose() source) never calls the
      // underlying THREE.WebGLRenderer's forceContextLoss() — the canvas can
      // be gone from the DOM while its WebGL context is still alive and
      // counted against the browser's small fixed per-page context limit.
      // Once enough quick navigations exhaust that limit, new Stages fail to
      // get a context and silently render nothing — again matching the
      // "works, then a mutation just doesn't show, and going back doesn't
      // bring back one that worked before" report. Force the real GPU
      // context to release, not just NGL's own cleanup.
      try {
        stage.viewer?.renderer?.forceContextLoss?.();
      } catch {
        // Best-effort — a missing/renamed internal shouldn't block teardown.
      }
      stage.dispose();
    };
  }, [mutation.pdb]);

  // Whole-chain base representation (ribbon / wire / balls / capsules) —
  // kept separate from the structure-loading effect above so switching
  // style just swaps this one representation instead of reloading the
  // whole PDB file. Applies to the ENTIRE primary chain, per Doron's
  // correction after the first version only let this control the mutation
  // residue ("התכוונתי לשינוי בכל החלבון").
  //
  // A second, dimmed representation draws every OTHER chain in the
  // deposited entry (e.g. 2WTK's STRADalpha/MO25alpha partners alongside
  // STK11/LKB1 itself) — added after Doron asked directly whether the
  // viewer could show all three chains of a multi-chain complex, and why
  // fpocket's own pockets on those other chains looked like they were
  // floating in empty space (they were real cavities on protein that
  // simply wasn't being drawn at all). Dimmed/desaturated rather than full
  // color so the primary chain (the actual mutation's own gene product)
  // stays visually unambiguous as the main subject.
  useEffect(() => {
    const component = structComponentRef.current;
    if (!component) return;
    if (baseReprRef.current) {
      component.removeRepresentation(baseReprRef.current);
      baseReprRef.current = null;
    }
    if (otherChainsReprRef.current) {
      component.removeRepresentation(otherChainsReprRef.current);
      otherChainsReprRef.current = null;
    }
    const style = STRUCT_STYLES.find(s => s.key === structStyle) || STRUCT_STYLES[0];
    baseReprRef.current = component.addRepresentation(style.repr, {
      sele: `:${mutation.chain}`,
      ...style.params,
    });
    otherChainsReprRef.current = component.addRepresentation(style.repr, {
      sele: `not :${mutation.chain}`,
      ...style.params,
      colorScheme: 'uniform',
      color: '#5a6a85',
      opacity: Math.min(style.params.opacity ?? 0.9, 0.9) * 0.4,
    });
    return () => {
      if (baseReprRef.current) {
        component.removeRepresentation(baseReprRef.current);
        baseReprRef.current = null;
      }
      if (otherChainsReprRef.current) {
        component.removeRepresentation(otherChainsReprRef.current);
        otherChainsReprRef.current = null;
      }
    };
  }, [structReady, structStyle, mutation.chain]);

  // Mutation-site highlight — a fixed, bold color (NOT the per-target accent
  // color, which can blend into the cartoon's own rainbow or a nearby pocket
  // sphere) so it never gets lost. Doron reported it "disappearing" when
  // zoomed in via GO TO MUTATION — two causes, both addressed: (1) the old
  // per-target color wasn't guaranteed to contrast against everything else
  // on screen, now a single unmistakable lime green used for nothing else
  // in this viewer; (2) licorice cylinders are impostor-rendered by default,
  // same underlying lighting issue the pocket spheres had up close (see that
  // fix's own comment) — side:'double' + diffuseInterior:true applied here
  // too. A spacefill ball is layered on top for visibility at any zoom
  // level, not just when close enough to resolve individual bonds.
  useEffect(() => {
    const component = structComponentRef.current;
    if (!component || !mutation.highlightRes || mutation.highlightRes.length === 0) return;
    if (mutSiteReprRef.current) {
      for (const r of mutSiteReprRef.current) component.removeRepresentation(r);
      mutSiteReprRef.current = null;
    }
    const sele = mutation.highlightRes.map(r => `${r}:${mutation.chain}`).join(' or ');
    const MUT_HIGHLIGHT_COLOR = '#39ff14'; // bold lime green — used nowhere else in this viewer
    const licoriceRepr = component.addRepresentation('licorice', {
      sele, color: MUT_HIGHLIGHT_COLOR, opacity: 1.0, radiusScale: 1.6,
      side: 'double', diffuseInterior: true,
    });
    const ballRepr = component.addRepresentation('spacefill', {
      sele, color: MUT_HIGHLIGHT_COLOR, opacity: 0.55, radiusScale: 0.55,
      side: 'double', diffuseInterior: true,
    });
    mutSiteReprRef.current = [licoriceRepr, ballRepr];
    return () => {
      if (mutSiteReprRef.current) {
        for (const r of mutSiteReprRef.current) component.removeRepresentation(r);
        mutSiteReprRef.current = null;
      }
    };
  }, [structReady, mutation.highlightRes, mutation.chain]);

  function toggleSpin() {
    const stage = stageRef.current;
    if (!stage) return;
    if (spinning) {
      stage.setSpin(null, 0);
    } else {
      stage.setSpin([0, 1, 0], 0.006);
    }
    setSpinning(s => !s);
  }

  const cc = hexToNGLColor(mutation.color);

  return (
    <div style={{
      position: 'absolute', inset: 0, zIndex: 60,
      background: '#020d1f', display: 'flex', flexDirection: 'column',
    }}>
      {/* Header bar */}
      <div style={{
        display: 'flex', alignItems: 'center', justifyContent: 'space-between',
        padding: '10px 18px', background: 'rgba(0,8,30,0.95)',
        borderBottom: `1px solid ${cc}44`, flexShrink: 0,
      }}>
        <div style={{ display: 'flex', alignItems: 'center', gap: 10 }}>
          <button
            onClick={onBack}
            style={{
              background: 'rgba(100,140,255,.12)', border: '1px solid rgba(100,140,255,.4)',
              color: 'rgba(160,200,255,0.9)', borderRadius: 8, padding: '5px 13px',
              cursor: 'pointer', fontSize: 12, letterSpacing: 1,
            }}
          >
            ← BACK
          </button>
          {(onPrev || onNext) && (
            <div style={{ display: 'flex', alignItems: 'center', gap: 6 }}>
              <button
                onClick={onPrev}
                disabled={!onPrev}
                style={{
                  background: 'rgba(100,140,255,.12)', border: '1px solid rgba(100,140,255,.4)',
                  color: onPrev ? 'rgba(160,200,255,.9)' : 'rgba(160,200,255,.3)',
                  borderRadius: 8, padding: '5px 10px',
                  cursor: onPrev ? 'pointer' : 'default', fontSize: 12,
                }}
              >
                ◀
              </button>
              {navPosition && (
                <span style={{ color: 'rgba(160,200,255,0.8)', fontSize: 13, letterSpacing: 1 }}>{navPosition}</span>
              )}
              <button
                onClick={onNext}
                disabled={!onNext}
                style={{
                  background: 'rgba(100,140,255,.12)', border: '1px solid rgba(100,140,255,.4)',
                  color: onNext ? 'rgba(160,200,255,.9)' : 'rgba(160,200,255,.3)',
                  borderRadius: 8, padding: '5px 10px',
                  cursor: onNext ? 'pointer' : 'default', fontSize: 12,
                }}
              >
                ▶
              </button>
            </div>
          )}
          <button
            onClick={toggleSpin}
            title={spinning ? 'Pause rotation' : 'Resume rotation'}
            style={{
              background: spinning ? 'rgba(255,180,50,.12)' : 'rgba(80,255,160,.12)',
              border: `1px solid ${spinning ? 'rgba(255,180,50,.5)' : 'rgba(80,255,160,.5)'}`,
              color: spinning ? 'rgba(255,200,80,.9)' : 'rgba(100,255,180,.9)',
              borderRadius: 8, padding: '5px 13px',
              cursor: 'pointer', fontSize: 12, letterSpacing: 1,
            }}
          >
            {spinning ? '⏸ PAUSE' : '▶ ROTATE'}
          </button>
          <button
            onClick={findPockets}
            disabled={pocketStatus === '◌ Running fpocket…'}
            title="Geometric pocket detection (fpocket) — the independent 3D-cavity axis, separate from the DMRG electronic-structure classification"
            style={{
              background: pocketResult ? 'rgba(255,210,30,.14)' : 'rgba(170,120,255,.12)',
              border: `1px solid ${pocketResult ? 'rgba(255,210,30,.55)' : 'rgba(170,120,255,.45)'}`,
              color: pocketResult ? 'rgba(255,220,100,.95)' : 'rgba(200,170,255,.9)',
              borderRadius: 8, padding: '5px 13px',
              cursor: pocketStatus === '◌ Running fpocket…' ? 'default' : 'pointer',
              fontSize: 12, letterSpacing: 1,
            }}
          >
            🔎 FIND POCKETS
          </button>
          {mutation.highlightRes && mutation.highlightRes.length > 0 && (
            <button
              onClick={goToMutation}
              title={`Zoom to the mutation residue (res ${mutation.highlightRes.join(', ')}, chain ${mutation.chain})`}
              style={{
                background: 'rgba(255,210,30,.12)', border: '1px solid rgba(255,210,30,.45)',
                color: 'rgba(255,220,100,.95)', borderRadius: 8, padding: '5px 13px',
                cursor: 'pointer', fontSize: 12, letterSpacing: 1,
              }}
            >
              📍 GO TO MUTATION
            </button>
          )}
          <div>
            <span style={{ color: cc, fontWeight: 'bold', fontSize: 16, letterSpacing: 3 }}>
              {mutation.id}
            </span>
            <span style={{ color: 'rgba(180,210,255,0.9)', fontSize: 14, letterSpacing: 2, marginLeft: 10 }}>
              {mutation.variant}
            </span>
          </div>
        </div>

        <div style={{ textAlign: 'right' }}>
          <div style={{ color: 'rgba(160,200,255,0.8)', fontSize: 12, letterSpacing: 2 }}>
            PDB CRYSTALLOGRAPHIC DATA
          </div>
          <div style={{ color: cc, fontSize: 13, fontWeight: 'bold', letterSpacing: 2 }}>
            {mutation.pdb.toUpperCase()}
          </div>
        </div>

        <div style={{ textAlign: 'right', maxWidth: 280 }}>
          <div style={{ color: 'rgba(160,200,255,0.8)', fontSize: 12, letterSpacing: 1.5, marginBottom: 2 }}>
            TARGETED THERAPY
          </div>
          <div style={{ color: 'rgba(220,235,255,0.9)', fontSize: 13 }}>
            {mutation.drug}
          </div>
          <div style={{ color: 'rgba(150,180,255,0.8)', fontSize: 12 }}>
            {mutation.phase}
          </div>
        </div>
      </div>

      {/* NGL canvas */}
      <div style={{ flex: 1, position: 'relative' }}>
        {/* NGL.Stage owns this div exclusively — its mount effect does
            `el.innerHTML = ''` on every structure change to defend against
            leaked canvases (see that effect's own comment). That call used
            to target THIS OUTER div, which also held the overlay cards
            below as React children — wiping it out silently deleted the
            overlay DOM nodes too, out from under React, every single
            reload. That's why the PDB-metadata card (and the mech card)
            intermittently vanished (reported live 2026-09-09). Giving NGL
            its own dedicated child div means clearing it can never touch
            the overlays, which live as siblings instead. */}
        <div ref={mountRef} style={{ position: 'absolute', inset: 0 }} />

        {hoverInfo && (
          <div style={{
            position: 'absolute', zIndex: 20, pointerEvents: 'none',
            left: hoverInfo.x + 12, top: hoverInfo.y + 12,
            background: 'rgba(2,6,18,.95)', border: `1px solid ${cc}66`, borderRadius: 6,
            padding: '4px 9px', fontSize: 12.5, color: 'rgba(230,240,255,0.98)',
            whiteSpace: 'nowrap',
          }}>
            chain <b style={{ color: cc }}>{hoverInfo.chain}</b>
            {' · '}{hoverInfo.resname}{hoverInfo.resno}
            {mutation.highlightRes?.includes(hoverInfo.resno) && hoverInfo.chain === mutation.chain && (
              <span style={{ color: '#39ff14' }}> ● mutation site</span>
            )}
          </div>
        )}

        {mutation.mech && (
          <div style={{
            position: 'absolute', top: 14, left: 14, zIndex: 10, maxWidth: 300,
            background: 'rgba(2,6,18,.88)', border: `1px solid ${cc}44`, borderRadius: 10,
            padding: '12px 16px', backdropFilter: 'blur(10px)',
          }}>
            {mutation.sub && (
              <div style={{ color: 'rgba(190,215,255,0.9)', fontSize: 13, marginBottom: 6 }}>{mutation.sub}</div>
            )}
            <div style={{ color: cc, fontSize: 12, letterSpacing: 2, marginBottom: 5 }}>● BINDING MECHANISM</div>
            <div style={{ color: 'rgba(220,235,255,0.95)', fontSize: 13.5, lineHeight: 1.6 }}>{mutation.mech}</div>
          </div>
        )}

        {(pocketStatus || pocketResult) && (
          <div style={{
            position: 'absolute', top: 14, right: 14, zIndex: 10, width: 290, maxHeight: 'calc(100% - 28px)',
            background: 'rgba(2,6,18,.93)', border: `1px solid ${cc}55`, borderRadius: 10,
            padding: '14px 16px', backdropFilter: 'blur(10px)',
            display: 'flex', flexDirection: 'column', gap: 10, overflow: 'hidden',
          }}>
            <div style={{ color: cc, fontSize: 14, letterSpacing: 2, fontWeight: 700 }}>● POCKET DETECTION — fpocket</div>

            {pocketStatus && (
              <div style={{ color: 'rgba(220,235,255,0.9)', fontSize: 13, lineHeight: 1.6 }}>{pocketStatus}</div>
            )}

            {pocketResult && (
              <>
                <div style={{ color: 'rgba(210,225,255,0.85)', fontSize: 12, lineHeight: 1.65 }}>
                  fpocket scans the <b>entire protein surface</b> for 3D cavities a drug-like
                  molecule could physically fit into — a question that's <b>independent</b> of
                  whether SOLANGE's DMRG classification found the electronic structure there
                  classically tractable. A target can be Class B with no pocket at all, or Class
                  A with a perfectly good one.
                </div>

                <div style={{ color: 'rgba(220,235,255,0.95)', fontSize: 12.5, lineHeight: 1.6 }}>
                  <b>{pocketResult.n_pockets_total}</b> candidate cavit{pocketResult.n_pockets_total === 1 ? 'y' : 'ies'} found
                  on the whole structure · <b>{pocketResult.n_single_chain_druggable}</b> clear the
                  conventional druggability bar (and aren't crystal-packing artifacts).
                </div>
                <div style={{ color: 'rgba(190,215,255,0.85)', fontSize: 11.5, lineHeight: 1.5 }}>
                  "BEST" below is ranked by fpocket's own combined <b>score</b>, not by
                  druggability alone — a pocket with lower druggability can still rank
                  higher if its overall score (shape, enclosure, etc.) is better. Both
                  numbers are shown on every row so you can see why.
                </div>

                {pocketResult.near_residue != null && (
                  atMutationSite.length > 0 ? (
                    <div style={{ color: '#fff', fontSize: 12.5, lineHeight: 1.6, background: 'rgba(255,255,255,.06)', borderRadius: 6, padding: '6px 9px' }}>
                      ✓ {atMutationSite.length} of these cavit{atMutationSite.length === 1 ? 'y overlaps' : 'ies overlap'} the
                      mutation residue itself ({pocketResult.near_residue}) — outlined in white below.
                    </div>
                  ) : (
                    <div style={{ color: 'rgba(190,215,255,0.95)', fontSize: 12.5, lineHeight: 1.6, background: 'rgba(255,255,255,.04)', borderRadius: 6, padding: '6px 9px' }}>
                      ✗ None of the candidate cavities overlap the mutation residue itself
                      ({pocketResult.near_residue}) — whatever pockets exist elsewhere on the
                      structure don't tell you the mutation site itself is druggable.
                    </div>
                  )
                )}

                {!pocketResult.best_single_chain_pocket && (
                  <div style={{ color: 'rgba(190,215,255,0.85)', fontSize: 12, lineHeight: 1.5 }}>
                    No single-chain candidate clears the druggability bar — consistent with a
                    genuinely non-druggable target by this criterion (not a tool failure).
                  </div>
                )}

                {shownPockets.length > 0 && (
                  <>
                    <div style={{ color: 'rgba(160,200,255,0.8)', fontSize: 13, letterSpacing: 1.5 }}>
                      SHOWING {shownPockets.length} OF {pocketResult.n_pockets_total} — CLICK TO FOCUS
                    </div>
                    <div style={{ overflowY: 'auto', display: 'flex', flexDirection: 'column', gap: 5, paddingRight: 2 }}>
                      {shownPockets
                        .slice()
                        .sort((a, b) => (b.score ?? -999) - (a.score ?? -999))
                        .map(p => {
                          const dotColor = !p.single_chain ? '#8c2626' : p.pocket_id === bestId ? '#ffd20a' : '#db40c7';
                          const selected = p.pocket_id === selectedPocketId;
                          // All chains in the deposited entry are drawn now (the primary
                          // one in full color, the rest dimmed grey) — so every pocket IS
                          // visible somewhere in the scene. This just flags which ones sit
                          // on a DIMMED (non-primary) chain, since "pocket 2 looks great
                          // but score low in my search" can still confuse without knowing
                          // it's on a secondary chain, not the mutation's own gene product.
                          const onPrimaryChain = p.chains.includes(mutation.chain);
                          return (
                            <button
                              key={p.pocket_id}
                              onClick={() => focusPocket(p)}
                              style={{
                                textAlign: 'left', cursor: 'pointer', fontFamily: 'inherit',
                                background: selected ? 'rgba(255,255,255,.12)' : 'rgba(255,255,255,.04)',
                                border: `1px solid ${selected ? 'rgba(255,255,255,.5)' : 'rgba(255,255,255,.1)'}`,
                                borderRadius: 7, padding: '6px 9px', color: 'rgba(220,235,255,0.95)', fontSize: 14.5,
                                display: 'flex', alignItems: 'center', gap: 8,
                              }}
                            >
                              <span style={{ width: 9, height: 9, borderRadius: '50%', background: dotColor, flexShrink: 0 }} />
                              <span style={{ flex: 1 }}>
                                pocket {p.pocket_id}
                                {' · '}druggability {(p.druggability_score ?? 0).toFixed(2)}
                                {' · '}<span title="fpocket's own combined score — ranks candidates here, NOT druggability alone (a higher druggability can still rank below a lower one; see module docstring on the TP53 2OCJ crystal-contact-artifact case this distinction was built to catch)">
                                  score {(p.score ?? 0).toFixed(2)}
                                </span>
                                {' · '}{(p.volume ?? 0).toFixed(0)} Å³
                                {' · '}chain{p.chains.length > 1 ? 's' : ''} {p.chains.join('/')}
                              </span>
                              {p.pocket_id === bestId && <span style={{ color: '#ffd20a', fontSize: 12, letterSpacing: 1 }}>BEST</span>}
                              {p.includes_target_residue && <span style={{ color: '#fff', fontSize: 12, letterSpacing: 1 }}>● SITE</span>}
                              {!onPrimaryChain && (
                                <span style={{ color: '#9aa8c4', fontSize: 11, letterSpacing: 0.5 }} title="This pocket is on a complex partner chain, dimmed grey in the 3D view — not the mutation's own chain">
                                  ◌ dimmed chain
                                </span>
                              )}
                            </button>
                          );
                        })}
                    </div>
                    {!showAllPockets && shownPockets.length < pocketResult.n_pockets_total && (
                      <button
                        onClick={() => setShowAllPockets(true)}
                        style={{
                          alignSelf: 'flex-start', background: 'transparent', border: '1px solid rgba(160,200,255,.3)',
                          color: 'rgba(160,200,255,0.8)', borderRadius: 6, padding: '4px 10px', fontSize: 13.5, cursor: 'pointer',
                        }}
                      >
                        show all {pocketResult.n_pockets_total} candidates (incl. low-score / crystal artifacts)
                      </button>
                    )}
                    {showAllPockets && (
                      <button
                        onClick={() => setShowAllPockets(false)}
                        style={{
                          alignSelf: 'flex-start', background: 'transparent', border: '1px solid rgba(160,200,255,.3)',
                          color: 'rgba(160,200,255,0.8)', borderRadius: 6, padding: '4px 10px', fontSize: 13.5, cursor: 'pointer',
                        }}
                      >
                        show top candidates only
                      </button>
                    )}
                  </>
                )}
              </>
            )}
          </div>
        )}

        {/* PDB structure data — RCSB-sourced, shown for every structure (the
            curated mech card above only exists for the five demo mutations;
            this is the real data a manually-searched gene/PDB ID otherwise
            had none of). */}
        <div style={{
          position: 'absolute', bottom: 14, left: 14, zIndex: 10, maxWidth: 260,
          background: 'rgba(2,6,18,.88)', border: `1px solid ${cc}44`, borderRadius: 10,
          padding: '10px 14px', backdropFilter: 'blur(10px)',
        }}>
          <div style={{ color: cc, fontSize: 12, letterSpacing: 2, marginBottom: 5 }}>● PDB STRUCTURE DATA — RCSB</div>
          {isAlphaFold ? (
            <div style={{ color: 'rgba(190,215,255,0.95)', fontSize: 13 }}>
              AlphaFold predicted model — no RCSB experimental record.
            </div>
          ) : pdbMeta ? (
            <>
              <div style={{ color: 'rgba(220,235,255,0.95)', fontSize: 13.5, lineHeight: 1.5 }}>{pdbMeta.title}</div>
              <div style={{ color: 'rgba(150,180,255,0.9)', fontSize: 12, marginTop: 4 }}>
                {pdbMeta.method || 'Method unknown'}
                {pdbMeta.resolution != null ? ` · ${pdbMeta.resolution.toFixed(2)} Å` : ''}
              </div>
            </>
          ) : (
            <div style={{ color: 'rgba(190,215,255,0.8)', fontSize: 13 }}>Loading structure data…</div>
          )}
          {chainInfo.length > 0 && (
            <div style={{ marginTop: 8, paddingTop: 8, borderTop: `1px solid ${cc}33` }}>
              <div style={{ color: cc, fontSize: 11, letterSpacing: 1.5, marginBottom: 3 }}>
                ● CHAINS IN THIS FILE {chainInfo.length > 1 ? '(dimmed in the 3D view unless it\'s the primary one)' : ''}
              </div>
              {chainInfo.map(c => (
                <div key={c.chainname} style={{ color: 'rgba(210,225,255,0.9)', fontSize: 12, lineHeight: 1.5 }}>
                  <b style={{ color: c.chainname === mutation.chain ? cc : 'rgba(210,225,255,0.9)' }}>
                    Chain {c.chainname}{c.chainname === mutation.chain ? ' (primary)' : ''}
                  </b>
                  {': '}{c.description}
                </div>
              ))}
            </div>
          )}
        </div>
      </div>

      {/* Footer legend */}
      <div style={{
        padding: '7px 18px', background: 'rgba(0,8,30,0.90)',
        borderTop: `1px solid ${cc}33`, flexShrink: 0,
        display: 'flex', gap: 24, alignItems: 'center',
        color: 'rgba(140,180,255,0.75)', fontSize: 12, letterSpacing: 1.5,
      }}>
        <div style={{ display: 'flex', alignItems: 'center', gap: 6 }}>
          <span>STRUCTURE STYLE:</span>
          <div style={{ display: 'flex', gap: 3 }}>
            {STRUCT_STYLES.map(s => (
              <button
                key={s.key}
                onClick={() => setStructStyle(s.key)}
                title={`Show the whole chain as ${s.label.toLowerCase()}`}
                style={{
                  background: structStyle === s.key ? `${cc}33` : 'rgba(255,255,255,.08)',
                  border: `1px solid ${structStyle === s.key ? cc : 'rgba(255,255,255,.2)'}`,
                  color: structStyle === s.key ? cc : 'rgba(210,225,255,.85)',
                  borderRadius: 5, padding: '3px 10px', fontSize: 12, letterSpacing: 0.5,
                  cursor: 'pointer', fontFamily: 'inherit',
                }}
              >
                {s.label}
              </button>
            ))}
          </div>
        </div>
        {mutation.highlightRes && mutation.highlightRes.length > 0 && (
          <span style={{ color: '#39ff14' }}>● MUTATION SITE — res {mutation.highlightRes.join(', ')}</span>
        )}
        <span style={{ color: '#aaffdd' }}>● Zn²⁺ ion (if present)</span>
        {pocketResult && pocketResult.all_pockets.length > 0 && (
          <>
            <span style={{ color: '#ffd20a' }}>◉ best druggable pocket</span>
            <span style={{ color: '#db40c7' }}>◉ other single-chain pocket</span>
            <span style={{ color: '#8c2626' }}>◉ multi-chain (crystal-contact artifact)</span>
          </>
        )}
        <span style={{ marginLeft: 'auto' }}>
          Source: RCSB PDB · Drag to rotate · Scroll to zoom · ⏸ to pause
        </span>
      </div>
    </div>
  );
}
