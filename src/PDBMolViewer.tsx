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

// Map hex color → NGL color string
function hexToNGLColor(hex: number): string {
  return '#' + hex.toString(16).padStart(6, '0');
}

export default function PDBMolViewer({ mutation, onBack, onPrev, onNext, navPosition }: Props) {
  const mountRef = useRef<HTMLDivElement>(null);
  const stageRef = useRef<any>(null);
  const [spinning, setSpinning] = useState(true);

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

    const stage = new NGL.Stage(el, {
      backgroundColor: '#020d1f',
      quality: 'high',
      antialias: true,
      impostor: true,
    });
    stageRef.current = stage;

    const mutColor = hexToNGLColor(mutation.color);
    const pdbUrl = mutation.url || `https://files.rcsb.org/download/${mutation.pdb}.pdb`;

    stage.loadFile(pdbUrl, { ext: 'pdb', defaultRepresentation: false }).then((component: any) => {
      const ch = mutation.chain;

      // Cartoon — restrict to one chain only
      component.addRepresentation('cartoon', {
        sele: `:${ch}`,
        colorScheme: 'residueindex',
        smoothSheet: true,
        opacity: 0.92,
      });

      if (mutation.highlightRes && mutation.highlightRes.length > 0) {
        // Mutation site spacefill — chain-qualified
        const sele = mutation.highlightRes.map(r => `${r}:${ch}`).join(' or ');
        component.addRepresentation('spacefill', {
          sele,
          color: mutColor,
          opacity: 1.0,
          radiusScale: 1.4,
        });
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
      console.error('NGL load error:', err);
    });

    stage.setSpin([0, 1, 0], 0.006);

    const onResize = () => stage.handleResize();
    window.addEventListener('resize', onResize);

    return () => {
      window.removeEventListener('resize', onResize);
      stageRef.current = null;
      stage.dispose();
    };
  }, [mutation.pdb]);

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
              color: 'rgba(160,200,255,.9)', borderRadius: 8, padding: '5px 13px',
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
                <span style={{ color: 'rgba(160,200,255,.6)', fontSize: 10, letterSpacing: 1 }}>{navPosition}</span>
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
          <div>
            <span style={{ color: cc, fontWeight: 'bold', fontSize: 16, letterSpacing: 3 }}>
              {mutation.id}
            </span>
            <span style={{ color: 'rgba(180,210,255,.7)', fontSize: 11, letterSpacing: 2, marginLeft: 10 }}>
              {mutation.variant}
            </span>
          </div>
        </div>

        <div style={{ textAlign: 'right' }}>
          <div style={{ color: 'rgba(160,200,255,.6)', fontSize: 9, letterSpacing: 2 }}>
            PDB CRYSTALLOGRAPHIC DATA
          </div>
          <div style={{ color: cc, fontSize: 13, fontWeight: 'bold', letterSpacing: 2 }}>
            {mutation.pdb.toUpperCase()}
          </div>
        </div>

        <div style={{ textAlign: 'right', maxWidth: 280 }}>
          <div style={{ color: 'rgba(160,200,255,.6)', fontSize: 9, letterSpacing: 1.5, marginBottom: 2 }}>
            TARGETED THERAPY
          </div>
          <div style={{ color: 'rgba(220,235,255,.9)', fontSize: 10 }}>
            {mutation.drug}
          </div>
          <div style={{ color: 'rgba(150,180,255,.6)', fontSize: 9 }}>
            {mutation.phase}
          </div>
        </div>
      </div>

      {/* NGL canvas */}
      <div ref={mountRef} style={{ flex: 1, position: 'relative' }}>
        {mutation.mech && (
          <div style={{
            position: 'absolute', top: 14, left: 14, zIndex: 10, maxWidth: 300,
            background: 'rgba(2,6,18,.88)', border: `1px solid ${cc}44`, borderRadius: 10,
            padding: '12px 16px', backdropFilter: 'blur(10px)',
          }}>
            {mutation.sub && (
              <div style={{ color: 'rgba(190,215,255,.9)', fontSize: 10, marginBottom: 6 }}>{mutation.sub}</div>
            )}
            <div style={{ color: cc, fontSize: 9, letterSpacing: 2, marginBottom: 5 }}>● BINDING MECHANISM</div>
            <div style={{ color: 'rgba(220,235,255,.95)', fontSize: 10.5, lineHeight: 1.6 }}>{mutation.mech}</div>
          </div>
        )}
      </div>

      {/* Footer legend */}
      <div style={{
        padding: '7px 18px', background: 'rgba(0,8,30,0.90)',
        borderTop: `1px solid ${cc}33`, flexShrink: 0,
        display: 'flex', gap: 24, alignItems: 'center',
        color: 'rgba(140,180,255,.55)', fontSize: 9, letterSpacing: 1.5,
      }}>
        <span>■ CARTOON — secondary structure</span>
        {mutation.highlightRes && mutation.highlightRes.length > 0 && (
          <span style={{ color: cc }}>● MUTATION SITE — res {mutation.highlightRes.join(', ')}</span>
        )}
        <span style={{ color: '#aaffdd' }}>● Zn²⁺ ion (if present)</span>
        <span style={{ marginLeft: 'auto' }}>
          Source: RCSB PDB · Drag to rotate · Scroll to zoom · ⏸ to pause
        </span>
      </div>
    </div>
  );
}
