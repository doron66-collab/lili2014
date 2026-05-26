import { useEffect, useState } from 'react';

interface Props { onEnter: () => void; }

export default function IntroScreen({ onEnter }: Props) {
  const [visible, setVisible] = useState(false);

  useEffect(() => {
    const t = setTimeout(() => setVisible(true), 80);
    return () => clearTimeout(t);
  }, []);

  return (
    <div style={{
      position: 'fixed', inset: 0, zIndex: 200,
      background: 'radial-gradient(ellipse at 60% 40%, rgba(11,79,108,0.45) 0%, #020610 65%)',
      display: 'flex', flexDirection: 'column', alignItems: 'center', justifyContent: 'center',
      fontFamily: 'Courier New, monospace',
      opacity: visible ? 1 : 0,
      transition: 'opacity 0.6s ease',
    }}>
      {/* Orbital decoration */}
      <svg width="220" height="220" viewBox="0 0 220 220" style={{ marginBottom: 24, opacity: 0.55 }}>
        <ellipse cx="110" cy="110" rx="100" ry="38" fill="none" stroke="#1d8fb0" strokeWidth="0.8" transform="rotate(-30 110 110)" />
        <ellipse cx="110" cy="110" rx="100" ry="38" fill="none" stroke="#1d8fb0" strokeWidth="0.8" transform="rotate(30 110 110)" />
        <ellipse cx="110" cy="110" rx="100" ry="38" fill="none" stroke="#0b4f6c" strokeWidth="0.5" />
        <circle cx="110" cy="110" r="7" fill="#b8441f" />
        <circle cx="207" cy="96" r="3.5" fill="#1d8fb0" />
        <circle cx="62"  cy="176" r="3.5" fill="#1d8fb0" />
      </svg>

      <div style={{ fontSize: 9, letterSpacing: 4, color: 'rgba(29,143,176,0.75)', marginBottom: 18, textTransform: 'uppercase' }}>
        IBM sqDRIFT · NSCLC · CGU 2026
      </div>

      <div style={{ fontSize: 28, fontWeight: 700, color: '#fff', letterSpacing: 3, marginBottom: 6 }}>
        SOLANGE
      </div>
      <div style={{ fontSize: 13, color: 'rgba(150,200,255,0.85)', letterSpacing: 2, marginBottom: 6 }}>
        3D QUANTUM SIMULATION PLATFORM
      </div>
      <div style={{ fontSize: 10, color: 'rgba(100,150,200,0.6)', letterSpacing: 1.5, marginBottom: 48 }}>
        Non-Druggable NSCLC Mutation Targets · TP53 · KEAP1 · STK11 · CDKN2A
      </div>

      <button
        onClick={onEnter}
        style={{
          background: 'transparent',
          border: '1.5px solid rgba(29,143,176,0.7)',
          color: 'rgba(100,200,255,0.95)',
          padding: '12px 38px',
          borderRadius: 30,
          fontSize: 12,
          letterSpacing: 3,
          cursor: 'pointer',
          fontFamily: 'Courier New, monospace',
          textTransform: 'uppercase',
          transition: 'all 0.2s',
          boxShadow: '0 0 24px rgba(29,143,176,0.15)',
        }}
        onMouseOver={e => {
          (e.target as HTMLButtonElement).style.background = 'rgba(29,143,176,0.15)';
          (e.target as HTMLButtonElement).style.boxShadow = '0 0 32px rgba(29,143,176,0.35)';
        }}
        onMouseOut={e => {
          (e.target as HTMLButtonElement).style.background = 'transparent';
          (e.target as HTMLButtonElement).style.boxShadow = '0 0 24px rgba(29,143,176,0.15)';
        }}
      >
        ⚛ Enter Simulation
      </button>

      <div style={{ position: 'absolute', bottom: 28, fontSize: 9, color: 'rgba(100,130,180,0.45)', letterSpacing: 2 }}>
        DORON COHEN · IST 697 · CLAREMONT GRADUATE UNIVERSITY
      </div>
    </div>
  );
}
