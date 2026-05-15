import { useEffect, useRef } from 'react';
import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';

interface Props { onBack: () => void; }

export default function TP53LoopsViewer({ onBack }: Props) {
  const mountRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    const el = mountRef.current;
    if (!el) return;
    const W = el.clientWidth, H = el.clientHeight;

    // ── Renderer ──────────────────────────────────────────────────────────────
    const renderer = new THREE.WebGLRenderer({ antialias: true });
    renderer.setSize(W, H);
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    el.appendChild(renderer.domElement);

    const scene = new THREE.Scene();
    scene.background = new THREE.Color(0x0a1628);
    scene.fog = new THREE.FogExp2(0x0a1628, 0.009);

    const camera = new THREE.PerspectiveCamera(44, W / H, 0.1, 200);
    camera.position.set(0, 4, 46);

    const controls = new OrbitControls(camera, renderer.domElement);
    controls.enableDamping = true;
    controls.dampingFactor = 0.08;
    controls.autoRotate = true;
    controls.autoRotateSpeed = 0.55;
    controls.target.set(0, -2, 0);

    // ── Lights ────────────────────────────────────────────────────────────────
    scene.add(new THREE.AmbientLight(0x4466aa, 4.5));
    const sun = new THREE.DirectionalLight(0xffffff, 2.2);
    sun.position.set(15, 20, 15); scene.add(sun);
    const fill = new THREE.DirectionalLight(0x4488cc, 1.4);
    fill.position.set(-10, -5, -10); scene.add(fill);
    const rim = new THREE.DirectionalLight(0x88bbff, 1.0);
    rim.position.set(0, -10, 20); scene.add(rim);

    // ── Helpers ───────────────────────────────────────────────────────────────
    function makeTube(pts: number[][], r: number, color: number, op = 1, emissive = 0x000000, emissiveIntensity = 0.5) {
      const curve = new THREE.CatmullRomCurve3(pts.map(p => new THREE.Vector3(p[0], p[1], p[2])));
      const geo = new THREE.TubeGeometry(curve, 7 * pts.length, r, 8, false);
      const mat = new THREE.MeshPhongMaterial({ color, emissive, emissiveIntensity, transparent: op < 1, opacity: op, shininess: 90 });
      return new THREE.Mesh(geo, mat);
    }
    function makeSph(x: number, y: number, z: number, r: number, color: number, em = 0x000000) {
      const m = new THREE.Mesh(
        new THREE.SphereGeometry(r, 22, 22),
        new THREE.MeshPhongMaterial({ color, emissive: em, shininess: 90 })
      );
      m.position.set(x, y, z); return m;
    }
    function makeBox(x1: number, y1: number, z1: number, x2: number, y2: number, z2: number, color: number) {
      const a = new THREE.Vector3(x1, y1, z1), b = new THREE.Vector3(x2, y2, z2);
      const mid = a.clone().add(b).multiplyScalar(0.5);
      const len = a.distanceTo(b);
      const m = new THREE.Mesh(
        new THREE.BoxGeometry(len, 0.6, 2.4),
        new THREE.MeshPhongMaterial({ color, shininess: 35 })
      );
      m.position.copy(mid); m.lookAt(b); m.rotateY(Math.PI / 2); return m;
    }

    // ── Beta-sandwich upper sheet ─────────────────────────────────────────────
    [[-6, 5.5], [-2.5, 6.3], [1.5, 6.3], [5.5, 5.5]].forEach(([vx, vy]) => {
      scene.add(makeBox(vx, vy, -4.5, vx, vy, 4.5, 0x3a6ecc));
    });
    // Lower sheet
    [[-6, -2.5], [-2.5, -1.5], [1.5, -1.5], [5.5, -2.5]].forEach(([vx, vy]) => {
      scene.add(makeBox(vx, vy, -3.5, vx, vy, 5.0, 0x2a55aa));
    });

    // Core ellipsoid
    const cg = new THREE.SphereGeometry(1, 28, 20);
    cg.applyMatrix4(new THREE.Matrix4().makeScale(9.5, 7, 5));
    scene.add(new THREE.Mesh(cg, new THREE.MeshPhongMaterial({ color: 0x162e60, emissive: 0x0a1830, emissiveIntensity: 0.4, transparent: true, opacity: 0.60 })));

    // H2 helix
    const hpts: number[][] = [];
    for (let i = 0; i < 18; i++) {
      const t = i / 17;
      hpts.push([8.8 + Math.cos(t * Math.PI * 2.2) * 1.3, 5.5 - t * 9, Math.sin(t * Math.PI * 2.2) * 1.3]);
    }
    scene.add(makeTube(hpts, 0.7, 0x3366ee, 1, 0x112266, 0.4));

    // Zinc ion
    const zinc = makeSph(-3.5, 0.5, 0, 1.0, 0xffee22, 0xaa8800);
    scene.add(zinc);
    [[-1.5, 2, -1.5], [-6.5, -0.3, -0.5], [-5, 1, 2.8], [-3, -1.5, 2.2]].forEach(p => {
      scene.add(makeTube([[-3.5, 0.5, 0], p], 0.14, 0xffdd44, 0.9, 0x664400, 0.5));
    });

    // L1 Loop — Minor Groove contact (cyan)
    scene.add(makeTube([
      [-6, 5.5, -4.5], [-8.8, 3.5, -5.8], [-10.5, 0.5, -5.2], [-9.8, -2.5, -4.2], [-6.5, -3.5, -3.2]
    ], 0.85, 0x00ddff, 1, 0x004466, 0.55));

    // L2 Loop — Zinc + Major Groove (amber)
    scene.add(makeTube([
      [-2.5, 6, -4.5], [-5, 3, -6.5], [-4, 0.5, -2.5], [-3.5, -2, 0], [-2.5, -4.5, 0.8]
    ], 0.85, 0xffaa00, 1, 0x552200, 0.55));

    // L3 Loop — Major Groove DNA reader (gold) — C275 sits here
    scene.add(makeTube([
      [1.5, -1.5, 4.8], [0, -4, 6.5], [-1.5, -7.5, 5.8],
      [-2, -11.5, 4.2], [-2.2, -13.8, 3.5], [-0.5, -12, 2.5],
      [2, -9, 1.5], [4.5, -5.5, 0.5], [5.5, -2.5, -3.2]
    ], 0.95, 0xFFBB30, 1, 0x664400, 0.6));

    // C275 residue — mutation site (pulsing green)
    const c275 = makeSph(-2.2, -13.8, 3.5, 1.45, 0x00ff88, 0x00dd55);
    scene.add(c275);
    // H-bond line to DNA
    scene.add(makeTube([[-2.2, -13.8, 3.5], [-2.2, -19, 3.5]], 0.14, 0x00ff88, 1.0, 0x00aa44, 0.6));

    // DNA double helix
    const dnaG = new THREE.Group();
    dnaG.position.set(0, -23, 3);
    dnaG.rotation.x = Math.PI / 2;
    const pA: THREE.Vector3[] = [], pB: THREE.Vector3[] = [];
    for (let d = 0; d <= 90; d++) {
      const td = d / 90, ad = td * Math.PI * 2 * 1.6, zd = (td - 0.5) * 30;
      pA.push(new THREE.Vector3(Math.cos(ad) * 3.8, Math.sin(ad) * 3.8, zd));
      pB.push(new THREE.Vector3(Math.cos(ad + Math.PI) * 3.8, Math.sin(ad + Math.PI) * 3.8, zd));
    }
    dnaG.add(new THREE.Mesh(
      new THREE.TubeGeometry(new THREE.CatmullRomCurve3(pA), 180, 0.62, 8),
      new THREE.MeshPhongMaterial({ color: 0x4488ee, emissive: 0x112244, emissiveIntensity: 0.5, shininess: 80 })
    ));
    dnaG.add(new THREE.Mesh(
      new THREE.TubeGeometry(new THREE.CatmullRomCurve3(pB), 180, 0.62, 8),
      new THREE.MeshPhongMaterial({ color: 0x2266bb, emissive: 0x0a1a44, emissiveIntensity: 0.5, shininess: 80 })
    ));
    for (let bp = 0; bp <= 90; bp += 9) {
      const tbp = bp / 90, abp = tbp * Math.PI * 2 * 1.6, zbp = (tbp - 0.5) * 30;
      const bpC = new THREE.LineCurve3(
        new THREE.Vector3(Math.cos(abp) * 3.8, Math.sin(abp) * 3.8, zbp),
        new THREE.Vector3(Math.cos(abp + Math.PI) * 3.8, Math.sin(abp + Math.PI) * 3.8, zbp)
      );
      dnaG.add(new THREE.Mesh(
        new THREE.TubeGeometry(bpC, 2, 0.24, 6),
        new THREE.MeshPhongMaterial({ color: 0x4488ff, transparent: true, opacity: 0.65 })
      ));
    }
    // Major groove highlight
    dnaG.add(new THREE.Mesh(
      new THREE.TorusGeometry(3.8, 1.6, 8, 22, Math.PI * 0.65),
      new THREE.MeshPhongMaterial({ color: 0xE8A020, transparent: true, opacity: 0.18 })
    ));
    scene.add(dnaG);

    // Background particles
    const pg = new THREE.BufferGeometry();
    const pp = new Float32Array(500 * 3);
    for (let i = 0; i < 500; i++) {
      pp[i * 3]     = (Math.random() - 0.5) * 110;
      pp[i * 3 + 1] = (Math.random() - 0.5) * 85;
      pp[i * 3 + 2] = (Math.random() - 0.5) * 65 - 20;
    }
    pg.setAttribute('position', new THREE.BufferAttribute(pp, 3));
    scene.add(new THREE.Points(pg, new THREE.PointsMaterial({ color: 0x334466, size: 0.26, transparent: true, opacity: 0.5 })));

    // ── Resize ────────────────────────────────────────────────────────────────
    const onResize = () => {
      const w = el.clientWidth, h = el.clientHeight;
      camera.aspect = w / h; camera.updateProjectionMatrix();
      renderer.setSize(w, h);
    };
    window.addEventListener('resize', onResize);

    // ── Animate ───────────────────────────────────────────────────────────────
    let tick = 0, rafId: number;
    function animate() {
      rafId = requestAnimationFrame(animate);
      tick += 0.016;
      const p = 1 + Math.sin(tick * 2.2) * 0.13;
      c275.scale.set(p, p, p);
      (c275.material as THREE.MeshPhongMaterial).emissiveIntensity = 0.65 + Math.sin(tick * 2.2) * 0.3;
      (zinc.material as THREE.MeshPhongMaterial).emissiveIntensity = 0.4 + Math.sin(tick * 1.7) * 0.2;
      controls.update();
      renderer.render(scene, camera);
    }
    animate();

    return () => {
      cancelAnimationFrame(rafId);
      window.removeEventListener('resize', onResize);
      controls.dispose(); renderer.dispose();
      if (renderer.domElement.parentElement) renderer.domElement.parentElement.removeChild(renderer.domElement);
    };
  }, []);

  return (
    <div style={{ width: '100%', height: '100%', position: 'relative', background: '#07101f' }}>
      <div ref={mountRef} style={{ width: '100%', height: '100%' }} />

      {/* Legend */}
      <div style={{
        position: 'absolute', top: 16, left: 16,
        background: 'rgba(5,12,28,.92)', border: '1px solid rgba(80,120,200,.28)',
        borderRadius: 10, padding: '14px 18px', fontFamily: 'Courier New',
        fontSize: 11, lineHeight: 2, color: '#ccd6f0',
        backdropFilter: 'blur(10px)',
      }}>
        <div style={{ fontSize: 12, fontWeight: 'bold', color: '#fff', marginBottom: 6, letterSpacing: 2 }}>STRUCTURE KEY</div>
        {[
          ['#1e3f88', 'Beta-Sandwich Core (scaffold)'],
          ['#ffee22', 'Zinc ion Zn²⁺'],
          ['#00bbdd', 'L1 Loop — Minor Groove contact'],
          ['#cc8800', 'L2 Loop — Zinc + Major Groove'],
          ['#E8A020', 'L3 Loop — Major Groove (DNA reader)'],
          ['#00ff88', 'C275 — H-bond to DNA (mutation site)'],
          ['#2255aa', 'DNA Double Helix (strand A)'],
          ['#113366', 'DNA Double Helix (strand B)'],
        ].map(([color, label]) => (
          <div key={label} style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
            <div style={{ width: 10, height: 10, borderRadius: '50%', background: color, flexShrink: 0, boxShadow: `0 0 5px ${color}` }} />
            <span style={{ fontSize: 10, color: 'rgba(200,220,255,.9)' }}>{label}</span>
          </div>
        ))}
      </div>

      {/* Title panel */}
      <div style={{
        position: 'absolute', top: 16, right: 16,
        background: 'rgba(5,12,28,.92)', border: '1px solid rgba(80,120,200,.28)',
        borderRadius: 10, padding: '14px 18px', fontFamily: 'Courier New',
        textAlign: 'right', backdropFilter: 'blur(10px)',
      }}>
        <div style={{ fontSize: 15, fontWeight: 'bold', color: '#fff', letterSpacing: 2, marginBottom: 4 }}>TP53 DNA-BINDING DOMAIN</div>
        <div style={{ fontSize: 10, color: 'rgba(180,210,255,.85)', lineHeight: 1.8 }}>
          Loops L1 · L2 · L3 + DNA Double Helix<br />
          <span style={{ color: '#00ff88' }}>Mutation: p.Cys275Phe (C275F)</span><br />
          <span style={{ color: 'rgba(255,221,34,.8)' }}>Active space: ~44e · ~88 qubits</span><br />
          <span style={{ color: 'rgba(140,170,220,.65)', fontSize: 9 }}>IST 697 · Doron Cohen · CGU 2026</span>
        </div>
      </div>

      {/* BQP badge */}
      <div style={{
        position: 'absolute', bottom: 60, right: 16,
        background: 'rgba(255,221,34,.08)', border: '1px solid rgba(255,221,34,.3)',
        borderRadius: 8, padding: '8px 14px', fontFamily: 'Courier New',
        backdropFilter: 'blur(8px)',
      }}>
        <div style={{ color: '#ffdd22', fontSize: 9, letterSpacing: 2, marginBottom: 3 }}>◈ BQP QUANTUM TARGET</div>
        <div style={{ color: 'rgba(255,240,180,.88)', fontSize: 9, lineHeight: 1.7 }}>
          Phe275 π-electrons: 6e / 6 orbitals<br />
          Full pocket: ~44e · Jordan-Wigner: <span style={{ color: '#44ffee', fontWeight: 'bold' }}>~88 qubits</span><br />
          <span style={{ color: 'rgba(100,220,255,.8)' }}>IBM Heron r3 · sqDRIFT-tractable</span>
        </div>
      </div>

      {/* Hint */}
      <div style={{
        position: 'absolute', bottom: 20, left: '50%', transform: 'translateX(-50%)',
        color: 'rgba(255,255,255,.3)', fontSize: 11, fontFamily: 'Courier New', letterSpacing: 1,
      }}>
        DRAG TO ROTATE · SCROLL TO ZOOM
      </div>

      {/* Back button */}
      <button onClick={onBack} style={{
        position: 'absolute', bottom: 16, left: 16,
        background: 'rgba(8,18,50,.92)', border: '1px solid rgba(70,140,255,.4)',
        color: 'rgba(150,200,255,.9)', borderRadius: 18, padding: '8px 18px',
        cursor: 'pointer', fontSize: 11, fontFamily: 'Courier New',
        letterSpacing: 1.5, backdropFilter: 'blur(8px)',
      }}>
        ← BACK TO VIEWER
      </button>
    </div>
  );
}
