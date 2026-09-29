#!/usr/bin/env python3
"""
make_boring_controls.py  --  chemistry-free negative-control geometries for SOLANGE

Builds closed-shell, saturated-hydrocarbon control systems from EXPERIMENTAL
internal coordinates (no software-dependent optimisation), so the geometry is
reproducible byte-for-byte on any machine and contains no stretched bonds,
no metal, no pi system and no near-degeneracy.

Two families:
  * n-alkane C_nH_(2n+2), all-anti  -- extended / quasi-1D orbital topology
  * adamantane C10H16               -- compact / 3-D orbital topology
    (built on the diamond lattice, which IS its experimental structure)

Experimental internals used:
  alkane      C-C 1.531 A, C-C-C 113.3 deg, C-H 1.096 A, H-C-H 106.5 deg
              (gas-phase electron diffraction / polyethylene crystal values)
  adamantane  C-C 1.540 A, C-C-C 109.47 deg (diamond lattice), C-H 1.099 A

Every structure is verified after construction: bond lengths, bond angles,
and the shortest non-bonded contact.
"""
import math, json, sys

# ---------------------------------------------------------------- vector maths
def sub(a, b): return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
def add(a, b): return (a[0]+b[0], a[1]+b[1], a[2]+b[2])
def mul(a, s): return (a[0]*s, a[1]*s, a[2]*s)
def dot(a, b): return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
def norm(a):   return math.sqrt(dot(a, a))
def unit(a):
    n = norm(a)
    return (a[0]/n, a[1]/n, a[2]/n)
def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
def dist(a, b): return norm(sub(a, b))
def angle(a, b, c):
    """angle a-b-c in degrees"""
    u, v = unit(sub(a, b)), unit(sub(c, b))
    return math.degrees(math.acos(max(-1.0, min(1.0, dot(u, v)))))

# --------------------------------------------------------------- H placement
def h_on_ch2(ri, ra, rb, d_ch, hch):
    """two H on a carbon with exactly two carbon neighbours a,b"""
    u, v = unit(sub(ra, ri)), unit(sub(rb, ri))
    bis  = unit(add(u, v))
    w    = unit(cross(u, v))
    beta = math.radians(hch / 2.0)
    d1 = add(mul(bis, -math.cos(beta)), mul(w,  math.sin(beta)))
    d2 = add(mul(bis, -math.cos(beta)), mul(w, -math.sin(beta)))
    return [add(ri, mul(unit(d1), d_ch)), add(ri, mul(unit(d2), d_ch))]

def h_on_ch3(ri, ra, r_ref, d_ch, hcc):
    """three H on a terminal carbon i bonded to a; r_ref (the atom beyond a)
    fixes the rotamer so that one H is anti-periplanar to it (staggered)."""
    ax  = unit(sub(ri, ra))
    ref = sub(r_ref, ra)
    perp = sub(ref, mul(ax, dot(ref, ax)))
    p = unit(perp)
    q = unit(cross(ax, p))
    tilt = math.radians(180.0 - hcc)          # angle of C-H away from ax
    out = []
    for k in range(3):
        phi = math.radians(120.0 * k)
        # start anti to r_ref  ->  in-plane component is -p at k=0
        rad = add(mul(p, -math.cos(phi)), mul(q, math.sin(phi)))
        d   = add(mul(ax, math.cos(tilt)), mul(rad, math.sin(tilt)))
        out.append(add(ri, mul(unit(d), d_ch)))
    return out

def h_on_ch(ri, neigh, d_ch):
    """one H on a carbon with three carbon neighbours"""
    s = (0.0, 0.0, 0.0)
    for rn in neigh:
        s = add(s, unit(sub(rn, ri)))
    return [add(ri, mul(unit(mul(s, -1.0)), d_ch))]

# ------------------------------------------------------------------- builders
CC_A, CCC_A, CH_A, HCH_A, HCC_A = 1.531, 113.3, 1.096, 106.5, 111.0

def build_alkane(n):
    """all-anti n-alkane, carbons zig-zagging in the xy plane"""
    assert n >= 4
    alpha = math.radians((180.0 - CCC_A) / 2.0)
    dx, dy = CC_A*math.cos(alpha), CC_A*math.sin(alpha)
    C = [(i*dx, dy if i % 2 else 0.0, 0.0) for i in range(n)]
    atoms = [("C", c) for c in C]
    H = []
    for i in range(n):
        if i == 0:
            H += h_on_ch3(C[0], C[1], C[2], CH_A, HCC_A)
        elif i == n-1:
            H += h_on_ch3(C[n-1], C[n-2], C[n-3], CH_A, HCC_A)
        else:
            H += h_on_ch2(C[i], C[i-1], C[i+1], CH_A, HCH_A)
    atoms += [("H", h) for h in H]
    bonds = [(i, i+1) for i in range(n-1)]
    return atoms, bonds, n

CC_AD, CH_AD = 1.540, 1.099

def build_adamantane():
    """C10H16 on the diamond lattice: 4 CH (tertiary) + 6 CH2 (bridging)"""
    a4 = CC_AD * 4.0 / math.sqrt(3.0) / 4.0      # (cubic a)/4
    ch_sites  = [(1,1,1), (1,3,3), (3,1,3), (3,3,1)]
    ch2_sites = [(0,2,2), (2,0,2), (2,2,0), (2,2,4), (2,4,2), (4,2,2)]
    ch  = [mul(s, a4) for s in ch_sites]
    ch2 = [mul(s, a4) for s in ch2_sites]
    C = ch + ch2                                  # indices 0-3 CH, 4-9 CH2
    # connectivity from geometry
    bonds = []
    for i in range(10):
        for j in range(i+1, 10):
            if dist(C[i], C[j]) < 1.75:
                bonds.append((i, j))
    atoms = [("C", c) for c in C]
    H = []
    for i in range(10):
        nb = [C[j] for (p, j) in [(b[0], b[1]) for b in bonds] if p == i]
        nb += [C[j] for (j, p) in [(b[0], b[1]) for b in bonds] if p == i]
        if len(nb) == 3:
            H += h_on_ch(C[i], nb, CH_AD)
        elif len(nb) == 2:
            H += h_on_ch2(C[i], nb[0], nb[1], CH_AD, HCH_A)
        else:
            raise RuntimeError("adamantane connectivity wrong at %d: %d" % (i, len(nb)))
    atoms += [("H", h) for h in H]
    return atoms, bonds, 10

# ---------------------------------------------------------------- verification
def verify(atoms, cc_bonds, label):
    pos = [a[1] for a in atoms]
    el  = [a[0] for a in atoms]
    nC  = el.count("C"); nH = el.count("H")
    cc  = [dist(pos[i], pos[j]) for i, j in cc_bonds]
    # C-H bonds by proximity
    chd = []
    for i, e in enumerate(el):
        if e != "H":
            continue
        best = min((dist(pos[i], pos[j]) for j, f in enumerate(el) if f == "C"))
        chd.append(best)
    # C-C-C angles
    nbr = {}
    for i, j in cc_bonds:
        nbr.setdefault(i, []).append(j); nbr.setdefault(j, []).append(i)
    ccc = []
    for b, ns in nbr.items():
        for x in range(len(ns)):
            for y in range(x+1, len(ns)):
                ccc.append(angle(pos[ns[x]], pos[b], pos[ns[y]]))
    # shortest non-bonded contact (any pair > 2 bonds apart, approximated by
    # excluding pairs closer than 1.3 A which are bonds)
    short = min(dist(pos[i], pos[j])
                for i in range(len(pos)) for j in range(i+1, len(pos))
                if dist(pos[i], pos[j]) > 1.30)
    n_sigma = len(cc_bonds) + nH
    return {
        "label": label, "n_C": nC, "n_H": nH,
        "formula": "C%dH%d" % (nC, nH),
        "n_sigma_bonds": n_sigma,
        "valence_electrons": 4*nC + nH,
        "avas_AO_ceiling_C2p_H1s": 3*nC + nH,
        "cc_min": round(min(cc), 4), "cc_max": round(max(cc), 4),
        "ch_min": round(min(chd), 4), "ch_max": round(max(chd), 4),
        "ccc_min": round(min(ccc), 2), "ccc_max": round(max(ccc), 2),
        "shortest_nonbonded": round(short, 3),
    }

def write_xyz(atoms, path, comment):
    with open(path, "w") as fh:
        fh.write("%d\n%s\n" % (len(atoms), comment))
        for e, (x, y, z) in atoms:
            fh.write("%-2s %14.8f %14.8f %14.8f\n" % (e, x, y, z))

# --------------------------------------------------------------------- ladder
def smallest_alkane_for(n_occ_act, n_vir_act):
    """3n+1 sigma bonds must supply both the active occupied and active
    virtual counts (sigma and sigma* come in pairs)."""
    need = max(n_occ_act, n_vir_act)
    n = 4
    while 3*n + 1 < need:
        n += 1
    return n

def main():
    out = {"structures": [], "ladder": []}

    # -- the two size-matched controls for FE2S2_PROXY, CAS(46,32)
    ad_atoms, ad_bonds, _ = build_adamantane()
    v = verify(ad_atoms, ad_bonds, "BORING-CAGE-32 (adamantane)")
    write_xyz(ad_atoms, "controls/adamantane_C10H16.xyz",
              "adamantane C10H16, diamond-lattice experimental geometry, "
              "C-C 1.540 A C-H 1.099 A; BORING-CAGE control for CAS(46,32)")
    out["structures"].append(v)

    oc_atoms, oc_bonds, _ = build_alkane(8)
    v = verify(oc_atoms, oc_bonds, "BORING-CHAIN-32 (n-octane)")
    write_xyz(oc_atoms, "controls/n_octane_C8H18.xyz",
              "n-octane C8H18 all-anti, experimental internals "
              "C-C 1.531 A C-C-C 113.3 deg C-H 1.096 A; "
              "BORING-CHAIN control for CAS(46,32)")
    out["structures"].append(v)

    # -- the size ladder: one boring alkane per target active-space size
    targets = [
        ("SDHB_2FE2S_site",          46,  28),
        ("FE2S2_PROXY",              46,  32),
        ("TP53_sweep_small",         54,  38),
        ("TP53_R175H_CAS96_54",      96,  54),
        ("NEGCTRL_BORING_replicate", 66,  55),
        ("TP53_R175H_CAS120_65",    120,  65),
        ("TP53_sweep_mid",          114,  80),
        ("TP53_sweep_large",        172, 120),
        ("TP53_sweep_max",          338, 235),
    ]
    made = {}
    for name, ne, no in targets:
        nocc, nvir = ne // 2, no - ne // 2
        n = smallest_alkane_for(nocc, nvir)
        if n not in made:
            at, bd, _ = build_alkane(n)
            write_xyz(at, "controls/n_alkane_C%02dH%02d.xyz" % (n, 2*n+2),
                      "n-alkane C%dH%d all-anti, experimental internals; "
                      "boring size-ladder member" % (n, 2*n+2))
            made[n] = verify(at, bd, "ladder n=%d" % n)
        dets = math.comb(no, nocc) ** 2
        out["ladder"].append({
            "matches_target": name, "n_elec": ne, "n_orb": no,
            "n_occ_act": nocc, "n_vir_act": nvir,
            "control_alkane_n": n,
            "control_formula": "C%dH%d" % (n, 2*n+2),
            "sigma_bonds_available": 3*n + 1,
            "qubits_JW": 2 * no,
            "determinants_closed_shell": dets,
        })
    out["ladder_members"] = made
    json.dump(out, open("controls/control_geometry_report.json", "w"), indent=2)

    for s in out["structures"]:
        print(json.dumps(s))
    print("---")
    for l in out["ladder"]:
        print("%-26s CAS(%3d,%3d) occ%3d vir%3d -> %-8s sigma=%3d dets=%.3e" % (
            l["matches_target"], l["n_elec"], l["n_orb"], l["n_occ_act"],
            l["n_vir_act"], l["control_formula"], l["sigma_bonds_available"],
            float(l["determinants_closed_shell"])))

if __name__ == "__main__":
    main()
