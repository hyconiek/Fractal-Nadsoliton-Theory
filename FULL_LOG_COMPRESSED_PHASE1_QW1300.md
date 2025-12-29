# SUPER EXTREME LOG PHASE 1
**Electromagnetism.**

## QW-1301
### S:QW-1301_Coulomb_Force.py
```python
REPORT_FILE = "RAPORT_QW1301_COULOMB_FORCE.md"
md = []
def log(msg):
    print(msg)
    md.append(msg)
log("=" * 80)
log("QW-1301: CZY SKYRMIONY ODDZIAŁUJĄ SIŁĄ COULOMBA?")
log("=" * 80)
def profile_f(r, R_size=1.0):
    return np.where(r < 1e-6, np.pi, 2.0 * np.arctan(R_size / r))
def energy_density_at_point(x, y, z, D):
    r1 = np.sqrt(x**2 + y**2 + z**2)
    r2 = np.sqrt(x**2 + y**2 + (z-D)**2)
    f1 = profile_f(r1)
    f2 = profile_f(r2)
    gf1 = 2.0 / (r1**2 + 1.0) 
    gf2 = 2.0 / (r2**2 + 1.0)
    return gf1 * gf2
def compute_interaction_energy(D):
    R_max = D + 10.0
    r_vals = np.linspace(0, R_max, 50)
    z_vals = np.linspace(-R_max, R_max, 100)
    rr, zz = np.meshgrid(r_vals, z_vals)
    dens = energy_density_at_point(rr, 0, zz, D)
    integrand = dens * rr 
    res = 2 * np.pi * simpson(simpson(integrand, r_vals), z_vals)
    return res
log("\n[1] SYMULACJA ZALEŻNOŚCI E(D)")
d_vals = np.linspace(2.0, 10.0, 9)
e_vals = []
log(f"{'Dystans D':<10} {'E_int':<15} {'E * D':<10} {'E * D^3':<10}")
log("-" * 50)
for d in d_vals:
    e = compute_interaction_energy(d)
    e_vals.append(e)
    col_test = e * d
    dip_test = e * d**3
    log(f"{d:<10.2f} {e:<15.6f} {col_test:<10.4f} {dip_test:<10.4f}")
log("\n[2] WNIOSKI")
log("-" * 80)
ed_vals = np.array(e_vals) * d_vals
slope = (ed_vals[-1] - ed_vals[0]) / len(ed_vals)
if abs(slope) < 0.1 * np.mean(ed_vals):
    log("WYNIK: E * D jest w przybliżeniu stałe.")
    log("Prawo oddziaływania: 1/D (COULOMB) ✅")
    log("Interpretacja: Teoria FIN generuje siły elektrostatyczne z geometrii!")
else:
    log("WYNIK: E * D nie jest stałe.")
    ed3_vals = np.array(e_vals) * d_vals**3
    if np.std(ed3_vals) < 0.2 * np.mean(ed3_vals):
        log("Prawo oddziaływania: 1/D^3 (DIPOL) ⚠️")
        log("Interpretacja: Skyrmiony zachowują się jak neutrane dipole, brakuje ładunku U(1).")
    else:
        log("Prawo oddziaływania: Złożone (Yukawa/Pośrednie).")
    f.write("
    f.write(f"**Data:** {datetime.now()}\n\n")
    f.write("```\n")
    f.write("\n".join(md))
    f.write("\n```\n")
print(f"Raport zapisano w {REPORT_FILE}")
```
### R:RAPORT_QW1301_COULOMB_FORCE.md
```markdown
# QW-1301: Emergentna Siła Coulomba
Prawo oddziaływania: 1/D (COULOMB) ✅
```
--------------------
## QW-1302
### S:QW-1302_Alpha_Constant.py
```python
REPORT_FILE = "RAPORT_QW1302_ALPHA_CALC.md"
md = []
def log(msg):
    print(msg)
    md.append(msg)
log("=" * 80)
log("QW-1302: POSZUKIWANIE GEOMETRYCZNEJ ALFY")
log("=" * 80)
def profile_deriv_sq(r, R_size=1.0):
    return (2.0 * R_size / (r**2 + R_size**2))**2
def compute_self_energy(R_size=1.0):
    r = np.linspace(0, 1000, 100000)
    integrand = profile_deriv_sq(r, R_size) * r**2
    return 4 * np.pi * simpson(integrand, r)
log("\n[1] OBLICZENIA ENERGII")
E_rest = compute_self_energy()
C_coulomb = 85.655 
log(f"Energia Spoczynkowa E_rest: {E_rest:.4f}")
log(f"Stała Coulomba C (E*D):     {C_coulomb:.4f}")
log("\n[2] TESTOWANIE HIPOTEZ NA ALPHĘ")
log("-" * 60)
alpha_1 = 1.0 / C_coulomb
log(f"H1: alpha = 1/C             = {alpha_1:.6f} (Exp: 0.007297)")
log(f"    Różnica: {abs(alpha_1 - 1/137.036)/(1/137.036)*100:.2f}%")
alpha_2 = 4 * np.pi / C_coulomb
log(f"H2: alpha = 4pi/C           = {alpha_2:.6f}")
alpha_3 = E_rest / C_coulomb
log(f"H3: alpha = E_rest/C        = {alpha_3:.6f}")
alpha_4 = 1.0 / E_rest
log(f"H4: alpha = 1/E_rest        = {alpha_4:.6f}")
alpha_5 = 1.0 / (C_coulomb + E_rest)
log(f"H5: alpha = 1/(C+E)         = {alpha_5:.6f}")
log("\n[3] ANALIZA WYNIKÓW")
log("-" * 80)
candidates = [alpha_1, alpha_2, alpha_3, alpha_4, alpha_5]
best = min(candidates, key=lambda x: abs(x - 1/137.036))
log(f"Najlepszy kandydat: {best:.6f}")
if abs(best - 1/137.036) < 0.001:
    log("SUKCES! Znaleziono korelację geometryczną.")
else:
    log("PORAŻKA: Żadna prosta kombinacja nie daje 1/137.")
    log("Wniosek: Alpha wymaga ekranowania kwantowego (polaryzacji próżni),")
    log("sama klasyczna geometria daje zbyt silne sprzężenie.")
    f.write("
    f.write(f"**Data:** {datetime.now()}\n\n")
    f.write("```\n")
    f.write("\n".join(md))
    f.write("\n```\n")
print(f"Raport zapisano w {REPORT_FILE}")
```
### R:RAPORT_QW1302_ALPHA_CALC.md
```markdown
# QW-1302: Obliczenie Stałej Alpha
```
--------------------
## QW-1303
### S:QW-1303_Vector_Coulomb.py
```python
REPORT_FILE = "RAPORT_QW1303_VECTOR_COULOMB.md"
md = []
def log(msg):
    print(msg)
    md.append(msg)
log("=" * 80)
log("QW-1303: COULOMB W PEŁNYM MODELU SU(2)")
log("=" * 80)
tau = [
    np.array([[0, 1], [1, 0]], dtype=complex),
    np.array([[0, -1j], [1j, 0]], dtype=complex),
    np.array([[1, 0], [0, -1]], dtype=complex)
]
Id = np.eye(2, dtype=complex)
def get_U(x, y, z):
    r = np.sqrt(x**2 + y**2 + z**2) + 1e-9
    R = 1.0
    f = 2.0 * np.arctan(R/r)
    n = [x/r, y/r, z/r]
    cf = np.cos(f)
    sf = np.sin(f)
    U = cf * Id + 1j * sf * (n[0]*tau[0] + n[1]*tau[1] + n[2]*tau[2])
    return U
def mat_inv(M):
    return M.conj().T
def trace_L_sq(U, dx=0.1):
    pass
def compute_energy_on_grid(D, box_L=8.0, pts=30):
    x = np.linspace(-box_L, box_L, pts)
    y = np.linspace(-box_L, box_L, pts)
    z = np.linspace(-box_L, box_L, pts)
    dx_grid = x[1] - x[0]
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    log("  Generowanie pól U1 i U2...")
    def get_quat(xx, yy, zz):
        r = np.sqrt(xx**2 + yy**2 + zz**2) + 1e-9
        f = 2.0 * np.arctan(1.0/r)
        sf = np.sin(f)
        cf = np.cos(f) 
        nx = xx/r * sf
        ny = yy/r * sf
        nz = zz/r * sf
        return cf, nx, ny, nz
    u1_0, u1_1, u1_2, u1_3 = get_quat(X, Y, Z)
    u2_0, u2_1, u2_2, u2_3 = get_quat(X, Y, Z - D)
    s1, v1x, v1y, v1z = u1_0, u1_1, u1_2, u1_3
    s2, v2x, v2y, v2z = u2_0, u2_1, u2_2, u2_3
    log("  Mnożenie pól (Product Ansatz)...")
    st = s1*s2 - (v1x*v2x + v1y*v2y + v1z*v2z)
    vtx = s1*v2x + s2*v1x + (v1y*v2z - v1z*v2y)
    vty = s1*v2y + s2*v1y + (v1z*v2x - v1x*v2z)
    vtz = s1*v2z + s2*v1z + (v1x*v2y - v1y*v2x)
    log("  Obliczanie gradientów...")
    grad_st = np.gradient(st, dx_grid) 
    grad_vx = np.gradient(vtx, dx_grid)
    grad_vy = np.gradient(vty, dx_grid)
    grad_vz = np.gradient(vtz, dx_grid)
    def norm_sq(grad_list):
        return grad_list[0]**2 + grad_list[1]**2 + grad_list[2]**2
    E_dens = norm_sq(grad_st) + norm_sq(grad_vx) + norm_sq(grad_vy) + norm_sq(grad_vz)
    E_total = np.sum(E_dens) * (dx_grid**3)
    return E_total * 0.5 
log("\n[1] SYMULACJA PEŁNEGO POLA (SU(2)/Kwaterniony)")
def run_experiment(d_list):
    energies = []
    box_L = 12.0
    pts = 60 
    pass
def compute_interaction_direct(D, box_L=12.0, pts=60):
    x = np.linspace(-box_L, box_L, pts)
    z = np.linspace(-box_L + D/2, box_L + D/2, pts) 
    y = np.linspace(-box_L, box_L, pts)
    dx = x[1] - x[0]
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    def get_quat(xx, yy, zz):
        r = np.sqrt(xx**2 + yy**2 + zz**2) + 1e-9
        f = 2.0 * np.arctan(1.0/r)
        sf = np.sin(f)
        cf = np.cos(f)
        return cf, (xx/r)*sf, (yy/r)*sf, (zz/r)*sf
    s1, v1x, v1y, v1z = get_quat(X, Y, Z) 
    s2, v2x, v2y, v2z = get_quat(X, Y, Z - D) 
    def grad_q(s, vx, vy, vz):
        gs = np.gradient(s, dx)
        gvx = np.gradient(vx, dx)
        gvy = np.gradient(vy, dx)
        gvz = np.gradient(vz, dx)
        return gs, gvx, gvy, gvz 
    g1s, g1x, g1y, g1z = grad_q(s1, v1x, v1y, v1z)
    g2s, g2x, g2y, g2z = grad_q(s2, v2x, v2y, v2z)
    E_int_tot = 0
    for axis in range(3): 
        da_s1, da_x1, da_y1, da_z1 = g1s[axis], g1x[axis], g1y[axis], g1z[axis]
        da_s2, da_x2, da_y2, da_z2 = g2s[axis], g2x[axis], g2y[axis], g2z[axis]
        def quat_mul(sa, xa, ya, za, sb, xb, yb, zb):
            st = sa*sb - (xa*xb + ya*yb + za*zb)
            xt = sa*xb + sb*xa + (ya*zb - za*yb)
            yt = sa*yb + sb*ya + (za*xb - xa*zb)
            zt = sa*zb + sb*za + (xa*yb - ya*xb)
            return st, xt, yt, zt
        t1s, t1x, t1y, t1z = quat_mul(da_s1, da_x1, da_y1, da_z1, s2, v2x, v2y, v2z)
        t2s, t2x, t2y, t2z = quat_mul(s1, v1x, v1y, v1z, da_s2, da_x2, da_y2, da_z2)
        dot = t1s*t2s + t1x*t2x + t1y*t2y + t1z*t2z
        E_int_tot += np.sum(2.0 * dot) * (dx**3)
    return E_int_tot
d_vals = np.linspace(3.0, 8.0, 6)
e_vals = []
log(f"{'D':<10} {'E_int':<15} {'E*D':<10} {'E*D^3':<10}")
log("-" * 50)
for d in d_vals:
    e = compute_interaction_direct(d)
    e_vals.append(e)
    log(f"{d:<10.2f} {e:<15.6f} {e*d:<10.4f} {e*d**3:<10.4f}")
evals = np.array(e_vals)
prod_ED = evals * d_vals
cv = np.std(prod_ED) / np.mean(prod_ED)
log("\n[2] WNIOSKI")
log("-" * 80)
if cv < 0.1:
    log(f"E * D jest stałe (CV = {cv*100:.2f}%).")
    log("POTWIERDZONO PRAWO COULOMBA W MODELU SU(2)!")
else:
    log(f"E * D nie jest stałe (CV = {cv*100:.2f}%).")
    prod_ED3 = evals * d_vals**3
    cv3 = np.std(prod_ED3) / np.mean(prod_ED3)
    if cv3 < 0.1:
         log("Prawo dipolowe (1/r^3).")
    else:
         log("Prawo złożone.")
    f.write("
    f.write(f"**Data:** {datetime.now()}\n\n")
    f.write("```\n")
    f.write("\n".join(md))
    f.write("\n```\n")
print(f"Raport zapisano w {REPORT_FILE}")
```
### R:RAPORT_QW1303_VECTOR_COULOMB.md
```markdown
# QW-1303: Vector Skyrmion Coulomb Check
```
--------------------
