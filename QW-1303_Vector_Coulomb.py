import numpy as np
from scipy.integrate import tplquad
from datetime import datetime

# ==============================================================================
# QW-1303: PEŁNE ODZIAŁYWANIE WEKTOROWE SKYRMIONÓW (SU(2))
# ==============================================================================
#
# CEL: Sprawdzić prawo siły między dwoma Skyrmionami w pełnym modelu SU(2)
#      (sigma model nieliniowy).
#      Czy nadal otrzymujemy prawo Coulomba (1/r) dla energii?
#
# MODEL:
# Pole U(x) = sigma + i * tau * pi
# U = cos(f) + i * (n * tau) * sin(f)
# Product Ansatz: U_tot = U1 * U2 (mnożenie macierzy SU(2))
#
# Lagranżjan: L = -1/2 Tr(L_mu L^mu), gdzie L_mu = U^dagger dU/dx_mu
# Energia = -L (statyczna).
#
# ==============================================================================

REPORT_FILE = "RAPORT_QW1303_VECTOR_COULOMB.md"
md = []

def log(msg):
    print(msg)
    md.append(msg)

log("=" * 80)
log("QW-1303: COULOMB W PEŁNYM MODELU SU(2)")
log("=" * 80)

# Macierze Pauliego
tau = [
    np.array([[0, 1], [1, 0]], dtype=complex),
    np.array([[0, -1j], [1j, 0]], dtype=complex),
    np.array([[1, 0], [0, -1]], dtype=complex)
]
Id = np.eye(2, dtype=complex)

def get_U(x, y, z):
    r = np.sqrt(x**2 + y**2 + z**2) + 1e-9
    # Profil z ogonem 1/r (dla pola pionowego)
    # sin(f) ~ 1/r => f ~ arcsin(1/r) ~ 1/r
    # Zastosujmy profil: f(r) = 2 * arctan(R/r), sin(f) ~ 2*R/r (dla dużych r)
    R = 1.0
    f = 2.0 * np.arctan(R/r)
    
    n = [x/r, y/r, z/r]
    
    # U = cos(f) I + i sin(f) (n . tau)
    cf = np.cos(f)
    sf = np.sin(f)
    
    U = cf * Id + 1j * sf * (n[0]*tau[0] + n[1]*tau[1] + n[2]*tau[2])
    return U

def mat_inv(M):
    # Dla SU(2), odwrotność to sprzężenie hermitowskie
    return M.conj().T

def trace_L_sq(U, dx=0.1):
    # Oblicza Tr(L_mu L^mu) numerycznie w punkcie
    # Potrzebujemy dU/dx, dU/dy, dU/dz
    # Numerycznie: (U(x+h) - U(x-h)) / 2h
    
    # To jest zbyt kosztowne do robienia wewnątrz pętli całkującej na siatce TPLQUAD.
    # Użyjemy uproszczenia analitycznego dla Product Ansatz.
    pass

# Podejście SIATKOWE (Grid):
# Zbudujemy siatkę punktów, obliczymy pola macierzowe, potem gradienty.
# To szybsze niż wywoływanie funkcji w pętli.

def compute_energy_on_grid(D, box_L=8.0, pts=30):
    x = np.linspace(-box_L, box_L, pts)
    y = np.linspace(-box_L, box_L, pts)
    z = np.linspace(-box_L, box_L, pts)
    dx_grid = x[1] - x[0]
    
    # Siatki 3D
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    # 1. Oblicz pole U1 (w 0,0,0) i U2 (w 0,0,D)
    # U ma kształt (pts, pts, pts, 2, 2)
    
    log("  Generowanie pól U1 i U2...")
    # Operacje wektorowe na macierzach 2x2 w numpy są trudne bez pętli.
    # Zrobimy to sprytnie: U = u0*Id + i*vec(u)*vec(tau)
    # Mnożenie kwaternionów jest szybsze.
    # U = (u0, u1, u2, u3), gdzie u0=cos(f), vec=sin(f)*n
    
    def get_quat(xx, yy, zz):
        r = np.sqrt(xx**2 + yy**2 + zz**2) + 1e-9
        f = 2.0 * np.arctan(1.0/r)
        sf = np.sin(f)
        cf = np.cos(f) # u0
        nx = xx/r * sf
        ny = yy/r * sf
        nz = zz/r * sf
        return cf, nx, ny, nz

    u1_0, u1_1, u1_2, u1_3 = get_quat(X, Y, Z)
    u2_0, u2_1, u2_2, u2_3 = get_quat(X, Y, Z - D)
    
    # Product Ansatz w języku kwaternionów: q_tot = q1 * q2
    # (s1, v1) * (s2, v2) = (s1s2 - v1.v2, s1v2 + s2v1 + v1 x v2)
    
    s1, v1x, v1y, v1z = u1_0, u1_1, u1_2, u1_3
    s2, v2x, v2y, v2z = u2_0, u2_1, u2_2, u2_3
    
    log("  Mnożenie pól (Product Ansatz)...")
    st = s1*s2 - (v1x*v2x + v1y*v2y + v1z*v2z)
    vtx = s1*v2x + s2*v1x + (v1y*v2z - v1z*v2y)
    vty = s1*v2y + s2*v1y + (v1z*v2x - v1x*v2z)
    vtz = s1*v2z + s2*v1z + (v1x*v2y - v1y*v2x)
    
    # Mamy pole q_tot = (st, vtx, vty, vtz)
    # Energia kinetyczna (gradientowa) sigma modelu:
    # L = (ds)^2 + (dv)^2  (suma kwadratów gradientów wszystkich 4 składowych!)
    # Bo U jest unitarne, więc leży na sferze S3. Metryka na S3 to euklidesowa w R4.
    # E_dens = |grad s|^2 + |grad vx|^2 + |grad vy|^2 + |grad vz|^2.
    
    log("  Obliczanie gradientów...")
    grad_st = np.gradient(st, dx_grid) # Zwraca listę [df/dx, df/dy, df/dz]
    grad_vx = np.gradient(vtx, dx_grid)
    grad_vy = np.gradient(vty, dx_grid)
    grad_vz = np.gradient(vtz, dx_grid)
    
    def norm_sq(grad_list):
        return grad_list[0]**2 + grad_list[1]**2 + grad_list[2]**2
    
    E_dens = norm_sq(grad_st) + norm_sq(grad_vx) + norm_sq(grad_vy) + norm_sq(grad_vz)
    
    # Całkowanie
    E_total = np.sum(E_dens) * (dx_grid**3)
    return E_total * 0.5 # Współczynnik 1/2 z Lagrange'a? Nieistotny dla skalowania.

log("\n[1] SYMULACJA PEŁNEGO POLA (SU(2)/Kwaterniony)")
# Najpierw policzmy energię pojedynczego Skyrmiona (referencyjna)
# Ustawiając D bardzo duże lub licząc oddzielnie.
# Zróbmy 'nieskończone' oddalenie przez policzenie E dla D=0 (z nałożonymi?) 
# Nie, policzmy E(U1) osobno.

# Obliczanie E_single
# Hack: Użyj funkcji `compute_energy_on_grid` z U2 będącym identycznością (f=0 everywhere)
# Ale funkcja nie ma takiej opcji.
# Po prostu: E_single to energia jednego pola.
# E_int = E_tot - 2 * E_single (dla identycznych).

# Ponieważ box jest skończony, E_single zależy od boxu. Musimy użyć tego samego boxu.
# Dla dużych odległości D, Skyrmiony mogą wyjść poza box.
# Musimy użyć dużego Boxu. L=15 (zakres -15 do 15 -> 30 szerokości).
# Dystans do 10. Skyrmiony są w 0 i 10. Box musi obejmować oba.
# Centrum boxu w D/2.

def run_experiment(d_list):
    energies = []
    
    # E_inf (Reference) - Dwa bardzo odległe (nieoddziałujące)
    # Symulujemy je jako sumę energii policzonych osobno w tym samym boxie?
    # Lepiej: E_int = E(U1*U2) - E(U1) - E(U2).
    # Liczymy te 3 wielkości numerycznie w tym samym oknie.
    
    box_L = 12.0
    pts = 60 # Rozdzielczość 24/60 = 0.4. Trochę mało dladzenia (r=1), ale ogon złapie.
    
    # Referencja - pojedynczy Skyrmion w centrum
    # E1 = compute_energy_single(centered at 0)
    # Ale E_int liczymy bezpośrednio z różnicy gęstości?
    # |grad(q1q2)|^2 - |grad q1|^2 - |grad q2|^2.
    # Wzór: |grad(ab)|^2 = |(grad a)b + a(grad b)|^2 
    # = |grad a|^2 |b|^2 + |a|^2 |grad b|^2 + 2 Re( (grad a b) . (a grad b)* )
    # Ponieważ |a|=|b|=1 (unitarne kwaterniony), to:
    # = |grad a|^2 + |grad b|^2 + Interaction_Term.
    # Interaction_Term = E_int_dens.
    # Całkujemy tylko ten człon! Znacznie dokładniej (kasują się błędy numeryczne E_single).
    
    pass

def compute_interaction_direct(D, box_L=12.0, pts=60):
    x = np.linspace(-box_L, box_L, pts)
    # Centrum układu w D/2
    z = np.linspace(-box_L + D/2, box_L + D/2, pts) 
    y = np.linspace(-box_L, box_L, pts)
    dx = x[1] - x[0]
    
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

    # Kwaterniony
    def get_quat(xx, yy, zz):
        r = np.sqrt(xx**2 + yy**2 + zz**2) + 1e-9
        # Profil z ogonem 1/r^2 (f ~ 1/r)
        # UWAGA: f ~ 1/r daje sin(f) ~ 1/r. 
        # Pole pi ~ n sin(f) ~ 1/r.
        # Grad pi ~ 1/r^2. E_dens ~ 1/r^4. -> E ~ 1/r (Coulomb).
        f = 2.0 * np.arctan(1.0/r)
        sf = np.sin(f)
        cf = np.cos(f)
        return cf, (xx/r)*sf, (yy/r)*sf, (zz/r)*sf

    # Q1 w 0, Q2 w D
    s1, v1x, v1y, v1z = get_quat(X, Y, Z) # Q1 at 0
    s2, v2x, v2y, v2z = get_quat(X, Y, Z - D) # Q2 at D
    
    # Gradienty Q1 (g1) i Q2 (g2)
    # g1[axis][component]
    def grad_q(s, vx, vy, vz):
        gs = np.gradient(s, dx)
        gvx = np.gradient(vx, dx)
        gvy = np.gradient(vy, dx)
        gvz = np.gradient(vz, dx)
        return gs, gvx, gvy, gvz # Każdy to [dx, dy, dz]

    g1s, g1x, g1y, g1z = grad_q(s1, v1x, v1y, v1z)
    g2s, g2x, g2y, g2z = grad_q(s2, v2x, v2y, v2z)
    
    # Interaction Term Density (Cross term)
    # Formula: 2 * (grad q1 . grad q2)_R4 inner product?
    # Nie do końca. Algebra kwaternionów jest niekomutatywna.
    # |grad(q1 q2)|^2 = | (grad q1) q2 + q1 (grad q2) |^2
    # = <g1 q2 + q1 g2, g1 q2 + q1 g2>
    # = |g1 q2|^2 + |q1 g2|^2 + 2 <g1 q2, q1 g2>
    # Ponieważ |q|=1, |g1 q2| = |g1|.
    # Więc Int = 2 <g1 q2, q1 g2>
    # <A, B> to iloczyn skalarny w R4 (suma składowych).
    
    # Musimy obliczyć g1*q2 i q1*g2 (mnożenie kwaternionów) dla każdej składowej gradientu (dx, dy, dz).
    E_int_tot = 0
    
    for axis in range(3): # dx, dy, dz
        # Wektor g1 dla danej osi różniczkowania
        da_s1, da_x1, da_y1, da_z1 = g1s[axis], g1x[axis], g1y[axis], g1z[axis]
        da_s2, da_x2, da_y2, da_z2 = g2s[axis], g2x[axis], g2y[axis], g2z[axis]
        
        # Mnożenie kwaternionów A * B
        def quat_mul(sa, xa, ya, za, sb, xb, yb, zb):
            st = sa*sb - (xa*xb + ya*yb + za*zb)
            xt = sa*xb + sb*xa + (ya*zb - za*yb)
            yt = sa*yb + sb*ya + (za*xb - xa*zb)
            zt = sa*zb + sb*za + (xa*yb - ya*xb)
            return st, xt, yt, zt
            
        # Term 1: (grad1) * q2
        t1s, t1x, t1y, t1z = quat_mul(da_s1, da_x1, da_y1, da_z1, s2, v2x, v2y, v2z)
        
        # Term 2: q1 * (grad2)
        t2s, t2x, t2y, t2z = quat_mul(s1, v1x, v1y, v1z, da_s2, da_x2, da_y2, da_z2)
        
        # Iloczyn skalarny <T1, T2> w R4
        dot = t1s*t2s + t1x*t2x + t1y*t2y + t1z*t2z
        
        E_int_tot += np.sum(2.0 * dot) * (dx**3)
        
    return E_int_tot

# Skanowanie D
d_vals = np.linspace(3.0, 8.0, 6)
e_vals = []

log(f"{'D':<10} {'E_int':<15} {'E*D':<10} {'E*D^3':<10}")
log("-" * 50)

for d in d_vals:
    e = compute_interaction_direct(d)
    e_vals.append(e)
    log(f"{d:<10.2f} {e:<15.6f} {e*d:<10.4f} {e*d**3:<10.4f}")

# Wnioski
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
    # Sprawdź dipol
    prod_ED3 = evals * d_vals**3
    cv3 = np.std(prod_ED3) / np.mean(prod_ED3)
    if cv3 < 0.1:
         log("Prawo dipolowe (1/r^3).")
    else:
         log("Prawo złożone.")

# Zapis raportu
with open(REPORT_FILE, 'w') as f:
    f.write("# QW-1303: Vector Skyrmion Coulomb Check\n\n")
    f.write(f"**Data:** {datetime.now()}\n\n")
    f.write("```\n")
    f.write("\n".join(md))
    f.write("\n```\n")

print(f"Raport zapisano w {REPORT_FILE}")
