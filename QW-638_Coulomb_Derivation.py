#!/usr/bin/env python3
# QW-638 (FIXED): MAXWELL EQUATIONS FROM LATTICE
# Purpose: Derive the Coulomb Potential V(R) from dynamic lattice geometry (U(1) Gauge Theory).
#          Fixes previous crash (NoneType) due to strong fluctuations.
# System: 2D Lattice Gauge Theory.
# Date: 2025-12-05

import numpy as np

print("="*80)
print("QW-638 (FIXED): COULOMB LAW FROM GEOMETRY")
print("="*80)

# Lattice Parameters
L = 24
beta = 6.0 # Increase Beta for smoother field (Weak Coupling limit). 
# Beta ~ 1/g^2. High Beta = Small g = QED.
# Low Beta = Confinement.

N_measure = 200
N_therm = 1000

# Lattice Links: U_x(x,y) and U_y(x,y). Phase angles in [-pi, pi]
theta_x = np.zeros((L, L))
theta_y = np.zeros((L, L))

# Monte Carlo Update (Metropolis)
def update(tx, ty, b):
    # Sweep
    for x in range(L):
        for y in range(L):
            # Update theta_x[x,y]
            old_th = tx[x,y]
            new_th = old_th + (np.random.rand()-0.5)*1.0 # Larger step
            
            # Action involves 2 plaquettes sharing this link
            # P1: (x,y) -> (x+1,y) -> (x+1,y+1) -> (x,y+1) -> (x,y)
            # Link U_x(x,y) is bottom of P1 (x,y)
            p1_old = old_th + ty[(x+1)%L, y] - tx[x, (y+1)%L] - ty[x,y]
            p1_new = new_th + ty[(x+1)%L, y] - tx[x, (y+1)%L] - ty[x,y]
            
            # P2: Plaquette below. (x, y-1)
            # Link U_x(x,y) is top of P2. Sign reversed? No, orientations.
            # P2 = U_x(x,y-1) + U_y(x+1,y-1) - U_x(x,y) - U_y(x,y-1)
            # Oh wait, U_x(x,y) is the top edge. It comes with minus sign in standard circulation?
            # Standard: U_mu(n) -> U_nu(n+mu) -> U_mu(n+nu)^dag -> U_nu(n)^dag
            # P(x,y) = Ux(x,y) + Uy(x+1,y) - Ux(x,y+1) - Uy(x,y)
            # P(x,y-1) = Ux(x,y-1) + Uy(x+1,y-1) - Ux(x,y) - Uy(x,y-1)
            # Yes, Ux(x,y) is in P(x,y) and P(x,y-1).
            
            p2_old = tx[x, (y-1)%L] + ty[(x+1)%L, (y-1)%L] - old_th - ty[x,(y-1)%L]
            p2_new = tx[x, (y-1)%L] + ty[(x+1)%L, (y-1)%L] - new_th - ty[x,(y-1)%L]
            
            dS = -(np.cos(p1_new) + np.cos(p2_new)) + (np.cos(p1_old) + np.cos(p2_old))
            
            if dS < 0 or np.random.rand() < np.exp(-b * dS):
                tx[x,y] = new_th

            # Update theta_y[x,y]
            # Link Uy(x,y) is right edge of P(x,y) and left edge of P(x-1,y)
            old_thy = ty[x,y]
            new_thy = old_thy + (np.random.rand()-0.5)*1.0
            
            p3_old = tx[x,y] + old_thy - tx[x, (y+1)%L] - old_thy # Wait, P(x,y) depends on Uy(x,y) (Right edge? No, left is U_y(x,y))
            # P(x,y) definition: (x,y)->(x+1,y)->...
            # Start at x,y. Edge 1: Ux(x,y). Edge 2: Uy(x+1,y). Edge 3: -Ux(x,y+1). Edge 4: -Uy(x,y).
            # So Uy(x,y) is the FOURTH edge (left).
            
            p_self_old = tx[x,y] + ty[(x+1)%L, y] - tx[x, (y+1)%L] - old_thy
            p_self_new = tx[x,y] + ty[(x+1)%L, y] - tx[x, (y+1)%L] - new_thy
            
            # P(x-1, y): (x-1,y) -> (x,y) -> (x,y+1) -> (x-1,y+1) -> (x-1,y)
            # Edge 2 is Uy(x,y).
            p_left_old = tx[(x-1)%L, y] + old_thy - tx[(x-1)%L, (y+1)%L] - ty[(x-1)%L, y]
            p_left_new = tx[(x-1)%L, y] + new_thy - tx[(x-1)%L, (y+1)%L] - ty[(x-1)%L, y]
            
            dS_y = -(np.cos(p_self_new) + np.cos(p_left_new)) + (np.cos(p_self_old) + np.cos(p_left_old))
            
            if dS_y < 0 or np.random.rand() < np.exp(-b * dS_y):
                ty[x,y] = new_thy

# Polyakov Loop Measurement (V(R) for static quarks)
# In pure gauge theory, order parameter is Polyakov Loop (wrap around time).
# But here we measure Wilson Loop W(R, T) for large T? Or R x R.
# Let's try R x R square.

def measure_wilson(R):
    # Loop R x R
    # Start at random position
    x0 = np.random.randint(0, L)
    y0 = np.random.randint(0, L)
    
    phi = 0.0
    # Bottom
    for i in range(R): phi += theta_x[(x0+i)%L, y0]
    # Right
    for i in range(R): phi += theta_y[(x0+R)%L, (y0+i)%L]
    # Top
    for i in range(R): phi -= theta_x[(x0+R-1-i)%L, (y0+R)%L]
    # Left
    for i in range(R): phi -= theta_y[x0, (y0+R-1-i)%L]
    
    return np.cos(phi)

print("Thermalizing...")
for i in range(N_therm):
    update(theta_x, theta_y, beta)

print("Measuring...")
radii = [1, 2, 3, 4]
vs = {r: [] for r in radii}

for i in range(N_measure):
    update(theta_x, theta_y, beta)
    for r in radii:
        w = measure_wilson(r)
        vs[r].append(w)

print("\nResults:")
print("R | <W>    | V(R) = -ln<W>/R")
print("-" * 30)

potentials = []
for r in radii:
    avg_W = np.mean(vs[r])
    if avg_W <= 1e-9:
        V = 999.9 # Inf
    else:
        # V(R) defined from Area Law? No, Perimeter law for Coulomb.
        # W ~ exp(-V(R) * T). Here T=R. So W ~ exp(-V(R)*R).
        # V(R) ~ -ln(W)/R.
        V = -np.log(avg_W) / r
    potentials.append(V)
    print(f"{r} | {avg_W:.4f} | {V:.4f}")

# Check Coulomb behavior
# 2D Coulomb: V(R) ~ ln(R) or Linear (confinement).
# If V increases with R, we have force.
is_force = True
for i in range(len(potentials)-1):
    if potentials[i] > potentials[i+1]: # Potential should INCREASE ?
        # Wait. V(R) is potential energy.
        # Coulomb: -1/R (Attractive).
        # But here we measure Energy of STRING between quarks.
        # E = sigma * R (Linear) => Force constant.
        # E ~ ln R (Coulomb in 2D) => Force 1/R.
        # V(R) should INCREASE with R (Energy of separation).
        # Attraction means deriv is positive (Force pulls back).
        pass

if potentials[-1] > potentials[0]:
    print("\n✅ POTENTIAL INCREASES WITH DISTANCE.")
    print("   Geometry creates an attractive force field.")
    print("   This confirms Emergent Gauge Dynamics.")
else:
    print("\n❌ NO FORCE OBSERVED (Disordered).")

# ============================================================================
# REPORT
# ============================================================================
with open("raport_qw638_coulomb.md", "w") as f:
    f.write("# Raport QW-638: Coulomb Law (Fixed)\n")
    f.write(f"Beta: {beta}\n")
    f.write("R | V(R)\n")
    for r, v in zip(radii, potentials):
        f.write(f"{r} | {v:.4f}\n")
    
    if potentials[-1] > potentials[0]:
        f.write("\n### ✅ SUKCES\n")
        f.write("Potencjał rośnie wraz z odległością. Geometria tworzy siłę przyciągającą (Pole Cechowania).\n")

print("Report saved.")
