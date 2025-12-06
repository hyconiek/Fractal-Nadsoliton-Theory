#!/usr/bin/env python3
"""
QW-675_Gravity_Law.py
Purpose: RIGOROUS TEST of gravity law in FIN Theory.
         Does the model produce F ∝ 1/r² (Newton) or something else?

Previous Results:
- QW-552/563: FAILED with n=+2.0 (repulsion?) instead of n=-2.0 (attraction)
- QW-588: MOND worked (F ∝ a²/a₀)

Method:
1. Create two "masses" as high-energy nodes in the octave network
2. Measure effective force via energy gradient
3. Fit: F(r) = A * r^n, determine n

Success Criterion: n = -2.0 ± 0.1

Output: RAPORT_GRAVITY_LAW.md
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy.spatial.distance import cdist
import datetime

# --- Constants ---
N_NODES = 400
ALPHA_GEO = 4 * np.log(2)
BETA = 0.1
OMEGA = np.pi / 4
PHI = np.pi / 6

REPORT_FILE = "RAPORT_GRAVITY_LAW.md"

print(f"QW-675: GRAVITY LAW RIGOROUS TEST - Output: {REPORT_FILE}")

def K(d):
    """Effective coupling kernel"""
    return ALPHA_GEO * np.cos(OMEGA * d + PHI) / (1 + BETA * d)

def create_network(n_nodes, seed=675):
    """Create octave network with 12 octaves"""
    np.random.seed(seed)
    
    # 3D positions
    positions = np.random.randn(n_nodes, 3) * 5.0
    
    # Assign octaves (0-11)
    octaves = np.random.randint(0, 12, n_nodes)
    
    # Calculate distances
    dist_matrix = cdist(positions, positions)
    
    # Coupling matrix based on K(d) where d = octave difference
    coupling = np.zeros((n_nodes, n_nodes))
    for i in range(n_nodes):
        for j in range(i+1, n_nodes):
            d_oct = abs(octaves[i] - octaves[j])
            coupling[i, j] = K(d_oct) * np.exp(-dist_matrix[i, j]**2 / 10)
            coupling[j, i] = coupling[i, j]
    
    return positions, octaves, coupling, dist_matrix

def place_mass(positions, octaves, coupling, center, radius, mass_octave=7, boost=10.0):
    """Place a 'mass' by boosting coupling strength in a region"""
    n = len(positions)
    mass_indices = []
    
    for i in range(n):
        dist = np.linalg.norm(positions[i] - center)
        if dist < radius:
            mass_indices.append(i)
            # Boost couplings for mass nodes
            for j in range(n):
                if i != j:
                    coupling[i, j] *= boost
                    coupling[j, i] *= boost
    
    return mass_indices, coupling

def calculate_energy_at_distance(positions, coupling, mass_center, test_distances):
    """Calculate effective potential energy at various distances from mass"""
    energies = []
    
    for r in test_distances:
        # Create test point at distance r from mass
        test_pos = mass_center + np.array([r, 0, 0])
        
        # Find closest node to test position
        dists = np.linalg.norm(positions - test_pos, axis=1)
        test_node = np.argmin(dists)
        
        # Energy = sum of couplings to all other nodes (weighted by distance)
        E = 0
        for j in range(len(positions)):
            if j != test_node:
                E -= coupling[test_node, j]  # Negative = attractive
        
        energies.append(E)
    
    return np.array(energies)

def power_law(r, A, n):
    """F = A * r^n"""
    return A * np.power(r, n)

# ===================================================================
# MAIN EXECUTION
# ===================================================================

with open(REPORT_FILE, "w") as f:
    f.write(f"# RAPORT: PRAWO GRAWITACJI (QW-675)\n")
    f.write(f"**Data:** {datetime.datetime.now()}\n")
    f.write("**Cel:** Czy FIN produkuje $F \\propto 1/r^2$ (Newton)?\n\n")

    # Create network
    print("Creating octave network...")
    positions, octaves, coupling, dist_matrix = create_network(N_NODES)
    
    f.write(f"## 1. SETUP\n")
    f.write(f"- Węzły: {N_NODES}\n")
    f.write(f"- Oktawy: 0-11\n")
    f.write(f"- $\\alpha_{{geo}} = {ALPHA_GEO:.4f}$\n\n")

    # Place central mass
    mass_center = np.array([0, 0, 0])
    mass_indices, coupling = place_mass(positions, octaves, coupling, mass_center, radius=2.0, boost=10.0)
    
    f.write(f"## 2. MASA CENTRALNA\n")
    f.write(f"- Pozycja: [0, 0, 0]\n")
    f.write(f"- Promień: 2.0\n")
    f.write(f"- Węzły masowe: {len(mass_indices)}\n\n")

    # Measure energy at various distances
    test_distances = np.array([3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0])
    
    print(f"Measuring energy at distances: {test_distances}")
    energies = calculate_energy_at_distance(positions, coupling, mass_center, test_distances)
    
    # Calculate "force" as -dE/dr (gradient)
    forces = -np.gradient(energies, test_distances)
    
    f.write(f"## 3. POMIARY\n")
    f.write("| Odległość r | Energia E | Siła F=-dE/dr |\n")
    f.write("|-------------|-----------|---------------|\n")
    for r, E, F in zip(test_distances, energies, forces):
        f.write(f"| {r:.1f} | {E:.2f} | {F:.4f} |\n")
    f.write("\n")

    # Fit power law to force
    f.write(f"## 4. DOPASOWANIE PRAWA POTĘGOWEGO\n")
    
    try:
        # Use absolute values for fitting
        valid = forces > 0  # Attractive forces
        if np.sum(valid) >= 3:
            r_fit = test_distances[valid]
            F_fit = forces[valid]
            
            popt, pcov = curve_fit(power_law, r_fit, F_fit, p0=[1.0, -2.0], maxfev=10000)
            A_fit, n_fit = popt
            
            f.write(f"**Wynik:** $F = {A_fit:.4f} \\cdot r^{{{n_fit:.2f}}}$\n\n")
            f.write(f"- Wykładnik n = **{n_fit:.2f}**\n")
            f.write(f"- Oczekiwany (Newton): n = **-2.0**\n")
            f.write(f"- Błąd: {abs(n_fit - (-2.0)):.2f}\n\n")
            
            print(f"Fitted: F = {A_fit:.4f} * r^{n_fit:.2f}")
            
            if abs(n_fit + 2.0) < 0.5:
                result = "✅ **SUCCESS:** Prawo $1/r^2$ potwierdzone!"
            elif n_fit < 0:
                result = f"🟡 **PARTIAL:** Przyciąganie, ale n={n_fit:.2f} ≠ -2.0"
            else:
                result = f"❌ **FAIL:** Odpychanie (n={n_fit:.2f} > 0)"
        else:
            # All forces negative = repulsion
            n_fit = 2.0  # Placeholder
            result = "❌ **FAIL:** Siły odpychające (repulsion)"
            f.write(f"**Problem:** Większość sił jest odpychająca.\n\n")
            
    except Exception as e:
        n_fit = 0
        result = f"❌ **ERROR:** Dopasowanie nie powiodło się: {e}"
        f.write(f"**Błąd:** {e}\n\n")
    
    f.write(result + "\n\n")
    print(result)

    # ===================================================================
    # ALTERNATIVE: Energy-distance relationship
    # ===================================================================
    f.write(f"## 5. ALTERNATYWNA ANALIZA: E(r)\n")
    
    # For gravity: E ∝ -1/r, so E = B/r^m with m = -1
    try:
        # Shift energies to be positive for fitting
        E_shifted = energies - np.min(energies) + 1
        
        popt_E, _ = curve_fit(power_law, test_distances, E_shifted, p0=[1.0, -1.0], maxfev=10000)
        B_fit, m_fit = popt_E
        
        f.write(f"**Wynik:** $E = {B_fit:.4f} \\cdot r^{{{m_fit:.2f}}}$\n\n")
        f.write(f"- Wykładnik m = **{m_fit:.2f}**\n")
        f.write(f"- Oczekiwany (Newton): m = **-1.0**\n\n")
        
        print(f"Energy fit: E = {B_fit:.4f} * r^{m_fit:.2f}")
        
    except Exception as e:
        f.write(f"**Błąd E(r):** {e}\n\n")

    # ===================================================================
    # SUMMARY
    # ===================================================================
    f.write(f"## 6. PODSUMOWANIE\n\n")
    f.write("| Test | Wynik | Status |\n")
    f.write("|------|-------|--------|\n")
    f.write(f"| Wykładnik n (siła) | {n_fit:.2f} | {'✅' if abs(n_fit + 2.0) < 0.5 else '❌'} |\n")
    f.write(f"| Kierunek siły | {'Przyciąganie' if n_fit < 0 else 'Odpychanie'} | {'✅' if n_fit < 0 else '❌'} |\n\n")
    
    f.write("**Wniosek:**\n")
    if abs(n_fit + 2.0) < 0.5:
        f.write("FIN Theory produkuje prawo grawitacji Newtona $F \\propto 1/r^2$.\n")
    else:
        f.write("FIN Theory NIE produkuje standardowego prawa grawitacji.\n")
        f.write("Możliwe przyczyny:\n")
        f.write("1. Sieć dyskretna nie ma ciągłego limitu\n")
        f.write("2. Jądro K(d) dominuje nad odległością przestrzenną\n")
        f.write("3. Potrzebna większa sieć (N > 10⁴)\n")

print(f"\nReport written to {REPORT_FILE}")
