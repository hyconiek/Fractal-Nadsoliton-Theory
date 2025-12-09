#!/usr/bin/env python3
"""
QW-697: K(d) MAXIMA VS STABLE ORBITS (EMERGENT OBSERVER)
=========================================================
Purpose: Verify that K(d) maxima match the previously derived stable orbits
         d₁=1.33, d₂=9.33, d₃=17.33 (FINAL_THEORY_REPORT.md)

Background:
  - QW-696: Forces correlate with K(d) (r=0.97)
  - Previous work: Stable orbits at d₁=1.33, d₂=9.33, d₃=17.33
  - Koide: 0.03% error with M ∝ d^α

Question: Do K(d) maxima (where K > 0, strongest attraction) coincide
          with the stable orbit positions?

Emergent Observer Perspective:
  - Observer at scale d sees physics at that scale
  - Particles exist where K(d) > 0 (attraction/resonance)
  - The specific d values define what observer "sees" as a particle
"""

import numpy as np
import datetime

print("="*80)
print("QW-697: K(d) MAXIMA VS STABLE ORBITS (EMERGENT OBSERVER)")
print("="*80)

# --- K(d) Parameters ---
ALPHA = 4 * np.log(2)  # 2.77
OMEGA = np.pi / 4
PHI = np.pi / 6
BETA = 0.1

def K(d):
    """The fundamental coupling kernel."""
    return ALPHA * np.cos(OMEGA * d + PHI) / (1 + BETA * d)

# --- Known stable orbits from FINAL_THEORY_REPORT ---
D_ORBITS_THEORY = [1.33, 9.33, 17.33]
print(f"\nTheoretical Stable Orbits: {D_ORBITS_THEORY}")

# --- Step 1: Find K(d) Maxima ---
print("\n[1] Finding K(d) Maxima...")

d_range = np.linspace(0.1, 30, 1000)
K_values = [K(d) for d in d_range]

# Find local maxima
from scipy.signal import find_peaks

peaks_idx, _ = find_peaks(K_values)
d_maxima = d_range[peaks_idx]
K_maxima = [K_values[i] for i in peaks_idx]

print(f"K(d) Maxima found at d = {d_maxima}")

# --- Step 2: Compare with Stable Orbits ---
print("\n[2] Comparing K(d) Maxima with Stable Orbits...")

print("| Stable Orbit d | Nearest K(d) Maximum | Difference |")
print("|----------------|----------------------|------------|")

matches = []
for d_orbit in D_ORBITS_THEORY:
    # Find nearest K(d) maximum
    if len(d_maxima) > 0:
        nearest_idx = np.argmin(np.abs(d_maxima - d_orbit))
        d_nearest = d_maxima[nearest_idx]
        diff = abs(d_orbit - d_nearest)
        matches.append((d_orbit, d_nearest, diff))
        print(f"| {d_orbit:.2f} | {d_nearest:.2f} | {diff:.2f} |")
    else:
        print(f"| {d_orbit:.2f} | N/A | N/A |")

# --- Step 3: Calculate K(d) at Orbit Positions ---
print("\n[3] K(d) Values at Stable Orbits...")

print("| Orbit | d | K(d) | Interpretation |")
print("|-------|---|------|----------------|")
for i, d_orbit in enumerate(D_ORBITS_THEORY):
    K_val = K(d_orbit)
    interp = "Attraction ✅" if K_val > 0 else "Repulsion ❌"
    particle = ["Electron", "Muon", "Tau"][i]
    print(f"| {particle} | {d_orbit:.2f} | {K_val:+.4f} | {interp} |")

# --- Step 4: Mass Scaling Test ---
print("\n[4] Mass Scaling: M ∝ d^α...")

ALPHA_GEO = 4 * np.log(2)
M_0 = 0.511  # Electron mass in MeV (reference)

masses_theory = [M_0 * (d / D_ORBITS_THEORY[0])**ALPHA_GEO for d in D_ORBITS_THEORY]
masses_exp = [0.511, 105.7, 1777]  # e, μ, τ in MeV

print("| Lepton | d | M_theory (MeV) | M_exp (MeV) | Error |")
print("|--------|---|----------------|-------------|-------|")
for i, name in enumerate(["e", "μ", "τ"]):
    err = abs(masses_theory[i] - masses_exp[i]) / masses_exp[i] * 100
    print(f"| {name} | {D_ORBITS_THEORY[i]:.2f} | {masses_theory[i]:.1f} | {masses_exp[i]:.1f} | {err:.1f}% |")

# --- Step 5: Emergent Observer Interpretation ---
print("\n[5] Emergent Observer Interpretation...")

print("""
In the Emergent Observer framework (QW-684/687):
- An observer at scale d_obs "sees" particles as local excitations of K(d).
- The observer WITHIN a particle (d_obs ≈ d_particle) sees quantum behavior (S > 0).
- The observer OUTSIDE (d_obs >> d_particle) sees classical behavior (S → 0).

For the K(d) force:
- Particles exist at K(d) MAXIMA (attraction).
- The K(d) value at d determines coupling strength.
- Larger d = larger particle = weaker K(d) = heavier mass.
""")

# --- Step 6: Verdict ---
print("="*80)
print("VERDICT")
print("="*80)

# Check if K(d) > 0 at all orbits
all_attractive = all(K(d) > 0 for d in D_ORBITS_THEORY)
mean_diff = np.mean([m[2] for m in matches]) if matches else float('inf')

if all_attractive and mean_diff < 2:
    verdict = "✅ STRONG MATCH"
    explanation = f"K(d) > 0 at ALL orbits. Mean deviation: {mean_diff:.2f}"
elif all_attractive:
    verdict = "🟡 PARTIAL MATCH"
    explanation = f"K(d) > 0 at all orbits, but positions differ by {mean_diff:.2f}"
else:
    verdict = "❌ MISMATCH"
    explanation = "Some orbits have K(d) < 0 (repulsion)."

print(f"\n{verdict}")
print(f"Explanation: {explanation}")

# --- Report ---
report_file = "raport_qw697_Kd_maxima_orbits.md"
print(f"\nSaving report to {report_file}...")

with open(report_file, "w") as f:
    f.write("# RAPORT QW-697: K(d) MAXIMA VS STABLE ORBITS\n")
    f.write(f"**Date:** {datetime.datetime.now()}\n\n")
    
    f.write("## 1. Question\n")
    f.write("Do K(d) maxima match the stable orbits d₁=1.33, d₂=9.33, d₃=17.33?\n\n")
    
    f.write("## 2. K(d) Maxima\n")
    f.write(f"Found at: {[f'{d:.2f}' for d in d_maxima]}\n\n")
    
    f.write("## 3. Comparison\n")
    f.write("| Stable Orbit | K(d) Maximum | Difference |\n")
    f.write("|--------------|--------------|------------|\n")
    for d_orbit, d_near, diff in matches:
        f.write(f"| {d_orbit:.2f} | {d_near:.2f} | {diff:.2f} |\n")
    f.write("\n")
    
    f.write("## 4. K(d) at Orbits\n")
    f.write("| Lepton | d | K(d) | Sign |\n")
    f.write("|--------|---|------|------|\n")
    for i, name in enumerate(["Electron", "Muon", "Tau"]):
        K_val = K(D_ORBITS_THEORY[i])
        sign = "+" if K_val > 0 else "-"
        f.write(f"| {name} | {D_ORBITS_THEORY[i]:.2f} | {K_val:+.4f} | {sign} |\n")
    f.write("\n")
    
    f.write("## 5. Verdict\n")
    f.write(f"### {verdict}\n")
    f.write(f"{explanation}\n\n")
    
    f.write("## 6. Emergent Observer Insight\n")
    f.write("Particles exist at K(d) maxima. The Emergent Observer at that scale\n")
    f.write("sees the particle as a stable resonance. Observers at other scales\n")
    f.write("see the particle as classical or as part of a larger structure.\n")

print("Done.")
