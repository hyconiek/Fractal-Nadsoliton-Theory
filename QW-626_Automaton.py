#!/usr/bin/env python3
# QW-626: 4-BIT AUTOMATON GENESIS
# Purpose: Test if simple 4-bit Cellular Automaton rules can generate
#          ordered structures (like 12-octave lattice) from random start.
# Hypothesis: The 12-octave structure is an attractor of 4-bit dynamics.
# Model: 1D Chain of 4-bit cells (Values 0-15).
# Rule: Update state based on neighbors (e.g. sum mod 16, XOR, etc.)
# Date: 2025-12-05

import numpy as np
import matplotlib.pyplot as plt

print("="*80)
print("QW-626: 4-BIT AUTOMATON GENESIS")
print("="*80)
print("Test: Czy dynamika 4-bitowa generuje struktury 12-elementowe?")
print("Hypothesis: Emergence of periodic structures length 12.")
print("="*80)

# Parameters
width = 100
steps = 200
states = 16 # 4 bits

# Initialize random state
grid = np.random.randint(0, states, size=(steps, width))

# Define Rules to test
# Rule A: Diffusion/Sum (Linear) - Expect spreading
# Rule B: XOR (Non-linear) - Expect fractals (Sierpinski)
# Rule C: Feedback (Difference) - Expect waves?

def run_automaton(rule_func, name):
    print(f"\nRunning Rule: {name}")
    g = np.zeros((steps, width), dtype=int)
    g[0] = np.random.randint(0, states, size=width)
    
    for t in range(steps-1):
        for x in range(width):
            # Neighbors (pbc)
            L = g[t, (x-1)%width]
            C = g[t, x]
            R = g[t, (x+1)%width]
            
            g[t+1, x] = rule_func(L, C, R) % states
            
    # Analyze periodicity in space
    # Look for domain sizes.
    final_state = g[-1]
    
    # Simple FFT to find dominant spatial frequencies
    fft_vals = np.abs(np.fft.fft(final_state))
    freqs = np.fft.fftfreq(width)
    
    peak_idx = np.argmax(fft_vals[1:width//2]) + 1
    dominant_period = 1.0 / freqs[peak_idx]
    
    print(f"  Dominant Spatial Period: {dominant_period:.2f}")
    if abs(dominant_period - 12.0) < 1.0:
        print("  ✅ MATCHES 12-OCTAVE STRUCTURE!")
    else:
        print(f"  Result: Period {dominant_period:.2f} (Target: 12)")
        
    return dominant_period

# Rules
def rule_diff(L, C, R):
    # Reaction-Diffusion like?
    return (L + R - C) 

def rule_xor(L, C, R):
    # Chaotic/Fractal
    return L ^ R

def rule_complex(L, C, R):
    # Trying to generate 12 using 4-bit logic
    # 4 bits = 0..15.
    # Maybe interaction involves mod 12? No, must be mod 16 naturally.
    # Let's try interaction that favors 3 or 4 neighbors.
    # Try: (L + R) * C ? nolinear
    return (L + R) + (C * 2) 

# Run tests
p1 = run_automaton(rule_diff, "Difference (Wave)")
p2 = run_automaton(rule_xor, "XOR (Fractal)")
p3 = run_automaton(rule_complex, "Weighted Sum")

print("\nConclusion:")
if any(abs(p - 12.0) < 1.5 for p in [p1, p2, p3]):
    print("✅ 4-BIT DYNAMICS CAN GENERATE 12-ELEMENT CYCLE.")
    print("   This supports the idea that 12-octave lattice emerges from 4-bit logic.")
else:
    print("❌ NO SPONTANEOUS 12-STRUCTURE from simple rules.")
    print("   The 12-structure might require specific topology (loops), not just CA rules.")

# ============================================================================
# REPORT
# ============================================================================
with open("raport_qw626_automaton.md", "w") as f:
    f.write("# Raport QW-626: 4-Bit Automaton\n")
    f.write("**Data:** 2025-12-05\n\n")
    f.write("## Wyniki Symulacji\n")
    f.write(f"- Rule Difference: Period {p1:.2f}\n")
    f.write(f"- Rule XOR: Period {p2:.2f}\n")
    f.write(f"- Rule Weighted: Period {p3:.2f}\n\n")
    
    if any(abs(p - 12.0) < 1.0 for p in [p1, p2, p3]):
        f.write("### ✅ SUKCES\n")
        f.write("Znaleziono regułę generującą okres ~12. Struktura 12-oktawowa może być atraktorem 4-bitowym.\n")
    else:
        f.write("### ❌ PORAŻKA\n")
        f.write("Nie uzyskano okresu 12 z prostych reguł CA.\n")

print("Report saved.")
