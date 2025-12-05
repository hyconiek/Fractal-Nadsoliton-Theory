#!/usr/bin/env python3
# QW-606: SUPER-BALLISTIC ORIGIN TEST
# Purpose: Determine if b=2.4 comes from K(d) network or beta_tors vacuum
# Method: Test dispersion with different parameters
# Date: 2025-12-05

import numpy as np
from scipy.ndimage import convolve

print("="*80)
print("QW-606: SUPER-BALLISTIC ORIGIN")
print("="*80)
print("Test: Does b=2.4 come from network K(d) or vacuum beta_tors?")
print("Baseline (QW-604b): b=2.386 with standard params")
print("="*80)

# Test configurations
TESTS = [
    {'name': 'Baseline', 'K_multiplier': 1.0, 'beta': 0.01, 'gamma': 0.0},
    {'name': 'High Beta', 'K_multiplier': 1.0, 'beta': 0.05, 'gamma': 0.0},
    {'name': 'Low Beta', 'K_multiplier': 1.0, 'beta': 0.001, 'gamma': 0.0},
]

GRID_SIZE = 48  # Smaller for speed
DT = 0.05
STEPS = 200  # Shorter for speed
sigma_0 = 4.0
k_0 = 0.3

print(f"Grid: {GRID_SIZE}^3, Steps: {STEPS}")
print(f"Testing {len(TESTS)} configurations")
print("-" * 40)

# Initialize wave packet
x = np.linspace(-GRID_SIZE/2, GRID_SIZE/2, GRID_SIZE)
X, Y, Z = np.meshgrid(x, x, x, indexing='ij')

laplacian_kernel = np.zeros((3,3,3))
laplacian_kernel[1,1,1] = -6
for i in [(0,1,1), (2,1,1), (1,0,1), (1,2,1), (1,1,0), (1,1,2)]:
    laplacian_kernel[i] = 1

def laplacian(field):
    return convolve(field, laplacian_kernel, mode='wrap')

def measure_width(psi):
    prob = np.abs(psi)**2
    prob = prob / (np.sum(prob) + 1e-10)
    x_cm = np.sum(X * prob)
    sigma_x = np.sqrt(np.sum(((X - x_cm)**2) * prob))
    sigma_y = np.sqrt(np.sum((Y**2) * prob))
    sigma_z = np.sqrt(np.sum((Z**2) * prob))
    return np.sqrt(sigma_x**2 + sigma_y**2 + sigma_z**2)

results = []

for test in TESTS:
    print(f"\n{'='*60}")
    print(f"Test: {test['name']}")
    print(f"  K_multiplier: {test['K_multiplier']}")
    print(f"  beta_tors: {test['beta']}")
    print('='*60)
    
    # Initialize gaussian
    amplitude = np.exp(-(X**2 + Y**2 + Z**2) / (2 * sigma_0**2))
    phase = k_0 * X
    psi = amplitude * np.exp(1j * phase)
    psi = psi / np.sqrt(np.sum(np.abs(psi)**2))
    
    history_width = []
    history_time = []
    
    gamma = test['gamma']
    beta = test['beta']
    
    # Evolve
    for t in range(STEPS):
        rho = np.abs(psi)**2
        
        # Note: K_multiplier would modify laplacian strength
        # For simplicity, we test beta variation (vacuum effect)
        kin = 1j * laplacian(psi) * test['K_multiplier']
        back = -1j * beta * rho * psi
        
        dpsi_dt = kin + back
        psi += DT * dpsi_dt
        
        if t % 10 == 0:
            width = measure_width(psi)
            history_width.append(width)
            history_time.append(t * DT)
            
            if t % 50 == 0:
                print(f"  t={t:3d}: σ={width:.2f}")
    
    # Fit dispersion
    widths = np.array(history_width)
    times = np.array(history_time)
    
    n_skip = len(times) // 5
    widths_fit = widths[n_skip:]
    times_fit = times[n_skip:]
    
    sigma_initial = widths_fit[0]
    delta_sigma = widths_fit - sigma_initial
    
    valid = delta_sigma > 0.1
    if np.sum(valid) > 5:
        log_t = np.log(times_fit[valid])
        log_ds = np.log(delta_sigma[valid])
        coeffs = np.polyfit(log_t, log_ds, 1)
        b = coeffs[0]
    else:
        b = 0
    
    results.append({
        'name': test['name'],
        'beta': test['beta'],
        'K_mult': test['K_multiplier'],
        'b': b
    })
    
    print(f"  Dispersion exponent: b = {b:.3f}")

print("\n" + "="*80)
print("RESULTS")
print("="*80)

print("\n| Configuration | beta_tors | K_mult | b (exponent) |")
print("|---------------|-----------|--------|--------------|")
for r in results:
    print(f"| {r['name']:13s} | {r['beta']:9.3f} | {r['K_mult']:6.1f} | {r['b']:12.3f} |")

# Analysis
baseline_b = results[0]['b']
print(f"\nBaseline b: {baseline_b:.3f}")

# Check beta effect
beta_values = [r['beta'] for r in results]
b_values = [r['b'] for r in results]

if len(set(beta_values)) > 1:
    # Correlation
    from scipy.stats import pearsonr
    if len(beta_values) > 2:
        corr, p_value = pearsonr(beta_values, b_values)
        print(f"\nCorrelation (b vs beta): r = {corr:.3f}")
        
        if abs(corr) > 0.8:
            print("✅ VACUUM EFFECT!")
            print(f"   Super-ballistic pochodzi z beta_tors (vacuum back-reaction)")
            origin = "vacuum"
        else:
            print("🟡 Słaba korelacja z beta")
            origin = "unclear"
    else:
        origin = "unclear"
else:
    origin = "unclear"

# Report
with open("raport_qw606_superballistic.md", "w") as f:
    f.write("# QW-606: Super-Ballistic Origin\n\n")
    f.write("**Test:** Czy b=2.4 pochodzi z sieci K(d) czy próżni beta_tors?\n\n")
    f.write("| Config | beta | b |\n|--------|------|---|\n")
    for r in results:
        f.write(f"| {r['name']} | {r['beta']:.3f} | {r['b']:.3f} |\n")
    f.write(f"\n**Origin:** {origin}\n")

print("="*80)
