
import numpy as np
import matplotlib.pyplot as plt

print('=' * 80)
print('QW-1160: THEORY VERIFICATION & SENSITIVITY ANALYSIS')
print('=' * 80)

# 1. SETUP DATA
particles = {
    'Top':      172760,
    'Bottom':   4180,
    'Tau':      1777,
    'Charm':    1275,
    'Muon':     105.7,
    'Strange':  95,
    'Down':     4.7,
    'Up':       2.2,
    'Electron': 0.511
}

M_REF = particles['Top']

def get_Q_errors(gamma, particles_dict):
    """
    Calculate how close the Q values are to integers for a given Gamma.
    Metric: Mean Squared Deviation from nearest integer.
    """
    total_sq_diff = 0
    count = 0
    
    # Q = -4 * log4(M/M_ref) / gamma
    #   = -4 * (ln(M/M_ref)/ln(4)) / gamma
    
    for name, m in particles_dict.items():
        if name == 'Top': continue # Always 0
        
        # Calculate Q
        val = -4 * np.log(m/M_REF) / (np.log(4) * gamma)
        
        # Distance to nearest integer
        nearest = round(val)
        diff = val - nearest
        total_sq_diff += diff**2
        count += 1
        
    return np.sqrt(total_sq_diff / count)

# 2. SENSITIVITY SCAN (GAMMA)
print('\n[1] SCANNING GAMMA PARAMETER (Standard: 1.52)')
gammas = np.linspace(1.0, 2.0, 1000)
errors = []

for g in gammas:
    errors.append(get_Q_errors(g, particles))

min_idx = np.argmin(errors)
best_gamma = gammas[min_idx]
min_error = errors[min_idx]

print(f'  Optimal Gamma found: {best_gamma:.4f}')
print(f'  Standard Gamma (1.52): {get_Q_errors(1.52, particles):.4f}')
print(f'  Improvement: {100*(get_Q_errors(1.52, particles)-min_error)/get_Q_errors(1.52, particles):.1f}%')

# 3. ROBUSTNESS CHECK
print('\n[2] CHECKING INTEGER LOCKING AT OPTIMAL GAMMA')
print(f'  Particle analysis using Gamma = {best_gamma:.4f}')
print(f'  {"Particle":<10} | {"M_exp":<10} | {"Q_calc":<8} | {"Int":<4} | {"Diff":<6}')
print('-' * 60)

qs = {}
for p, m in particles.items():
    q_val = -4 * np.log(m/M_REF) / (np.log(4) * best_gamma)
    qs[p] = q_val
    print(f'  {p:<10} | {m:<10.1f} | {q_val:<8.3f} | {round(q_val):<4} | {q_val-round(q_val):+.3f}')

# 4. RUNNING COUPLING HYPOTHESIS
print('\n[3] RUNNING COUPLING TEST')
print('  Hypothesis: Gamma is not constant, but scales with Energy (Mass).')
print('  Gamma(E) = Gamma_0 + b * log(E/M_top)')

# Simple fit for Gamma(E) that forces Q to be EXACTLY integers?
# Q_int = -4 * log4(M/M_top) / Gamma
# Gamma_req = -4 * log4(M/M_top) / Q_int

print(f'  {"Particle":<10} | {"Q_target":<8} | {"Gamma_req":<8} | {"Scale In(M/M_top)"}')
print('-' * 70)

gamma_reqs = []
scales = []
valid_particles = []

target_map = {
    'Top': 0, 'Bottom': 7, 'Tau': 9, 'Charm': 10, 
    'Muon': 14, 'Strange': 14, 'Down': 20, 'Up': 21, 'Electron': 24
}

for p, m in particles.items():
    if p == 'Top': continue
    q_target = target_map[p]
    
    # Calculate required gamma
    # Gamma * Q = -4 * log4(M/M_top)
    term = -4 * np.log(m/M_REF) / np.log(4)
    g_req = term / q_target
    
    scale = np.log(m/M_REF)
    
    gamma_reqs.append(g_req)
    scales.append(scale)
    valid_particles.append(p)
    
    print(f'  {p:<10} | {q_target:<8} | {g_req:<8.4f} | {scale:<8.2f}')


# Linear regression on Gamma vs Log(Scale)
# Gamma(scale) = slope * scale + intercept
coeffs = np.polyfit(scales, gamma_reqs, 1)
slope = coeffs[0]
intercept = coeffs[1]

print('\n  FIT RESULTS:')
print(f'  Gamma(E) = {intercept:.4f} + ({slope:.5f}) * ln(M/M_top)')
print(f'  Base Gamma (at M_top): {intercept:.4f} (Compare to 1.52)')

# 5. VERDICT
print('\n[4] VERIFICATION VERDICT')
print('-' * 80)
if abs(best_gamma - 1.52) < 0.1:
    print('  [PASS] Gamma factor 1.52 is close to optimal global fit.')
else:
    print(f'  [WARN] Optimal Gamma {best_gamma:.4f} deviates from 1.52.')

if abs(slope) < 0.05:
    print('  [PASS] ' + f'Running of Gamma is weak (slope={slope:.4f}), constant approximation holds.')
else:
    print('  [NOTE] ' + f'Significant running of Gamma detected (slope={slope:.4f}). High energy physics needed.')

print('=' * 80)
