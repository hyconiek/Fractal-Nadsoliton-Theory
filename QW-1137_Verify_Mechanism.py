
import numpy as np

print('VERIFICATION OF EMERGENCE MECHANISM')
print('=' * 60)

# 1. Kernel Parameters
ALPHA = 4 * np.log(2)
OMEGA = np.pi / 4
print(f'Kernel Parameters:')
print(f'  alpha = {ALPHA:.5f} (4*ln2)')
print(f'  omega = {OMEGA:.5f} (pi/4)')

# 2. Gravitational Exponent n
# Theoretical: n = 2 + omega/pi
n_theo = 2 + OMEGA / np.pi
# Experimental (QW-722): n = 2.26
n_exp = 2.26

print(f'\nGravitational Exponent n:')
print(f'  Theoretical (n = 2 + omega/pi): {n_theo:.5f}') 
print(f'  Simulation QW-722: {n_exp:.5f}')

# 3. Mass Exponent gamma
# gamma = -(3 - 2n) = 2n - 3
gamma_theo = abs(3 - 2 * n_theo)
gamma_exp = abs(3 - 2 * n_exp)

print(f'\nMass Exponent gamma (M ~ r^-gamma):')
print(f'  Theoretical (gamma = 2n - 3): {gamma_theo:.5f}')
print(f'  Simulation (gamma = 2n - 3):   {gamma_exp:.5f}')

# 4. Mass Predictions
# M(d) = M_ref * 4^(-gamma * d)
particles = {
    'Top': 0.0,
    'Bottom': 1.75,
    'Tau': 2.25,
    'Muon': 3.5,
    'Electron': 6.0
}
m_ref = 172760 # MeV (Top)
m_actual = {
    'Top': 172760,
    'Bottom': 4180,
    'Tau': 1777,
    'Muon': 105.7,
    'Electron': 0.511
}

print(f'\nMass Predictions (M_ref = {m_ref} MeV):')
print(f'Particle   | d    | Actual     | Theo (g=1.5) | Sim (g=1.52) | Err(Theo) | Err(Sim)')
print('-' * 90)

err_theo_sum = 0
err_sim_sum = 0

for p, d in particles.items():
    m_act = m_actual[p]
    
    # Theoretical (gamma = 1.5)
    m_theo = m_ref * 4 ** (-gamma_theo * d)
    err_theo = abs(m_theo - m_act) / m_act * 100
    err_theo_sum += err_theo
    
    # Simulation (gamma = 1.52)
    m_sim = m_ref * 4 ** (-gamma_exp * d)
    err_sim = abs(m_sim - m_act) / m_act * 100
    err_sim_sum += err_sim
    
    print(f'{p:<10} | {d:<4} | {m_act:<10.1f} | {m_theo:<12.1f} | {m_sim:<12.1f} | {err_theo:6.1f}% | {err_sim:6.1f}%')

print('-' * 90)
print(f'AVERAGE ERROR:                                   | {err_theo_sum/5:6.1f}% | {err_sim_sum/5:6.1f}%')

if err_sim_sum < err_theo_sum:
     print('\nCONCLUSION: Empirical gamma=1.52 works better.')
else:
     print('\nCONCLUSION: Theoretical gamma=1.50 works better.')
