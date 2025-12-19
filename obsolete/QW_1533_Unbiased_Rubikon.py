# OBSOLETE - Superceded by QW_1533_Rubikon_Final_Audit.py (Scientific Audit Round 3)
import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.interpolate import interp1d

# ============================================================
# QW-1533: UNBIASED RUBIKON TEST (MONTE CARLO SELECTION KERNEL)
# ============================================================

np.random.seed(42)

# ------------------------
# CONFIGURATION
# ------------------------
N_EVENTS = 25        # Increased for better statistics
TRUE_N = 1.00        # Set to 1.00 for GR Sanity Check
D_CALIB = 40.0
D_MIN, D_MAX = 40.0, 1500.0  # Larger volume for FIN reach

SIGMA_GW_REL = 0.12
SIGMA_CAL = 0.04
SIGMA_EM_REL = 0.05

SNR_THRESHOLD = 10.0 # Standard detection threshold
SNR_CALIB = 20.0     # Realistic SNR at D_CALIB (e.g. GW150914-like)

N_MCMC = 50000
BURNIN = 10000
STEP = 0.012

# ------------------------
# SNR KERNEL (PHYSICAL)
# ------------------------
def get_orientation_factor():
    """Simplified Finn's Theta factor for isotropic orientations."""
    cos_i = np.random.uniform(-1, 1)
    # Approximation of antenna pattern response
    Fp = 0.5 * (1 + cos_i**2)
    Fx = cos_i
    return np.sqrt(Fp**2 + Fx**2)

def get_chirp_mass():
    """Sample realistic chirp mass from log-normal distribution (approx BBH)."""
    # Mc mean ~ 30 M_sun, sigma ~ 0.3
    return np.random.lognormal(mean=np.log(30.0), sigma=0.3)

def calculate_snr(D, n, orientation, Mc):
    # SNR ~ (Mc/30)^(5/6) * (D_calib / D)^n * orientation
    # Reference: SNR_CALIB is defined at D_CALIB for Mc=30
    return SNR_CALIB * (Mc/30.0)**(5/6.0) * (D_CALIB / D)**n * orientation

# ------------------------
# SELECTION PROBABILITY P_det(n)
# ------------------------
print("Pre-computing Monte Carlo Selection Kernel (including masses)...")
N_SAMP_MC = 300000
mc_n_grid = np.linspace(0.4, 1.4, 40)
mc_p_det = []

for n_val in mc_n_grid:
    # Sample population from uniform volume
    d_samp = (np.random.uniform(0, 1, N_SAMP_MC) * (D_MAX**3 - D_MIN**3) + D_MIN**3)**(1/3)
    # Sample orientations and masses
    theta_samp = np.array([get_orientation_factor() for _ in range(N_SAMP_MC)])
    mc_samp = np.array([get_chirp_mass() for _ in range(N_SAMP_MC)])
    
    snrs = calculate_snr(d_samp, n_val, theta_samp, mc_samp)
    detected = snrs > SNR_THRESHOLD
    mc_p_det.append(np.mean(detected))

# Interpolator for fast lookup during MCMC
log_p_det_interp = interp1d(mc_n_grid, np.log(mc_p_det), kind='cubic', fill_value="extrapolate")

# ------------------------
# SIMULATION
# ------------------------
def simulate_unbiased_catalog(n_events, true_n):
    D_EM, D_GW, SIGMA = [], [], []
    
    total = 0
    while len(D_EM) < n_events:
        total += 1
        D_true = (np.random.uniform(0, 1) * (D_MAX**3 - D_MIN**3) + D_MIN**3)**(1/3)
        orientation = get_orientation_factor()
        Mc = get_chirp_mass()
        
        snr = calculate_snr(D_true, true_n, orientation, Mc)
        if snr < SNR_THRESHOLD:
            continue
            
        # Observables
        D_em = D_true * (1 + np.random.normal(0, SIGMA_EM_REL))
        D_gw = D_CALIB * (D_true / D_CALIB)**true_n
        
        # D_gw noise effectively reflects SNR
        D_gw *= (1 + np.random.normal(0, SIGMA_GW_REL))
        D_gw *= (1 + np.random.normal(0, SIGMA_CAL))
        
        sigma_tot = D_gw * np.sqrt(SIGMA_GW_REL**2 + SIGMA_CAL**2)
        
        D_EM.append(D_em)
        D_GW.append(D_gw)
        SIGMA.append(sigma_tot)
        
    efficiency = n_events / total
    print(f"Simulation efficiency: {efficiency:.3%}")
    return {"D_EM": np.array(D_EM), "D_GW": np.array(D_GW), "SIGMA": np.array(SIGMA), "efficiency": efficiency}

data = simulate_unbiased_catalog(N_EVENTS, TRUE_N)

# ------------------------
# HIERARCHICAL MODEL
# ------------------------
def log_prior(n):
    return 0.0 if 0.4 < n < 1.3 else -np.inf

def log_likelihood_events(n, data):
    D_pred = D_CALIB * (data["D_EM"]/D_CALIB)**n
    return -0.5 * np.sum(((data["D_GW"] - D_pred)/data["SIGMA"])**2)

def log_posterior(n, data):
    lp = log_prior(n)
    if not np.isfinite(lp): return -np.inf
    
    ll = log_likelihood_events(n, data)
    # Selection correction using MC Kernel
    ln_norm = N_EVENTS * log_p_det_interp(n)
    
    return lp + ll - ln_norm

# ------------------------
# MCMC
# ------------------------
print(f"Running MCMC (Injection n={TRUE_N})...")
chain = np.zeros(N_MCMC)
n_curr = 1.0 # Start at GR
logp_curr = log_posterior(n_curr, data)

accepted = 0
for i in range(N_MCMC):
    n_prop = n_curr + np.random.normal(0, STEP)
    logp_prop = log_posterior(n_prop, data)
    
    if np.log(np.random.rand()) < logp_prop - logp_curr:
        n_curr = n_prop
        logp_curr = logp_prop
        accepted += 1
    chain[i] = n_curr

print(f"Acceptance rate: {accepted/N_MCMC:.2f}")
chain = chain[BURNIN:]

# ------------------------
# RESULTS
# ------------------------
n_mean = np.mean(chain)
n_std = np.std(chain)
z_gr = abs(n_mean - 1.0) / n_std

print("\n=== QW-1533 REFINED RESULT ===")
print(f"Recovered n = {n_mean:.4f} ± {n_std:.4f}")
print(f"True n      = {TRUE_N}")
print(f"Bias        = {n_mean - TRUE_N:+.4f}")
print(f"Tension GR  = {z_gr:.2f} sigma")

# ------------------------
# SAVE REPORT
# ------------------------
timestamp = int(time.time())
scenario = "FIN" if TRUE_N < 0.9 else "GR"
report_name = f"QW_1533_Refined_Report_{scenario}_{timestamp}.md"
with open(report_name, "w") as f:
    f.write(f"# QW-1533 Refined Rubikon Report ({scenario})\n\n")
    f.write(f"**Method:** Monte Carlo Selection Kernel (Mass-Integrated SNR)\n")
    f.write(f"**Status:** PHYSICALLY MOTIVATED UNBIASED RECOVERY (SIMPLIFIED)\n\n")
    f.write(f"- True n: {TRUE_N}\n")
    f.write(f"- Recovered n: {n_mean:.4f} ± {n_std:.4f}\n")
    f.write(f"- Bias: {n_mean - TRUE_N:+.4f}\n")
    f.write(f"- Tension with GR: {z_gr:.2f} sigma\n")
    f.write(f"- Efficiency: {data['efficiency']:.3%}\n")

print(f"Report saved to {report_name}")

# ------------------------
# PLOT
# ------------------------
plt.figure(figsize=(10, 6))
plt.hist(chain, bins=60, density=True, color='lightgreen', edgecolor='black', alpha=0.7)
plt.axvline(1.0, color='red', linestyle='--', label='GR (n=1.0)')
plt.axvline(TRUE_N, color='black', linewidth=2, label=f'True n ({TRUE_N})')
plt.axvline(n_mean, color='blue', label=f'Mean recovered ({n_mean:.4f})')
plt.legend()
plt.title(f"QW-1533 Unbiased Population Inference ({scenario})")
plt.xlabel("Propagation Exponent n")
plt.ylabel("Posterior Density")
plot_name = f"QW_1533_Plot_{scenario}_{timestamp}.png"
plt.savefig(plot_name)
print(f"Plot saved to {plot_name}")
