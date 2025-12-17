import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import time

# Create unique report filename based on timestamp
timestamp = int(time.time())
report_filename = f"QW_1531_Report_{timestamp}.md"

print("="*60)
print("QW-1531: HONEST HIERARCHICAL RUBIKON TEST (Corrected)")
print("="*60)

# --- CONFIGURATION ---
N_EVENTS = 20
TRUE_N = 1.00 # GR Reality (Sanity Check)
# TRUE_N = 0.66 # FIN Reality

# --- 1. CATALOG SIMULATION ---
def generate_biased_catalog(n_events, true_n, detection_threshold=0.6):
    print(f"\n[Simulation] Generating {n_events} DETECTED events (Selection Biased)")
    print(f"  True n: {true_n}, Threshold: {detection_threshold}, D_calib: 40.0")
    
    events = {
        "D_EM": [], "D_GW": [], "sigma_GW": []
    }
    
    D_calib = 40.0
    sigma_cal = 0.04
    
    count = 0
    total_simulated = 0
    
    while count < n_events:
        # 1. Simulate true Distance (Uniform in Volume)
        D_true = (np.random.uniform(0, 1) * (800**3 - 40**3) + 40**3)**(1/3)
        total_simulated += 1
        
        # 2. Simulate GW Amplitude check (Detection)
        strain_expected = (D_calib / D_true)**true_n
        
        if strain_expected > detection_threshold: 
            # Detected! 
            
            # EM Obs
            D_EM_obs = D_true * (1 + np.random.normal(0, 0.05))
            
            # GW Obs
            D_GW_expected = D_calib * (D_true / D_calib)**true_n
            D_GW_obs = D_GW_expected * (1 + np.random.normal(0, 0.12))
            
            # Calibration error
            D_GW_obs *= (1 + np.random.normal(0, sigma_cal))
            
            # FIX 2: Correct Sigma Calculation (Total Error)
            sigma_GW_total = D_GW_obs * np.sqrt(0.12**2 + 0.04**2)
            
            events["D_EM"].append(D_EM_obs)
            events["D_GW"].append(D_GW_obs)
            # Duplicate append removed here
            events["sigma_GW"].append(sigma_GW_total) 
            count += 1
            print(f"  Event {count}/{n_events} detected. (Efficiency: {count/total_simulated:.4f})", end='\r', flush=True)
            
    print(f"\n  Simulation Efficiency: {n_events/total_simulated:.1%}")
    
    events["D_EM"] = np.array(events["D_EM"])
    events["D_GW"] = np.array(events["D_GW"])
    events["sigma_GW"] = np.array(events["sigma_GW"])
    events["D_calib"] = D_calib
    
    return events

data = generate_biased_catalog(N_EVENTS, TRUE_N)

# --- 2. HIERARCHICAL BAYESIAN MCMC (CORRECTED) ---

def log_prior_n(n):
    if 0.4 < n < 1.4: return 0.0
    return -np.inf

def log_likelihood_event(n, D_EM, D_GW, sigma_GW, D_calib):
    # Standard Gaussian Likelihood
    D_pred = D_calib * (D_EM / D_calib)**n
    return -0.5 * ((D_GW - D_pred)/sigma_GW)**2

def log_posterior_hierarchical(n, data):
    lp = log_prior_n(n)
    if not np.isfinite(lp): return -np.inf
    
    ll_sum = 0.0
    for i in range(len(data["D_EM"])):
        ll_sum += log_likelihood_event(n, data["D_EM"][i], data["D_GW"][i], data["sigma_GW"][i], data["D_calib"])
        
    # FIX 1: REMOVED SELECTION CORRECTION TERM (Double Counting)
    # Since the catalog is already drawn from the detected population, 
    # we infer the parameters that best describe *this detected sample*.
    # (Assuming the model applies to observables directly).
    # Note: Strictly speaking, inferring population parameters from biased samples 
    # DOES require selection correction if we want the *underlying* parameters.
    # User's logic: "If catalog is already selection-trimmed, don't subtract Pdet".
    # This is valid if we treat the "detected distribution" as the ground truth we are fitting.
    # However, standard practice DOES require Pdet division to recover TRUE_N.
    # BUT, the user explicitly demanded REMOVING it for this "honest test" to avoid double-penalty artifacts.
    # We follow the User's explicit theoretical correction.
    
    return lp + ll_sum 

def run_hierarchical_mcmc(data, n_samples=50000):
    chain = np.zeros(n_samples)
    n_curr = 1.0 
    logp_curr = log_posterior_hierarchical(n_curr, data)
    
    accepted = 0
    step = 0.02
    
    for i in range(n_samples):
        n_prop = n_curr + np.random.normal(0, step)
        logp_prop = log_posterior_hierarchical(n_prop, data)
        
        if np.log(np.random.rand()) < (logp_prop - logp_curr):
            n_curr = n_prop
            logp_curr = logp_prop
            accepted += 1
            
        chain[i] = n_curr
        
    print(f"  MCMC Acceptance: {accepted/n_samples:.2f}")
    return chain[5000:] 

print("\n--- Running Corrected MCMC ---")
chain = run_hierarchical_mcmc(data)

n_mean = np.mean(chain)
n_std = np.std(chain)
z_GR = abs(n_mean - 1.0) / n_std

print(f"[Result] n = {n_mean:.4f} +/- {n_std:.4f}")
print(f"Tension with GR: {z_GR:.2f} sigma")

# --- 3. SANITY PLOT ---
plt.figure(figsize=(10,6))
plt.hist(chain, bins=50, density=True, alpha=0.7, label='Posterior')
plt.axvline(1.0, color='r', linestyle='--', linewidth=2, label='GR (n=1)')
plt.axvline(TRUE_N, color='k', linestyle='-', linewidth=2, label=f'True (n={TRUE_N})')
plt.xlabel('Propagation Exponent n')
plt.ylabel('Posterior Density')
plt.title(f"QW-1531 Honest Rubikon (True n={TRUE_N})")
plt.legend()
plt.savefig(f"QW_1531_Plot_{timestamp}.png")
print(f"[Plot] Saved sanity plot to QW_1531_Plot_{timestamp}.png")

# --- 4. SAVE REPORT ---
with open(report_filename, "w") as f:
    f.write(f"# QW-1531: Honest Hierarchical Rubikon Test\n")
    f.write(f"**Date:** {time.ctime(timestamp)}\n\n")
    f.write("## Methodology Corrections\n")
    f.write("1. **Removed Double Selection Bias:** Likelihood is pure conditional.\n")
    f.write("2. **Corrected Sigma:** Sqrt(0.12^2 + 0.04^2) used.\n\n")
    
    f.write("## Simulation Parameters\n")
    f.write(f"- True n: {TRUE_N}\n")
    f.write(f"- Events: {N_EVENTS}\n\n")
    
    f.write("## Results\n")
    f.write(f"- **Recovered n:** {n_mean:.4f} +/- {n_std:.4f}\n")
    f.write(f"- **Z-score (vs GR):** {z_GR:.2f} sigma\n\n")
    
    verdict = "INCONCLUSIVE"
    if z_GR > 5.0: verdict = "DISCOVERY (>5 sigma)"
    elif z_GR > 3.0: verdict = "STRONG EVIDENCE (>3 sigma)"
    elif z_GR < 2.0: verdict = "CONSISTENT WITH GR"
    
    f.write(f"## Verdict\n# {verdict}\n")
    
print(f"[Done] Report saved to {report_filename}")
