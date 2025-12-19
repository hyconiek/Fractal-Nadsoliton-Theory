# OBSOLETE - Superceded by QW_1530_Rubikon_Audit.py (Scientific Audit Round 3)
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

print("="*60)
print("QW-1530: SELECTION-AWARE HIERARCHICAL RUBIKON TEST")
print("="*60)

# --- CONFIGURATION ---
# Simulate a future catalog of BNS/NSBH events with EM counterparts
N_EVENTS = 20
TRUE_N = 1.00 # GR Reality (Sanity Check)
# TRUE_N = 0.66 # FIN Reality

# --- 1. SELECTION-BIASED CATALOG SIMULATION ---
def generate_biased_catalog(n_events, true_n, detection_threshold=1.0):
    print(f"\n[Simulation] Generating {n_events} DETECTED events (Selection Biased)")
    
    events = {
        "D_EM": [], "D_GW": [], "sigma_GW": []
    }
    
    # We need to simulate a population and reject undetectable events
    # To get N_EVENTS detected, we simulate batch by batch
    
    D_calib = 40.0
    # Calibration error
    sigma_cal = 0.04
    
    count = 0
    total_simulated = 0
    
    while count < n_events:
        # 1. Simulate true Distance (Uniform in Volume)
        D_true = (np.random.uniform(0, 1) * (800**3 - 40**3) + 40**3)**(1/3)
        total_simulated += 1
        
        # 2. Simulate GW Amplitude (Signal)
        # Model: h ~ (D_calib / D)^n
        # Let's say h_calib = 10 (SNR units). So threshold=1 means SNR>1 is not realistic, 
        # usually SNR > 8. Let's work in "Distance Equivalent" units.
        # Detectable if D_GW < Horizon.
        # Horizon D_H(n) = D_calib * (SNR_calib / SNR_thresh)^(1/n)
        # Let's simplify: Detection Prob P(det) is step function on strain.
        
        # Let's calculate expected strain (arbitrary units, threshold=1)
        strain_expected = (D_calib / D_true)**true_n
        
        # Add noise to strain? No, selection is on observed strain usually.
        # Let's keep it simple as User requested:
        # P_det(n) propto D_max(n)^3.
        # This implies we just define a Horizon for the simulation.
        # Horizon for TRUE physics: h_true > h_th
        
        if strain_expected > detection_threshold: 
            # Detected! Now generate observables.
            
            # EM Obs
            D_EM_obs = D_true * (1 + np.random.normal(0, 0.05))
            
            # GW Obs (Inferred Distance assuming GR n=1 usually, but let's stick to raw distance D_GW)
            # D_GW in user's model is D_calib * (D_EM/D_calib)^n
            D_GW_expected = D_calib * (D_true / D_calib)**true_n
            D_GW_obs = D_GW_expected * (1 + np.random.normal(0, 0.12))
            
            # Calibration error
            D_GW_obs *= (1 + np.random.normal(0, sigma_cal))
            
            events["D_EM"].append(D_EM_obs)
            events["D_GW"].append(D_GW_obs)
            events["sigma_GW"].append(D_GW_obs * 0.12) # Approximate
            count += 1
            
    print(f"  Simulation Efficiency: {n_events/total_simulated:.1%}")
    
    events["D_EM"] = np.array(events["D_EM"])
    events["D_GW"] = np.array(events["D_GW"])
    events["sigma_GW"] = np.array(events["sigma_GW"])
    events["D_calib"] = D_calib
    events["h_thresh"] = detection_threshold
    
    return events

data = generate_biased_catalog(N_EVENTS, TRUE_N)


# --- 2. SELECTION-AWARE STACKING (MLE) ---

def detection_probability(n, D_calib, h_thresh=1.0):
    # D_max(n) is distance where (D_calib/D)^n = h_thresh
    # => D_max = D_calib * h_thresh^(-1/n)
    D_max = D_calib * h_thresh**(-1.0/n)
    # P_det proportional to Volume ~ D_max^3
    # We strictly need density integral, but for uniform distribution, P_det ~ V_max.
    # Note: If D_max > Universe Boundary, this saturates. Assuming D_max < 800 roughly.
    return D_max**3

def neg_log_likelihood_selection(n_arr, data):
    n = n_arr[0]
    # 1. Standard Fit Term (Chi2)
    D_model = data["D_calib"] * (data["D_EM"] / data["D_calib"])**n
    resid = data["D_GW"] - D_model
    chi2 = np.sum((resid / data["sigma_GW"])**2)
    
    # 2. Selection Correction Term (Normalization)
    # The likelihood for each event is P(Data|n) / P(Det|n)
    # ln L = Sum [ ln P(Data|n) - ln P(Det|n) ]
    # ln P(Det|n) approx ln(D_max(n)^3) (normalization constant)
    
    # Correction: We must normalize properly. 
    # If we assume uniform prior in volume, P(Det|n) ~ D_max(n)^3.
    # ln P_det = 3 * ln D_max(n)
    
    D_max = data["D_calib"] * data["h_thresh"]**(-1.0/n)
    ln_P_det = 3 * np.log(D_max)
    
    # Total Log Likelihood
    lnL = -0.5 * chi2 - len(data["D_EM"]) * ln_P_det
    
    return -lnL

print("\n--- STEP 1: Selection-Aware Optimization (MLE) ---")
res = minimize(
    neg_log_likelihood_selection,
    x0=[1.0],
    args=(data,),
    method="Nelder-Mead"
)

n_sel = res.x[0]
print(f"Selection-Aware MLE: n = {n_sel:.4f}")


# --- 3. HIERARCHICAL BAYESIAN MCMC ---

def log_prior_n(n):
    if 0.4 < n < 1.4: return 0.0
    return -np.inf

def log_likelihood_event(n, D_EM, D_GW, sigma_GW, D_calib):
    D_pred = D_calib * (D_EM / D_calib)**n
    return -0.5 * ((D_GW - D_pred)/sigma_GW)**2

# Selection term for Bayes
def log_selection_correction(n, D_calib, h_thresh):
    # P(Det|n) ~ D_max^3
    D_max = D_calib * h_thresh**(-1.0/n)
    return 3 * np.log(D_max)

def log_posterior_hierarchical(n, data):
    lp = log_prior_n(n)
    if not np.isfinite(lp): return -np.inf
    
    ll_sum = 0.0
    # Sum over events
    for i in range(len(data["D_EM"])):
        ll_sum += log_likelihood_event(n, data["D_EM"][i], data["D_GW"][i], data["sigma_GW"][i], data["D_calib"])
        
    # Subtract Selection Bias (Normalization) for EACH event?
    # Yes, the likelihood is Product_i ( P(Di|n) / P(Det|n) )
    # So we subtract N * ln P(Det|n)
    ln_norm = len(data["D_EM"]) * log_selection_correction(n, data["D_calib"], data["h_thresh"])
    
    return lp + ll_sum - ln_norm

def run_hierarchical_mcmc(data, n_samples=30000):
    chain = np.zeros(n_samples)
    n_curr = 1.0 # Start at GR
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
        
    print(f"Bayesian MCMC Acceptance rate: {accepted/n_samples:.2f}")
    return chain[5000:] # Burn-in

print("\n--- STEP 2: Hierarchical Bayesian MCMC ---")
chain = run_hierarchical_mcmc(data)

n_mean = np.mean(chain)
n_std = np.std(chain)
z_GR = abs(n_mean - 1.0) / n_std

print(f"[Hierarchical Result] n = {n_mean:.4f} +/- {n_std:.4f}")
print(f"Tension with GR: {z_GR:.2f} sigma")


# --- 4. THE VERDICT ---
print("\n--- FINAL VERDICT ---")
verdict = "UNDECIDED"
if z_GR > 5.0 and abs(n_mean - 1.0) > 0.05:
    verdict = "🔴 GR INCOMPLETE (FIN CONFIRMED)"
    details = "Robust selection-aware hierarchical evidence > 5 sigma."
elif z_GR < 2.0:
    verdict = "⚠️ CONSISTENT WITH GR"
    details = "No statistically significant deviation found."
else:
    verdict = "❓ INCONCLUSIVE"
    details = f"Tension {z_GR:.2f} sigma is suggestive but < 5 sigma."
    
print(f"Verdict: {verdict}")
print(f"Details: {details}")

# Save Report
with open("QW_1530_Final_Verdict.md", "w") as f:
    f.write("# QW-1530: The Final Rubikon Verdict\n\n")
    f.write("## 1. Methodology\n")
    f.write("- **Corrected for Selection Bias:** (Malmquist)\n")
    f.write("- **Hierarchical Bayesian Inference:** (MCMC)\n")
    f.write("- **Catalog:** 20 Simulated Events (BNS/NSBH)\n\n")
    
    f.write("## 2. Results\n")
    f.write(f"- **Recovered Exponent:** n = {n_mean:.4f} +/- {n_std:.4f}\n")
    f.write(f"- **Tension with GR:** {z_GR:.2f} sigma\n\n")
    
    f.write("## 3. Verdict\n")
    f.write(f"# {verdict}\n")
    f.write(f"{details}\n")
    
print("\n[Done] Saved to QW_1530_Final_Verdict.md")
