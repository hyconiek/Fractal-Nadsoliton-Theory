import numpy as np
import matplotlib.pyplot as plt

print("="*60)
print("QW-1528b: Relative Propagation Test (Anchor Method)")
print("="*60)

# --- 1. Model: Relative Amplitude ---
# We compare two events:
# Event 1 (Anchor, e.g., GW170817): measurements h1, Mc1, DL1_prior
# Event 2 (Target, e.g., GW150914): measurements h2, Mc2
#
# Ratio R = h2 / h1
# Model R = (Mc2/Mc1)^(5/3) * (DL1/DL2)^n
#
# Parameters to estimate: [DL1, DL2, n]
# Note: DL1 is tightly constrained by EM prior!

def model_ratio(Mc1, DL1, Mc2, DL2, n):
    return (Mc2 / Mc1)**(5/3) * (DL1 / DL2)**n

def log_prior(DL1, DL2, n):
    # DL1 (Anchor): Gaussian Prior from EM (passed in data, but here we hardcode checks)
    if DL1 <= 0 or DL2 <= 0: return -np.inf
    
    # DL2 (Target): Volumetric Prior ~ DL^2
    # But wait, we are estimating parameters. 
    # Let's handle the specific priors in log_posterior
    
    # n: Uniform [0.5, 1.2]
    if n < 0.5 or n > 1.2: return -np.inf
    
    return 2 * np.log(DL2) # Volumetric for Target only

def log_posterior(params, data):
    # params = [DL1, DL2, n]
    # data = [R_obs, sigma_R, Mc1, Mc2, DL1_mean, DL1_std]
    
    DL1, DL2, n = params
    R_obs, sigma_R, Mc1, Mc2, DL1_mean, DL1_std = data
    
    # 1. Priors
    # Anchor Distance (Gaussian)
    lp_DL1 = -0.5 * ((DL1 - DL1_mean) / DL1_std)**2
    
    # Target Distance & n (General)
    lp_general = log_prior(DL1, DL2, n)
    
    if not np.isfinite(lp_general): return -np.inf
    
    loss_prior = lp_DL1 + lp_general
    
    # 2. Likelihood (Gaussian on Ratio)
    R_pred = model_ratio(Mc1, DL1, Mc2, DL2, n)
    ll = -0.5 * ((R_obs - R_pred) / sigma_R)**2
    
    return loss_prior + ll

# --- 2. MCMC Sampler ---
def run_mcmc_relative(data, n_samples=50000):
    # Unpack for initialization
    R_obs, sigma_R, Mc1, Mc2, DL1_mean, DL1_std = data
    
    # Init params: DL1 near mean, DL2 near expected from ratio, n=1
    # Rough guess for DL2: R ~ (Mc2/Mc1)^(5/3) * (DL1/DL2)^1 => DL2 ~ DL1 * (Mc2/Mc1)^(5/3) / R
    mc_ratio_factor = (Mc2/Mc1)**(5/3)
    dl2_guess = DL1_mean * mc_ratio_factor / R_obs
    
    current_params = np.array([DL1_mean, dl2_guess, 1.0])
    current_log = log_posterior(current_params, data)
    
    chain = np.zeros((n_samples, 3))
    accepted = 0
    
    # Proposal widths
    # DL1: constrained by EM (sigma ~ 3)
    # DL2: broad
    # n: small
    props = np.array([1.0, 20.0, 0.02])
    
    for i in range(n_samples):
        proposal = current_params + np.random.normal(0, props)
        prop_log = log_posterior(proposal, data)
        
        if np.log(np.random.rand()) < (prop_log - current_log):
            current_params = proposal
            current_log = prop_log
            accepted += 1
            
        chain[i] = current_params
        
    print(f"  Acceptance Rate: {accepted/n_samples:.2f}")
    return chain[10000:]

# --- 3. Scenario Runner ---
def run_scenario(name, anchor_event, target_event):
    print(f"\nComparing: {anchor_event['name']} (Anchor) vs {target_event['name']}")
    
    # 1. Generate Synthetic Data (Injecting n = 0.66)
    TRUE_N = 0.66
    
    # Amplitudes (arbitrary units, but consistent)
    A_ref = 1e-21
    h1_true = A_ref * (anchor_event['Mc']/30)**(5/3) / (anchor_event['DL']/400)**TRUE_N
    h2_true = A_ref * (target_event['Mc']/30)**(5/3) / (target_event['DL']/400)**TRUE_N
    
    # Ratio
    R_true = h2_true / h1_true
    
    # Noise limits
    # Relative error add in quadrature: sigma_R/R ~ sqrt( (s1/h1)^2 + (s2/h2)^2 )
    # s/h = 1/SNR
    rel_err_1 = 1.0 / anchor_event['SNR']
    rel_err_2 = 1.0 / target_event['SNR']
    rel_err_R = np.sqrt(rel_err_1**2 + rel_err_2**2)
    
    sigma_R = R_true * rel_err_R
    R_obs = R_true + np.random.normal(0, sigma_R)
    
    print(f"  True Ratio: {R_true:.3f}")
    print(f"  Obs Ratio:  {R_obs:.3f} +/- {sigma_R:.3f} ({100*rel_err_R:.1f}%)")
    
    # 2. Prepare Data Package
    # Anchor EM info: DL ~ N(True, 3Mpc)
    DL1_mean = anchor_event['DL'] + np.random.normal(0, 3.0) 
    DL1_std = 3.0
    
    data = [R_obs, sigma_R, anchor_event['Mc'], target_event['Mc'], DL1_mean, DL1_std]
    
    # 3. Run
    chain = run_mcmc_relative(data)
    
    # 4. Results
    n_mean = np.mean(chain[:, 2])
    n_std = np.std(chain[:, 2])
    
    print(f"  Recovered n: {n_mean:.3f} +/- {n_std:.3f}")
    
    z_score = abs(1.0 - n_mean) / n_std
    print(f"  Tension with GR (n=1): {z_score:.2f} sigma")
    
    # Bayes Factor vs GR (Density Ratio)
    hist, edges = np.histogram(chain[:, 2], bins=50, density=True, range=(0.4, 1.2))
    centers = 0.5*(edges[:-1]+edges[1:])
    try:
        dens_1 = hist[np.argmin(np.abs(centers - 1.0))]
        dens_066 = hist[np.argmin(np.abs(centers - 0.66))]
        bayes_fac = dens_066 / (dens_1 + 1e-6)
        print(f"  Bayes Factor (FIN/GR): {bayes_fac:.2f}")
        return n_mean, n_std, bayes_fac
    except:
        return n_mean, n_std, 0.0

# --- 4. Main Execution ---

# Event Definitions
GW170817 = {
    "name": "GW170817",
    "Mc": 1.188,
    "DL": 40.0,
    "SNR": 32.4
}

GW150914 = {
    "name": "GW150914",
    "Mc": 28.0,
    "DL": 410.0,
    "SNR": 24.0
}

BBH_HighSNR = {
    "name": "Sim-BBH-HighSNR",
    "Mc": 30.0,
    "DL": 800.0,
    "SNR": 50.0
}

results = []

# Pair 1: BNS Anchor vs Tyical BBH
res1 = run_scenario("Pair 1", GW170817, GW150914)
results.append({"pair": "GW170817 vs GW150914", "res": res1})

# Pair 2: BNS Anchor vs High SNR BBH (Simulated)
res2 = run_scenario("Pair 2", GW170817, BBH_HighSNR)
results.append({"pair": "GW170817 vs Sim-BBH", "res": res2})

# Save Output
with open("QW_1528b_Relative_Report.md", "w") as f:
    f.write("# QW-1528b: Relative Propagation Test (Anchor Method)\n\n")
    f.write("| Pair | Recovered n | Sigma vs GR | Bayes Factor |\n")
    f.write("|---|---|---|---|\n")
    for row in results:
        param = row['res']
        sigma = abs(1.0 - param[0]) / param[1]
        f.write(f"| {row['pair']} | {param[0]:.3f} +/- {param[1]:.3f} | {sigma:.2f} | {param[2]:.2f} |\n")
        
print("\n[Done] saved to QW_1528b_Relative_Report.md")
