# OBSOLETE - Superceded by QW_1531_1532_Sanity_Audit.py (Scientific Audit Round 3)
import numpy as np
import matplotlib.pyplot as plt
import time

# ============================================================
# QW-1532: FULL HIERARCHICAL POPULATION INFERENCE (WITH SELECTION)
# ============================================================

np.random.seed(42)

# ------------------------
# CONFIGURATION
# ------------------------
N_EVENTS = 20
TRUE_N = 1.00        # change to 1.00 for GR sanity check
D_CALIB = 40.0
D_MIN, D_MAX = 40.0, 800.0

SIGMA_GW_REL = 0.12
SIGMA_CAL = 0.04
SIGMA_EM_REL = 0.05

DETECTION_THRESHOLD = 0.6
N_MCMC = 60000
BURNIN = 10000
STEP = 0.015

# ------------------------
# POPULATION SIMULATION
# ------------------------
def simulate_detected_catalog(n_events, true_n):
    D_EM, D_GW, SIGMA = [], [], []

    total = 0
    while len(D_EM) < n_events:
        total += 1

        # True distance ~ uniform in volume
        D_true = (np.random.uniform(0,1)*(D_MAX**3-D_MIN**3)+D_MIN**3)**(1/3)

        # Detection condition (selection)
        strain = (D_CALIB / D_true)**true_n
        if strain < DETECTION_THRESHOLD:
            continue

        # Observables
        D_em = D_true * (1 + np.random.normal(0, SIGMA_EM_REL))
        D_gw = D_CALIB * (D_true / D_CALIB)**true_n
        D_gw *= (1 + np.random.normal(0, SIGMA_GW_REL))
        D_gw *= (1 + np.random.normal(0, SIGMA_CAL))

        sigma_tot = D_gw * np.sqrt(SIGMA_GW_REL**2 + SIGMA_CAL**2)

        D_EM.append(D_em)
        D_GW.append(D_gw)
        SIGMA.append(sigma_tot)

    efficiency = n_events / total
    print(f"Detection efficiency: {efficiency:.3%}")

    return {
        "D_EM": np.array(D_EM),
        "D_GW": np.array(D_GW),
        "SIGMA": np.array(SIGMA),
        "efficiency": efficiency
    }

data_dict = simulate_detected_catalog(N_EVENTS, TRUE_N)
data = data_dict
efficiency = data_dict["efficiency"]

# ------------------------
# HIERARCHICAL MODEL
# ------------------------
def log_prior(n):
    return 0.0 if 0.3 < n < 1.4 else -np.inf

def log_likelihood_events(n, data):
    D_pred = D_CALIB * (data["D_EM"]/D_CALIB)**n
    return -0.5 * np.sum(((data["D_GW"] - D_pred)/data["SIGMA"])**2)

def log_selection_normalization(n):
    # P_det(n) ∝ V_max(n)
    # D_max(n): (D_calib/D)^n = threshold
    D_max = D_CALIB * DETECTION_THRESHOLD**(-1.0/n)
    D_max = min(D_max, D_MAX)
    return 3 * np.log(D_max)

def log_posterior(n, data):
    lp = log_prior(n)
    if not np.isfinite(lp):
        return -np.inf

    ll = log_likelihood_events(n, data)
    ln_norm = N_EVENTS * log_selection_normalization(n)

    return lp + ll - ln_norm

# ------------------------
# MCMC
# ------------------------
chain = np.zeros(N_MCMC)
n_curr = 1.0
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

print(f"MCMC acceptance rate: {accepted/N_MCMC:.2f}")

chain = chain[BURNIN:]

# ------------------------
# RESULTS
# ------------------------
n_mean = np.mean(chain)
n_std = np.std(chain)
z_GR = abs(n_mean - 1.0) / n_std

print("\n=== QW-1532 RESULT ===")
print(f"Recovered n = {n_mean:.4f} ± {n_std:.4f}")
print(f"Tension with GR: {z_GR:.2f} sigma")

# ------------------------
# SAVE REPORT
# ------------------------
timestamp = int(time.time())
scenario = "FIN" if TRUE_N < 0.9 else "GR"
report_name = f"QW_1532_Partial_Report_{scenario}_{timestamp}.md"
with open(report_name, "w") as f:
    f.write(f"# QW-1532 Partial Result ({scenario})\n\n")
    f.write(f"- True n: {TRUE_N}\n")
    f.write(f"- Recovered n: {n_mean:.4f} ± {n_std:.4f}\n")
    f.write(f"- Tension with GR: {z_GR:.2f} sigma\n")
    f.write(f"- Detection efficiency: {efficiency:.3%}\n")
    f.write(f"- MCMC Acceptance: {accepted/N_MCMC:.2f}\n")

print(f"Report saved to {report_name}")

# ------------------------
# PLOT
# ------------------------
plt.figure(figsize=(10, 6))
plt.hist(chain, bins=50, density=True, color='skyblue', alpha=0.7)
plt.axvline(1.0, color='red', linestyle='--', label='General Relativity (n=1.0)')
plt.axvline(TRUE_N, color='black', linewidth=2, label=f'True n ({TRUE_N})')
plt.axvline(n_mean, color='green', label=f'Recovered n ({n_mean:.4f})')
plt.legend()
plt.xlabel("Propagation exponent n")
plt.ylabel("Posterior density")
plt.title(f"QW-1532 Population Inference ({scenario})")
plot_name = f"QW_1532_Plot_{scenario}_{timestamp}.png"
plt.savefig(plot_name)
print(f"Plot saved to {plot_name}")
# plt.show() # Disabled for headless execution
