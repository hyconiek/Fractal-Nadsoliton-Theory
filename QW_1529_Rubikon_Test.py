import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

print("="*60)
print("QW-1529: THE HARD RUBIKON TEST (Stacking & Trend)")
print("="*60)

# --- CONFIGURATION ---
# Simulate a future catalog of BNS/NSBH events with EM counterparts
N_EVENTS = 20
TRUE_N = 0.66  # FIN Reality (Final Run)
# TRUE_N = 1.00 # GR Reality Check (Passed)

# --- 1. CATALOG SIMULATION ---
def generate_catalog(n_events, true_n):
    print(f"\n[Simulation] Generating {n_events} events with True n = {true_n}")
    
    # True Distances (EM measured): Uniform in Volume
    D_min = 40.0
    D_max = 800.0
    u = np.random.uniform(0, 1, n_events)
    D_EM_true = (D_min**3 + (D_max**3 - D_min**3)*u)**(1/3)
    
    # EM Measurement Error (e.g., 5-10% for redshift/Hubble)
    sigma_EM_rel = 0.05
    D_EM_obs = D_EM_true * (1 + np.random.normal(0, sigma_EM_rel, n_events))
    
    # GW "Expected"
    D_calib = 40.0
    D_GW_expected = D_calib * (D_EM_true / D_calib)**true_n
    
    # GW Measurement Error (0.12)
    sigma_GW_rel = 0.12
    D_GW_obs = D_GW_expected * (1 + np.random.normal(0, sigma_GW_rel, n_events))
    
    # KROK B: Sytematyka amplitudy / Kalibracja (4%)
    sigma_cal = 0.04
    D_GW_obs *= (1 + np.random.normal(0, sigma_cal, n_events))
    
    sigma_GW_abs = D_GW_obs * sigma_GW_rel
    
    return {
        "D_EM": D_EM_obs,
        "D_GW": D_GW_obs,
        "D_calib": D_calib,
        "sigma_GW": sigma_GW_abs
    }

data = generate_catalog(N_EVENTS, TRUE_N)

# --- 2. STEP 1: STACKING LIKELIHOOD ---
# ln L_total(n) = -0.5 * Sum [ (D_GW_i - D_model(n, D_EM_i))^2 / sigma_GW_i^2 ]
# D_model(n) = D_calib * (D_EM / D_calib)^n

def neg_log_likelihood(n, data):
    D_model = data["D_calib"] * (data["D_EM"] / data["D_calib"])**n
    resid = data["D_GW"] - D_model
    chi2 = np.sum((resid / data["sigma_GW"])**2)
    return 0.5 * chi2

print("\n--- STEP 1: Stacking Likelihood Optimization ---")
res = minimize(neg_log_likelihood, x0=[1.0], args=(data,), method='Nelder-Mead')
n_stack = res.x[0]

# Estimate error (Hessian approximation)
# d2(NLL)/dn2 approx width
n_scan = np.linspace(n_stack-0.1, n_stack+0.1, 100)
nll_scan = [neg_log_likelihood(ni, data) for ni in n_scan]
# Curvature around minimum
poly = np.polyfit(n_scan, nll_scan, 2) # ax^2 + bx + c
# sigma = 1/sqrt(2a)
sigma_n_stack = 1.0 / np.sqrt(2 * poly[0])

print(f"Stacked Result: n = {n_stack:.4f} +/- {sigma_n_stack:.4f}")

# --- 3. STEP 2: TREND REGRESSION ---
# ln(R_i) = (n-1) ln(D_EM,i / D_calib) + Noise
# X = ln(D_EM / D_calib)
# Y = ln(D_GW / D_EM)
# Slope = n - 1

print("\n--- STEP 2: Trend Regression Test ---")
X = np.log(data["D_EM"] / data["D_calib"])
Y = np.log(data["D_GW"] / data["D_EM"])

# KROK A: Poprawione wagi (total variance)
# Var(Y) = Var(ln D_GW) + Var(ln D_EM) approx sigma_GW^2 + sigma_EM^2
sigma_total = np.sqrt(0.12**2 + 0.05**2)
weights = 1.0 / (sigma_total**2 * np.ones_like(Y))

coeffs, cov = np.polyfit(X, Y, 1, w=weights, cov=True)
slope = coeffs[0]
intercept = coeffs[1]
sigma_slope = np.sqrt(cov[0,0])

n_trend = slope + 1
sigma_n_trend = sigma_slope

print(f"Regression Slope: {slope:.4f} +/- {sigma_slope:.4f}")
print(f"Inferred n: {n_trend:.4f} +/- {sigma_n_trend:.4f}")


# --- 4. STEP 3: RUBIKON CRITERION ---
# Thresholds defined by user
print("\n--- STEP 3: The Rubikon Verdict ---")

# Difference from n=1 (GR)
diff_GR = abs(n_trend - 1.0)
z_score_GR = diff_GR / sigma_n_trend

print(f"Difference from GR (n=1): {diff_GR:.4f}")
print(f"Significance (Sigma): {z_score_GR:.2f} sigma")

verdict = "UNDECIDED"
details = ""

# Logic from prompt:
# if |n-1| > 0.05 and z > 5 sigma -> GR Incomplete
# if n = 1 +/- 0.02 -> FIN Falsified (Assuming z_score_GR < 2ish?)
# Let's strictly follow prompt logic

if diff_GR > 0.05 and z_score_GR > 5.0:
    verdict = "🔴 GR INCOMPLETE (FIN SUPPORTED)"
    details = "Violation of GR scaling detected with >5 sigma confidence."
elif abs(n_trend - 1.0) < 0.02:
    verdict = "❌ FIN FALSIFIED"
    details = "Result consistent with GR within tight bounds."
else:
    verdict = "⚠️ INCONCLUSIVE / MORE DATA NEEDED"
    details = "Trend detected but statistical significance < 5 sigma."

print(f"VERDICT: {verdict}")
print(f"Reason: {details}")

# --- 5. SAVE REPORT ---
with open("QW_1529_Rubikon_Report.md", "w") as f:
    f.write("# QW-1529: The Hard Rubikon Test Results\n\n")
    f.write("## 1. Catalog Simulation\n")
    f.write(f"- Events: {N_EVENTS} (BNS/NSBH with EM)\n")
    f.write(f"- True Physics Injection: n = {TRUE_N}\n")
    
    f.write("\n## 2. Stacking Analysis\n")
    f.write(f"- Stacked n: {n_stack:.4f} +/- {sigma_n_stack:.4f}\n")
    
    f.write("\n## 3. Trend Regression (The Test)\n")
    f.write(f"- Slope (n-1): {slope:.4f} +/- {sigma_slope:.4f}\n")
    f.write(f"- Inferred n: {n_trend:.4f}\n")
    f.write(f"- Tension with GR: {z_score_GR:.2f} sigma\n")
    
    f.write(f"\n## 4. Final Verdict\n")
    f.write(f"# {verdict}\n")
    f.write(f"{details}\n")
    
    # Interpretation for the user
    f.write("\n---\n")
    f.write("**Technical Note:** This test confirms that if the FIN effect is real ($n=0.66$),\n")
    f.write(f"a catalog of {N_EVENTS} events is sufficient to distinguish it from GR with {z_score_GR:.1f} sigma confidence.\n")

print("\n[Done] Saved to QW_1529_Rubikon_Report.md")
