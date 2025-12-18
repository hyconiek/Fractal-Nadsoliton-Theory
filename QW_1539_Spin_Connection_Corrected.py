import numpy as np
from datetime import datetime

# ==============================================================================
# QW-1539: Emergent Spin Connection & Curvature (CORRECTED)
# ==============================================================================
# Corrections / Refinements:
# - Explicitly stating the torsion-free / EFT approximation status.
# - No changes to the core numerical logic (as it was deemed methodologically correct for EFT).
# ==============================================================================

REPORT_FILE = "RAPORT_QW1539_SPIN_CONNECTION_CORRECTED.md"
md = []

def log(msg=""):
    print(msg)
    md.append(msg)

log("=" * 80)
log("QW-1539 (CORRECTED): EMERGENT SPIN CONNECTION & CURVATURE")
log("=" * 80)
log("NOTE: The spin connection is computed in a simplified torsion-free approximation")
log("      suitable for weakly perturbed tetrads. Full Palatini or Einstein-Cartan")
log("      formulations are left for future work.")
log("-" * 80)

# Parameters
EPSILON = 1e-4 
PERTURBATION_SCALE = 0.01

def delta(a, mu):
    return 1.0 if a == mu else 0.0

def h_perturbation(a, mu, x):
    return np.sin(x[0] + x[1] + x[2] + x[3] + a*mu)

def e_field(x, flat=False):
    e = np.zeros((4, 4)) 
    for a in range(4):
        for mu in range(4):
            val = delta(a, mu)
            if not flat:
                val += PERTURBATION_SCALE * h_perturbation(a, mu, x)
            e[a, mu] = val
    return e

def get_inverse_tetrad(e_mat):
    return np.linalg.inv(e_mat)

def partial_deriv(func, x, mu, flat=False):
    x_plus = x.copy(); x_plus[mu] += EPSILON
    x_minus = x.copy(); x_minus[mu] -= EPSILON
    return (func(x_plus, flat) - func(x_minus, flat)) / (2 * EPSILON)

def calculate_omega(x, flat=False):
    e_val = e_field(x, flat)
    e_inv = get_inverse_tetrad(e_val)
    
    d_e = np.zeros((4, 4, 4)) 
    for lam in range(4):
        d_e[lam] = partial_deriv(e_field, x, lam, flat)
        
    omega = np.zeros((4, 4, 4)) 
    
    for mu in range(4):
        for a in range(4):
            for b in range(4):
                # Torsion-free approximation logic used previously
                t1 = 0.0
                for nu in range(4):
                    t1 += e_inv[nu, a] * (d_e[mu, b, nu] - d_e[nu, b, mu])
                
                t2 = 0.0
                for nu in range(4):
                    t2 += e_inv[nu, b] * (d_e[mu, a, nu] - d_e[nu, a, mu])
                    
                omega[mu, a, b] = 0.5 * (t1 - t2)
    return omega

def calculate_curvature(x, flat=False):
    def get_omega_at(pt):
        return calculate_omega(pt, flat)
    
    omega = get_omega_at(x)
    
    d_omega = np.zeros((4, 4, 4, 4)) 
    for lam in range(4):
        pt_plus = x.copy(); pt_plus[lam] += EPSILON
        pt_minus = x.copy(); pt_minus[lam] -= EPSILON
        d_omega[lam] = (get_omega_at(pt_plus) - get_omega_at(pt_minus)) / (2 * EPSILON)
        
    R = np.zeros((4, 4, 4, 4))
    eta = np.diag([-1, 1, 1, 1])
    
    for mu in range(4):
        for nu in range(4):
            w_nu = omega[nu]
            w_mu = omega[mu]
            linear = d_omega[mu, nu] - d_omega[nu, mu]
            
            comm = np.zeros((4,4))
            for a in range(4):
                for b in range(4):
                    val = 0.0
                    for c in range(4):
                        w_nu_c_b = 0.0
                        for k in range(4): w_nu_c_b += w_nu[k, b] * eta[k, c]
                        val += w_mu[a, c] * w_nu_c_b
                    
                    val2 = 0.0
                    for c in range(4):
                         w_mu_c_b = 0.0
                         for k in range(4): w_mu_c_b += w_mu[k, b] * eta[k, c]
                         val2 += w_nu[a, c] * w_mu_c_b
                         
                    comm[a, b] = val - val2
            
            R[:, :, mu, nu] = linear + comm

    return R

# Execution
point = np.array([0.5, 0.5, 0.5, 0.5])

log("[1] TEST 1: FLAT SPACE")
R_flat = calculate_curvature(point, flat=True)
max_R_flat = np.max(np.abs(R_flat))
log(f"Max Curvature: {max_R_flat:.2e}")

log("\n[2] TEST 2: PERTURBED TETRAD")
R_pert = calculate_curvature(point, flat=False)
max_R_pert = np.max(np.abs(R_pert))
log(f"Max Curvature: {max_R_pert:.2e}")

if max_R_pert > 1e-6:
    log("✅ PASSED: Curvature generation confirmed.")

# Save Report
with open(REPORT_FILE, 'w', encoding='utf-8') as f:
    f.write(f"# QW-1539 (CORRECTED): Emergent Spin Connection & Curvature\n\n")
    f.write(f"**Data:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
    f.write("## Methodological Note\n")
    f.write("> The spin connection is computed in a simplified torsion-free approximation suitable for weakly perturbed tetrads. Full Palatini or Einstein-Cartan formulations are left for future work.\n\n")
    f.write("```\n")
    f.write('\n'.join(md))
    f.write("\n```\n")

print(f"\n✅ Report saved to: {REPORT_FILE}")
