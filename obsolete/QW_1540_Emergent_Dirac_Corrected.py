# OBSOLETE - Superceded by QW_1540_Dirac_Audit.py (Scientific Audit Round 3)
import numpy as np
from datetime import datetime

# ==============================================================================
# QW-1540: Emergent Dirac Equation (CORRECTED)
# ==============================================================================
# Corrections:
# 1. Added explicit imaginary factor '1j' to the Dirac operator D.
# 2. Corrected Spin Connection term to use the commutator [gamma_b, gamma_c] / 4.
#    Actually, standard is Sigma_ab = i/4 [gamma_a, gamma_b]? Or 1/4 [gamma, gamma].
#    The reviewer specified: 1/4 * omega * [gamma, gamma].
# ==============================================================================

REPORT_FILE = "RAPORT_QW1540_DIRAC_EQ_CORRECTED.md"
md = []

def log(msg=""):
    print(msg)
    md.append(msg)

log("=" * 80)
log("QW-1540 (CORRECTED): EMERGENT DIRAC EQUATION")
log("=" * 80)

# Dirac Matrices
I2 = np.eye(2, dtype=complex); Z2 = np.zeros((2, 2), dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
g0 = np.block([[I2, Z2], [Z2, -I2]])
g1 = np.block([[Z2, sx], [-sx, Z2]])
g2 = np.block([[Z2, sy], [-sy, Z2]])
g3 = np.block([[Z2, sz], [-sz, Z2]])
gammas = [g0, g1, g2, g3]

EPSILON = 1e-4

def get_flat_tetrad(x): return np.eye(4)

def get_perturbed_tetrad(x):
    e = np.eye(4)
    val = 0.01 * np.sin(np.sum(x))
    e[0, 0] += val
    return e

def get_inverse_tetrad(e): return np.linalg.inv(e)

def get_spin_connection(x, e_func):
    # Simplified calculation reuse
    e_val = e_func(x)
    e_inv = get_inverse_tetrad(e_val)
    d_e = np.zeros((4, 4, 4))
    for lam in range(4):
        x_p = x.copy(); x_p[lam] += EPSILON
        x_m = x.copy(); x_m[lam] -= EPSILON
        d_e[lam] = (e_func(x_p) - e_func(x_m)) / (2*EPSILON)
    
    omega = np.zeros((4, 4, 4))
    for mu in range(4):
        for a in range(4):
            for b in range(4):
                t1 = 0.0; t2 = 0.0
                for nu in range(4): t1 += e_inv[nu, a] * (d_e[mu, b, nu] - d_e[nu, b, mu])
                for nu in range(4): t2 += e_inv[nu, b] * (d_e[mu, a, nu] - d_e[nu, a, mu])
                omega[mu, a, b] = 0.5 * (t1 - t2)
    return omega

def apply_dirac_op(x, psi_func, e_func):
    e_val = e_func(x)
    e_inv = get_inverse_tetrad(e_val)
    omega = get_spin_connection(x, e_func)
    psi = psi_func(x)
    
    # 1. d_mu psi
    d_psi = np.zeros((4, 4), dtype=complex)
    for mu in range(4):
        x_p = x.copy(); x_p[mu] += EPSILON
        x_m = x.copy(); x_m[mu] -= EPSILON
        d_psi[mu] = (psi_func(x_p) - psi_func(x_m)) / (2*EPSILON)
        
    term_deriv = np.zeros(4, dtype=complex)
    for a in range(4):
        sum_deriv = np.zeros(4, dtype=complex)
        for mu in range(4):
            sum_deriv += e_inv[mu, a] * d_psi[mu]
        term_deriv += np.dot(gammas[a], sum_deriv)
        
    # 2. Spin connection term: correction to commutators
    # 1/4 * omega_mu^bc * sigma_bc
    # sigma_bc = 0.5 * [gamma_b, gamma_c] (geometric convention often)
    # The prompt explicitly suggested: 1/4 omega_bc [gamma^b, gamma^c].
    
    term_spin = np.zeros(4, dtype=complex)
    for a in range(4):
        # e_a^mu factor for the whole D operator
        # D = gamma^a e_a^mu ( ... )
        for mu in range(4):
            factor = e_inv[mu, a]
            if abs(factor) < 1e-9: continue
            
            spin_conn_mat = np.zeros((4, 4), dtype=complex)
            for b in range(4):
                for c in range(4):
                    w = omega[mu, b, c]
                    if abs(w) < 1e-9: continue
                    
                    # Correction: Commutator
                    # [gamma_b, gamma_c]
                    comm = gammas[b] @ gammas[c] - gammas[c] @ gammas[b]
                    
                    # Factor 1/4 from definition, maybe another 1/2 from sigma?
                    # Reviewer: "1/4 omega [gamma, gamma]"
                    spin_conn_mat += w * comm
            
            spin_conn_mat *= 0.25 
            
            term_spin += factor * np.dot(gammas[a], np.dot(spin_conn_mat, psi))

    # Correction: Add '1j' factor to the whole operator D = i * gamma ...
    # Standard: i gamma^mu D_mu
    result = 1j * (term_deriv + term_spin)
    
    return result

def psi_plane(x):
    p = np.array([1.0, 0.0, 0.0, 0.0])
    phase = np.exp(-1j * np.dot(p, x))
    u = np.array([1, 0, 0, 0], dtype=complex)
    return u * phase

x0 = np.array([0., 0., 0., 0.])

log("[1] FLAT LIMIT TEST")
res_flat = apply_dirac_op(x0, psi_plane, get_flat_tetrad)
# Expected: i * gamma^0 * (-i m) u = i * (-i) * 1 * u = 1 * u
# Wait. D psi = i gamma^mu d_mu psi.
# d_0 psi = -i m psi.
# i gamma^0 (-i m) = i(-i) m gamma^0 = m gamma^0.
# gamma^0 u = u.
# So expected is m * u = 1.0 * u.
expected_flat = 1.0 * np.array([1, 0, 0, 0], dtype=complex)
diff_flat = np.linalg.norm(res_flat - expected_flat)

log(f"Result Flat: {res_flat}")
log(f"Expected:    {expected_flat}")
log(f"Diff:        {diff_flat:.2e}")

res_curved = apply_dirac_op(x0, psi_plane, get_perturbed_tetrad)
log(f"\n[2] CURVED TEST: {res_curved}")

with open(REPORT_FILE, 'w', encoding='utf-8') as f:
    f.write(f"# QW-1540 (CORRECTED): Emergent Dirac Equation\n\n")
    f.write(f"**Data:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.write("```\n")
    f.write('\n'.join(md))
    f.write("\n```\n")

print(f"\n✅ Report saved to: {REPORT_FILE}")
