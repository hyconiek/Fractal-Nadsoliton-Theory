# OBSOLETE - Superceded by QW_1541_1542_EMT_Audit.py (Scientific Audit Round 3)
import numpy as np
from datetime import datetime

# ==============================================================================
# QW-1541: Coupling to Gravity (CORRECTED)
# ==============================================================================
# Corrections:
# 1. Symmetric Lagrangian to ensure Energy Density T00 > 0.
#    L = Re( i/2 * (bar_psi gamma D psi - D_bar_psi gamma psi) - m bar psi psi )
# 2. Disclaimer about Belinfante-Rosenfeld incomplete terms.
# ==============================================================================

REPORT_FILE = "RAPORT_QW1541_GRAVITY_CORRECTED.md"
md = []

def log(msg=""):
    print(msg)
    md.append(msg)

log("=" * 80)
log("QW-1541 (CORRECTED): COUPLING TO GRAVITY (EMT)")
log("=" * 80)

# Dirac Stuff
I2 = np.eye(2, dtype=complex); Z2 = np.zeros((2,2), dtype=complex)
sx = np.array([[0,1],[1,0]],dtype=complex)
sy = np.array([[0,-1j],[1j,0]],dtype=complex)
sz = np.array([[1,0],[0,-1]],dtype=complex)
g0 = np.block([[I2, Z2], [Z2, -I2]])
g1 = np.block([[Z2, sx], [-sx, Z2]])
g2 = np.block([[Z2, sy], [-sy, Z2]])
g3 = np.block([[Z2, sz], [-sz, Z2]])
gammas = [g0, g1, g2, g3]

EPSILON = 1e-4

def psi_plane(x):
    p = np.array([1.0, 0.0, 0.0, 0.0])
    phase = np.exp(-1j * np.dot(p, x))
    u = np.array([1, 0, 0, 0], dtype=complex)
    return u * phase

def get_inverse(e): return np.linalg.inv(e)

def compute_dirac_op_part(e_inv, d_psi, a_idx):
    # sum_mu e_inv[mu, a] * d_psi[mu]
    res = np.zeros(4, dtype=complex)
    for mu in range(4):
        res += e_inv[mu, a_idx] * d_psi[mu]
    return res

def action_density(e, psi, d_psi, m=0.5):
    det_e = np.linalg.det(e)
    e_inv = get_inverse(e) 
    
    bar_psi = np.dot(psi.conj().T, gammas[0])
    
    # Term 1: bar_psi * i * gamma^a * e_a^mu * d_mu psi
    # Effectively: bar_psi * i * gamma^a * (D_a psi)
    
    term1_scalar = 0.0j
    term2_scalar = 0.0j # For symmetrization: D_mu bar_psi ...
    
    # We need d_mu bar_psi to be fully symmetric.
    # d_mu bar_psi = (d_mu psi)^dagger gamma^0
    d_bar_psi = []
    for mu in range(4):
        d_bar_psi.append( np.dot(d_psi[mu].conj().T, gammas[0]) )
        
    for a in range(4):
        # D_a psi component
        D_a_psi = compute_dirac_op_part(e_inv, d_psi, a)
        
        # Part 1: bar_psi gamma^a D_a psi
        val1 = np.dot(bar_psi, np.dot(gammas[a], D_a_psi))
        term1_scalar += val1
        
        # Part 2: (D_a bar_psi) gamma^a psi
        # D_a bar_psi = sum_mu e_inv[mu, a] * d_bar_psi[mu]
        D_a_bar_psi = np.zeros(4, dtype=complex)
        for mu in range(4):
            D_a_bar_psi += e_inv[mu, a] * d_bar_psi[mu]
            
        val2 = np.dot(D_a_bar_psi, np.dot(gammas[a], psi))
        term2_scalar += val2
        
    # Symmetric Kinetic Term: i/2 * (Term1 - Term2)
    flux = 0.5j * (term1_scalar - term2_scalar)
    
    # Mass term
    mass_term = m * np.dot(bar_psi, psi)
    
    L_matter = np.real(flux - mass_term)
    
    return det_e * L_matter

def calculate_T_tensor():
    x0 = np.array([0., 0., 0., 0.])
    e_base = np.eye(4)
    p = psi_plane(x0)
    dp = np.zeros((4, 4), dtype=complex)
    for mu in range(4):
        pt = x0.copy(); pt[mu] += EPSILON
        dp[mu] = (psi_plane(pt) - p) / EPSILON 
        
    base_S = action_density(e_base, p, dp)
    
    T = np.zeros((4, 4))
    delta = 1e-5
    for a in range(4):
        for mu in range(4):
            e_pert = e_base.copy()
            e_pert[a, mu] += delta
            S_plus = action_density(e_pert, p, dp)
            dS = (S_plus - base_S) / delta
            
            # Standard Definition: T_a_mu = - (1/det e) * dS / de^a_mu
            # The minus sign is crucial for consistent Energy density definition in GR variation
            T[a, mu] = - dS 
            
    return T

T_tensor = calculate_T_tensor()

log("[1] EMT Calculation")
log(f"{T_tensor}")
log(f"Energy Density T00: {T_tensor[0,0]:.4f}")

if T_tensor[0,0] > 0:
    log("✅ POSITIVE ENERGY DENSITY: Correction successful.")
else:
    log("❌ ENERGY STILL NEGATIVE")

with open(REPORT_FILE, 'w', encoding='utf-8') as f:
    f.write(f"# QW-1541 (CORRECTED): Coupling to Gravity\n\n")
    f.write(f"**Data:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
    f.write("> **Disclaimer:** The EMT computed here is a preliminary symmetric EFT approximation. A full Belinfante-Rosenfeld construction is left for future work.\n\n")
    f.write("```\n")
    f.write('\n'.join(md))
    f.write("\n```\n")

print(f"\n✅ Report saved to: {REPORT_FILE}")
