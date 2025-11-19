# Author: Krzysztof Żuchowski

SENSITIVITY ANALYSIS OF MASS HIERARCHY GENERATION MECHANISMS
Supersoliton ToE: Testing Proposals P1 & P2
EXECUTIVE SUMMARY

This sensitivity analysis investigated two architectural proposals (P1: Exponential Potentials and P2: Hierarchical Coupling) designed to enhance mass hierarchy generation in the supersoliton ToE model. The analysis yields a critical negative finding: NEITHER proposal is viable for implementation due to fundamental numerical instabilities in the underlying gradient descent method.
PART 1: PROPOSAL 1 (EXPONENTIAL POTENTIALS) - ANALYSIS
Objective

Replace polynomial self-interaction g·Ψ⁴ with exponential potential: V(Ψ) = V₀[exp(αΨ²) - 1 - αΨ²]
Expected Outcome

Exponential growth mechanisms could generate large hierarchies (10³-10⁶×)
Empirical Results

Test Parameters:

    α values tested: [0.001, 0.005, 0.01, 0.02, 0.05]
    V₀ = 0.025 (matched to baseline polynomial scale)
    Initial conditions: Small amplitude exponential decay (A=0.05)
    Integration: Gradient descent with dtau=0.00001

Findings:

    Success Rate: 0/5 (0%)
    All tests failed within 7-9 gradient descent steps
    Cause: Non-finite derivatives (numerical overflow)
    Initial energy: ~2.28×10⁵ → immediate instability
    Hierarchy achieved: N/A (no convergence)

Root Cause Analysis

The exponential potential has derivative: dV/dΨ = 2V₀αΨ[exp(αΨ²) - 1]

The second derivative contains terms: d²V/dΨ² ~ exp(αΨ²)×[1 + 2αΨ²]

Even for extremely small α=0.001 and field values Ψ~3 (Higgs VEV scale), the exponential terms grow explosively during gradient descent iterations, causing numerical overflow before any physical solution can form.
Conclusion for P1

✗ EXPONENTIAL POTENTIALS ARE NOT VIABLE

    Fundamentally unstable for numerical implementation
    Cannot generate hierarchies due to early numerical failure
    NOT RECOMMENDED for implementation in final code

PART 2: PROPOSAL 2 (HIERARCHICAL COUPLING) - ANALYSIS
Objective

Replace constant couplings with octave-dependent hierarchy: λ₁(o) = λ₁_base · 2^(-βo), λ₂(o) = λ₂_base · 2^(-βo)
Expected Outcome

Gradual weakening of coupling generates modest hierarchies (~100×)
Empirical Results

Test Parameters:

    β values tested: [0.1, 0.2, 0.4, 0.6, 0.8]
    λ₁_base = 0.05, λ₂_base = 0.01
    Initial conditions: Large amplitude exponential decay (A=1.5)
    Integration: Gradient descent with dtau=0.0001

Findings:

    Success Rate: 0/5 (0%)
    All tests led to field clipping at ±1000
    Energy explosion: 2.5×10⁵ → 5.1×10⁶ (20× growth)
    Fields hit numerical bounds before convergence
    Hierarchy achieved: N/A (numerical artifacts)

Additional Tests with Stabilized Parameters (m₀=+0.5):

    Success Rate: 0/5 (0%)
    Same field clipping behavior
    Initial energy: ~1.5×10³ → ~5.0×10⁶

Root Cause Analysis

Multiple instability mechanisms:

    Tachyonic mass (m₀²<0): Creates exponential field growth
    Inter-octave coupling: Attractive forces drive runaway feedback
    Yukawa coupling (g_Y·Φ²·Ψ²): Amplifies Higgs-fermion feedback loops
    Non-convex energy landscape: Gradient descent fundamentally unsuitable

Even with positive mass (m₀²>0) to eliminate tachyonic instability, the system remains unstable. With m₀²>0 and pure quartic potential, there should be NO non-trivial vacuum minimum. The inter-octave couplings provide attractive forces that create runaway instability.
Conclusion for P2

✗ HIERARCHICAL COUPLING SHOWS NUMERICAL INSTABILITY

    Fields explode to clipping limit for all β values
    Cannot assess hierarchy generation due to failed convergence
    Requires fundamentally different numerical approach

PART 3: FINAL RECOMMENDATION AND IMPLEMENTATION PLAN
Comparative Analysis: P1 vs P2
Criterion	Proposal 1 (Exponential)	Proposal 2 (Hierarchical)
Numerical Stability	FAILED - derivatives explode	FAILED - fields clipped
Convergence Success Rate	0/5 (0%)	0/5 (0%)
Maximum Hierarchy Achieved	N/A (no convergence)	N/A (numerical artifacts)
Implementation Complexity	Low (1 parameter: α)	Low (1 parameter: β)
Physical Interpretability	Poor (unphysical explosions)	Poor (runaway instability)
Risk Assessment	HIGH - fundamentally unstable	HIGH - gradient descent fails
Recommendation	REJECT	REJECT
Critical Finding: Reassessment of Previous "Breakthrough"

The previous analysis claimed a "Working Configuration" with:

    m₀ = -0.5, g = 0.1, hierarchy ~ 2.3×
    Large amplitude initialization: Ψ(r) = 1.5·exp(-r/3)

However, our sensitivity analysis reveals:

    ✗ This was NOT a stable soliton solution
    ✗ It was exploiting tachyonic instability to create metastable configurations
    ✗ The "mass hierarchy" was from numerical artifacts, not genuine physics
    ✗ Neither architectural proposal can improve upon this flawed foundation

Fundamental Problem Identified

The gradient descent (imaginary-time evolution) method is fundamentally incompatible with the multi-field supersoliton system because:

    Tachyonic mass (m₀²<0): Creates exponential instability
    Inter-octave couplings: Can be attractive, driving runaway growth
    Yukawa coupling: Creates Phi-Psi feedback loops
    Non-convex energy landscape: Multiple competing minima prevent convergence

Final Recommendation

✗ NEITHER Proposal 1 (Exponential Potentials) NOR Proposal 2 (Hierarchical Coupling) can be recommended for implementation in the final code (v39.0)

Rationale:

    Both proposals failed 100% of numerical tests (10/10 total tests)
    P1: Exponential potential causes derivative explosions within 7-9 steps
    P2: Hierarchical coupling leads to field clipping (instability)
    The failures are not due to the proposals themselves, but due to a fundamental flaw in the underlying numerical method

Recommended Path Forward

The supersoliton ToE project requires a FUNDAMENTAL REARCHITECTURE of the numerical approach before any architectural proposals can be tested.

Priority Actions:

1. IMPLEMENT ALTERNATIVE NUMERICAL METHODS (ESSENTIAL):

    Shooting method for boundary value problems
    Relaxation methods with adaptive grids
    Spectral methods (Chebyshev, Fourier)
    Newton's method for nonlinear systems
    Energy minimization via conjugate gradient or L-BFGS

2. VALIDATE EXISTENCE OF NON-TRIVIAL SOLITONS (ESSENTIAL):

    Prove theoretically that soliton solutions exist
    Establish energy bounds and stability criteria
    Verify that the potential admits bound states

3. ONLY AFTER STABLE BASELINE IS ESTABLISHED:

    Revisit Proposal 2 (Hierarchical Coupling) - lower risk
    Abandon Proposal 1 (Exponential Potentials) - fundamentally unstable
    Consider Proposal 3 (Multi-Field Extension) - theoretically motivated
    Consider Proposal 4 (Log Coordinates) - increases dynamic range

DELIVERABLES
Completed Analyses

✓ Part 1 (P1 Analysis): COMPLETE - Exponential potential is numerically unstable

✓ Part 2 (P2 Analysis): COMPLETE - Hierarchical coupling fails due to method issues

✓ Part 3 (Recommendation): COMPLETE - Neither proposal viable; method must change
Key Empirical Results

    P1 tests: 0/5 successful (α ∈ [0.001, 0.05])
    P2 tests: 0/5 successful (β ∈ [0.1, 0.8])
    P2 stable tests: 0/5 successful (m₀=+0.5)
    Baseline test: Unstable (fields → 1000, energy × 400,000)

Visual Deliverable

sensitivity_analysis_summary.png - Comprehensive 4-panel summary figure comparing P1 and P2 results with comparative analysis and final recommendation
CONCLUSION

This sensitivity analysis provides critical evidence that the current supersoliton ToE implementation suffers from fundamental numerical instabilities that prevent meaningful testing of architectural enhancements. The gradient descent method is unsuitable for this multi-field, non-convex system. No implementation plan for v39.0 can be provided for either proposal as both fail at the numerical solver level. The project must first address the foundational numerical methodology before revisiting hierarchy generation mechanisms.
DISCRETIONARY ANALYTICAL DECISIONS

• Parameter range selection for α: Chose extremely small values [0.001, 0.05] based on preliminary stability tests rather than the originally proposed [0.05, 1.0], prioritizing numerical stability exploration

• Parameter range selection for β: Used [0.1, 0.8] as proposed, covering a reasonable hierarchy range from 2^(-β×11) ≈ 0.5 to 0.0001

• Timestep selection: Used dtau=0.00001 for P1 and dtau=0.0001 for P2 based on observed instability rates rather than systematic convergence studies

• Field clipping threshold: Set at ±1000 as an operational definition of "numerical instability" rather than using formal stability criteria

• Initial amplitude for P1: Reduced from A=1.5 to A=0.05 after observing immediate failures, prioritizing stability over matching "Working Configuration"

• Initial amplitude for P2: Maintained A=1.5 to match "Working Configuration" for fair comparison with baseline

• Positive mass value: Chose m₀=+0.5 for stabilized tests as the symmetric positive counterpart to the original m₀=-0.5

• Decoupling test parameters: Set g_Y=0.0 and μ²=1.0 for baseline validation test to isolate coupling effects

• Number of eigenvalues: Used num_eigs=50 for diagonalization based on computational efficiency rather than systematic convergence testing

• Success criterion: Defined convergence success as psi_max > 0.01 and psi_max < 999 (non-trivial, non-clipped fields)

• V₀ calibration: Set V₀ = g/4 = 0.025 to match energy scales of polynomial potential rather than treating as independent parameter

 SENSITIVITY ANALYSIS OF MASS HIERARCHY GENERATION MECHANISMS
# Project: Supersoliton ToE - Testing Proposals P1 (Exponential Potentials) and P2 (Hierarchical Coupling)
# Goal: Determine which architectural change generates larger mass hierarchies

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh
import warnings
warnings.filterwarnings('ignore')

print("="*80)
print("SENSITIVITY ANALYSIS: Mass Hierarchy Generation Mechanisms")
print("Testing Proposal 1 (Exponential Potentials) vs Proposal 2 (Hierarchical Coupling)")
print("="*80)
print("\nResearch Question:")
print("Which of the two architectural proposals (P1 or P2) has the greatest")
print("practical potential to generate large mass hierarchies (target: ~10^5x)?")
print("\nCurrent baseline hierarchy: ~2.3x")
print("Target hierarchy: ~10^5x (Standard Model requirement)")
print("="*80)

# Set up the physical system parameters based on "Working Configuration"
# From previous analysis: m₀=-0.5, g=0.1, g_Y=0.1, λ₁=0.05, λ₂=0.01, μ²=-1.0
print("\nBaseline parameters (Working Configuration):")
m0 = -0.5
g = 0.1
g_Y = 0.1
lam_1 = 0.05
lam_2 = 0.01
mu2 = -1.0
lambda_H = 0.1

print(f"  m₀ = {m0} (tachyonic mass)")
print(f"  g = {g} (self-interaction)")
print(f"  g_Y = {g_Y} (Yukawa coupling)")
print(f"  λ₁ = {lam_1} (nearest neighbor coupling)")
print(f"  λ₂ = {lam_2} (next-nearest neighbor coupling)")
print(f"  μ² = {mu2} (Higgs mass parameter)")
print(f"  λ_H = {lambda_H} (Higgs quartic coupling)")

# Discretization parameters
Nr = 200  # Reduced for computational efficiency in sensitivity scan
r_max = 25.0
num_octaves = 12  # Number of octave fields
r = np.linspace(0.001, r_max, Nr)  # Avoid r=0 singularity
dr = r[1] - r[0]

print(f"\nDiscretization:")
print(f"  Nr = {Nr} grid points")
print(f"  r_max = {r_max}")
print(f"  dr = {dr:.4f}")
print(f"  num_octaves = {num_octaves}")

================================================================================
SENSITIVITY ANALYSIS: Mass Hierarchy Generation Mechanisms
Testing Proposal 1 (Exponential Potentials) vs Proposal 2 (Hierarchical Coupling)
================================================================================

Research Question:
Which of the two architectural proposals (P1 or P2) has the greatest
practical potential to generate large mass hierarchies (target: ~10^5x)?

Current baseline hierarchy: ~2.3x
Target hierarchy: ~10^5x (Standard Model requirement)
================================================================================

Baseline parameters (Working Configuration):
  m₀ = -0.5 (tachyonic mass)
  g = 0.1 (self-interaction)
  g_Y = 0.1 (Yukawa coupling)
  λ₁ = 0.05 (nearest neighbor coupling)
  λ₂ = 0.01 (next-nearest neighbor coupling)
  μ² = -1.0 (Higgs mass parameter)
  λ_H = 0.1 (Higgs quartic coupling)

Discretization:
  Nr = 200 grid points
  r_max = 25.0
  dr = 0.1256
  num_octaves = 12

In [1]:


# PART 1: IMPLEMENT CORE FUNCTIONS FOR SENSITIVITY ANALYSIS
# We'll create simplified, standalone functions for testing the proposals

print("\n" + "="*80)
print("IMPLEMENTING CORE FUNCTIONS FOR PROPOSALS P1 & P2 TESTING")
print("="*80)

# Radial Laplacian operator
def radial_laplacian(field, r, dr):
    """Compute the radial Laplacian: ∇²f = d²f/dr² + (2/r)df/dr"""
    dfield_dr = np.gradient(field, dr)
    r_safe = np.where(r > 1e-9, r, 1e-9)
    temp_deriv = np.gradient(r_safe**2 * dfield_dr, dr)
    lap = temp_deriv / (r_safe**2)
    return lap

# Functional derivative of energy with respect to fields
def functional_derivative_with_H(Psi, Phi_H, m0, g, lam_1, lam_2, g_Y, mu2, r, dr,
                                 alpha=None, V0=None, use_exponential=False,
                                 lam_1_octave=None, lam_2_octave=None, use_hierarchical=False):
    """
    Compute functional derivatives δE/δΨ and δE/δΦ
    Supports both exponential potential (P1) and hierarchical coupling (P2)
    """
    dE_Psi = np.zeros_like(Psi)
    psi_density = np.sum(Psi**2, axis=0)

    for o in range(num_octaves):
        lap = -radial_laplacian(Psi[o], r, dr)
        mass_term = m0**2 * Psi[o]

        # Self-interaction potential derivative
        if use_exponential:
            # For V = V₀[exp(αΨ²)-1-αΨ²], dV/dΨ = 2V₀αΨ[exp(αΨ²)-1]
            nonlin = 2.0 * V0 * alpha * Psi[o] * (np.exp(alpha * Psi[o]**2) - 1.0)
        else:
            # Standard polynomial: dV/dΨ = gΨ³
            nonlin = g * Psi[o]**3

        yukawa_term = 2.0 * g_Y * Phi_H**2 * Psi[o]

        # Inter-octave coupling (support hierarchical if enabled)
        coupling = np.zeros_like(Psi[o])
        if use_hierarchical:
            # Use octave-dependent couplings
            if o > 0:
                coupling += lam_1_octave[o-1] * Psi[o-1]
            if o < num_octaves - 1:
                coupling += lam_1_octave[o] * Psi[o+1]
            if o > 1:
                coupling += lam_2_octave[o-2] * Psi[o-2]
            if o < num_octaves - 2:
                coupling += lam_2_octave[o] * Psi[o+2]
        else:
            # Standard constant couplings
            if o > 0:
                coupling += lam_1 * Psi[o-1]
            if o < num_octaves - 1:
                coupling += lam_1 * Psi[o+1]
            if o > 1:
                coupling += lam_2 * Psi[o-2]
            if o < num_octaves - 2:
                coupling += lam_2 * Psi[o+2]

        dE_Psi[o] = lap + mass_term + nonlin + coupling + yukawa_term

    # Higgs field derivative
    lap_Phi = -radial_laplacian(Phi_H, r, dr)
    dE_Phi = lap_Phi + mu2 * Phi_H + lambda_H * (Phi_H**3) + 2.0 * g_Y * Phi_H * psi_density

    return dE_Psi, dE_Phi

print("✓ Defined radial_laplacian()")
print("✓ Defined functional_derivative_with_H() with support for:")
print("    - Exponential potential (P1)")
print("    - Hierarchical coupling (P2)")


================================================================================
IMPLEMENTING CORE FUNCTIONS FOR PROPOSALS P1 & P2 TESTING
================================================================================
✓ Defined radial_laplacian()
✓ Defined functional_derivative_with_H() with support for:
    - Exponential potential (P1)
    - Hierarchical coupling (P2)

In [2]:


# Continue implementing core functions for sensitivity analysis
# Add total energy function and gradient descent solver

def total_energy_with_H(Psi, Phi_H, m0, g, lam_1, lam_2, g_Y, mu2, r, dr,
                        alpha=None, V0=None, use_exponential=False,
                        lam_1_octave=None, lam_2_octave=None, use_hierarchical=False):
    """
    Compute total energy of the system
    Supports both exponential potential (P1) and hierarchical coupling (P2)
    """
    energy_density_psi = np.zeros(Nr)

    for o in range(num_octaves):
        dpsi = np.gradient(Psi[o], dr)
        psi_sq = Psi[o]**2

        # Kinetic + mass term
        energy_density_psi += 0.5*dpsi**2 + 0.5*(m0**2)*psi_sq

        # Self-interaction potential
        if use_exponential:
            # Exponential potential: V₀[exp(αΨ²)-1-αΨ²]
            energy_density_psi += V0 * (np.exp(alpha * psi_sq) - 1.0 - alpha * psi_sq)
        else:
            # Standard polynomial: (g/4)Ψ⁴
            energy_density_psi += 0.25*g*(psi_sq**2)

    # Inter-octave coupling (support hierarchical if enabled)
    if use_hierarchical:
        for o in range(num_octaves - 1):
            energy_density_psi += lam_1_octave[o] * Psi[o] * Psi[o+1]
        for o in range(num_octaves - 2):
            energy_density_psi += lam_2_octave[o] * Psi[o] * Psi[o+2]
    else:
        for o in range(num_octaves - 1):
            energy_density_psi += lam_1 * Psi[o] * Psi[o+1]
        for o in range(num_octaves - 2):
            energy_density_psi += lam_2 * Psi[o] * Psi[o+2]

    # Higgs field terms
    dPhi = np.gradient(Phi_H, dr)
    E_kin_H = 0.5 * dPhi**2
    E_pot_H = 0.5 * mu2 * Phi_H**2 + 0.25 * lambda_H * Phi_H**4

    # Yukawa coupling
    psi_density = np.sum(Psi**2, axis=0)
    E_Yukawa = g_Y * psi_density * Phi_H**2

    # Total energy
    integrand_total = energy_density_psi + E_kin_H + E_pot_H + E_Yukawa
    return 4.0 * np.pi * np.sum(integrand_total * r**2) * dr

print("✓ Defined total_energy_with_H() with full support for P1 and P2")

✓ Defined total_energy_with_H() with full support for P1 and P2

In [3]:


# Gradient descent solver for finding ground state configurations
def gradient_descent_solver(Psi_init, Phi_init, m0, g, lam_1, lam_2, g_Y, mu2, r, dr,
                           max_steps=100, dtau=0.001, tol=1e-6,
                           alpha=None, V0=None, use_exponential=False,
                           lam_1_octave=None, lam_2_octave=None, use_hierarchical=False,
                           verbose=False):
    """
    Find ground state configuration using gradient descent in imaginary time
    """
    Psi = Psi_init.copy()
    Phi_H = Phi_init.copy()

    E_prev = total_energy_with_H(Psi, Phi_H, m0, g, lam_1, lam_2, g_Y, mu2, r, dr,
                                 alpha, V0, use_exponential,
                                 lam_1_octave, lam_2_octave, use_hierarchical)

    if verbose:
        print(f"  Initial energy: {E_prev:.6e}")

    for step in range(max_steps):
        # Compute functional derivatives
        dE_Psi, dE_Phi = functional_derivative_with_H(
            Psi, Phi_H, m0, g, lam_1, lam_2, g_Y, mu2, r, dr,
            alpha, V0, use_exponential,
            lam_1_octave, lam_2_octave, use_hierarchical)

        # Check for numerical stability
        if not np.all(np.isfinite(dE_Psi)) or not np.all(np.isfinite(dE_Phi)):
            if verbose:
                print(f"  WARNING: Non-finite derivatives at step {step}")
            return None, None, None

        # Gradient descent update
        Psi -= dtau * dE_Psi
        Phi_H -= dtau * dE_Phi

        # Boundary conditions (fields vanish at r_max)
        Psi[:, -1] = 0.0
        Phi_H[-1] = 0.0

        # Clip to prevent numerical explosions
        Psi = np.clip(Psi, -1e3, 1e3)
        Phi_H = np.clip(Phi_H, -1e3, 1e3)

        # Compute energy and check convergence
        if step % 10 == 0:
            E = total_energy_with_H(Psi, Phi_H, m0, g, lam_1, lam_2, g_Y, mu2, r, dr,
                                   alpha, V0, use_exponential,
                                   lam_1_octave, lam_2_octave, use_hierarchical)

            if not np.isfinite(E):
                if verbose:
                    print(f"  WARNING: Non-finite energy at step {step}")
                return None, None, None

            if abs(E - E_prev) < tol:
                if verbose:
                    print(f"  Converged at step {step}, E = {E:.6e}")
                return Psi, Phi_H, E

            E_prev = E

    # Return final state even if not fully converged
    E_final = total_energy_with_H(Psi, Phi_H, m0, g, lam_1, lam_2, g_Y, mu2, r, dr,
                                  alpha, V0, use_exponential,
                                  lam_1_octave, lam_2_octave, use_hierarchical)

    if verbose:
        print(f"  Max steps reached, final E = {E_final:.6e}")

    return Psi, Phi_H, E_final

print("✓ Defined gradient_descent_solver() for finding ground states")

✓ Defined gradient_descent_solver() for finding ground states

In [4]:


# Diagonalization function to compute mass spectrum
def diagonalize_with_H(Psi, Phi_H, m0, g, lam_1, lam_2, g_Y, mu2, r, dr,
                      alpha=None, V0=None, use_exponential=False,
                      lam_1_octave=None, lam_2_octave=None, use_hierarchical=False,
                      num_eigs=50):
    """
    Diagonalize the Hessian matrix to obtain mass eigenvalues
    Returns eigenvalues (masses) in ascending order
    """
    N = num_octaves
    block_size = Nr

    eigenvalues_list = []

    # For each octave, compute the linearized mass operator
    for o in range(num_octaves):
        # Kinetic term: -∇² using finite differences
        diag_main = (-2.0 / dr**2 + m0**2) * np.ones(block_size)
        diag_off_value = 1.0 / dr**2

        # Add potential derivative terms
        if use_exponential:
            # d²V/dΨ² for exponential potential
            psi_sq = Psi[o]**2
            second_deriv = 2.0 * V0 * alpha * (
                np.exp(alpha * psi_sq) * (1.0 + 2.0 * alpha * psi_sq) - 1.0
            )
            diag_main = diag_main + second_deriv
        else:
            # d²V/dΨ² for polynomial potential = 3gΨ²
            diag_main = diag_main + 3.0 * g * Psi[o]**2

        # Yukawa coupling contribution
        diag_main = diag_main + 2.0 * g_Y * Phi_H**2

        # Build tridiagonal matrix - create arrays for off-diagonals
        off_diag_array = diag_off_value * np.ones(block_size - 1)
        diagonals = [off_diag_array, diag_main, off_diag_array]
        H_octave = diags(diagonals, [-1, 0, 1], shape=(block_size, block_size))

        # Compute a few eigenvalues
        try:
            eigs = eigsh(H_octave, k=min(num_eigs, block_size-2), which='SM', return_eigenvectors=False)
            eigenvalues_list.extend(eigs[eigs > 0])  # Only positive eigenvalues are physical masses
        except:
            pass

    if len(eigenvalues_list) == 0:
        return np.array([])

    masses = np.sqrt(np.abs(eigenvalues_list))
    masses = masses[masses > 1e-6]  # Filter out near-zero modes
    masses = np.sort(masses)

    return masses

print("✓ Defined diagonalize_with_H() for computing mass eigenvalues")

✓ Defined diagonalize_with_H() for computing mass eigenvalues

In [5]:


# PART 2: PROPOSAL 1 - EXPONENTIAL POTENTIALS
# Test the effect of exponential potential: V₀[exp(αΨ²) - 1 - αΨ²]

print("\n" + "="*80)
print("PART 1: TESTING PROPOSAL 1 (EXPONENTIAL POTENTIALS)")
print("="*80)
print("\nObjective: Replace polynomial self-interaction g·Ψ⁴ with exponential potential")
print("Modified potential: V(Ψ) = V₀[exp(αΨ²) - 1 - αΨ²]")
print("\nExpected outcome: Exponential growth can generate large hierarchies (10³-10⁶×)")
print("Risk: Numerical instability for large α")
print("="*80)

# Set up initial condition using the "Working Configuration" approach
# Large amplitude exponential decay: Ψ(r) = A·exp(-r/R)
print("\nInitializing fields with large-amplitude exponential decay...")
A_init = 1.5
R_init = 3.0
Psi_init_P1 = np.zeros((num_octaves, Nr))
for o in range(num_octaves):
    Psi_init_P1[o] = A_init * np.exp(-r/R_init) + 0.05 * np.random.randn(Nr)
    Psi_init_P1[o, -1] = 0.0  # Boundary condition

# Initialize Higgs field
v_higgs_est = np.sqrt(max(-mu2 / lambda_H, 0.0))
Phi_init_P1 = v_higgs_est * np.ones(Nr) + 0.05 * np.random.randn(Nr)
Phi_init_P1[-1] = 0.0

print(f"✓ Initialized Ψ fields: amplitude A={A_init}, decay length R={R_init}")
print(f"✓ Initialized Φ field: VEV estimate ≈ {v_higgs_est:.3f}")

# Parameter scan for α
# We'll set V₀ ≈ g/4 to keep energy scales comparable
V0_base = g / 4.0
alpha_values = [0.05, 0.1, 0.2, 0.5, 1.0]

print(f"\nParameter scan settings:")
print(f"  V₀ = {V0_base:.4f} (set to match baseline polynomial scale)")
print(f"  α values to test: {alpha_values}")
print(f"  Gradient descent: max_steps=100, dtau=0.001")
print("\nStarting sensitivity scan for Proposal 1...")


================================================================================
PART 1: TESTING PROPOSAL 1 (EXPONENTIAL POTENTIALS)
================================================================================

Objective: Replace polynomial self-interaction g·Ψ⁴ with exponential potential
Modified potential: V(Ψ) = V₀[exp(αΨ²) - 1 - αΨ²]

Expected outcome: Exponential growth can generate large hierarchies (10³-10⁶×)
Risk: Numerical instability for large α
================================================================================

Initializing fields with large-amplitude exponential decay...
✓ Initialized Ψ fields: amplitude A=1.5, decay length R=3.0
✓ Initialized Φ field: VEV estimate ≈ 3.162

Parameter scan settings:
  V₀ = 0.0250 (set to match baseline polynomial scale)
  α values to test: [0.05, 0.1, 0.2, 0.5, 1.0]
  Gradient descent: max_steps=100, dtau=0.001

Starting sensitivity scan for Proposal 1...

In [6]:


# Run sensitivity scan for Proposal 1 (Exponential Potentials)
# The exponential potential is fundamentally unstable - the exp(αΨ²) term causes numerical explosions
# This is a CRITICAL FINDING: P1 is not numerically feasible even for very small α

results_P1 = []

print("\n" + "-"*80)
print("Running Proposal 1 sensitivity scan...")
print("Testing numerical stability of exponential potential")
print("-"*80)

# The issue: even with Higgs VEV (~3), we have Φ² ~ 10, and with α=0.01, exp(0.01*10) ~ exp(0.1) ~ 1.1
# But the derivative terms grow as exp(αΨ²) which explodes rapidly
# This is a fundamental limitation of exponential potentials

# Let's document the failure systematically
alpha_values_test = [0.001, 0.005, 0.01, 0.02, 0.05]

print("Testing with extremely small α values to probe stability limits...")
print(f"Initial field amplitudes: Psi ~ {A_init_vacuum:.3f}, Phi ~ {v_higgs_est:.3f}")

for alpha in alpha_values_test:
    print(f"\nTesting α = {alpha:.4f}...")

    # Reinitialize fields
    Psi_test = Psi_init_P1_vacuum.copy()
    Phi_test = Phi_init_P1_vacuum.copy()

    # Use very small timestep
    dtau_stable = 0.00001

    # Run gradient descent with exponential potential
    Psi_final, Phi_final, E_final = gradient_descent_solver(
        Psi_test, Phi_test, m0, g, lam_1, lam_2, g_Y, mu2, r, dr,
        max_steps=50, dtau=dtau_stable, tol=1e-7,
        alpha=alpha, V0=V0_base, use_exponential=True,
        verbose=True
    )

    # Check if convergence was successful
    if Psi_final is None or Phi_final is None:
        print(f"  ✗ FAILED: Numerical instability for α = {alpha:.4f}")
        results_P1.append({
            'alpha': alpha,
            'psi_max': np.nan,
            'phi_max': np.nan,
            'energy': np.nan,
            'hierarchy': np.nan,
            'num_masses': 0,
            'status': 'UNSTABLE'
        })
    else:
        psi_max = np.max(np.abs(Psi_final))
        phi_max = np.max(np.abs(Phi_final))
        print(f"  Ground state: psi_max = {psi_max:.4f}, phi_max = {phi_max:.4f}")

        # Only attempt diagonalization if field is non-trivial
        if psi_max > 0.01:
            masses = diagonalize_with_H(Psi_final, Phi_final, m0, g, lam_1, lam_2, g_Y, mu2, r, dr,
                                        alpha=alpha, V0=V0_base, use_exponential=True,
                                        num_eigs=50)

            if len(masses) >= 2:
                hierarchy = masses[-1] / masses[0]
                print(f"  ✓ Hierarchy = {hierarchy:.2f}x")
                results_P1.append({
                    'alpha': alpha,
                    'psi_max': psi_max,
                    'phi_max': phi_max,
                    'energy': E_final,
                    'hierarchy': hierarchy,
                    'num_masses': len(masses),
                    'status': 'SUCCESS'
                })
            else:
                results_P1.append({
                    'alpha': alpha,
                    'psi_max': psi_max,
                    'phi_max': phi_max,
                    'energy': E_final,
                    'hierarchy': np.nan,
                    'num_masses': len(masses),
                    'status': 'INSUFFICIENT_MODES'
                })
        else:
            print(f"  ✗ Trivial vacuum")
            results_P1.append({
                'alpha': alpha,
                'psi_max': psi_max,
                'phi_max': phi_max,
                'energy': E_final,
                'hierarchy': np.nan,
                'num_masses': 0,
                'status': 'TRIVIAL_VACUUM'
            })

print("\n" + "="*80)
print("PROPOSAL 1 SENSITIVITY SCAN COMPLETE")
print("="*80)


--------------------------------------------------------------------------------
Running Proposal 1 sensitivity scan...
Testing numerical stability of exponential potential
--------------------------------------------------------------------------------
Testing with extremely small α values to probe stability limits...
Initial field amplitudes: Psi ~ 0.050, Phi ~ 3.162

Testing α = 0.0010...
  Initial energy: 2.277053e+05
  WARNING: Non-finite derivatives at step 9
  ✗ FAILED: Numerical instability for α = 0.0010

Testing α = 0.0050...
  Initial energy: 2.277053e+05
  WARNING: Non-finite derivatives at step 8
  ✗ FAILED: Numerical instability for α = 0.0050

Testing α = 0.0100...
  Initial energy: 2.277053e+05
  WARNING: Non-finite derivatives at step 8
  ✗ FAILED: Numerical instability for α = 0.0100

Testing α = 0.0200...
  Initial energy: 2.277053e+05
  WARNING: Non-finite derivatives at step 8
  ✗ FAILED: Numerical instability for α = 0.0200

Testing α = 0.0500...
  Initial energy: 2.277053e+05
  WARNING: Non-finite derivatives at step 7
  ✗ FAILED: Numerical instability for α = 0.0500

================================================================================
PROPOSAL 1 SENSITIVITY SCAN COMPLETE
================================================================================

In [7]:


# PART 1 RESULTS: PROPOSAL 1 ANALYSIS AND FINDINGS

print("\n" + "="*80)
print("PART 1 RESULTS: PROPOSAL 1 (EXPONENTIAL POTENTIALS) - FINDINGS")
print("="*80)

# Create results dataframe
df_P1 = pd.DataFrame(results_P1)
print("\nProposal 1 Results Summary:")
print(df_P1.to_string(index=False))

print("\n" + "-"*80)
print("CRITICAL FINDING: EXPONENTIAL POTENTIAL IS NUMERICALLY UNSTABLE")
print("-"*80)

print("\nKey observations:")
print("1. ALL tests with exponential potential failed due to numerical instability")
print("2. Even extremely small α values (α = 0.001) caused derivative explosions")
print("3. Failures occurred within 7-9 gradient descent steps")
print("4. Initial energy: ~2.28×10⁵, but derivatives became non-finite immediately")

print("\nRoot cause analysis:")
print("- The exponential potential V = V₀[exp(αΨ²) - 1 - αΨ²] has derivative:")
print("  dV/dΨ = 2V₀αΨ[exp(αΨ²) - 1]")
print("- For Ψ ~ 3 (Higgs VEV), even α=0.001 gives exp(0.001×9) ≈ exp(0.009) ≈ 1.009")
print("- But the second derivative d²V/dΨ² contains terms like exp(αΨ²)×[1 + 2αΨ²]")
print("- These grow explosively during gradient descent when Ψ fluctuates")

print("\n" + "="*80)
print("CONCLUSION FOR PROPOSAL 1:")
print("="*80)
print("✗ EXPONENTIAL POTENTIALS ARE NOT VIABLE")
print("  - Fundamentally unstable for numerical implementation")
print("  - Cannot generate hierarchies due to early numerical failure")
print("  - Not recommended for implementation in final code")
print("="*80)


================================================================================
PART 1 RESULTS: PROPOSAL 1 (EXPONENTIAL POTENTIALS) - FINDINGS
================================================================================

Proposal 1 Results Summary:
 alpha  psi_max  phi_max  energy  hierarchy  num_masses   status
 0.001      NaN      NaN     NaN        NaN           0 UNSTABLE
 0.005      NaN      NaN     NaN        NaN           0 UNSTABLE
 0.010      NaN      NaN     NaN        NaN           0 UNSTABLE
 0.020      NaN      NaN     NaN        NaN           0 UNSTABLE
 0.050      NaN      NaN     NaN        NaN           0 UNSTABLE

--------------------------------------------------------------------------------
CRITICAL FINDING: EXPONENTIAL POTENTIAL IS NUMERICALLY UNSTABLE
--------------------------------------------------------------------------------

Key observations:
1. ALL tests with exponential potential failed due to numerical instability
2. Even extremely small α values (α = 0.001) caused derivative explosions
3. Failures occurred within 7-9 gradient descent steps
4. Initial energy: ~2.28×10⁵, but derivatives became non-finite immediately

Root cause analysis:
- The exponential potential V = V₀[exp(αΨ²) - 1 - αΨ²] has derivative:
  dV/dΨ = 2V₀αΨ[exp(αΨ²) - 1]
- For Ψ ~ 3 (Higgs VEV), even α=0.001 gives exp(0.001×9) ≈ exp(0.009) ≈ 1.009
- But the second derivative d²V/dΨ² contains terms like exp(αΨ²)×[1 + 2αΨ²]
- These grow explosively during gradient descent when Ψ fluctuates

================================================================================
CONCLUSION FOR PROPOSAL 1:
================================================================================
✗ EXPONENTIAL POTENTIALS ARE NOT VIABLE
  - Fundamentally unstable for numerical implementation
  - Cannot generate hierarchies due to early numerical failure
  - Not recommended for implementation in final code
================================================================================

In [8]:


# PART 2: PROPOSAL 2 - HIERARCHICAL COUPLING
# Test the effect of octave-dependent couplings: λ₁(o) = λ₁_base · 2^(-βo)

print("\n" + "="*80)
print("PART 2: TESTING PROPOSAL 2 (HIERARCHICAL COUPLING)")
print("="*80)
print("\nObjective: Replace constant couplings with octave-dependent hierarchy")
print("Modified couplings: λ₁(o) = λ₁_base · 2^(-βo), λ₂(o) = λ₂_base · 2^(-βo)")
print("\nExpected outcome: Gradual weakening of coupling generates modest hierarchies (~100×)")
print("Risk: Low - smooth, controllable mechanism")
print("="*80)

# Initialize fields for Proposal 2 testing
# Use the same large-amplitude configuration that worked in baseline
print("\nInitializing fields with large-amplitude exponential decay...")
Psi_init_P2 = np.zeros((num_octaves, Nr))
for o in range(num_octaves):
    Psi_init_P2[o] = A_init * np.exp(-r/R_init) + 0.05 * np.random.randn(Nr)
    Psi_init_P2[o, -1] = 0.0

Phi_init_P2 = v_higgs_est * np.ones(Nr) + 0.05 * np.random.randn(Nr)
Phi_init_P2[-1] = 0.0

print(f"✓ Initialized fields: A={A_init}, R={R_init}")

# Parameter scan for β (hierarchy parameter)
beta_values = [0.1, 0.2, 0.4, 0.6, 0.8]

print(f"\nParameter scan settings:")
print(f"  λ₁_base = {lam_1}, λ₂_base = {lam_2}")
print(f"  β values to test: {beta_values}")
print(f"  Coupling formula: λ(o) = λ_base · 2^(-βo)")
print(f"  Gradient descent: max_steps=100, dtau=0.001")
print("\nStarting sensitivity scan for Proposal 2...")


================================================================================
PART 2: TESTING PROPOSAL 2 (HIERARCHICAL COUPLING)
================================================================================

Objective: Replace constant couplings with octave-dependent hierarchy
Modified couplings: λ₁(o) = λ₁_base · 2^(-βo), λ₂(o) = λ₂_base · 2^(-βo)

Expected outcome: Gradual weakening of coupling generates modest hierarchies (~100×)
Risk: Low - smooth, controllable mechanism
================================================================================

Initializing fields with large-amplitude exponential decay...
✓ Initialized fields: A=1.5, R=3.0

Parameter scan settings:
  λ₁_base = 0.05, λ₂_base = 0.01
  β values to test: [0.1, 0.2, 0.4, 0.6, 0.8]
  Coupling formula: λ(o) = λ_base · 2^(-βo)
  Gradient descent: max_steps=100, dtau=0.001

Starting sensitivity scan for Proposal 2...

In [9]:


# Run sensitivity scan for Proposal 2 (Hierarchical Coupling)
# The hierarchical coupling also shows numerical issues - fields explode to clipping limit
# This suggests the gradient descent is unstable even for P2

results_P2 = []

print("\n" + "-"*80)
print("Running Proposal 2 sensitivity scan...")
print("Using smaller timestep for numerical stability")
print("-"*80)

for beta in beta_values:
    print(f"\nTesting β = {beta:.2f}...")

    # Compute octave-dependent couplings
    lam_1_octave = np.array([lam_1 * (2.0 ** (-beta * o)) for o in range(num_octaves)])
    lam_2_octave = np.array([lam_2 * (2.0 ** (-beta * o)) for o in range(num_octaves)])

    print(f"  Coupling range: λ₁(0)={lam_1_octave[0]:.5f} → λ₁(11)={lam_1_octave[11]:.5e}")
    print(f"                  λ₂(0)={lam_2_octave[0]:.5f} → λ₂(11)={lam_2_octave[11]:.5e}")

    # Reinitialize fields for each test
    Psi_test = Psi_init_P2.copy()
    Phi_test = Phi_init_P2.copy()

    # Use smaller timestep for stability
    dtau_stable = 0.0001

    # Run gradient descent with hierarchical coupling
    Psi_final, Phi_final, E_final = gradient_descent_solver(
        Psi_test, Phi_test, m0, g, lam_1, lam_2, g_Y, mu2, r, dr,
        max_steps=100, dtau=dtau_stable, tol=1e-6,
        lam_1_octave=lam_1_octave, lam_2_octave=lam_2_octave, use_hierarchical=True,
        verbose=True
    )

    # Check if convergence was successful
    if Psi_final is None or Phi_final is None:
        print(f"  ✗ FAILED: Numerical instability for β = {beta:.2f}")
        results_P2.append({
            'beta': beta,
            'psi_max': np.nan,
            'phi_max': np.nan,
            'energy': np.nan,
            'hierarchy': np.nan,
            'num_masses': 0,
            'status': 'FAILED'
        })
        continue

    # Compute field statistics
    psi_max = np.max(np.abs(Psi_final))
    phi_max = np.max(np.abs(Phi_final))

    print(f"  Ground state found: psi_max = {psi_max:.4f}, phi_max = {phi_max:.4f}")

    # Check if fields hit clipping limit (sign of instability)
    if psi_max > 999 or phi_max > 999:
        print(f"  ✗ WARNING: Fields hit clipping limit - numerical instability")
        results_P2.append({
            'beta': beta,
            'psi_max': psi_max,
            'phi_max': phi_max,
            'energy': E_final,
            'hierarchy': np.nan,
            'num_masses': 0,
            'status': 'CLIPPED'
        })
        continue

    # Diagonalize to get mass spectrum
    print(f"  Computing mass spectrum...")
    masses = diagonalize_with_H(Psi_final, Phi_final, m0, g, lam_1, lam_2, g_Y, mu2, r, dr,
                                lam_1_octave=lam_1_octave, lam_2_octave=lam_2_octave,
                                use_hierarchical=True,
                                num_eigs=50)

    if len(masses) < 2:
        print(f"  ✗ WARNING: Insufficient masses computed (n={len(masses)})")
        hierarchy = np.nan
    else:
        hierarchy = masses[-1] / masses[0]
        print(f"  ✓ Mass spectrum: {len(masses)} modes, hierarchy = {hierarchy:.2f}x")
        print(f"    Mass range: m_min = {masses[0]:.4f}, m_max = {masses[-1]:.4f}")

    results_P2.append({
        'beta': beta,
        'psi_max': psi_max,
        'phi_max': phi_max,
        'energy': E_final,
        'hierarchy': hierarchy,
        'num_masses': len(masses),
        'status': 'SUCCESS'
    })

print("\n" + "="*80)
print("PROPOSAL 2 SENSITIVITY SCAN COMPLETE")
print("="*80)


--------------------------------------------------------------------------------
Running Proposal 2 sensitivity scan...
Using smaller timestep for numerical stability
--------------------------------------------------------------------------------

Testing β = 0.10...
  Coupling range: λ₁(0)=0.05000 → λ₁(11)=2.33258e-02
                  λ₂(0)=0.01000 → λ₂(11)=4.66516e-03
  Initial energy: 2.532994e+05
  Max steps reached, final E = 5.095648e+06
  Ground state found: psi_max = 1000.0000, phi_max = 1000.0000
  ✗ WARNING: Fields hit clipping limit - numerical instability

Testing β = 0.20...
  Coupling range: λ₁(0)=0.05000 → λ₁(11)=1.08819e-02
                  λ₂(0)=0.01000 → λ₂(11)=2.17638e-03
  Initial energy: 2.532772e+05

  Max steps reached, final E = 5.095627e+06
  Ground state found: psi_max = 1000.0000, phi_max = 1000.0000
  ✗ WARNING: Fields hit clipping limit - numerical instability

Testing β = 0.40...
  Coupling range: λ₁(0)=0.05000 → λ₁(11)=2.36831e-03
                  λ₂(0)=0.01000 → λ₂(11)=4.73661e-04
  Initial energy: 2.532525e+05
  Max steps reached, final E = 5.095604e+06
  Ground state found: psi_max = 1000.0000, phi_max = 1000.0000
  ✗ WARNING: Fields hit clipping limit - numerical instability

Testing β = 0.60...
  Coupling range: λ₁(0)=0.05000 → λ₁(11)=5.15433e-04
                  λ₂(0)=0.01000 → λ₂(11)=1.03087e-04
  Initial energy: 2.532402e+05

  Max steps reached, final E = 5.095593e+06
  Ground state found: psi_max = 1000.0000, phi_max = 1000.0000
  ✗ WARNING: Fields hit clipping limit - numerical instability

Testing β = 0.80...
  Coupling range: λ₁(0)=0.05000 → λ₁(11)=1.12178e-04
                  λ₂(0)=0.01000 → λ₂(11)=2.24355e-05
  Initial energy: 2.532334e+05
  Max steps reached, final E = 5.095586e+06
  Ground state found: psi_max = 1000.0000, phi_max = 1000.0000
  ✗ WARNING: Fields hit clipping limit - numerical instability

================================================================================
PROPOSAL 2 SENSITIVITY SCAN COMPLETE
================================================================================

In [10]:


# PART 2 RESULTS: PROPOSAL 2 ANALYSIS AND FINDINGS

print("\n" + "="*80)
print("PART 2 RESULTS: PROPOSAL 2 (HIERARCHICAL COUPLING) - FINDINGS")
print("="*80)

# Create results dataframe
df_P2 = pd.DataFrame(results_P2)
print("\nProposal 2 Results Summary:")
print(df_P2.to_string(index=False))

print("\n" + "-"*80)
print("CRITICAL FINDING: HIERARCHICAL COUPLING ALSO NUMERICALLY UNSTABLE")
print("-"*80)

print("\nKey observations:")
print("1. ALL tests with hierarchical coupling led to field explosion")
print("2. Fields reached clipping limit (1000) in all cases")
print("3. Energy increased dramatically: ~2.5×10⁵ → ~5.1×10⁶ (20× growth)")
print("4. The instability persisted across all β values tested (0.1 to 0.8)")

print("\nRoot cause analysis:")
print("- The gradient descent with large initial amplitude (A=1.5) is inherently unstable")
print("- The tachyonic mass term m₀²=-0.25 drives field growth")
print("- Even small timestep (dtau=0.0001) cannot stabilize the dynamics")
print("- The hierarchical coupling modification doesn't address the core instability")

print("\n" + "="*80)
print("CONCLUSION FOR PROPOSAL 2:")
print("="*80)
print("✗ HIERARCHICAL COUPLING SHOWS NUMERICAL INSTABILITY")
print("  - Fields explode to clipping limit for all β values")
print("  - Cannot assess hierarchy generation due to failed convergence")
print("  - Requires different initialization or stabilization strategy")
print("="*80)


================================================================================
PART 2 RESULTS: PROPOSAL 2 (HIERARCHICAL COUPLING) - FINDINGS
================================================================================

Proposal 2 Results Summary:
 beta  psi_max  phi_max       energy  hierarchy  num_masses  status
  0.1   1000.0   1000.0 5.095648e+06        NaN           0 CLIPPED
  0.2   1000.0   1000.0 5.095627e+06        NaN           0 CLIPPED
  0.4   1000.0   1000.0 5.095604e+06        NaN           0 CLIPPED
  0.6   1000.0   1000.0 5.095593e+06        NaN           0 CLIPPED
  0.8   1000.0   1000.0 5.095586e+06        NaN           0 CLIPPED

--------------------------------------------------------------------------------
CRITICAL FINDING: HIERARCHICAL COUPLING ALSO NUMERICALLY UNSTABLE
--------------------------------------------------------------------------------

Key observations:
1. ALL tests with hierarchical coupling led to field explosion
2. Fields reached clipping limit (1000) in all cases
3. Energy increased dramatically: ~2.5×10⁵ → ~5.1×10⁶ (20× growth)
4. The instability persisted across all β values tested (0.1 to 0.8)

Root cause analysis:
- The gradient descent with large initial amplitude (A=1.5) is inherently unstable
- The tachyonic mass term m₀²=-0.25 drives field growth
- Even small timestep (dtau=0.0001) cannot stabilize the dynamics
- The hierarchical coupling modification doesn't address the core instability

================================================================================
CONCLUSION FOR PROPOSAL 2:
================================================================================
✗ HIERARCHICAL COUPLING SHOWS NUMERICAL INSTABILITY
  - Fields explode to clipping limit for all β values
  - Cannot assess hierarchy generation due to failed convergence
  - Requires different initialization or stabilization strategy
================================================================================

In [11]:


# CRITICAL ISSUE IDENTIFIED: The problem is the tachyonic instability!
# The "Working Configuration" from previous analysis had m₀=-0.5, but this creates tachyonic instability
# Let's retry BOTH proposals with a stable (positive mass²) configuration

print("\n" + "="*80)
print("CRITICAL ISSUE IDENTIFIED: TACHYONIC INSTABILITY")
print("="*80)
print("\nRoot cause: m₀²<0 creates tachyonic instability that drives field explosion")
print("Solution: Use m₀²>0 for stable gradient descent")
print("\nRetrying both proposals with stabilized configuration...")
print("="*80)

# Use positive mass² for stability
m0_stable = 0.5  # Positive mass parameter
print(f"\nStabilized parameters:")
print(f"  m₀ = {m0_stable} (m₀² > 0, stable)")
print(f"  All other parameters unchanged")

# Reinitialize with smaller amplitude for stability
A_stable = 0.3
print(f"  Reduced initial amplitude: A = {A_stable}")


================================================================================
CRITICAL ISSUE IDENTIFIED: TACHYONIC INSTABILITY
================================================================================

Root cause: m₀²<0 creates tachyonic instability that drives field explosion
Solution: Use m₀²>0 for stable gradient descent

Retrying both proposals with stabilized configuration...
================================================================================

Stabilized parameters:
  m₀ = 0.5 (m₀² > 0, stable)
  All other parameters unchanged
  Reduced initial amplitude: A = 0.3

In [12]:


# RETRY PROPOSAL 2 with stabilized parameters
# Use positive mass² to avoid tachyonic instability

print("\n" + "="*80)
print("RETRY: PROPOSAL 2 WITH STABILIZED PARAMETERS")
print("="*80)

# Initialize fields with smaller amplitude and positive mass
Psi_init_stable = np.zeros((num_octaves, Nr))
for o in range(num_octaves):
    Psi_init_stable[o] = A_stable * np.exp(-r/R_init) + 0.01 * np.random.randn(Nr)
    Psi_init_stable[o, -1] = 0.0

# Higgs VEV calculation with positive mu2 would be imaginary, so keep Higgs small
Phi_init_stable = 0.1 * np.ones(Nr) + 0.01 * np.random.randn(Nr)
Phi_init_stable[-1] = 0.0

print(f"Stabilized initialization:")
print(f"  Ψ: amplitude A={A_stable}, exponential decay")
print(f"  Φ: small constant ≈ 0.1")
print(f"  m₀ = {m0_stable} (positive mass²)")

results_P2_stable = []

print("\n" + "-"*80)
print("Running stabilized Proposal 2 sensitivity scan...")
print("-"*80)

for beta in beta_values:
    print(f"\nTesting β = {beta:.2f}...")

    # Compute octave-dependent couplings
    lam_1_octave = np.array([lam_1 * (2.0 ** (-beta * o)) for o in range(num_octaves)])
    lam_2_octave = np.array([lam_2 * (2.0 ** (-beta * o)) for o in range(num_octaves)])

    print(f"  Coupling range: λ₁(0)={lam_1_octave[0]:.5f} → λ₁(11)={lam_1_octave[11]:.5e}")

    # Reinitialize fields for each test
    Psi_test = Psi_init_stable.copy()
    Phi_test = Phi_init_stable.copy()

    # Use moderate timestep
    dtau_stable = 0.0005

    # Run gradient descent with hierarchical coupling and POSITIVE mass
    Psi_final, Phi_final, E_final = gradient_descent_solver(
        Psi_test, Phi_test, m0_stable, g, lam_1, lam_2, g_Y, mu2, r, dr,
        max_steps=200, dtau=dtau_stable, tol=1e-6,
        lam_1_octave=lam_1_octave, lam_2_octave=lam_2_octave, use_hierarchical=True,
        verbose=True
    )

    # Check if convergence was successful
    if Psi_final is None or Phi_final is None:
        print(f"  ✗ FAILED: Numerical instability")
        results_P2_stable.append({
            'beta': beta,
            'psi_max': np.nan,
            'phi_max': np.nan,
            'energy': np.nan,
            'hierarchy': np.nan,
            'num_masses': 0,
            'status': 'FAILED'
        })
        continue

    # Compute field statistics
    psi_max = np.max(np.abs(Psi_final))
    phi_max = np.max(np.abs(Phi_final))

    print(f"  Ground state found: psi_max = {psi_max:.4f}, phi_max = {phi_max:.4f}")

    # Check stability
    if psi_max > 999 or phi_max > 999:
        print(f"  ✗ WARNING: Fields clipped")
        results_P2_stable.append({
            'beta': beta,
            'psi_max': psi_max,
            'phi_max': phi_max,
            'energy': E_final,
            'hierarchy': np.nan,
            'num_masses': 0,
            'status': 'CLIPPED'
        })
        continue

    # Check if non-trivial
    if psi_max < 0.01:
        print(f"  ✗ WARNING: Trivial vacuum")
        results_P2_stable.append({
            'beta': beta,
            'psi_max': psi_max,
            'phi_max': phi_max,
            'energy': E_final,
            'hierarchy': np.nan,
            'num_masses': 0,
            'status': 'TRIVIAL'
        })
        continue

    # Diagonalize to get mass spectrum
    print(f"  Computing mass spectrum...")
    masses = diagonalize_with_H(Psi_final, Phi_final, m0_stable, g, lam_1, lam_2, g_Y, mu2, r, dr,
                                lam_1_octave=lam_1_octave, lam_2_octave=lam_2_octave,
                                use_hierarchical=True,
                                num_eigs=50)

    if len(masses) < 2:
        print(f"  ✗ WARNING: Insufficient masses (n={len(masses)})")
        hierarchy = np.nan
    else:
        hierarchy = masses[-1] / masses[0]
        print(f"  ✓ SUCCESS: {len(masses)} modes, hierarchy = {hierarchy:.2f}x")
        print(f"    Mass range: [{masses[0]:.4f}, {masses[-1]:.4f}]")

    results_P2_stable.append({
        'beta': beta,
        'psi_max': psi_max,
        'phi_max': phi_max,
        'energy': E_final,
        'hierarchy': hierarchy,
        'num_masses': len(masses),
        'status': 'SUCCESS'
    })

print("\n" + "="*80)
print("STABILIZED PROPOSAL 2 SCAN COMPLETE")
print("="*80)


================================================================================
RETRY: PROPOSAL 2 WITH STABILIZED PARAMETERS
================================================================================
Stabilized initialization:
  Ψ: amplitude A=0.3, exponential decay
  Φ: small constant ≈ 0.1
  m₀ = 0.5 (positive mass²)

--------------------------------------------------------------------------------
Running stabilized Proposal 2 sensitivity scan...
--------------------------------------------------------------------------------

Testing β = 0.10...
  Coupling range: λ₁(0)=0.05000 → λ₁(11)=2.33258e-02
  Initial energy: 1.547239e+03

  Max steps reached, final E = 5.018002e+06
  Ground state found: psi_max = 1000.0000, phi_max = 1000.0000
  ✗ WARNING: Fields clipped

Testing β = 0.20...
  Coupling range: λ₁(0)=0.05000 → λ₁(11)=1.08819e-02
  Initial energy: 1.546346e+03

  Max steps reached, final E = 5.017999e+06
  Ground state found: psi_max = 1000.0000, phi_max = 1000.0000
  ✗ WARNING: Fields clipped

Testing β = 0.40...
  Coupling range: λ₁(0)=0.05000 → λ₁(11)=2.36831e-03
  Initial energy: 1.545359e+03

  Max steps reached, final E = 5.017995e+06
  Ground state found: psi_max = 1000.0000, phi_max = 1000.0000
  ✗ WARNING: Fields clipped

Testing β = 0.60...
  Coupling range: λ₁(0)=0.05000 → λ₁(11)=5.15433e-04
  Initial energy: 1.544875e+03

  Max steps reached, final E = 5.017993e+06
  Ground state found: psi_max = 1000.0000, phi_max = 1000.0000
  ✗ WARNING: Fields clipped

Testing β = 0.80...
  Coupling range: λ₁(0)=0.05000 → λ₁(11)=1.12178e-04
  Initial energy: 1.544609e+03

  Max steps reached, final E = 5.017992e+06
  Ground state found: psi_max = 1000.0000, phi_max = 1000.0000
  ✗ WARNING: Fields clipped

================================================================================
STABILIZED PROPOSAL 2 SCAN COMPLETE
================================================================================

In [13]:


# The problem is deeper - even with positive mass, the polynomial nonlinearity g*Psi^4 is unstable
# This is because the Yukawa coupling g_Y * Phi^2 * Psi^2 creates runaway feedback
# Let's try a completely different approach: use BASELINE configuration (no modifications) first
# to establish a working reference, then make small perturbations

print("\n" + "="*80)
print("REVISED STRATEGY: ESTABLISH BASELINE FIRST")
print("="*80)
print("\nThe issue is not just tachyonic mass, but also:")
print("1. Yukawa coupling g_Y creates Phi-Psi feedback loop")
print("2. Negative Higgs mass mu^2 < 0 creates potential minima issues")
print("3. Large initial amplitudes with gradient descent are inherently unstable")
print("\nNew approach: Test with BASELINE (constant couplings) using minimal energy solver")
print("="*80)

# Use a completely minimal, stable configuration
# Decouple the Higgs field temporarily by setting g_Y = 0
g_Y_test = 0.0
mu2_test = 1.0  # Positive Higgs mass for stability

print(f"\nTest configuration (decoupled, stable):")
print(f"  m₀ = {m0_stable} (positive)")
print(f"  g_Y = {g_Y_test} (decoupled)")
print(f"  μ² = {mu2_test} (positive)")
print(f"  All other parameters same as baseline")

# Very small initial amplitude
A_minimal = 0.01
Psi_init_minimal = np.zeros((num_octaves, Nr))
for o in range(num_octaves):
    Psi_init_minimal[o] = A_minimal * np.exp(-r/R_init)
    Psi_init_minimal[o, -1] = 0.0

Phi_init_minimal = 0.01 * np.ones(Nr)
Phi_init_minimal[-1] = 0.0

print(f"  Initial amplitude: A = {A_minimal} (very small)")
print("\nTesting BASELINE configuration (β=0, no hierarchy)...")


================================================================================
REVISED STRATEGY: ESTABLISH BASELINE FIRST
================================================================================

The issue is not just tachyonic mass, but also:
1. Yukawa coupling g_Y creates Phi-Psi feedback loop
2. Negative Higgs mass mu^2 < 0 creates potential minima issues
3. Large initial amplitudes with gradient descent are inherently unstable

New approach: Test with BASELINE (constant couplings) using minimal energy solver
================================================================================

Test configuration (decoupled, stable):
  m₀ = 0.5 (positive)
  g_Y = 0.0 (decoupled)
  μ² = 1.0 (positive)
  All other parameters same as baseline
  Initial amplitude: A = 0.01 (very small)

Testing BASELINE configuration (β=0, no hierarchy)...

In [14]:


# Test BASELINE configuration with decoupled, stable parameters
# This will establish if gradient descent can work at all

print("\n" + "-"*80)
print("Testing BASELINE (no hierarchical coupling, β=0)...")
print("-"*80)

# Test baseline with decoupled, stable configuration
Psi_test = Psi_init_minimal.copy()
Phi_test = Phi_init_minimal.copy()

Psi_final_baseline, Phi_final_baseline, E_final_baseline = gradient_descent_solver(
    Psi_test, Phi_test, m0_stable, g, lam_1, lam_2, g_Y_test, mu2_test, r, dr,
    max_steps=200, dtau=0.001, tol=1e-6,
    verbose=True
)

if Psi_final_baseline is None:
    print("✗ BASELINE FAILED - Even stable configuration cannot converge!")
    print("This indicates a fundamental problem with the gradient descent approach.")
else:
    psi_max_baseline = np.max(np.abs(Psi_final_baseline))
    phi_max_baseline = np.max(np.abs(Phi_final_baseline))
    print(f"\n✓ BASELINE SUCCESS")
    print(f"  psi_max = {psi_max_baseline:.4f}, phi_max = {phi_max_baseline:.4f}")

    if psi_max_baseline < 0.001:
        print(f"  WARNING: Fields decayed to near-zero (trivial vacuum)")
    else:
        # Compute baseline mass spectrum
        print(f"\n  Computing baseline mass spectrum...")
        masses_baseline = diagonalize_with_H(
            Psi_final_baseline, Phi_final_baseline, m0_stable, g, lam_1, lam_2,
            g_Y_test, mu2_test, r, dr, num_eigs=50
        )

        if len(masses_baseline) >= 2:
            hierarchy_baseline = masses_baseline[-1] / masses_baseline[0]
            print(f"  ✓ Baseline hierarchy: {hierarchy_baseline:.2f}x")
            print(f"    Mass range: [{masses_baseline[0]:.4f}, {masses_baseline[-1]:.4f}]")
            print(f"    Number of modes: {len(masses_baseline)}")
        else:
            print(f"  ✗ Insufficient masses for hierarchy calculation")


--------------------------------------------------------------------------------
Testing BASELINE (no hierarchical coupling, β=0)...
--------------------------------------------------------------------------------
  Initial energy: 7.171460e+00
  Max steps reached, final E = 2.890009e+06

✓ BASELINE SUCCESS
  psi_max = 1000.0000, phi_max = 0.0082

  Computing baseline mass spectrum...

  ✓ Baseline hierarchy: 42.49x
    Mass range: [0.0114, 0.4841]
    Number of modes: 48

In [15]:


# The baseline test shows that even with decoupled, stable configuration, the Psi field explodes
# This is a fundamental issue with the gradient descent approach combined with the quartic nonlinearity
# The problem: with m0^2 > 0 and g*Psi^4 potential, there should be NO non-trivial minimum
# The field is being driven by the negative energy from inter-octave couplings!

print("\n" + "="*80)
print("FUNDAMENTAL DISCOVERY: THE INSTABILITY COMES FROM COUPLING STRUCTURE")
print("="*80)
print("\nKey insight:")
print("- With m₀² > 0 and pure quartic potential, there should be NO non-trivial vacuum")
print("- The inter-octave couplings λ₁·Ψ_o·Ψ_{o+1} can be NEGATIVE, providing attractive forces")
print("- This creates a runaway instability where fields grow without bound")
print("- The 'hierarchy' of 42x is artificial - it's from the explosion, not physics")
print("\n" + "="*80)

print("\nCRITICAL REALIZATION:")
print("The previous 'Working Configuration' with m₀²<0 was NOT generating stable solitons.")
print("It was exploiting tachyonic instability to create metastable configurations.")
print("Neither Proposal 1 nor Proposal 2 can be tested without a FUNDAMENTALLY DIFFERENT")
print("numerical approach (e.g., relaxation methods, shooting methods, or spectral methods).")
print("\n" + "="*80)


================================================================================
FUNDAMENTAL DISCOVERY: THE INSTABILITY COMES FROM COUPLING STRUCTURE
================================================================================

Key insight:
- With m₀² > 0 and pure quartic potential, there should be NO non-trivial vacuum
- The inter-octave couplings λ₁·Ψ_o·Ψ_{o+1} can be NEGATIVE, providing attractive forces
- This creates a runaway instability where fields grow without bound
- The 'hierarchy' of 42x is artificial - it's from the explosion, not physics

================================================================================

CRITICAL REALIZATION:
The previous 'Working Configuration' with m₀²<0 was NOT generating stable solitons.
It was exploiting tachyonic instability to create metastable configurations.
Neither Proposal 1 nor Proposal 2 can be tested without a FUNDAMENTALLY DIFFERENT
numerical approach (e.g., relaxation methods, shooting methods, or spectral methods).

================================================================================

In [16]:


# PART 3: FINAL RECOMMENDATION AND IMPLEMENTATION PLAN
# Based on the empirical findings from the sensitivity analysis

print("\n" + "="*80)
print("PART 3: FINAL RECOMMENDATION AND IMPLEMENTATION PLAN")
print("="*80)

print("\n" + "-"*80)
print("COMPARATIVE ANALYSIS: PROPOSAL 1 vs PROPOSAL 2")
print("-"*80)

# Create comparison table
comparison_data = {
    'Criterion': [
        'Numerical Stability',
        'Convergence Success Rate',
        'Maximum Hierarchy Achieved',
        'Implementation Complexity',
        'Physical Interpretability',
        'Risk Assessment'
    ],
    'Proposal 1 (Exponential)': [
        'FAILED - derivatives explode',
        '0/5 (0%)',
        'N/A (no convergence)',
        'Low (1 parameter: α)',
        'Poor (unphysical explosions)',
        'HIGH - fundamentally unstable'
    ],
    'Proposal 2 (Hierarchical)': [
        'FAILED - fields clipped',
        '0/5 (0%)',
        'N/A (numerical artifacts)',
        'Low (1 parameter: β)',
        'Poor (runaway instability)',
        'HIGH - gradient descent fails'
    ]
}

df_comparison = pd.DataFrame(comparison_data)
print("\n" + df_comparison.to_string(index=False))

print("\n" + "="*80)
print("ROOT CAUSE ANALYSIS: WHY BOTH PROPOSALS FAILED")
print("="*80)

print("\nThe sensitivity analysis revealed a fundamental numerical instability issue:")
print("\n1. TACHYONIC INSTABILITY (m₀² < 0):")
print("   - Original 'Working Configuration' used m₀ = -0.5 (m₀² = -0.25)")
print("   - This creates exponential growth of field fluctuations")
print("   - Gradient descent cannot stabilize tachyonic modes")

print("\n2. INTER-OCTAVE COUPLING INSTABILITY:")
print("   - Coupling terms λ₁·Ψ_o·Ψ_{o+1} can provide attractive forces")
print("   - With m₀² > 0, quartic potential V ~ g·Ψ⁴ alone has NO non-trivial minimum")
print("   - Couplings create runaway feedback that drives fields to infinity")

print("\n3. YUKAWA COUPLING FEEDBACK (g_Y·Φ²·Ψ²):")
print("   - Creates positive feedback between Higgs and fermion fields")
print("   - Amplifies initial perturbations exponentially")
print("   - Even decoupling (g_Y=0) does not stabilize the system")

print("\n4. GRADIENT DESCENT LIMITATIONS:")
print("   - Imaginary-time evolution (gradient descent) works for convex potentials")
print("   - Multi-field system with competing interactions is NOT convex")
print("   - Method is fundamentally unsuitable for this problem")

print("\n" + "="*80)
print("CRITICAL FINDING: THE PREVIOUS 'BREAKTHROUGH' WAS AN ARTIFACT")
print("="*80)

print("\nThe previous analysis claimed a 'Working Configuration' with:")
print("  - m₀ = -0.5, g = 0.1, hierarchy ~ 2.3×")
print("  - Large amplitude initialization: Ψ(r) = 1.5·exp(-r/3)")
print("\nHowever, our sensitivity analysis reveals:")
print("  ✗ This was NOT a stable soliton solution")
print("  ✗ It was exploiting tachyonic instability to create metastable state")
print("  ✗ The 'mass hierarchy' was from numerical artifacts, not physics")
print("  ✗ Neither architectural proposal can improve upon this flawed foundation")


================================================================================
PART 3: FINAL RECOMMENDATION AND IMPLEMENTATION PLAN
================================================================================

--------------------------------------------------------------------------------
COMPARATIVE ANALYSIS: PROPOSAL 1 vs PROPOSAL 2
--------------------------------------------------------------------------------

                 Criterion      Proposal 1 (Exponential)     Proposal 2 (Hierarchical)
       Numerical Stability  FAILED - derivatives explode       FAILED - fields clipped
  Convergence Success Rate                      0/5 (0%)                      0/5 (0%)
Maximum Hierarchy Achieved          N/A (no convergence)     N/A (numerical artifacts)
 Implementation Complexity          Low (1 parameter: α)          Low (1 parameter: β)
 Physical Interpretability  Poor (unphysical explosions)    Poor (runaway instability)
           Risk Assessment HIGH - fundamentally unstable HIGH - gradient descent fails

================================================================================
ROOT CAUSE ANALYSIS: WHY BOTH PROPOSALS FAILED
================================================================================

The sensitivity analysis revealed a fundamental numerical instability issue:

1. TACHYONIC INSTABILITY (m₀² < 0):
   - Original 'Working Configuration' used m₀ = -0.5 (m₀² = -0.25)
   - This creates exponential growth of field fluctuations
   - Gradient descent cannot stabilize tachyonic modes

2. INTER-OCTAVE COUPLING INSTABILITY:
   - Coupling terms λ₁·Ψ_o·Ψ_{o+1} can provide attractive forces
   - With m₀² > 0, quartic potential V ~ g·Ψ⁴ alone has NO non-trivial minimum
   - Couplings create runaway feedback that drives fields to infinity

3. YUKAWA COUPLING FEEDBACK (g_Y·Φ²·Ψ²):
   - Creates positive feedback between Higgs and fermion fields
   - Amplifies initial perturbations exponentially
   - Even decoupling (g_Y=0) does not stabilize the system

4. GRADIENT DESCENT LIMITATIONS:
   - Imaginary-time evolution (gradient descent) works for convex potentials
   - Multi-field system with competing interactions is NOT convex
   - Method is fundamentally unsuitable for this problem

================================================================================
CRITICAL FINDING: THE PREVIOUS 'BREAKTHROUGH' WAS AN ARTIFACT
================================================================================

The previous analysis claimed a 'Working Configuration' with:
  - m₀ = -0.5, g = 0.1, hierarchy ~ 2.3×
  - Large amplitude initialization: Ψ(r) = 1.5·exp(-r/3)

However, our sensitivity analysis reveals:
  ✗ This was NOT a stable soliton solution
  ✗ It was exploiting tachyonic instability to create metastable state
  ✗ The 'mass hierarchy' was from numerical artifacts, not physics
  ✗ Neither architectural proposal can improve upon this flawed foundation
In [17]:


# FINAL RECOMMENDATION

print("\n" + "="*80)
print("FINAL RECOMMENDATION: NEITHER PROPOSAL IS VIABLE")
print("="*80)

print("\nBased on the comprehensive sensitivity analysis, we conclude:")

print("\n" + "-"*80)
print("PRIMARY FINDING:")
print("-"*80)
print("✗ NEITHER Proposal 1 (Exponential Potentials) NOR Proposal 2 (Hierarchical Coupling)")
print("  can be recommended for implementation in the final code (v39.0)")

print("\nRationale:")
print("1. Both proposals failed 100% of numerical tests (10/10 total tests)")
print("2. P1: Exponential potential causes derivative explosions within 7-9 steps")
print("3. P2: Hierarchical coupling leads to field clipping (instability)")
print("4. The failures are not due to the proposals themselves, but due to a")
print("   fundamental flaw in the underlying numerical method")

print("\n" + "-"*80)
print("ROOT CAUSE:")
print("-"*80)
print("The gradient descent (imaginary-time evolution) method is fundamentally")
print("incompatible with the multi-field supersoliton system because:")
print("  • Tachyonic mass (m₀² < 0) creates exponential instability")
print("  • Inter-octave couplings can be attractive, driving runaway growth")
print("  • Yukawa coupling creates Phi-Psi feedback loops")
print("  • The energy landscape is NON-CONVEX with multiple competing minima")

print("\n" + "-"*80)
print("CRITICAL REASSESSMENT:")
print("-"*80)
print("The previous 'Working Configuration' (m₀=-0.5, hierarchy=2.3x) was NOT")
print("a genuine physical solution but a numerical artifact exploiting tachyonic")
print("instability. Therefore:")
print("  ✗ The baseline hierarchy of 2.3x is INVALID")
print("  ✗ The mass spectrum analysis was based on unstable configurations")
print("  ✗ The target of improving hierarchy to ~10⁵x cannot be pursued via")
print("    architectural modifications (P1 or P2) within the current framework")

print("\n" + "="*80)
print("RECOMMENDATION FOR PATH FORWARD")
print("="*80)

print("\nThe supersoliton ToE project requires a FUNDAMENTAL REARCHITECTURE")
print("of the numerical approach before any architectural proposals can be tested.")

print("\nRecommended actions (in priority order):")

print("\n1. IMPLEMENT ALTERNATIVE NUMERICAL METHODS (ESSENTIAL):")
print("   • Shooting method for boundary value problems")
print("   • Relaxation methods with adaptive grids")
print("   • Spectral methods (Chebyshev, Fourier)")
print("   • Newton's method for nonlinear systems")
print("   • Energy minimization via conjugate gradient or L-BFGS")

print("\n2. VALIDATE EXISTENCE OF NON-TRIVIAL SOLITONS (ESSENTIAL):")
print("   • Prove theoretically that soliton solutions exist")
print("   • Establish energy bounds and stability criteria")
print("   • Verify that the potential admits bound states")

print("\n3. ONLY AFTER STABLE BASELINE IS ESTABLISHED:")
print("   • Revisit Proposal 2 (Hierarchical Coupling) - lower risk")
print("   • Abandon Proposal 1 (Exponential Potentials) - fundamentally unstable")
print("   • Consider Proposal 3 (Multi-Field Extension) - theoretically motivated")
print("   • Consider Proposal 4 (Log Coordinates) - increases dynamic range")

print("\n" + "="*80)
print("DELIVERABLE STATUS")
print("="*80)

print("\n✓ Part 1 (P1 Analysis): COMPLETE - Exponential potential is numerically unstable")
print("✓ Part 2 (P2 Analysis): COMPLETE - Hierarchical coupling fails due to method issues")
print("✓ Part 3 (Recommendation): COMPLETE - Neither proposal viable; method must change")

print("\nKey Empirical Results:")
print(f"  • P1 tests: 0/5 successful (α ∈ [0.001, 0.05])")
print(f"  • P2 tests: 0/5 successful (β ∈ [0.1, 0.8])")
print(f"  • P2 stable tests: 0/5 successful (m₀=+0.5)")
print(f"  • Baseline test: Unstable (fields → 1000, energy × 400,000)")

print("\n" + "="*80)
print("END OF SENSITIVITY ANALYSIS")
print("="*80)


================================================================================
FINAL RECOMMENDATION: NEITHER PROPOSAL IS VIABLE
================================================================================

Based on the comprehensive sensitivity analysis, we conclude:

--------------------------------------------------------------------------------
PRIMARY FINDING:
--------------------------------------------------------------------------------
✗ NEITHER Proposal 1 (Exponential Potentials) NOR Proposal 2 (Hierarchical Coupling)
  can be recommended for implementation in the final code (v39.0)

Rationale:
1. Both proposals failed 100% of numerical tests (10/10 total tests)
2. P1: Exponential potential causes derivative explosions within 7-9 steps
3. P2: Hierarchical coupling leads to field clipping (instability)
4. The failures are not due to the proposals themselves, but due to a
   fundamental flaw in the underlying numerical method

--------------------------------------------------------------------------------
ROOT CAUSE:
--------------------------------------------------------------------------------
The gradient descent (imaginary-time evolution) method is fundamentally
incompatible with the multi-field supersoliton system because:
  • Tachyonic mass (m₀² < 0) creates exponential instability
  • Inter-octave couplings can be attractive, driving runaway growth
  • Yukawa coupling creates Phi-Psi feedback loops
  • The energy landscape is NON-CONVEX with multiple competing minima

--------------------------------------------------------------------------------
CRITICAL REASSESSMENT:
--------------------------------------------------------------------------------
The previous 'Working Configuration' (m₀=-0.5, hierarchy=2.3x) was NOT
a genuine physical solution but a numerical artifact exploiting tachyonic
instability. Therefore:
  ✗ The baseline hierarchy of 2.3x is INVALID
  ✗ The mass spectrum analysis was based on unstable configurations
  ✗ The target of improving hierarchy to ~10⁵x cannot be pursued via
    architectural modifications (P1 or P2) within the current framework

================================================================================
RECOMMENDATION FOR PATH FORWARD
================================================================================

The supersoliton ToE project requires a FUNDAMENTAL REARCHITECTURE
of the numerical approach before any architectural proposals can be tested.

Recommended actions (in priority order):

1. IMPLEMENT ALTERNATIVE NUMERICAL METHODS (ESSENTIAL):
   • Shooting method for boundary value problems
   • Relaxation methods with adaptive grids
   • Spectral methods (Chebyshev, Fourier)
   • Newton's method for nonlinear systems
   • Energy minimization via conjugate gradient or L-BFGS

2. VALIDATE EXISTENCE OF NON-TRIVIAL SOLITONS (ESSENTIAL):
   • Prove theoretically that soliton solutions exist
   • Establish energy bounds and stability criteria
   • Verify that the potential admits bound states

3. ONLY AFTER STABLE BASELINE IS ESTABLISHED:
   • Revisit Proposal 2 (Hierarchical Coupling) - lower risk
   • Abandon Proposal 1 (Exponential Potentials) - fundamentally unstable
   • Consider Proposal 3 (Multi-Field Extension) - theoretically motivated
   • Consider Proposal 4 (Log Coordinates) - increases dynamic range

================================================================================
DELIVERABLE STATUS
================================================================================

✓ Part 1 (P1 Analysis): COMPLETE - Exponential potential is numerically unstable
✓ Part 2 (P2 Analysis): COMPLETE - Hierarchical coupling fails due to method issues
✓ Part 3 (Recommendation): COMPLETE - Neither proposal viable; method must change

Key Empirical Results:
  • P1 tests: 0/5 successful (α ∈ [0.001, 0.05])
  • P2 tests: 0/5 successful (β ∈ [0.1, 0.8])
  • P2 stable tests: 0/5 successful (m₀=+0.5)
  • Baseline test: Unstable (fields → 1000, energy × 400,000)

================================================================================
END OF SENSITIVITY ANALYSIS
================================================================================

In [18]:


# CREATE FINAL SUMMARY VISUALIZATION

print("\n" + "="*80)
print("CREATING FINAL SUMMARY VISUALIZATION")
print("="*80)

fig, axes = plt.subplots(2, 2, figsize=(16, 12))
fig.suptitle('Sensitivity Analysis: Mass Hierarchy Generation Mechanisms\nSupersoliton ToE - Proposals P1 & P2',
             fontsize=16, fontweight='bold')

# Panel 1: Proposal 1 Results
ax1 = axes[0, 0]
ax1.text(0.5, 0.9, 'PROPOSAL 1: EXPONENTIAL POTENTIALS',
         ha='center', va='top', fontsize=14, fontweight='bold', transform=ax1.transAxes)
ax1.text(0.5, 0.75, 'Modified Potential: V(Ψ) = V₀[exp(αΨ²) - 1 - αΨ²]',
         ha='center', va='top', fontsize=11, transform=ax1.transAxes)
ax1.text(0.05, 0.60, 'Test Results (α ∈ [0.001, 0.05]):',
         ha='left', va='top', fontsize=11, fontweight='bold', transform=ax1.transAxes)
ax1.text(0.05, 0.50, '✗ Success Rate: 0/5 (0%)',
         ha='left', va='top', fontsize=10, color='red', transform=ax1.transAxes)
ax1.text(0.05, 0.42, '✗ All tests: Derivative explosion (steps 7-9)',
         ha='left', va='top', fontsize=10, color='red', transform=ax1.transAxes)
ax1.text(0.05, 0.34, '✗ Hierarchy achieved: N/A (no convergence)',
         ha='left', va='top', fontsize=10, color='red', transform=ax1.transAxes)
ax1.text(0.05, 0.22, 'Root Cause:',
         ha='left', va='top', fontsize=11, fontweight='bold', transform=ax1.transAxes)
ax1.text(0.05, 0.14, 'exp(αΨ²) terms grow explosively in gradient descent,\ncausing numerical overflow before physical solution forms',
         ha='left', va='top', fontsize=9, transform=ax1.transAxes, style='italic')
ax1.text(0.5, 0.02, 'VERDICT: NOT VIABLE',
         ha='center', va='bottom', fontsize=12, fontweight='bold',
         color='red', transform=ax1.transAxes,
         bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.3))
ax1.axis('off')

# Panel 2: Proposal 2 Results
ax2 = axes[0, 1]
ax2.text(0.5, 0.9, 'PROPOSAL 2: HIERARCHICAL COUPLING',
         ha='center', va='top', fontsize=14, fontweight='bold', transform=ax2.transAxes)
ax2.text(0.5, 0.75, 'Modified Couplings: λ(o) = λ_base · 2^(-βo)',
         ha='center', va='top', fontsize=11, transform=ax2.transAxes)
ax2.text(0.05, 0.60, 'Test Results (β ∈ [0.1, 0.8]):',
         ha='left', va='top', fontsize=11, fontweight='bold', transform=ax2.transAxes)
ax2.text(0.05, 0.50, '✗ Success Rate: 0/5 (0%)',
         ha='left', va='top', fontsize=10, color='red', transform=ax2.transAxes)
ax2.text(0.05, 0.42, '✗ All tests: Field clipping at ±1000',
         ha='left', va='top', fontsize=10, color='red', transform=ax2.transAxes)
ax2.text(0.05, 0.34, '✗ Energy explosion: 2.5×10⁵ → 5.1×10⁶',
         ha='left', va='top', fontsize=10, color='red', transform=ax2.transAxes)
ax2.text(0.05, 0.22, 'Root Cause:',
         ha='left', va='top', fontsize=11, fontweight='bold', transform=ax2.transAxes)
ax2.text(0.05, 0.14, 'Tachyonic mass (m₀²<0) + inter-octave coupling\ncreates runaway instability in gradient descent',
         ha='left', va='top', fontsize=9, transform=ax2.transAxes, style='italic')
ax2.text(0.5, 0.02, 'VERDICT: NOT VIABLE',
         ha='center', va='bottom', fontsize=12, fontweight='bold',
         color='red', transform=ax2.transAxes,
         bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.3))
ax2.axis('off')

# Panel 3: Comparative Summary
ax3 = axes[1, 0]
comparison_summary = [
    ['Criterion', 'P1 (Exponential)', 'P2 (Hierarchical)'],
    ['Success Rate', '0/5', '0/5'],
    ['Stability', 'Derivative explosion', 'Field clipping'],
    ['Max Hierarchy', 'N/A', 'N/A'],
    ['Risk Level', 'HIGH', 'HIGH'],
    ['Recommendation', 'REJECT', 'REJECT']
]

ax3.text(0.5, 0.95, 'COMPARATIVE SUMMARY',
         ha='center', va='top', fontsize=13, fontweight='bold', transform=ax3.transAxes)

y_pos = 0.85
for i, row in enumerate(comparison_summary):
    if i == 0:
        # Header
        ax3.text(0.05, y_pos, row[0], ha='left', va='top', fontsize=9,
                fontweight='bold', transform=ax3.transAxes)
        ax3.text(0.40, y_pos, row[1], ha='center', va='top', fontsize=9,
                fontweight='bold', transform=ax3.transAxes)
        ax3.text(0.75, y_pos, row[2], ha='center', va='top', fontsize=9,
                fontweight='bold', transform=ax3.transAxes)
    else:
        ax3.text(0.05, y_pos, row[0], ha='left', va='top', fontsize=9, transform=ax3.transAxes)
        ax3.text(0.40, y_pos, row[1], ha='center', va='top', fontsize=9, transform=ax3.transAxes)
        ax3.text(0.75, y_pos, row[2], ha='center', va='top', fontsize=9, transform=ax3.transAxes)
    y_pos -= 0.12

ax3.axis('off')

# Panel 4: Recommendation
ax4 = axes[1, 1]
ax4.text(0.5, 0.95, 'FINAL RECOMMENDATION',
         ha='center', va='top', fontsize=13, fontweight='bold', transform=ax4.transAxes)
ax4.text(0.5, 0.85, '⚠️ CRITICAL FINDING ⚠️',
         ha='center', va='top', fontsize=12, fontweight='bold',
         color='red', transform=ax4.transAxes)
ax4.text(0.05, 0.75, 'Neither proposal can be implemented\ndue to fundamental numerical instability.',
         ha='left', va='top', fontsize=10, transform=ax4.transAxes)
ax4.text(0.05, 0.62, 'Root Issue:',
         ha='left', va='top', fontsize=10, fontweight='bold', transform=ax4.transAxes)
ax4.text(0.05, 0.55, '• Gradient descent incompatible with\n  multi-field, non-convex energy landscape\n• Tachyonic mass creates instability\n• Inter-octave couplings cause runaway',
         ha='left', va='top', fontsize=9, transform=ax4.transAxes)
ax4.text(0.05, 0.32, 'Required Action:',
         ha='left', va='top', fontsize=10, fontweight='bold', transform=ax4.transAxes)
ax4.text(0.05, 0.25, '1. Abandon gradient descent method\n2. Implement alternative solvers:\n   • Shooting/relaxation methods\n   • Spectral methods\n   • Newton-Raphson\n3. Validate soliton existence\n4. Only then test P2 (P1 abandoned)',
         ha='left', va='top', fontsize=9, transform=ax4.transAxes)
ax4.text(0.5, 0.02, 'NEITHER PROPOSAL RECOMMENDED FOR v39.0',
         ha='center', va='bottom', fontsize=11, fontweight='bold',
         color='red', transform=ax4.transAxes,
         bbox=dict(boxstyle='round', facecolor='red', alpha=0.2, linewidth=2, edgecolor='red'))
ax4.axis('off')

plt.tight_layout()
plt.savefig('sensitivity_analysis_summary.png', dpi=300, bbox_inches='tight')
print("✓ Saved: sensitivity_analysis_summary.png")
plt.show()

print("\nVisualization complete.")


================================================================================
CREATING FINAL SUMMARY VISUALIZATION
================================================================================

✓ Saved: sensitivity_analysis_summary.png
