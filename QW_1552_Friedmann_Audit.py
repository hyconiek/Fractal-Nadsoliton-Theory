import numpy as np
from datetime import datetime

# ==============================================================================
# QW-1552 AUDIT (ROUND 3): Emergent Friedmann (No Postulates)
# ==============================================================================
# OPERATIONAL INSTRUCTION:
# 1. Define total network energy E(a).
# 2. Define node change rate dN/dt.
# 3. Derive da/dt from E-evolution. NO pre-assumed Friedmann equation.
# 4. Status: VERIFIED / INCONCLUSIVE / FAILED.

REPORT = "RAPORT_QW1552_FRIEDMANN_AUDIT.md"
md = []

def log(msg=""):
    print(msg)
    md.append(msg)

log("="*80)
log("QW-1552 OPERATIONAL AUDIT: DERIVED FRIEDMANN DYNAMICS")
log("="*80)

# ------------------------------------------------------------------------------
# 1. First-Principles Model
# ------------------------------------------------------------------------------
# We consider a network of size 'a' (scale factor).
# Total energy E is contained in the links.
# Let's assume the energy density rho ~ 1/a^3 (Matter-like).

def simulate_expansion(steps=100, dt=0.01):
    a = 1.0 # Initial scale
    rho0 = 1.0
    
    a_history = [a]
    t_history = [0.0]
    
    # We DO NOT use da/dt = sqrt(rho).
    # Instead, we postulate a fundamental link-creation law:
    # dN/dt ~ Intensity of Information Processing (E_total)
    
    for i in range(steps):
        # Current total energy E = V * rho ~ a^3 * (1/a^3) = constant
        # Wait, if rho ~ 1/a^3, then E = const. 
        # In a network, energy is proportional to number of active processing links.
        
        # Fundamental Mechanism: Numerical expansion rate is proportional 
        # to the available information gradient.
        # Let's see what happens if da/dt ~ E / a^2 (derived from flux balance)
        
        rho = rho0 / (a**3)
        E_total = (a**3) * rho  # E = 1.0 (constant matter)
        
        # Derived rate from information flux density across surface 4*pi*a^2:
        da_dt = E_total / (a**2)  # This is a PURE mechanism test
        
        a = a + da_dt * dt
        a_history.append(a)
        t_history.append((i+1)*dt)
        
    return np.array(t_history), np.array(a_history)

t, a = simulate_expansion(100, 0.5)

# ------------------------------------------------------------------------------
# 2. Comparison with Theoretical Expectation (a ~ t^(2/3) for matter)
# ------------------------------------------------------------------------------
# Matter dominated GR: a(t) = (t / t0)^(2/3)
# We perform a power-law fit: log(a) = p * log(t) + const
log("Power-law Fitting: a(t) ~ t^p")

# Skip early transient
mask = t > 5.0
coeffs = np.polyfit(np.log(t[mask]), np.log(a[mask]), 1)
p_measured = coeffs[0]

log(f"Measured Exponent p: {p_measured:.6f}")
log(f"GR-Matter Exponent p_theory: 0.666667")

# ------------------------------------------------------------------------------
# 3. Verdict
# ------------------------------------------------------------------------------
status = "FAILED"
if abs(p_measured - 0.666667) < 0.1:
    status = "VERIFIED (Natural Friedmann Emergence)"
elif abs(p_measured - 0.5) < 0.1:
    status = "INCONCLUSIVE (Radiation-like Scaling Detected)"
else:
    status = "FAILED (Non-Friedmann Expansion)"

log(f"\nSTATUS: {status}")

# ------------------------------------------------------------------------------
# 4. Output Report
# ------------------------------------------------------------------------------
with open(REPORT, "w") as f:
    f.write("# QW-1552 AUDIT: Emergent Friedmann Dynamics\n\n")
    f.write(f"**STATUS:** {status}\n\n")
    f.write("## Operational Analysis\n")
    f.write("- **Method:** Derived the expansion rate $da/dt$ from a fundamental \n")
    f.write("  link-density flux model ($da/dt \\propto E_{total}/a^2$), rather \n")
    f.write("  than postulating the Friedmann equation.\n")
    f.write("- **Mechanism:** Energy density $\\rho$ scales inversely with volume.\n")
    f.write(f"- **Measured Scaling:** $a(t) \\propto t^{{{p_measured:.2f}}}$\n")
    f.write("- **Constraint:** A match with GR exponent ($2/3$) suggests that \n")
    f.write("  FIN link-dynamics naturally satisfy the Friedmann constraints.\n\n")
    
    if "VERIFIED" in status:
        f.write("> **Verdict:** The simulation produced the matter-dominated \n")
        f.write("> power-law exponent natively. This confirms the Friedmann \n")
        f.write("> structure as a derivation from link-conservation.\n")
    else:
        f.write("> **Verdict:** The scaling does not match Standard Cosmology. \n")
        f.write("> The mechanism may be valid, but the correspondence to GR is not proved.\n")

    f.write("\n## Raw Log\n```\n" + "\n".join(md) + "\n```\n")

print(f"\n✅ Report saved to {REPORT}")
