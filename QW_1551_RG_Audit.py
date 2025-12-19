import numpy as np
from datetime import datetime

# ==============================================================================
# QW-1551 AUDIT (ROUND 3): RG Flow Scaling Consistency
# ==============================================================================
# OPERATIONAL INSTRUCTION:
# 1. Define TWO different mass measures: sigma_width and E_grad.
# 2. Check if both flow to the same fixed point.
# 3. Status: VERIFIED / INCONCLUSIVE / FAILED.

REPORT = "RAPORT_QW1551_RG_AUDIT.md"
md = []

def log(msg=""):
    print(msg)
    md.append(msg)

log("="*80)
log("QW-1551 OPERATIONAL AUDIT: RG FLOW CONSISTENCY")
log("="*80)

# ------------------------------------------------------------------------------
# 1. Coarse-Graining Simulation
# ------------------------------------------------------------------------------
def simulate_flow(steps=5):
    # Initial "particle" as a localized field in 1D
    N = 256
    phi = np.exp(-np.linspace(-10, 10, N)**2 / 2.0)
    
    results = []
    
    for s in range(steps):
        # Measure A: Sigma Width (1/std)
        # We define mass as the inverse of the effective volume
        std = np.sqrt(np.sum(phi * np.linspace(-10, 10, N)**2) / np.sum(phi))
        m_width = 1.0 / (std + 1e-9)
        
        # Measure B: Gradient Energy (Integrated (grad phi)^2)
        grad_phi = np.gradient(phi)
        m_grad = np.sum(grad_phi**2)
        
        results.append((m_width, m_grad))
        
        # Coarse-graining: Simple smoothing (Heat equation step)
        phi = phi + 0.1 * (np.roll(phi, 1) - 2*phi + np.roll(phi, -1))
        
    return results

log(f"{'Step':<5} | {'m_width':<15} | {'m_grad':<15}")
log("-" * 40)

flow_data = simulate_flow(8)
for i, (mw, mg) in enumerate(flow_data):
    log(f"{i:<5d} | {mw:<15.6f} | {mg:<15.6f}")

# ------------------------------------------------------------------------------
# 2. Consistency Analysis
# ------------------------------------------------------------------------------
# Normalize to initial values to check scaling correlation
mw_norm = [x[0]/flow_data[0][0] for x in flow_data]
mg_norm = [x[1]/flow_data[0][1] for x in flow_data]

correlation = np.corrcoef(mw_norm, mg_norm)[0, 1]
log(f"\nScaling Correlation (m_width vs m_grad): {correlation:.6f}")

# ------------------------------------------------------------------------------
# 3. Verdict
# ------------------------------------------------------------------------------
# If they scale perfectly together, the mechanism is robust.
status = "FAILED"
if correlation > 0.99:
    status = "VERIFIED (Mechanism Robust)"
elif correlation > 0.9:
    status = "INCONCLUSIVE (Partial Scaling Correlation)"

log(f"\nSTATUS: {status}")

# ------------------------------------------------------------------------------
# 4. Output Report
# ------------------------------------------------------------------------------
with open(REPORT, "w") as f:
    f.write("# QW-1551 AUDIT: RG Flow Scaling Consistency\n\n")
    f.write(f"**STATUS:** {status}\n\n")
    f.write("## Operational Analysis\n")
    f.write("- **Method:** Monitored the flow of two independent mass measures \n")
    f.write("  ($m_{width}$ and $m_{grad}$) during numerical coarse-graining.\n")
    f.write("- **Constraint:** A true physical flow must be independent of the \n")
    f.write("  specific definition of the parameter (Measurement Invariance).\n")
    f.write(f"- **Correlation Result:** {correlation:.6f}\n\n")
    
    if "VERIFIED" in status:
        f.write("> **Verdict:** Both measures exhibit synchronized scaling, confirming \n")
        f.write("> that the RG flow identifies a coherent physical property of the soliton.\n")
    else:
        f.write("> **Verdict:** Divergent scaling detected. The 'mass flow' depends on \n")
        f.write("> the definition, indicating an ad-hoc construct, not a physical law.\n")

    f.write("\n## Raw Log\n```\n" + "\n".join(md) + "\n```\n")

print(f"\n✅ Report saved to {REPORT}")
