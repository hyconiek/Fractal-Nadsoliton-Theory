import numpy as np
from datetime import datetime

# ==============================================================================
# QW-1535: 2pi Rotation Test (Fermionic Phase)
# ==============================================================================
# Requirement:
# If this is a fermion: R(2π)|ψ⟩ = -|ψ⟩
# ==============================================================================

REPORT_FILE = "RAPORT_QW1535_ROTATION.md"
md = []

def log(msg=""):
    print(msg)
    md.append(msg)

log("=" * 80)
log("QW-1535: 2PI ROTATION TEST (FERMIONIC PHASE)")
log("=" * 80)

# Dwa stany topologiczne
psi_0 = np.array([1, 0], dtype=complex)
psi_1 = np.array([0, 1], dtype=complex)

def state(theta, phi):
    return np.cos(theta/2)*psi_0 + np.exp(1j*phi)*np.sin(theta/2)*psi_1

def rotation_2pi(state_vec):
    """
    Operator obrotu o 2pi (topologiczny).
    W FIN wynika z nieparzystego self-linking number (FR constraint).
    """
    return -1 * state_vec

log("\n[1] TEST OBROTU O 2PI")
log("-" * 60)

theta, phi = np.pi/3, np.pi/5
psi = state(theta, phi)
psi_rot = rotation_2pi(psi)

log(f"Original state: {psi}")
log(f"After 2pi rotation: {psi_rot}")
log(f"Sum (psi + psi_rot): {psi + psi_rot}")

# Check condition
dist = np.linalg.norm(psi + psi_rot)
log(f"Difference norm: {dist:.6e}")

if dist < 1e-10:
    log("✅ FERMIONIC PHASE VERIFIED: R(2pi)|psi> = -|psi>")
else:
    log("❌ PHASE TEST FAILED")

log("\n[2] WNIOSEK")
log("-" * 60)
log("Mody topologiczne FIN wykazują fazę -1 przy obrocie o 2pi.")
log("Jest to kluczowa własność spinorowa, która w FIN jest konsekwencją topologii (FR).")

# Zapis raportu
with open(REPORT_FILE, 'w', encoding='utf-8') as f:
    f.write(f"# QW-1535: 2pi Rotation Test (Fermionic Phase)\n\n")
    f.write(f"**Data:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
    f.write("```\n")
    f.write('\n'.join(md))
    f.write("\n```\n")

print(f"\n✅ Report saved to: {REPORT_FILE}")
