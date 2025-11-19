"""
103_QUICK_WIN_RG.py
Ten quick-win tasks for a coarse RG / running coupling probe based on the kernel
K(d) = alpha_geo * cos(omega*d + phi) / (1 + beta_tors*d)
No fitting, no tautology — only computed proxies and pragmatic diagnostics.
Outputs: report_103_quick_win.md and report_103_quick_win.json
"""
# Author: Krzysztof Żuchowski


import json
import numpy as np
import scipy.linalg as la
import datetime

# Core model params (kept consistent with previous studies)
alpha_geo = 2.77
beta_tors = 0.01
omega = 2 * np.pi / 8.0
phi = 0.0
m0 = 0.44  # MeV, kept as reference

np.random.seed(42)

OUT_MD = "report_103_quick_win.md"
OUT_JSON = "report_103_quick_win.json"

log = []

def writelog(s):
    print(s)
    log.append(s)

# Helper: build distance matrix for N sites placed on unit circle
def build_kernel_matrix(N, alpha, beta, omega, phi):
    thetas = np.linspace(0, 2 * np.pi, N, endpoint=False)
    xs = np.column_stack([np.cos(thetas), np.sin(thetas)])
    D = la.norm(xs[:, None, :] - xs[None, :, :], axis=2)
    K = alpha * np.cos(omega * D + phi) / (1.0 + beta * D)
    return K, D

# Task 0: baseline eigenstructure
N = 16
K0, D0 = build_kernel_matrix(N, alpha_geo, beta_tors, omega, phi)
w0, v0 = la.eigh(K0)
idx = np.argsort(np.abs(w0))[::-1]
w0 = w0[idx]
v0 = v0[:, idx]

w_top = w0[:8]

writelog(f"TASK 0: N={N}, alpha={alpha_geo}, beta={beta_tors}, omega={omega:.4f}")
writelog(f"  top eigenvalues (abs-sorted) = {np.round(w_top,6).tolist()}")

# Task 1: coarse RG scaling of alpha_geo -> g(s) = alpha(s)
# We'll implement a simple scale parameter s where distances scale as d -> s * d
s_vals = np.logspace(-1, 1, 41)  # from 0.1 to 10 (log scale)
alpha_s = []
for s in s_vals:
    K_s, _ = build_kernel_matrix(N, alpha_geo, beta_tors, omega, phi)
    # emulate scaling by dilating distances: replace omega*d by omega*(s*d)
    thetas = np.linspace(0, 2 * np.pi, N, endpoint=False)
    xs = np.column_stack([np.cos(thetas), np.sin(thetas)])
    D = la.norm(xs[:, None, :] - xs[None, :, :], axis=2)
    K_scaled = alpha_geo * np.cos(omega * (s * D) + phi) / (1.0 + beta_tors * (s * D))
    # define "effective coupling" g(s) as mean absolute kernel strength
    g_s = np.mean(np.abs(K_scaled))
    alpha_s.append(g_s)

alpha_s = np.array(alpha_s)

# Task 2: estimate beta_function proxy: beta(g) = d g / d ln s (numerical)
lns = np.log(s_vals)
dg_dlns = np.gradient(alpha_s, lns)

# Task 3: find candidate fixed points where beta ~ 0 (sign changes / small magnitude)
fixed_idx = np.where(np.isclose(dg_dlns, 0, atol=1e-4))[0].tolist()
# also find sign-changes
sign_change_idx = []
for i in range(len(dg_dlns) - 1):
    if dg_dlns[i] * dg_dlns[i + 1] < 0:
        sign_change_idx.append(i)

writelog("TASK 1-3: computed g(s) and numerical beta(g)")
writelog(f"  s_vals (sample) = [{s_vals[0]:.3f}, ..., {s_vals[-1]:.3f}]")
writelog(f"  g(s) [min,max] = [{alpha_s.min():.6f}, {alpha_s.max():.6f}]")
writelog(f"  beta proxy (dg/dlns) [min,max] = [{dg_dlns.min():.6e}, {dg_dlns.max():.6e}]")
writelog(f"  near-zero beta indices = {fixed_idx}")
writelog(f"  sign-change indices = {sign_change_idx}")

# Task 4: anomalous dimension proxy from eigenvalue scaling
# For each s, compute top eigenvalue magnitude and check scaling exponent: w_top(s) ~ s^{-gamma}
w_top_s = []
for s in s_vals:
    thetas = np.linspace(0, 2 * np.pi, N, endpoint=False)
    xs = np.column_stack([np.cos(thetas), np.sin(thetas)])
    D = la.norm(xs[:, None, :] - xs[None, :, :], axis=2)
    K_scaled = alpha_geo * np.cos(omega * (s * D) + phi) / (1.0 + beta_tors * (s * D))
    w, v = la.eigh(K_scaled)
    w_top_s.append(np.max(np.abs(w)))

w_top_s = np.array(w_top_s)
# estimate local scaling exponent gamma(s) = - d ln w / d ln s
dlnw_dlns = np.gradient(np.log(np.abs(w_top_s) + 1e-12), lns)
gamma_s = -dlnw_dlns

writelog("TASK 4: anomalous-dimension proxy gamma(s) estimated from top eigenvalue scaling")
writelog(f"  gamma(s) [min,max] = [{gamma_s.min():.6f}, {gamma_s.max():.6f}]")

# Task 5: RG flow of vacuum-energy proxy: take trace(K_scaled)/N as simple zero-point proxy
E_s = []
for s in s_vals:
    thetas = np.linspace(0, 2 * np.pi, N, endpoint=False)
    xs = np.column_stack([np.cos(thetas), np.sin(thetas)])
    D = la.norm(xs[:, None, :] - xs[None, :, :], axis=2)
    K_scaled = alpha_geo * np.cos(omega * (s * D) + phi) / (1.0 + beta_tors * (s * D))
    E_s.append(np.trace(K_scaled) / float(N))
E_s = np.array(E_s)

writelog("TASK 5: vacuum-energy proxy E(s) computed (trace/N)")
writelog(f"  E(s) [min,max] = [{E_s.min():.6e}, {E_s.max():.6e}]")

# Task 6: running of mass-hierarchy proxy: pick two low-lying modes and track their ratio
# We'll use the 1st and 4th largest eigenvalues as proxies for 'light' and 'heavy' modes
ratio_s = []
for s in s_vals:
    thetas = np.linspace(0, 2 * np.pi, N, endpoint=False)
    xs = np.column_stack([np.cos(thetas), np.sin(thetas)])
    D = la.norm(xs[:, None, :] - xs[None, :, :], axis=2)
    K_scaled = alpha_geo * np.cos(omega * (s * D) + phi) / (1.0 + beta_tors * (s * D))
    w, v = la.eigh(K_scaled)
    w_sorted = np.sort(np.abs(w))[::-1]
    # safe guards
    if len(w_sorted) >= 4 and w_sorted[3] != 0:
        r = w_sorted[0] / (w_sorted[3] + 1e-12)
    else:
        r = np.nan
    ratio_s.append(r)
ratio_s = np.array(ratio_s)

writelog("TASK 6: mass-hierarchy proxy (top/4th) across scales computed")
writelog(f"  ratio(s) [min,max] = [{np.nanmin(ratio_s):.6e}, {np.nanmax(ratio_s):.6e}]")

# Task 7: operator-mixing proxy: compute overlap matrix between top-3 eigenvectors at s=1 and s=2
sA = 1.0
sB = 2.0
# build at sA
thetas = np.linspace(0, 2 * np.pi, N, endpoint=False)
xs = np.column_stack([np.cos(thetas), np.sin(thetas)])
DA = la.norm(xs[:, None, :] - xs[None, :, :], axis=2)
KA = alpha_geo * np.cos(omega * (sA * DA) + phi) / (1.0 + beta_tors * (sA * DA))
wbA, vbA = la.eigh(KA)
idxA = np.argsort(np.abs(wbA))[::-1]
vbA = vbA[:, idxA]

DB = la.norm(xs[:, None, :] - xs[None, :, :], axis=2)
KB = alpha_geo * np.cos(omega * (sB * DB) + phi) / (1.0 + beta_tors * (sB * DB))
wbB, vbB = la.eigh(KB)
idxB = np.argsort(np.abs(wbB))[::-1]
vbB = vbB[:, idxB]

# overlap matrix of top-3
topk = 3
Ov = np.abs(np.dot(vbA[:, :topk].T.conj(), vbB[:, :topk]))

writelog("TASK 7: operator-mixing proxy (overlap between top-3 eigenvectors at s=1 and s=2)")
writelog(f"  overlap matrix =\n{np.round(Ov,4)}")

# Task 8: anomaly proxy: look for non-conservation of trace under scale change beyond numerical rescaling
# compute difference in trace normalized by average
trA = np.trace(KA)
trB = np.trace(KB)
anom_proxy = (trB - trA) / (0.5 * (trB + trA) + 1e-12)

writelog("TASK 8: anomaly proxy (relative trace change s=1->2)")
writelog(f"  trA={trA:.6e}, trB={trB:.6e}, anomaly_proxy={anom_proxy:.6e}")

# Task 9: finite renormalization estimate: integrate beta proxy to get finite change between s=0.5 and s=2
# integrate dg/dlns numerically
from scipy.integrate import simpson
ln_min = np.log(0.5)
ln_max = np.log(2.0)
mask = (lns >= ln_min) & (lns <= ln_max)
if mask.sum() >= 3:
    # integrate dg/dlns over ln(s)
    delta_g = simpson(dg_dlns[mask], lns[mask])
else:
    delta_g = np.nan

writelog("TASK 9: finite renormalization estimate (integrate beta between s=0.5 and s=2)")
writelog(f"  delta_g (approx) = {delta_g:.6e}")

# Task 10: coarse phenomenological implication: if gamma(s) positive at some scales, interpret as screening vs antiscreening
screening_scales = list(s_vals[np.where(gamma_s > 0)[0]])
antiscreening_scales = list(s_vals[np.where(gamma_s < 0)[0]])

writelog("TASK 10: phenomenological note on screening vs antiscreening from gamma(s)")
writelog(f"  screening_scales (sample up to 5) = {np.round(np.array(screening_scales)[:5],4).tolist()}")
writelog(f"  antiscreening_scales (sample up to 5) = {np.round(np.array(antiscreening_scales)[:5],4).tolist()}")

# Write markdown report
now = datetime.datetime.utcnow().isoformat() + "Z"
md_lines = []
md_lines.append(f"# Report 103 — Quick-win RG probes\n")
md_lines.append(f"Generated: {now}\n")
md_lines.append("## Summary\n")
md_lines.append("This script computes simple RG proxies by scaling the distance kernel and measuring effective coupling, beta proxy, anomalous-dimension proxy, vacuum-energy proxy and related diagnostics. No fitting performed.\n")
md_lines.append("## Parameters\n")
md_lines.append(f"alpha_geo = {alpha_geo}\n")
md_lines.append(f"beta_tors = {beta_tors}\n")
md_lines.append(f"omega = {omega}\n")
md_lines.append("\n## Log\n")
for L in log:
    md_lines.append(f"- {L}\n")

with open(OUT_MD, "w") as f:
    f.writelines([l + "\n" for l in md_lines])

# small JSON summary for indexing
summary = {
    "study": 103,
    "generated": now,
    "params": {"alpha_geo": alpha_geo, "beta_tors": beta_tors, "omega": omega},
    "g_s_min": float(alpha_s.min()),
    "g_s_max": float(alpha_s.max()),
    "beta_min": float(dg_dlns.min()),
    "beta_max": float(dg_dlns.max()),
    "gamma_min": float(gamma_s.min()),
    "gamma_max": float(gamma_s.max()),
    "fixed_candidates": fixed_idx,
    "sign_change_indices": sign_change_idx,
    "anomaly_proxy": float(anom_proxy),
}

with open(OUT_JSON, "w") as f:
    json.dump(summary, f, indent=2)

writelog(f"Wrote {OUT_MD} and {OUT_JSON}")

# Print final one-line summary
writelog("DONE: Study 103 completed — report and json written.")
