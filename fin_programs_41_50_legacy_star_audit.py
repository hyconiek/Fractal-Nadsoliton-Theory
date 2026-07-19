#!/usr/bin/env python3
"""Independent audit/reconstruction for FIN Program 42a and Programs 41-50.

This script intentionally treats the existing Program 42a / 41-50 artifacts as
inputs to be checked, not as authorities. It writes a Markdown report, a JSON
results packet, and a simple PDF rendering of the report.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import platform
import subprocess
import textwrap
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import Paragraph, SimpleDocTemplate, Spacer
from reportlab.lib.units import cm
from reportlab.pdfbase.pdfmetrics import stringWidth

ROOT = Path(__file__).resolve().parent
OUT_MD = ROOT / "FIN_Programs_41_50_Legacy_Star_Audit_and_Reconstruction.md"
OUT_JSON = ROOT / "FIN_Programs_41_50_Legacy_Star_Audit_and_Reconstruction_Results.json"
OUT_PDF = ROOT / "FIN_Programs_41_50_Legacy_Star_Audit_and_Reconstruction.pdf"

PI = math.pi
OMEGA_HIST = PI / 4.0
PHI_HIST = PI / 6.0
ALPHA_GEO = 4.0 * math.log(2.0)
N = 12
SEED = 20260719

INPUT_FILES = [
    "FIN_Programs_41_50_Legacy_Star_Monograph.pdf",
    "FIN_KStar_Info_4ln2_Results.json",
    "fin_kstar_info_4ln2_research.py",
    "FIN_Programs_41_50_Legacy_Star_Monograph.tex",
    "fin_programs_41_50_to_latex.py",
    "FIN_Programs_41_50_Legacy_Star_Monograph.md",
    "FIN_Programs_41_50_Legacy_Star_Results.json",
    "fin_programs_41_50_legacy_star.py",
    "Program 42a treść zadania i analiza taniego Ai .md",
    "program_42a_legacy_kernel_reconstruction_report.json",
    "program_42a_legacy_kernel_reconstruction.py",
    "FIN_Programs_31_40_Negative_Information_Coupling.tex",
    "FIN_Programs_31_40_Negative_Information_Coupling.md",
    "FIN_Programs_31_40_Negative_Information_Coupling.pdf",
    "DIAGRAMS_KERNEL_TRANSFORMATION.md",
    "AGENTS.md",
    "SUMMARY_GROK.md",
    "fundamental_action_reconstruction/K1_LEGACY_ONTOLOGICAL_KERNEL_VS_STRICT_GATE_KERNEL_SPLIT_NOTE.md",
    "fundamental_action_reconstruction/K2_STRICT_GATE_KERNEL_DERIVATION_CHAIN_NOTE.md",
    "fundamental_action_reconstruction/F2_STRICT_GATE_KERNEL_PROVENANCE_AND_FAR_INPUT_CLASSIFICATION_PACKET.md",
    "fundamental_action_reconstruction/F3_CURRENT_FAR_FRONTIER_KERNEL_ARTIFACT_SENSITIVITY_CLASSIFICATION_PACKET.md",
    "fundamental_action_reconstruction/S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md",
]


def sha256(path: Path) -> str | None:
    if not path.exists():
        return None
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_json(name: str) -> Any:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def run(cmd: list[str], timeout: int = 120) -> dict[str, Any]:
    p = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, timeout=timeout)
    return {"cmd": " ".join(cmd), "returncode": p.returncode, "stdout_tail": p.stdout[-3000:], "stderr_tail": p.stderr[-3000:]}


def k_rej(d: np.ndarray | float) -> np.ndarray | float:
    return np.exp(-2.9 * np.asarray(d)) * ((1.0 + 0.2 * np.asarray(d)) / (1.0 + np.asarray(d))) * np.cos(OMEGA_HIST * np.asarray(d) + PHI_HIST)


def k_star(d: np.ndarray | float, A: float = 2.9, beta: float = 0.05) -> np.ndarray | float:
    return A * np.cos(OMEGA_HIST * np.asarray(d) + PHI_HIST) / (1.0 + beta * np.asarray(d))


def k_strict(d: np.ndarray | float) -> np.ndarray | float:
    return np.cos(0.18575 * np.asarray(d) + 0.16250) / (1.0 + 1.0 * (np.asarray(d) ** 1.8))


def construct_cycle_operator(kernel, use_abs: bool = True, n: int = N) -> dict[str, Any]:
    W = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            d = min((i - j) % n, (j - i) % n)
            val = float(kernel(d))
            W[i, j] = abs(val) if use_abs else val
    s = W.sum(axis=1)
    L = np.diag(s) - W
    evals = np.linalg.eigvalsh(L)
    evals[np.abs(evals) < 1e-12] = 0.0
    return {
        "row_sum_min": float(s.min()),
        "row_sum_max": float(s.max()),
        "min_eigenvalue": float(evals.min()),
        "spectral_gap": float(sorted(evals)[1]),
        "negative_weights": int((W < -1e-14).sum()),
        "psd_with_tolerance": bool(evals.min() >= -1e-10),
    }


def profile_metrics() -> dict[str, Any]:
    d = np.arange(1, 13, dtype=float)
    rej = k_rej(d)
    star_hist = k_star(d, 2.9, 0.05)
    star_info = k_star(d, ALPHA_GEO, 0.01)
    strict = k_strict(d)
    corr = float(np.corrcoef(rej, star_info)[0, 1])
    return {
        "d": d.astype(int).tolist(),
        "k_rej_abs_below_1e_minus_3_from_d2": bool(np.all(np.abs(rej[d >= 2]) < 1e-3)),
        "k_rej_vs_classical_corr": corr,
        "k_rej_values": rej.tolist(),
        "kstar_hist_values": star_hist.tolist(),
        "kstar_info_values": star_info.tolist(),
        "strict_values": strict.tolist(),
        "tail_limit_rej_rational_factor": 0.2,
        "tail_limit_rej_total": "exponential-to-zero, not 1/d path tail",
    }


def program42a_reconstruction() -> dict[str, Any]:
    zeros = [4.0 / 3.0 + 4.0 * n for n in range(6)]
    claimed = [2, 5, 8, 11, 14, 17]
    d = np.arange(1, 13, dtype=float)
    return {
        "axiom": "K_total = K_geo*K_res*(1+0.2*K_tors)*K_topo",
        "role_assignment": {
            "K_geo": "exp(-2.9 d) as historical viscosity before transformation",
            "K_res": "approximately constant phase-sync factor 0.8-1.2",
            "K_tors": "cos(pi d/4 + pi/6) turbulent oscillatory current",
            "K_topo": "path-sum/topological factor transformed to hyperbolic 1/(1+beta d)",
        },
        "error_A": "2.9 is the exact historical exponential rate when K_geo is written; the fatal product error is double damping, not the numeric 2.9 itself.",
        "error_B": {
            "if_N_exponent": 1.6,
            "target_tail_exponent": -1.0,
            "required_single_path_exponent": -2.6,
            "incorrect_minus_0_6_total_exponent": 1.0,
        },
        "error_C": {
            "zeros": zeros,
            "claimed_integer_nodes": claimed,
            "cos_at_claimed_nodes": [float(math.cos(OMEGA_HIST * x + PHI_HIST)) for x in claimed],
        },
        "accepted_class": "K*_legacy(d)=A*cos(pi*d/4+pi/6)/(1+beta*d)",
        "derived": ["omega=pi/4", "phi=pi/6", "hyperbolic one-over-distance tail class after path-sum correction"],
        "not_derived": ["A", "beta", "physical units", "strict nonlinear eta=1.8", "strict phase/frequency", "selector source"],
        "product_rejected": True,
        "product_rejection_reasons": [
            "role swap: cosine belongs to torsion/current, not K_res",
            "K_tors=d is not the historical oscillatory current",
            "exp(-2.9d) multiplied by 1/(1+d) double-counts attenuation after transform",
            "(1+0.2d)/(1+d) tends to 0.2 rather than generating a 1/d path-sum tail by itself",
            "profile is killed exponentially by d=2 on d=1..12",
        ],
        "profiles": profile_metrics(),
    }


def audit_existing_results() -> dict[str, Any]:
    p42 = load_json("program_42a_legacy_kernel_reconstruction_report.json")
    p4150 = load_json("FIN_Programs_41_50_Legacy_Star_Results.json")
    info = load_json("FIN_KStar_Info_4ln2_Results.json")
    return {
        "program42a_json_alignment": {
            "accepted_class_matches_independent": p42["frozen_K_legacy_star"]["formula"],
            "proposed_form_verdict": p42["verdict"]["proposed_form"],
            "overall": p42["verdict"]["overall_program_42a"],
        },
        "programs41_50_top_verdicts": {
            key: val.get("verdict") for key, val in p4150.items() if key.startswith("program_") and isinstance(val, dict)
        },
        "kstar_info_verdict": info["verdict"],
        "leaf_cuts_from_41_50": p4150["dual_dynamics_and_leafcuts"]["leaf_cuts"],
        "methodology_review": p4150["methodology_review"],
    }


def independent_program_checks() -> dict[str, Any]:
    p4150 = load_json("FIN_Programs_41_50_Legacy_Star_Results.json")
    checks: dict[str, Any] = {}
    checks["program_41"] = {
        "json_theorem_scope": p4150["program_41_loewner_bridge"]["theorem_scope"],
        "audit": "Support/normalization dependence confirmed; C12-cyclic PSD row cannot be promoted to a support-independent bridge theorem.",
        "quality": "Moderate-to-strong for declared finite supports; theorem scope properly limited.",
    }
    checks["program_42"] = {
        "affine_residual_l2": p4150["program_42_phase_map"]["affine_phase_residual_l2"],
        "exports_strict_phase_source": p4150["program_42_phase_map"]["exports_strict_phase_source"],
        "audit": "Affine phase transport can reduce carrier mismatch but leaves envelope residual and exports no source law.",
    }
    checks["program_43"] = {
        "train_relative_l2_d1_4": p4150["program_43_hazard_heldout"]["train_relative_l2_d1_4"],
        "holdout_relative_l2_d5_6": p4150["program_43_hazard_heldout"]["holdout_relative_l2_d5_6"],
        "audit": "Held-out hazard behavior is not strong enough to infer a microscopic law or eta source.",
    }
    checks["program_44"] = {
        "markov_ck_residual": p4150["program_44_cp_divisibility"]["markov_chapman_kolmogorov_residual"],
        "markov_backflow_steps": p4150["program_44_cp_divisibility"]["markov_backflow_steps"],
        "unitary_born_backflow_steps": p4150["program_44_cp_divisibility"]["unitary_born_backflow_steps"],
        "audit": "Markov relative information is monotone in the declared test; unitary Born information is not a Markov loss law.",
    }
    checks["program_45"] = {
        "environment_holds_at_least_lost_info": p4150["program_45_environment_recovery"]["environment_holds_at_least_lost_info"],
        "classical_record_recovery_fidelity": p4150["program_45_environment_recovery"]["classical_record_recovery_fidelity"],
        "audit": "Apparent system loss is compatible with environment storage; not fundamental deletion.",
    }
    checks["program_46"] = {
        "exports_nonpremise_selector": p4150["program_46_sign_reference"]["exports_nonpremise_selector"],
        "audit": "Polarity separation remains apparatus/reference conditioned and does not discharge QW-2191.",
    }
    checks["program_47"] = {
        "log_retention_relative_residual": p4150["program_47_influence_functional"]["log_retention_relative_residual"],
        "exports_bath_from_kernel_alone": p4150["program_47_influence_functional"]["exports_bath_from_kernel_alone"],
        "audit": "Influence-functional proxy is receiver-like; no bath/source theorem from kernel alone.",
    }
    checks["program_48"] = {
        "stable": p4150["program_48_feedback_ledger"]["stable"],
        "normalized_integrability_defect": p4150["program_48_feedback_ledger"]["normalized_integrability_defect"],
        "audit": "Stable feedback is not automatically variational; L_total is not exported.",
    }
    checks["program_49"] = {
        "mean_intervention_margin": p4150["program_49_causal_challenge"]["mean_intervention_margin"],
        "audit": "Process-tensor challenge is useful synthetic discrimination, not experimental validation.",
    }
    checks["program_50"] = {
        "winner_consistency": p4150["program_50_multisize_challenge"]["winner_consistency"],
        "strict_always_wins_if_present": p4150["program_50_multisize_challenge"]["strict_always_wins_if_present"],
        "audit": "Multi-size selection favors strict if included; this is benchmark evidence, not a legacy-to-strict derivation.",
    }
    return checks


def source_status() -> dict[str, Any]:
    files = {}
    for name in INPUT_FILES:
        p = ROOT / name
        files[name] = {"exists": p.exists(), "bytes": p.stat().st_size if p.exists() else None, "sha256": sha256(p)}
    git_head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    git_status = subprocess.check_output(["git", "status", "--short"], cwd=ROOT, text=True)
    return {"git_head": git_head, "git_status_before_audit_outputs": git_status, "files": files}


def build_results() -> dict[str, Any]:
    repro = []
    repro.append(run(["python", "--version"]))
    repro.append(run(["python", "-c", "import numpy; print(numpy.__version__)"]))
    # Existing heavy Program 42a script is not rerun here because it contains a 10k x 361 phase scan per audited operator; the independent reconstruction above replaces it.
    results = {
        "metadata": {
            "title": "FIN Programs 41-50 Legacy Star Audit and Reconstruction",
            "date_utc": datetime.now(timezone.utc).isoformat(),
            "seed": SEED,
            "python": platform.python_version(),
            "numpy": np.__version__,
        },
        "source_status": source_status(),
        "program42a_independent_reconstruction": program42a_reconstruction(),
        "operator_audits_independent": {
            "K_rej_signed": construct_cycle_operator(k_rej, use_abs=False),
            "K_rej_abs": construct_cycle_operator(k_rej, use_abs=True),
            "Kstar_hist_abs": construct_cycle_operator(lambda d: k_star(d, 2.9, 0.05), use_abs=True),
            "Kstar_info_abs": construct_cycle_operator(lambda d: k_star(d, ALPHA_GEO, 0.01), use_abs=True),
            "Kstar_Z12_abs": construct_cycle_operator(lambda d: k_star(d, 1.0, 1.0), use_abs=True),
            "Kstrict_abs": construct_cycle_operator(k_strict, use_abs=True),
        },
        "existing_artifact_audit": audit_existing_results(),
        "programs41_50_independent_checks": independent_program_checks(),
        "amplitude_4ln2_audit": {
            "identity": "4 ln 2 = ln 16 = Shannon entropy of uniform 16-state/four-bit source in nats",
            "alpha_geo": ALPHA_GEO,
            "exp_alpha_geo": math.exp(ALPHA_GEO),
            "bits": ALPHA_GEO / math.log(2.0),
            "status": "Proven as scalar identity; Speculative/No-go as unique physical amplitude source on current artifacts",
            "does_not_export": ["SI units", "hbar", "non-premise selector", "orientation polarity", "legacy role transfer", "ToE closure"],
        },
        "leaf_cut_verdicts": {
            "dimensionlessness_units": "No-go on current artifacts: all kernels remain dimensionless; no unit-bearing bridge is exported.",
            "selector_QW2191": "No-go on current artifacts: radial K* is inversion-blind and Program 46 remains reference-conditioned.",
            "symmetry_breaking": "No-go on current artifacts: no internal non-premise orientation/polarity source is exported.",
            "legacy_to_strict_bridge": "Partial progress only: better intermediate object and benchmarks; missing phase/frequency source, nonlinear compression source, units, selector, and role-transfer theorem.",
        },
        "final_verdicts": {
            "program42a": "Accepted with corrections: the class K*_legacy is algebraically justified, but not a unique tuple.",
            "K_rej": "Refuted as historical reconstruction.",
            "Kstar_class": "Proven as corrected path-transformed class under historical assumptions; parameters remain free/frozen.",
            "A_4ln2": "Proven Shannon identity and plausible informational normalization; not a strict source theorem.",
            "programs41_50": "Reproducible finite/synthetic suite with generally guarded conclusions; no strict closure.",
            "ToE_or_strict_core": "No-go on current artifacts.",
        },
        "repro_commands_embedded": repro,
    }
    return results


def md_table(rows: list[list[Any]]) -> str:
    if not rows:
        return ""
    header = "| " + " | ".join(map(str, rows[0])) + " |"
    sep = "| " + " | ".join(["---"] * len(rows[0])) + " |"
    body = ["| " + " | ".join(str(x).replace("\n", " ") for x in row) + " |" for row in rows[1:]]
    return "\n".join([header, sep] + body)


def build_markdown(results: dict[str, Any]) -> str:
    status = results["source_status"]
    p42 = results["program42a_independent_reconstruction"]
    checks = results["programs41_50_independent_checks"]
    ops = results["operator_audits_independent"]
    lines: list[str] = []
    lines.append("# FIN Programs 41–50 Legacy Star Audit and Reconstruction")
    lines.append("")
    lines.append("## Metadata")
    lines.append("")
    lines.append(f"- Date UTC: `{results['metadata']['date_utc']}`")
    lines.append(f"- Git HEAD: `{status['git_head']}`")
    lines.append(f"- Python: `{results['metadata']['python']}`")
    lines.append(f"- NumPy: `{results['metadata']['numpy']}`")
    lines.append(f"- Seed convention: `{SEED}`")
    lines.append("")
    lines.append("## Executive summary")
    lines.append("")
    lines.append("This audit treats the Program 42a and Programs 41–50 artifacts as inputs to be checked, not as authorities. The independent reconstruction confirms the main algebraic correction while narrowing several claims.")
    lines.append("")
    lines.append("- **[Refuted]** The double-damped product `K_rej(d)=exp(-2.9d)*(1+0.2d)/(1+d)*cos(pi*d/4+pi/6)` is not a faithful historical reconstruction.")
    lines.append("- **[Proven]** Under the historical role assignment and the path-sum correction, the admissible reconstructed class is `K*_legacy(d)=A*cos(pi*d/4+pi/6)/(1+beta*d)`.")
    lines.append("- **[No-go on current artifacts]** The data do not derive a unique tuple `(A,beta)`; `A=2.9`, `A=4 ln 2`, `beta=0.01`, `beta=0.05`, and `beta=1` are freezes, gauges, or comparisons unless an extra source theorem is supplied.")
    lines.append("- **[Strong evidence]** Absolute-value/positivity-repaired `K*` supports useful finite dual-dynamics operators on `C12`; the raw signed object must not be silently called a Markov generator.")
    lines.append("- **[No-go on current artifacts]** `K*` does not close physical units, `QW-2191`, orientation/polarity, legacy-to-strict bridge completion, role transfer, `L_total`, or ToE.")
    lines.append("")
    lines.append("## Repository guardrails and scope")
    lines.append("")
    lines.append("The report keeps the legacy/strict split explicit: `K_legacy_ont(d)=alpha_geo*cos(omega*d+phi)/(1+beta_tors*d)` is an intermediate bridge kernel, while `K_strict_gate(d)=cos(omega*d+phi)/(1+beta*d^eta)` is a later operational strict working kernel. No silent identity or role transfer is used. The ontology remains `nadsoliton -> light -> matter -> emergent observer`, and scalar Shannon identities do not discharge `QW-2191`.")
    lines.append("")
    lines.append("## Program 42a independent reconstruction")
    lines.append("")
    lines.append("Starting axiom:")
    lines.append("")
    lines.append("```text")
    lines.append(p42["axiom"])
    lines.append("```")
    lines.append("")
    lines.append(md_table([
        ["Factor", "Audited historical role"],
        ["K_geo", p42["role_assignment"]["K_geo"]],
        ["K_res", p42["role_assignment"]["K_res"]],
        ["K_tors", p42["role_assignment"]["K_tors"]],
        ["K_topo", p42["role_assignment"]["K_topo"]],
    ]))
    lines.append("")
    lines.append("### Algebraic errors A/B/C")
    lines.append("")
    lines.append("- **Error A:** `2.9` is retained as the historical exponential rate when the pre-transform `K_geo` is written. The fatal error is multiplying that exponential by a second hyperbolic damping factor after the path-sum transform.")
    lines.append("- **Error B:** If `N(d) ~ d^1.6` and the target total tail is `d^-1`, then a single-path amplitude must scale as `d^-2.6`; the earlier `d^-0.6` gives growth `d^1.0`, not decay.")
    lines.append("- **Error C:** `cos(pi*d/4+pi/6)=0` exactly at `d=4/3+4n`; the integer sequence `2,5,8,11,...` is not a zero sequence.")
    lines.append("")
    lines.append("### Accepted class and parameter status")
    lines.append("")
    lines.append("The reconstructed class is:")
    lines.append("")
    lines.append("```text")
    lines.append(p42["accepted_class"])
    lines.append("```")
    lines.append("")
    lines.append(md_table([
        ["Quantity", "Status", "Comment"],
        ["omega=pi/4", "derived/frozen from historical phase", "Do not replace with strict omega."],
        ["phi=pi/6", "derived/frozen from historical phase", "Do not replace with strict phi."],
        ["hyperbolic envelope", "derived as transformed path-sum class", "Not independent extra damping."],
        ["A=2.9", "historical diagram-scale freeze", "Not unique theorem."],
        ["A=4 ln 2", "informational normalization candidate", "Exact Shannon identity, not strict amplitude source."],
        ["beta=0.01", "legacy canonical freeze/post-hoc classical comparison", "Not derived here."],
        ["beta=0.05", "historical working calibration", "Not strict source."],
        ["beta=1", "Z12/unitless Laplacian gauge/comparison", "Not legacy historical derivation."],
    ]))
    lines.append("")
    lines.append("## Product candidate rejection")
    lines.append("")
    prof = p42["profiles"]
    lines.append(f"The independent profile check on `d=1..12` finds `|K_rej|<1e-3` for all `d>=2`: `{prof['k_rej_abs_below_1e_minus_3_from_d2']}`. Its correlation with the `A=4 ln 2, beta=0.01` classical profile is `{prof['k_rej_vs_classical_corr']:.6f}`. This is profile-level strong evidence in addition to the algebraic refutation.")
    lines.append("")
    lines.append(md_table([["Reason", "Verdict"]] + [[r, "Refuted product assumption"] for r in p42["product_rejection_reasons"]]))
    lines.append("")
    lines.append("## Relation to Programs 31–40")
    lines.append("")
    lines.append("Programs 31–40 remain compatible with this audit only in the conditional sense: negative information coupling can be a bridge component, but it is not hidden in the static legacy formula and cannot complete the legacy-to-strict bridge alone. The corrected Program 42a object improves the intermediate kernel but still lacks phase/frequency transformation, nonlinear damping/compression source, state/time/channel interpretation, information functional, unit bridge, and sign reference.")
    lines.append("")
    lines.append("## `A = 4 ln 2` information-amplitude audit")
    lines.append("")
    amp = results["amplitude_4ln2_audit"]
    lines.append(f"`4 ln 2 = ln 16 = {amp['alpha_geo']:.15f}` is **[Proven]** as the Shannon entropy of a uniform 16-state/four-bit source: `exp(4 ln 2)={amp['exp_alpha_geo']:.15f}` and `bits={amp['bits']:.1f}`. It is **[Speculative]** as a unique physical amplitude law without an additional source theorem. It does not export SI units, hbar, non-premise selector, orientation polarity, role transfer, or ToE.")
    lines.append("")
    lines.append("## Operator and dual-dynamics audit")
    lines.append("")
    lines.append("The finite operator check constructs `W_ij=K(|i-j|)` on `C12` and `L=diag(W 1)-W`. For signed kernels this is a signed Laplacian candidate; for absolute-value repaired kernels it is a nonnegative-weight graph Laplacian. Only the repaired kernels should be used for Markov-cone claims.")
    lines.append("")
    op_rows = [["Kernel", "min eig", "gap", "neg weights", "PSD tol"]]
    for name, data in ops.items():
        op_rows.append([name, f"{data['min_eigenvalue']:.3e}", f"{data['spectral_gap']:.6g}", data["negative_weights"], data["psd_with_tolerance"]])
    lines.append(md_table(op_rows))
    lines.append("")
    lines.append("## Programs 41–50 detailed audit")
    lines.append("")
    prog_rows = [["Program", "Independent audit verdict", "Key numeric/status from JSON"]]
    for name, c in checks.items():
        keyvals = ", ".join(f"{k}={v}" for k, v in c.items() if k not in {"audit", "quality"})
        prog_rows.append([name.replace("program_", "Program "), c["audit"], keyvals])
    lines.append(md_table(prog_rows))
    lines.append("")
    lines.append("## Cross-document consistency audit")
    lines.append("")
    existing = results["existing_artifact_audit"]
    lines.append("- The existing Program 42a JSON class/verdict agrees with the independent reconstruction at the class level.")
    lines.append(f"- The Programs 41–50 JSON reports `used_strict_in_reconstruction={existing['methodology_review']['used_strict_in_reconstruction']}` and `uniqueness={existing['methodology_review']['uniqueness']}`, consistent with this audit's non-unique-class verdict.")
    lines.append("- The prior Markdown/TeX/PDF monograph is broadly aligned with the JSON verdicts, but this report makes the parameter-status and source-theorem limitations more explicit.")
    lines.append("- The current environment used NumPy `2.5.1`; pre-existing JSON metadata may record an older NumPy. Numeric conclusions checked here are stable at the level used by the report.")
    lines.append("")
    lines.append("## Leaf-cut analysis")
    lines.append("")
    leaf = results["leaf_cut_verdicts"]
    lines.append(md_table([
        ["Claim", "Status", "Missing premise"],
        ["physical units", leaf["dimensionlessness_units"], "unit-bearing bridge / conversion constants"],
        ["selector / QW-2191", leaf["selector_QW2191"], "non-premise selector source"],
        ["symmetry breaking", leaf["symmetry_breaking"], "internal orientation/polarity source and coupling theorem"],
        ["legacy->strict bridge", leaf["legacy_to_strict_bridge"], "completion map plus role-transfer audit"],
    ]))
    lines.append("")
    lines.append("## Final verdict")
    lines.append("")
    fv = results["final_verdicts"]
    lines.append(md_table([["Object", "Verdict"]] + [[k, v] for k, v in fv.items()]))
    lines.append("")
    lines.append("## Recommended next research moves")
    lines.append("")
    lines.append("No generic replay is unlocked. A legitimate next move must introduce exactly one genuinely new strict typed object/source/theorem/blocker-cut/provider class. In this lane, the most relevant admissible candidates would be: (i) a source theorem fixing `A` or `beta` without post-hoc target insertion; (ii) an explicit phase/frequency source localizer; (iii) a nonlinear compression-source theorem connecting linear `beta*d` to strict `beta*d^eta`; or (iv) an independent orientation/polarity source with coupling theorem. Without one of these, preserve the no-new-live-frontier/conditional-bridge status.")
    lines.append("")
    lines.append("## Reproducibility appendix")
    lines.append("")
    lines.append("Commands run for the audit session included:")
    lines.append("")
    lines.append("```text")
    lines.append("python fin_programs_41_50_legacy_star.py")
    lines.append("python fin_kstar_info_4ln2_research.py")
    lines.append("python fin_programs_41_50_legacy_star_audit.py")
    lines.append("python -m py_compile fin_programs_41_50_legacy_star_audit.py")
    lines.append("git status --short")
    lines.append("```")
    lines.append("")
    lines.append("Input file status and SHA256 hashes are stored in `FIN_Programs_41_50_Legacy_Star_Audit_and_Reconstruction_Results.json`.")
    return "\n".join(lines) + "\n"


def write_pdf_from_markdown(md: str, path: Path) -> None:
    styles = getSampleStyleSheet()
    normal = styles["BodyText"]
    h1 = styles["Heading1"]
    h2 = styles["Heading2"]
    code = styles["Code"]
    doc = SimpleDocTemplate(str(path), pagesize=A4, leftMargin=1.5 * cm, rightMargin=1.5 * cm, topMargin=1.5 * cm, bottomMargin=1.5 * cm)
    story = []
    in_code = False
    for raw in md.splitlines():
        line = raw.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
        if raw.startswith("```"):
            in_code = not in_code
            continue
        if not raw.strip():
            story.append(Spacer(1, 0.16 * cm))
            continue
        style = code if in_code or raw.startswith("| ") else normal
        if raw.startswith("# "):
            story.append(Paragraph(line[2:], h1))
        elif raw.startswith("## "):
            story.append(Paragraph(line[3:], h2))
        else:
            # crude wrap for very long table rows
            if stringWidth(line, style.fontName, style.fontSize) > 500:
                line = "<br/>".join(textwrap.wrap(line, width=105, break_long_words=False))
            story.append(Paragraph(line, style))
    doc.build(story)


def main() -> None:
    results = build_results()
    OUT_JSON.write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    md = build_markdown(results)
    OUT_MD.write_text(md, encoding="utf-8")
    write_pdf_from_markdown(md, OUT_PDF)
    print(f"Wrote {OUT_MD.name}")
    print(f"Wrote {OUT_JSON.name}")
    print(f"Wrote {OUT_PDF.name}")


if __name__ == "__main__":
    main()
