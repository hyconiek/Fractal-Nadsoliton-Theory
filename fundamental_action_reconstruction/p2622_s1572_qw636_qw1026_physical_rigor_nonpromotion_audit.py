#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Any

import numpy as np

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2622_s1572_qw636_qw1026_physical_rigor_nonpromotion_audit.json"
MD = GEN / "p2622_s1572_qw636_qw1026_physical_rigor_nonpromotion_audit.md"

SOURCE_FILES = {
    "P2366_PHASE_REFERENCE_CHI11_AUDIT": GEN / "p2366_s1316_selector_candidate_phase_reference_chi11_audit_probe.json",
    "P2619_SELECTOR_LATTICE": GEN / "p2619_s1569_p2618_selector_source_obligation_lattice.json",
    "P2620_TWO_OBSTRUCTION_CUT": GEN / "p2620_s1570_p2618_p2619_bridge_two_obstruction_cut.json",
    "P2621_CHIRAL_HOPFION_CONDITIONAL_AUDIT": GEN / "p2621_s1571_qw636_qw1026_chiral_hopfion_selector_source_audit.json",
}

PRIOR_ART_FILES = {
    "QW636_FINAL_PARITY_CHECK": REPO / "TOE_Proven_Hypotheses_v3.1_FINAL" / "QW-636_Parity_Check.py",
    "QW637_GAUGE_CHECK": REPO / "TOE_Proven_Hypotheses_v3.1_FINAL" / "QW-637_Gauge_Check.py",
    "QW1026_SYNTHESIS_SUITE": REPO / "QW-1017_to_QW-1036_Synthesis_Suite.py",
    "RELEASE8_SIGN_GAUGE_GUARD": REPO / "TOE_FINAL_DOCUMENTATION_RELEASE_8_STRICT_FULL.tex",
}

NEGATIVE_EXPORT_FLAGS = [
    "qw636_promoted_to_strict_selector_source",
    "qw1026_promoted_to_strict_selector_source",
    "p2621_conditional_schema_promoted_to_unconditional_source",
    "orientation_odd_selector_source_exported_from_prior_art_alone",
    "p2620_bridge_source_cut_repaired",
    "nonlinear_damping_completion_source_exported",
    "full_bridge_revalidated",
    "role_transfer_revalidated",
    "role_bearing_ltotal_reenabled",
    "qw2191_discharged_by_this_packet",
    "toe_closure_claimed",
]

ALPHA_GEO = 4.0 * math.log(2.0)
OMEGA = math.pi / 4.0
PHI = math.pi / 6.0
BETA_TORS = 0.01
N_OCTAVES = 12


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else "MISSING"


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, ".",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex",
            "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**",
        ],
        cwd=REPO,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:160]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2622|S1572|physical rigor nonpromotion|QW-636.*QW-1026.*nonpromotion|selector prior-art nonpromotion",
        "qw636_qw637_prior_art": "QW-636|QW-637|Hopfion Phase|Complex Phase.*break Parity|GAUGE SYMMETRY|fixed links|Dynamic Geometry|Peierls|Flux",
        "qw1026_prior_art": "QW-1026|Chiral Anomaly|gamma5|γ₅|Tr\\(γ₅|K³|ANOMALY_FREE|alternating signs",
        "selector_guard_prior_art": "P2366|P2619|P2620|P2621|orientation_odd_selector_source|phase_origin_plus_chiral_bispectrum|sign-gauge|directed.*gauge|QW-2191",
        "methodology_risk_terms": "heuristic|candidate-not-theorem|basis-dependent|gauge/convention|not a strict-core source|no directed physical orientation datum|source theorem",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def q636_energy(sigma: int, k: np.ndarray, theta: float) -> np.ndarray:
    return 2.0 * np.cos(k + sigma * theta)


def q636_covariance_certificate() -> dict[str, Any]:
    theta = math.pi / 5.0
    k_grid = np.linspace(-math.pi, math.pi, 1001)
    e_plus = q636_energy(1, k_grid, theta)
    e_minus = q636_energy(-1, k_grid, theta)
    parity_covariance_defect = float(np.max(np.abs(e_plus - q636_energy(-1, -k_grid, theta))))
    spectral_set_defect = float(np.max(np.abs(np.sort(e_plus) - np.sort(e_minus))))

    # Open-chain Peierls phases are pure gauge: H_sigma(theta) = U H_0 U^†.
    n = 8
    h0 = np.zeros((n, n), dtype=complex)
    htheta = np.zeros((n, n), dtype=complex)
    for i in range(n - 1):
        h0[i, i + 1] = h0[i + 1, i] = 1.0
        htheta[i, i + 1] = np.exp(1j * theta)
        htheta[i + 1, i] = np.exp(-1j * theta)
    u = np.diag([np.exp(-1j * j * theta) for j in range(n)])
    gauge_equivalence_defect = float(np.max(np.abs(htheta - u @ h0 @ u.conj().T)))
    eigen_defect = float(np.max(np.abs(np.linalg.eigvalsh(htheta) - np.linalg.eigvalsh(h0))))

    return {
        "claim_checked": "QW-636 Bloch parity asymmetry is covariance under (sigma,k)->(-sigma,-k), not a physical sign selector unless a directed momentum/orientation source is already fixed.",
        "theta": theta,
        "pointwise_parity_covariance_identity": "E_sigma(k)=E_-sigma(-k)",
        "parity_covariance_max_defect": parity_covariance_defect,
        "sorted_spectrum_sigma_plus_vs_minus_max_defect": spectral_set_defect,
        "open_chain_peierls_phase_is_pure_gauge_identity": "H_theta = U H_0 U^dagger for U_j=exp(i*j*theta)",
        "open_chain_gauge_equivalence_max_defect": gauge_equivalence_defect,
        "open_chain_eigenvalue_max_defect": eigen_defect,
        "methodological_verdict": "QW636_PARITY_DELTA_IS_A_CONDITIONAL_DIAGNOSTIC_NOT_A_STRICT_SELECTOR_SOURCE",
        "failed_requirements": [
            "gauge-invariant directed sign observable",
            "strict source of the orientation of k or loop traversal",
            "proof that the sign is not a convention-layer relabeling",
        ],
    }


def k_complex(d: int) -> complex:
    return ALPHA_GEO * np.exp(1j * (OMEGA * d + PHI)) / (1.0 + BETA_TORS * abs(d))


def q1026_matrix() -> np.ndarray:
    k = np.zeros((N_OCTAVES, N_OCTAVES), dtype=complex)
    for i in range(N_OCTAVES):
        for j in range(N_OCTAVES):
            k[i, j] = k_complex(abs(i - j))
    return k


def q1026_anomaly_certificate() -> dict[str, Any]:
    k = q1026_matrix()
    k3 = k @ k @ k
    gamma5_even_start = np.diag([(-1) ** i for i in range(N_OCTAVES)])
    gamma5_odd_start = -gamma5_even_start
    anomaly_even = np.trace(gamma5_even_start @ k3)
    anomaly_odd = np.trace(gamma5_odd_start @ k3)
    shifted_labels = np.roll(np.arange(N_OCTAVES), 1)
    p = np.eye(N_OCTAVES)[shifted_labels]
    shifted_gamma5 = p @ gamma5_even_start @ p.T
    shifted_anomaly = np.trace(shifted_gamma5 @ k3)

    return {
        "claim_checked": "QW-1026 uses gamma5=diag((-1)^i), so the anomaly sign is tied to an alternating-label convention unless gamma5 is independently sourced.",
        "N_OCTAVES": N_OCTAVES,
        "anomaly_even_start_real": float(np.real(anomaly_even)),
        "anomaly_even_start_imag": float(np.imag(anomaly_even)),
        "anomaly_odd_start_real": float(np.real(anomaly_odd)),
        "anomaly_odd_start_imag": float(np.imag(anomaly_odd)),
        "opposite_gamma5_sign_sum_real": float(np.real(anomaly_even + anomaly_odd)),
        "opposite_gamma5_sign_sum_imag": float(np.imag(anomaly_even + anomaly_odd)),
        "shifted_gamma5_equals_minus_original": bool(np.allclose(shifted_gamma5, -gamma5_even_start)),
        "shifted_anomaly_real": float(np.real(shifted_anomaly)),
        "shifted_anomaly_imag": float(np.imag(shifted_anomaly)),
        "sign_flips_under_allowed_alternating_origin_change": bool(np.real(anomaly_even) * np.real(anomaly_odd) < 0),
        "methodological_verdict": "QW1026_ANOMALY_SIGN_IS_CONVENTION_DEPENDENT_WITHOUT_A_GAMMA5_SOURCE_THEOREM",
        "failed_requirements": [
            "basis-independent chirality operator",
            "orientation-sourced gamma5 convention",
            "index/anomaly theorem tying Re Tr(gamma5 K^3) to a physical topological charge",
        ],
    }


def strict_acceptance_table(q636: dict[str, Any], q1026: dict[str, Any]) -> list[dict[str, Any]]:
    candidates = [
        {
            "candidate": "QW636_hopfion_parity_delta_prior_art_alone",
            "nonzero_odd_diagnostic": True,
            "gauge_invariant_physical_sign": False,
            "basis_or_convention_independent": False,
            "strict_internal_source_theorem": False,
            "reason": q636["methodological_verdict"],
        },
        {
            "candidate": "QW1026_alternating_gamma5_trace_prior_art_alone",
            "nonzero_odd_diagnostic": bool(abs(q1026["anomaly_even_start_real"]) > 1e-12 or abs(q1026["anomaly_even_start_imag"]) > 1e-12),
            "gauge_invariant_physical_sign": False,
            "basis_or_convention_independent": False,
            "strict_internal_source_theorem": False,
            "reason": q1026["methodological_verdict"],
        },
        {
            "candidate": "P2621_typed_nonzero_chiral_source_schema",
            "nonzero_odd_diagnostic": True,
            "gauge_invariant_physical_sign": "premise_required",
            "basis_or_convention_independent": "premise_required",
            "strict_internal_source_theorem": "premise_required",
            "reason": "P2621 is a valid conditional schema only after a typed nonzero chiral source is admitted; QW-636/QW-1026 prior art alone does not supply that premise.",
        },
    ]
    for row in candidates:
        row["accepted_as_orientation_odd_selector_source_now"] = all(
            row[key] is True for key in [
                "nonzero_odd_diagnostic",
                "gauge_invariant_physical_sign",
                "basis_or_convention_independent",
                "strict_internal_source_theorem",
            ]
        )
    return candidates


def p2620_update(accepted_orientation_now: bool) -> dict[str, Any]:
    rows = []
    for damping in (False, True):
        rows.append({
            "assignment": {
                "nonlinear_damping_completion_source": damping,
                "orientation_odd_selector_source_from_qw636_qw1026_prior_art_alone": accepted_orientation_now,
            },
            "bridge_source_cut_repaired": damping and accepted_orientation_now,
            "missing": [name for name, value in {
                "nonlinear_damping_completion_source": damping,
                "orientation_odd_selector_source_from_qw636_qw1026_prior_art_alone": accepted_orientation_now,
            }.items() if not value],
        })
    return {
        "orientation_source_from_prior_art_alone": accepted_orientation_now,
        "rows": rows,
        "accepting_row_count": sum(1 for row in rows if row["bridge_source_cut_repaired"]),
        "verdict": "P2620 remains unrepaired by QW-636/QW-1026 prior art alone.",
    }


def prior_art_reading() -> dict[str, Any]:
    return {
        name: {"path": rel(path), "exists": path.exists(), "sha256": sha256_file(path)}
        for name, path in PRIOR_ART_FILES.items()
    }


def build_payload() -> dict[str, Any]:
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    q636 = q636_covariance_certificate()
    q1026 = q1026_anomaly_certificate()
    table = strict_acceptance_table(q636, q1026)
    accepted_now = any(row["accepted_as_orientation_odd_selector_source_now"] is True for row in table)

    theorem_export = {
        "theorem_name": "P2622 QW-636/QW-1026 prior-art nonpromotion theorem",
        "claim": (
            "QW-636 and QW-1026 contain useful parity/chirality diagnostics, but as written they do not export the missing strict "
            "orientation_odd_selector_source.  QW-636's phase asymmetry is gauge/covariance sensitive and needs an already-directed "
            "momentum or loop orientation.  QW-1026's sign is tied to an unsourced alternating gamma5 convention.  Therefore P2621 "
            "remains only a conditional schema, not an unconditional source theorem."
        ),
        "positive_content": [
            "QW-636 can be retained as a diagnostic for a typed chiral/flux source once a gauge-invariant orientation source is supplied.",
            "QW-1026 can be retained as a diagnostic for a typed chirality operator once gamma5 is physically sourced and basis-safe.",
            "P2621 remains a conditional implication: typed nonzero chiral source implies an odd selector.",
        ],
        "obstructions": [
            "QW-636: E_sigma(k)=E_-sigma(-k), sorted spectra for sigma=+/- match, and open-chain Peierls phases are pure gauge.",
            "QW-1026: gamma5=diag((-1)^i) changes sign under the alternating-origin choice, flipping Re Tr(gamma5 K^3).",
            "Neither script supplies the independent nonlinear damping completion atom required by P2620.",
        ],
        "next_admissible_research_targets": [
            "derive a gauge-invariant Wilson-loop/flux orientation source on a closed typed cycle with fixed orientation convention sourced internally",
            "derive gamma5/chirality as a basis-independent operator from strict nadsoliton field content rather than alternating labels",
            "continue nonlinear damping completion separately; selector progress alone cannot repair P2620",
        ],
        "not_licensed": [
            "promotion of QW-636/QW-1026 prior art alone to strict selector source",
            "unconditional P2621 source closure",
            "P2620 bridge-source cut repair",
            "role-bearing L_total",
            "QW-2191 discharge",
            "ToE closure",
        ],
    }

    payload = {
        "packet_id": "P2622",
        "slice_id": "S1572",
        "status": "P2622_QW636_QW1026_PRIOR_ART_NONPROMOTION_NO_SELECTOR_SOURCE_NO_BRIDGE_NO_LTOTAL_NO_QW2191_NO_TOE",
        "source_artifacts": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "prior_art_files": prior_art_reading(),
        "rg_antiduplication_audit": rg_audit(),
        "inherited_artifact_status": {name: artifact.get("status", artifact.get("packet_id", "UNKNOWN")) for name, artifact in artifacts.items()},
        "q636_covariance_and_gauge_certificate": q636,
        "q1026_gamma5_convention_certificate": q1026,
        "strict_acceptance_table": table,
        "accepted_orientation_source_from_qw636_qw1026_prior_art_alone": accepted_now,
        "p2620_prior_art_update": p2620_update(accepted_now),
        "theorem_export": theorem_export,
        "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    payload["certificate_hash"] = sha256_json({k: v for k, v in payload.items() if k != "certificate_hash"})
    return payload


def write_markdown(payload: dict[str, Any]) -> None:
    theorem = payload["theorem_export"]
    lines = [
        "# P2622/S1572 QW-636/QW-1026 physical-rigor nonpromotion audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Anti-duplication grep audit",
        "",
    ]
    for name, data in payload["rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits; samples retained in JSON certificate.")
    lines.extend([
        "",
        "## Theorem export",
        "",
        f"**Claim.** {theorem['claim']}",
        "",
        "Positive retained content:",
    ])
    for item in theorem["positive_content"]:
        lines.append(f"- {item}")
    lines.append("")
    lines.append("Obstructions:")
    for item in theorem["obstructions"]:
        lines.append(f"- {item}")
    q636 = payload["q636_covariance_and_gauge_certificate"]
    q1026 = payload["q1026_gamma5_convention_certificate"]
    lines.extend([
        "",
        "## Computational certificates",
        "",
        f"- QW-636 parity covariance max defect: `{q636['parity_covariance_max_defect']}`.",
        f"- QW-636 sorted sigma-spectrum defect: `{q636['sorted_spectrum_sigma_plus_vs_minus_max_defect']}`.",
        f"- QW-636 open-chain gauge-equivalence defect: `{q636['open_chain_gauge_equivalence_max_defect']}`.",
        f"- QW-1026 opposite gamma5 sign sum real/imag: `{q1026['opposite_gamma5_sign_sum_real']}`, `{q1026['opposite_gamma5_sign_sum_imag']}`.",
        f"- Accepted orientation source from QW-636/QW-1026 prior art alone: `{payload['accepted_orientation_source_from_qw636_qw1026_prior_art_alone']}`.",
        f"- P2620 accepting rows from prior art alone: `{payload['p2620_prior_art_update']['accepting_row_count']}`.",
        "",
        "## Next admissible targets",
    ])
    for item in theorem["next_admissible_research_targets"]:
        lines.append(f"- {item}")
    lines.extend(["", "Not licensed:"])
    for item in theorem["not_licensed"]:
        lines.append(f"- {item}")
    lines.extend(["", f"Certificate hash: `{payload['certificate_hash']}`", ""])
    MD.write_text("\n".join(lines), encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2622/S1572 QW-636/QW-1026 physical-rigor nonpromotion audit

`P2622/S1572` corrects the P2621 lane by auditing whether the QW-636/QW-1026 prior art itself supplies the missing strict orientation source.  It does not: QW-636's Hopfion phase asymmetry is parity-covariant and gauge/convention sensitive unless a directed momentum/loop orientation source is already supplied, while QW-1026's `Tr(gamma5 K^3)` sign flips with the unsourced alternating `gamma5` origin.  Thus P2621 remains a conditional schema only; QW-636/QW-1026 prior art alone does not export `orientation_odd_selector_source`, does not repair P2620, and does not reopen role-bearing `L_total`, QW-2191, or ToE closure.
"""
    lag_section = """
## P2622/S1572 QW-636/QW-1026 nonpromotion Ltotal guard

`P2622/S1572` keeps role-bearing `L_total` closed and sharpens the selector guard: the old QW-636 and QW-1026 scripts are retained as diagnostics for a future typed chiral/flux source, but are not themselves a strict source theorem because their signs remain gauge/convention dependent without additional physical orientation and `gamma5` premises.  The independent nonlinear damping completion atom also remains open.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2622/S1572 QW-636/QW-1026 physical-rigor nonpromotion audit", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2622/S1572 QW-636/QW-1026 nonpromotion Ltotal guard", lag_section)


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps({
        "packet_id": payload["packet_id"],
        "status": payload["status"],
        "accepted_orientation_source_from_prior_art_alone": payload["accepted_orientation_source_from_qw636_qw1026_prior_art_alone"],
        "q636_open_chain_gauge_defect": payload["q636_covariance_and_gauge_certificate"]["open_chain_gauge_equivalence_max_defect"],
        "q1026_gamma5_sign_flips": payload["q1026_gamma5_convention_certificate"]["sign_flips_under_allowed_alternating_origin_change"],
        "p2620_accepting_row_count": payload["p2620_prior_art_update"]["accepting_row_count"],
        "certificate_hash": payload["certificate_hash"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
