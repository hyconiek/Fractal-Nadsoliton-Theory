#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import itertools
import json
import math
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2621_s1571_qw636_qw1026_chiral_hopfion_selector_source_audit.json"
MD = GEN / "p2621_s1571_qw636_qw1026_chiral_hopfion_selector_source_audit.md"

SOURCE_FILES = {
    "P2618_ANALYTIC_OBSTRUCTION": GEN / "p2618_s1568_analytic_legacy_to_strict_completion_obstruction.json",
    "P2619_SELECTOR_LATTICE": GEN / "p2619_s1569_p2618_selector_source_obligation_lattice.json",
    "P2620_TWO_OBSTRUCTION_CUT": GEN / "p2620_s1570_p2618_p2619_bridge_two_obstruction_cut.json",
    "P2366_PHASE_REFERENCE_CHI11_AUDIT": GEN / "p2366_s1316_selector_candidate_phase_reference_chi11_audit_probe.json",
}

PRIOR_ART_FILES = {
    "QW636_FINAL_PARITY_CHECK": REPO / "TOE_Proven_Hypotheses_v3.1_FINAL" / "QW-636_Parity_Check.py",
    "QW636_ROOT_PARITY_CHECK": REPO / "QW-636_Parity_Check.py",
    "QW1026_SYNTHESIS_SUITE": REPO / "QW-1017_to_QW-1036_Synthesis_Suite.py",
}

NEGATIVE_EXPORT_FLAGS = [
    "strict_core_selector_source_exported_unconditionally",
    "p2620_bridge_source_cut_repaired",
    "nonlinear_damping_completion_source_exported",
    "full_bridge_revalidated",
    "role_transfer_revalidated",
    "role_bearing_ltotal_reenabled",
    "gf2_bridge_revalidated",
    "beta_tors_chi11_route_reopened",
    "qw2191_discharged_by_this_packet",
    "toe_closure_claimed",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else "MISSING"


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, ".",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "!fundamental_action_reconstruction/generated/**",
            "-g", "!.git/**",
        ],
        cwd=REPO,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:120]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2621|S1571|chiral hopfion selector|QW-636.*QW-1026|Hopfion.*anomaly.*selector",
        "qw636_prior_art": "QW-636|Hopfion Phase|Complex Phase.*break Parity|t\\(r\\).*exp\\(i.*theta|E\\(k\\).*E\\(-k\\)",
        "qw1026_prior_art": "QW-1026|Chiral Anomaly|Tr\\(γ₅|Tr\\[γ_5|gamma5|K³|K\\^3|anomaly_real",
        "selector_chirality_prior_art": "chiral.*selector|chirality.*selector|phase_origin_plus_chiral_bispectrum|orientation_odd_selector_source|QW-2191|chi11_selector_source",
        "bridge_guards": "P2618|P2619|P2620|nonlinear_damping_completion_source|role-bearing L_total|role-transfer theorem|ToE closure",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def hopfion_energy(sigma: int, k: float, theta: float, amplitudes: list[float]) -> float:
    return sum(2.0 * amp * math.cos((idx + 1) * k + sigma * theta) for idx, amp in enumerate(amplitudes))


def hopfion_selector_rows() -> list[dict[str, Any]]:
    theta = math.pi / 7.0
    k = math.pi / 5.0
    amplitudes = [1.0, 0.5, 1.0 / 3.0, 0.25]
    sine_sum = sum(amp * math.sin((idx + 1) * k) for idx, amp in enumerate(amplitudes))
    rows = []
    for sigma in (-1, 1):
        e_plus = hopfion_energy(sigma, k, theta, amplitudes)
        e_minus = hopfion_energy(sigma, -k, theta, amplitudes)
        delta = e_plus - e_minus
        recovered = -delta / (4.0 * math.sin(theta) * sine_sum)
        rows.append({
            "sigma": sigma,
            "theta": theta,
            "k_probe": k,
            "positive_amplitudes": amplitudes,
            "sine_probe_sum": sine_sum,
            "E_sigma_k": e_plus,
            "E_sigma_minus_k": e_minus,
            "parity_odd_energy_delta": delta,
            "exact_identity": "Delta_sigma(k) = -4*sigma*sin(theta)*sum_r a_r*sin(k*r)",
            "recovered_sigma_from_delta": round(recovered),
            "recovery_error_abs": abs(recovered - sigma),
            "selector_defined": abs(math.sin(theta) * sine_sum) > 1e-12,
        })
    return rows


def matmul(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    n = len(a)
    return [[sum(a[i][j] * b[j][k] for j in range(n)) for k in range(n)] for i in range(n)]


def trace_gamma5_k3(sigma: int, a: float = 2.0) -> float:
    # Minimal exact chiral block: K_sigma = diag(sigma*a, -sigma*a), gamma5 = diag(1, -1).
    # Then Tr(gamma5 K_sigma^3) = 2*sigma*a^3; it flips under orientation reversal.
    k = [[sigma * a, 0.0], [0.0, -sigma * a]]
    gamma5 = [[1.0, 0.0], [0.0, -1.0]]
    k3 = matmul(matmul(k, k), k)
    return sum(gamma5[i][i] * k3[i][i] for i in range(2))


def anomaly_selector_rows() -> list[dict[str, Any]]:
    rows = []
    for sigma in (-1, 1):
        anomaly_real = trace_gamma5_k3(sigma)
        rows.append({
            "sigma": sigma,
            "minimal_chiral_block": "K_sigma=diag(sigma*a,-sigma*a), gamma5=diag(1,-1), a=2",
            "Re_Tr_gamma5_K3": anomaly_real,
            "closed_form": "Re A_sigma = 2*sigma*a^3",
            "recovered_sigma_from_sign_Re_A": 1 if anomaly_real > 0 else -1,
            "orientation_odd": trace_gamma5_k3(-sigma) == -anomaly_real,
            "selector_defined": anomaly_real != 0.0,
        })
    return rows


def c2_selector_enumeration() -> dict[str, Any]:
    points = ["chiral_minus", "chiral_plus"]

    def act(g: int, point: str) -> str:
        if g == 0:
            return point
        return "chiral_plus" if point == "chiral_minus" else "chiral_minus"

    def sign_act(g: int, value: int) -> int:
        return value if g == 0 else -value

    equivariant = []
    for values in itertools.product([-1, 1], repeat=2):
        fmap = dict(zip(points, values))
        if all(fmap[act(g, x)] == sign_act(g, fmap[x]) for g in (0, 1) for x in points):
            equivariant.append(fmap)
    canonical = {"chiral_minus": -1, "chiral_plus": 1}
    return {
        "input_action": "C2 acts freely by chiral_plus <-> chiral_minus",
        "codomain_action": "C2 acts by sign negation on Sigma={-1,+1}",
        "candidate_function_count": 4,
        "equivariant_function_count": len(equivariant),
        "equivariant_functions": equivariant,
        "canonical_orientation_convention_selector": canonical,
        "canonical_selector_is_equivariant": canonical in equivariant,
    }


def p2620_conditional_update() -> dict[str, Any]:
    rows = []
    for nonlinear in (False, True):
        for chiral_source in (False, True):
            rows.append({
                "assignment": {
                    "nonlinear_damping_completion_source": nonlinear,
                    "orientation_odd_selector_source_from_chiral_hopfion_anomaly": chiral_source,
                },
                "p2620_bridge_source_cut_repaired": nonlinear and chiral_source,
                "orientation_obstruction_repaired_conditionally": chiral_source,
                "missing": [name for name, value in {
                    "nonlinear_damping_completion_source": nonlinear,
                    "orientation_odd_selector_source_from_chiral_hopfion_anomaly": chiral_source,
                }.items() if not value],
            })
    return {
        "interpretation": "QW-636/QW-1026 can supply only the orientation atom when promoted to a typed nonzero chiral source premise; P2620 still also requires nonlinear damping completion.",
        "rows": rows,
        "accepting_row_count": sum(1 for row in rows if row["p2620_bridge_source_cut_repaired"]),
    }


def prior_art_reading() -> dict[str, Any]:
    return {
        name: {
            "path": rel(path),
            "exists": path.exists(),
            "sha256": sha256_file(path),
        }
        for name, path in PRIOR_ART_FILES.items()
    }


def build_payload() -> dict[str, Any]:
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    hopfion_rows = hopfion_selector_rows()
    anomaly_rows = anomaly_selector_rows()
    c2 = c2_selector_enumeration()
    p2620_update = p2620_conditional_update()

    theorem_export = {
        "theorem_name": "P2621 conditional chiral-Hopfion orientation-odd selector theorem",
        "claim": (
            "If the strict kernel is supplied with a typed nonzero chiral Hopfion/anomaly source, "
            "then sigma=sign(Re Tr(gamma5 K^3)) (equivalently the normalized QW-636 parity-odd "
            "Hopfion energy delta) is a C2-odd topological phase selector.  This repairs the P2619/P2620 "
            "orientation atom only under that explicit source premise."
        ),
        "proof_steps": [
            "For QW-636-style directional Hermitian hopping t_sigma(r)=a_r exp(i sigma theta), E_sigma(k)-E_sigma(-k)=-4 sigma sin(theta) sum_r a_r sin(k r).",
            "When sin(theta) and the declared sine probe are nonzero, the normalized parity-odd energy delta recovers sigma exactly and flips under orientation reversal.",
            "For QW-1026-style anomaly A_sigma=Tr(gamma5 K_sigma^3), a typed chiral kernel satisfying Re A_-sigma=-Re A_sigma and Re A_sigma != 0 yields sigma=sign(Re A_sigma) after an orientation convention.",
            "The C2 enumeration verifies that a freely acting chiral two-point torsor admits equivariant sign maps, unlike the invariant scalar inputs blocked by P2619.",
        ],
        "professorial_assessment_of_prompt": (
            "Scientifically usable only after narrowing: QW-636 and QW-1026 are heuristic/numerical prior-art scripts, not strict proofs by themselves. "
            "They can be formalized into a conditional theorem exporting an orientation-odd source premise, but they do not unconditionally discharge QW-2191 or repair the full P2620 bridge."
        ),
        "orientation_atom_status": "CONDITIONALLY_REPAIRED_IF_TYPED_NONZERO_CHIRAL_SOURCE_IS_ADMITTED",
        "not_licensed": [
            "unconditional strict-core selector closure",
            "nonlinear damping completion source",
            "full P2620 bridge-source cut repair without the damping atom",
            "GF(2) bridge revalidation",
            "beta_tors -> chi11 route reopening",
            "role-bearing L_total re-enable",
            "legacy physical-role transfer",
            "QW-2191 discharge",
            "ToE closure",
        ],
    }

    payload = {
        "packet_id": "P2621",
        "slice_id": "S1571",
        "status": "P2621_CONDITIONAL_CHIRAL_HOPFION_SELECTOR_SOURCE_AUDIT_NO_FULL_BRIDGE_NO_LTOTAL_NO_QW2191_NO_TOE",
        "source_artifacts": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "prior_art_files": prior_art_reading(),
        "rg_antiduplication_audit": rg_audit(),
        "inherited_artifact_status": {
            name: artifact.get("status", artifact.get("packet_id", "UNKNOWN")) for name, artifact in artifacts.items()
        },
        "hopfion_energy_selector_certificate": {
            "analytic_identity": "For t_sigma(r)=a_r exp(i sigma theta), E_sigma(k)=2 sum_r a_r cos(k r + sigma theta), hence Delta_sigma(k)=-4 sigma sin(theta) sum_r a_r sin(k r).",
            "rows": hopfion_rows,
            "all_rows_recover_sigma": all(row["recovered_sigma_from_delta"] == row["sigma"] and row["recovery_error_abs"] < 1e-12 for row in hopfion_rows),
            "all_rows_selector_defined": all(row["selector_defined"] for row in hopfion_rows),
        },
        "chiral_anomaly_selector_certificate": {
            "analytic_identity": "Typed nonzero QW-1026 anomaly source: A_sigma=Tr(gamma5 K_sigma^3), Re A_-sigma=-Re A_sigma; minimal block gives Re A_sigma=2 sigma a^3.",
            "rows": anomaly_rows,
            "all_rows_recover_sigma": all(row["recovered_sigma_from_sign_Re_A"] == row["sigma"] for row in anomaly_rows),
            "all_rows_orientation_odd": all(row["orientation_odd"] for row in anomaly_rows),
            "all_rows_selector_defined": all(row["selector_defined"] for row in anomaly_rows),
        },
        "c2_equivariant_selector_enumeration": c2,
        "p2620_conditional_update": p2620_update,
        "theorem_export": theorem_export,
        "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    payload["certificate_hash"] = sha256_json({k: v for k, v in payload.items() if k != "certificate_hash"})
    return payload


def write_markdown(payload: dict[str, Any]) -> None:
    theorem = payload["theorem_export"]
    lines = [
        "# P2621/S1571 QW-636/QW-1026 chiral-Hopfion selector source audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Professorial verdict on the proposed next step",
        "",
        theorem["professorial_assessment_of_prompt"],
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
        "Proof skeleton:",
    ])
    for step in theorem["proof_steps"]:
        lines.append(f"- {step}")
    lines.extend([
        "",
        "## Computational certificates",
        "",
        f"- Hopfion energy selector recovers sigma on all rows: `{payload['hopfion_energy_selector_certificate']['all_rows_recover_sigma']}`.",
        f"- Chiral anomaly selector recovers sigma on all rows: `{payload['chiral_anomaly_selector_certificate']['all_rows_recover_sigma']}`.",
        f"- C2 equivariant sign maps from the chiral torsor: `{payload['c2_equivariant_selector_enumeration']['equivariant_function_count']}`.",
        f"- P2620 conditional accepting rows: `{payload['p2620_conditional_update']['accepting_row_count']}` of 4.",
        "",
        "## Guarded conclusion",
        "",
        f"Orientation atom status: `{theorem['orientation_atom_status']}`.",
        "",
        "This does **not** repair the nonlinear damping completion atom and therefore does **not** by itself reopen role-bearing `L_total`, role transfer, QW-2191 discharge, or ToE closure.",
        "",
        "Not licensed:",
    ])
    for item in theorem["not_licensed"]:
        lines.append(f"- {item}")
    lines.extend(["", f"Certificate hash: `{payload['certificate_hash']}`", ""])
    MD.write_text("\n".join(lines), encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2621/S1571 QW-636/QW-1026 chiral-Hopfion selector source audit

`P2621/S1571` audits the proposed use of the old QW-636 Hopfion parity phase and QW-1026 chiral anomaly as the missing P2620 orientation source.  The usable theorem is conditional: if a typed, nonzero chiral Hopfion/anomaly source is admitted, then the normalized parity-odd energy difference and `sign(Re Tr(gamma5 K^3))` define a `C2`-odd selector for the strict phase sign.  This conditionally supplies only `orientation_odd_selector_source`; it does not supply `nonlinear_damping_completion_source`, does not revalidate the full bridge, and does not reopen role-bearing `L_total`, role-transfer, QW-2191, or ToE closure.
"""
    lag_section = """
## P2621/S1571 chiral-Hopfion selector-source Ltotal guard

`P2621/S1571` keeps role-bearing `L_total` closed.  QW-636/QW-1026 prior art can be formalized as a conditional chiral source theorem for the orientation atom, but the source must be typed and nonzero; the old scripts are heuristic prior art rather than unconditional strict-core closure.  P2620 still requires the independent nonlinear damping completion atom before bridge validity or role transfer may be rerun.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2621/S1571 QW-636/QW-1026 chiral-Hopfion selector source audit", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2621/S1571 chiral-Hopfion selector-source Ltotal guard", lag_section)


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps({
        "packet_id": payload["packet_id"],
        "status": payload["status"],
        "orientation_atom_status": payload["theorem_export"]["orientation_atom_status"],
        "hopfion_rows_recover_sigma": payload["hopfion_energy_selector_certificate"]["all_rows_recover_sigma"],
        "anomaly_rows_recover_sigma": payload["chiral_anomaly_selector_certificate"]["all_rows_recover_sigma"],
        "p2620_accepting_row_count": payload["p2620_conditional_update"]["accepting_row_count"],
        "certificate_hash": payload["certificate_hash"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
