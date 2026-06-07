#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2555_s1505_legacy_to_strict_damping_denominator_nonrenormalization_certificate.json"
MD = GEN / "p2555_s1505_legacy_to_strict_damping_denominator_nonrenormalization_certificate.md"

SOURCE_FILES = {
    "P2554_LOCAL_EXHAUSTION_REORIENTATION": GEN / "p2554_s1504_strict_damping_local_exhaustion_bridge_reorientation_certificate.json",
    "F158_DAMPING_NONRENORMALIZATION_PACKET": ROOT / "F158_FIRST_ACTUAL_LEGACY_TO_STRICT_KERNEL_DAMPING_NONRENORMALIZATION_OBSTRUCTION_WITNESS_PACKET.md",
    "S2_STRATEGIC_REORIENTATION": ROOT / "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md",
}

DOMAIN = list(range(1, 12))
BETA_TORS_CANDIDATES = [0.01, 0.05]
STRICT_BETA = 1.0
STRICT_ETA = 9.0 / 5.0
NEGATIVE_EXPORT_FLAGS = [
    "damping_compression_bridge_component_ready", "legacy_to_strict_completion_bridge_exported", "full_bridge_theorem_exported",
    "role_transfer_theorem_exported", "legacy_role_transfer_claimed", "beta_tors_to_beta_eta_translation_exported",
    "strict_damping_beta_eta_source_exported", "source_obligation_discharge_exported", "qw2191_discharged_by_this_certificate",
    "role_bearing_ltotal_exported", "toe_closure_claimed",
]


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run([
        "rg", "-n", pattern, "fundamental_action_reconstruction", "material_dowodowy",
        "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "!fundamental_action_reconstruction/generated/**",
    ], cwd=REPO, check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:80]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2555|S1505|damping denominator nonrenormalization|denominator nonrenormalization|linear nonlinear damping bridge",
        "intended_research_nonduplication": "legacy.*strict.*damping.*denominator|beta_tors.*beta.*eta|linear.*nonlinear.*denominator|damping compression passage|nonrenormalization obstruction",
        "bridge_precursors": "F158|P2554|S1504|legacy -> strict completion bridge|K_legacy_ont|K_strict_gate",
        "guardrails": "role-transfer audit|QW-2191|ToE closure|selector guardrail|silent.*role transfer",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def denom_legacy(beta_tors: float, d: int) -> float:
    return 1.0 + beta_tors * d


def denom_strict(d: int) -> float:
    return 1.0 + STRICT_BETA * (d ** STRICT_ETA)


def second_differences(values: list[float]) -> list[float]:
    return [values[i + 2] - 2.0 * values[i + 1] + values[i] for i in range(len(values) - 2)]


def damping_rows() -> list[dict[str, Any]]:
    rows = []
    strict_values = [denom_strict(d) for d in DOMAIN]
    strict_second = second_differences(strict_values)
    for beta_tors in BETA_TORS_CANDIDATES:
        legacy_values = [denom_legacy(beta_tors, d) for d in DOMAIN]
        legacy_second = second_differences(legacy_values)
        ratios = [strict / legacy for strict, legacy in zip(strict_values, legacy_values)]
        inverse_ratios = [legacy / strict for strict, legacy in zip(strict_values, legacy_values)]
        rows.append({
            "beta_tors": beta_tors,
            "legacy_denominator": "1 + beta_tors*d",
            "strict_denominator": "1 + beta*d^eta with beta=1, eta=9/5",
            "legacy_values_d_1_to_11": legacy_values,
            "strict_values_d_1_to_11": strict_values,
            "legacy_second_difference_max_abs": max(abs(x) for x in legacy_second),
            "strict_second_difference_min": min(strict_second),
            "strict_second_differences_all_positive": all(x > 0 for x in strict_second),
            "raw_denominator_equality_on_domain": all(abs(a - b) < 1e-12 for a, b in zip(legacy_values, strict_values)),
            "constant_proportionality_ratio_spread": max(ratios) - min(ratios),
            "inverse_ratio_spread": max(inverse_ratios) - min(inverse_ratios),
            "constant_amplitude_can_absorb_denominator_mismatch": max(ratios) - min(ratios) < 1e-12,
            "endpoint_ratio_d1": ratios[0],
            "endpoint_ratio_d11": ratios[-1],
        })
    return rows


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2554_payload = load_json(SOURCE_FILES["P2554_LOCAL_EXHAUSTION_REORIENTATION"])
    p2554 = theorem(p2554_payload, "strict_damping_local_exhaustion_bridge_reorientation_certificate")
    text_sources = {name: path.read_bytes() for name, path in SOURCE_FILES.items() if path.suffix == ".md"}
    rows = damping_rows()
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2555_T1_legacy_to_strict_damping_denominator_nonrenormalization_certificate",
        "audited_chain": ["F158", "P2554/S1504", "S2"],
        "bridge_component_under_attack": "damping_compression_passage_from_legacy_linear_torsion_denominator_to_strict_nonlinear_compression_denominator",
        "p2554_bridge_reorientation_inherited": p2554.get("local_route_exhaustion_exported_as_reorientation_not_closure") is True,
        "domain": DOMAIN,
        "strict_parameters_used_for_obstruction": {"beta": STRICT_BETA, "eta": STRICT_ETA},
        "beta_tors_candidates_audited": BETA_TORS_CANDIDATES,
        "damping_denominator_nonrenormalization_rows": rows,
        "legacy_linear_second_difference_zero_for_all_candidates": all(row["legacy_second_difference_max_abs"] < 1e-12 for row in rows),
        "strict_nonlinear_second_difference_positive": all(row["strict_second_differences_all_positive"] for row in rows),
        "raw_denominator_identity_refuted_on_domain": all(not row["raw_denominator_equality_on_domain"] for row in rows),
        "constant_amplitude_absorption_refuted_on_domain": all(not row["constant_amplitude_can_absorb_denominator_mismatch"] for row in rows),
        "beta_tors_to_beta_eta_scalar_renormalization_refuted_for_audited_class": True,
        "required_bridge_addition": "a non-scalar nonlinear compression/source map that supplies beta=1 and eta=9/5 rather than a silent beta_tors renormalization",
        "recommended_next_honest_step": (
            "Build the explicit legacy->strict damping/compression completion map, if it exists, by specifying the non-scalar nonlinear source that changes the linear beta_tors*d denominator into beta*d^eta. Do not transfer legacy EM/Weinberg/gravity roles until that bridge and a separate role-transfer theorem exist."
        ),
        "not_licensed": [
            "No beta_tors -> (beta, eta) translation is exported.",
            "No raw identity or constant-amplitude renormalization from K_legacy_ont to K_strict_gate is exported.",
            "No bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing L_total, or ToE closure is exported.",
        ],
    }
    for flag in NEGATIVE_EXPORT_FLAGS:
        theorem_export[flag] = False
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "p2554_reorientation_inherited": theorem_export["p2554_bridge_reorientation_inherited"],
        "linear_vs_nonlinear_curvature_separation": theorem_export["legacy_linear_second_difference_zero_for_all_candidates"] and theorem_export["strict_nonlinear_second_difference_positive"],
        "raw_identity_refuted": theorem_export["raw_denominator_identity_refuted_on_domain"],
        "constant_amplitude_absorption_refuted": theorem_export["constant_amplitude_absorption_refuted_on_domain"],
        "no_false_bridge_or_role_transfer": theorem_export["full_bridge_theorem_exported"] is False and theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_or_toe_claim": theorem_export["qw2191_discharged_by_this_certificate"] is False and theorem_export["toe_closure_claimed"] is False,
    }
    source_fingerprints = {"P2554_LOCAL_EXHAUSTION_REORIENTATION": sha256_json(p2554_payload)}
    source_fingerprints.update({name: sha256_bytes(data) for name, data in text_sources.items()})
    return {
        "packet_id": "P2555",
        "stage_id": "S1505",
        "status": "LEGACY_TO_STRICT_DAMPING_DENOMINATOR_NONRENORMALIZATION_CERTIFICATE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "legacy_to_strict_damping_denominator_nonrenormalization_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": source_fingerprints,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["legacy_to_strict_damping_denominator_nonrenormalization_certificate"]["theorem_export"]
    lines = [
        "# P2555/S1505 legacy-to-strict damping denominator nonrenormalization certificate", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Bridge component under attack: `{t['bridge_component_under_attack']}`.",
        f"- Audited beta_tors candidates: `{t['beta_tors_candidates_audited']}`.",
        f"- Legacy linear second differences vanish: `{t['legacy_linear_second_difference_zero_for_all_candidates']}`.",
        f"- Strict nonlinear second differences are positive: `{t['strict_nonlinear_second_difference_positive']}`.",
        f"- Raw denominator identity refuted on domain: `{t['raw_denominator_identity_refuted_on_domain']}`.",
        f"- Constant-amplitude absorption refuted on domain: `{t['constant_amplitude_absorption_refuted_on_domain']}`.", "", "## Interpretation", "",
        "The legacy damping denominator `1+beta_tors*d` is linear in `d`, while the strict denominator `1+beta*d^eta` with `beta=1, eta=9/5` is strictly convex on the audited domain.  Therefore the strict nonlinear compression cannot be obtained by a raw identity or scalar amplitude renormalization of the legacy linear torsion damping.",
        "", "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No `beta_tors -> (beta, eta)` translation, bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, legacy-role transfer, or ToE closure is exported.",
        "", "## Fingerprint", "", f"`{payload['legacy_to_strict_damping_denominator_nonrenormalization_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
`P2555/S1505` begins the recommended `legacy -> strict` bridge audit at the damping denominator.  On `d=1..11`, the legacy denominator `1+beta_tors*d` has zero second finite differences, while the strict denominator `1+d^(9/5)` has positive second finite differences.  For audited legacy candidates `beta_tors=0.01` and `0.05`, neither raw denominator identity nor constant-amplitude absorption can map the linear legacy torsion damping into the nonlinear strict compression.  A valid bridge therefore needs a real nonlinear compression/source map, not a scalar `beta_tors -> (beta,eta)` renormalization.
""".strip()
    lag_section = """
`P2555/S1505` blocks promotion of legacy linear torsion damping into the strict `L_total` compression term by scalar renormalization.  The damping/compression bridge still requires an explicit nonlinear source map, followed only later by a separate role-transfer theorem.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2555/S1505 legacy-to-strict damping denominator nonrenormalization guard", "## P2555/S1505 legacy-to-strict damping denominator nonrenormalization guard\n\n" + eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2555/S1505 damping denominator nonrenormalization Ltotal guard", "## P2555/S1505 damping denominator nonrenormalization Ltotal guard\n\n" + lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
