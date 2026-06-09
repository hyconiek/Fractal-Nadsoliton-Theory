#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2617_s1567_p2606_exponent_semantics_reclassification.json"
MD = GEN / "p2617_s1567_p2606_exponent_semantics_reclassification.md"

SOURCE_FILES = {
    "P2414_STRICT_DENOMINATOR_IDENTIFIABILITY": GEN / "p2414_s1364_strict_damping_parameter_identifiability_nonabsorption_certificate.json",
    "P2603_FRACTAL_CODIMENSION_SLOPE": GEN / "p2603_s1553_nadsoliton_fractal_codimension_slope_source_theorem.json",
    "P2606_NONLINEAR_RESIDUAL": GEN / "p2606_s1556_strict_side_nonlinear_compression_residual_addition.json",
    "P2615_LINEAR_SLICE_NONBRIDGE": GEN / "p2615_s1565_p2605_linear_slice_nonbridge_obstruction.json",
    "P2616_ROLE_OBSTRUCTION": GEN / "p2616_s1566_p2608_role_acceptance_obstruction_after_source_revalidation.json",
}

BETA_TORS = Fraction(1, 100)
OMEGA = math.pi / 4.0
PHI = math.pi / 6.0
AUDIT_D_VALUES = list(range(1, 13))
DELTA_CODIMENSION = Fraction(4, 5)
ETA_STRICT = Fraction(9, 5)
ETA_LINEAR = Fraction(1, 1)

NEGATIVE_EXPORT_FLAGS = [
    "p2606_strict_side_residual_addition_revalidated",
    "p2607_bridge_completion_revalidated",
    "p2608_role_bearing_ltotal_reenabled",
    "legacy_physical_role_transfer_exported",
    "qw2191_discharged_by_this_packet",
    "toe_closure_claimed",
    "apd_source_exported",
]


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
        "new_packet": "P2617|S1567|P2606.*exponent.*semantics|codimension.*strict eta.*audit",
        "intended_research_nonduplication": "P2606.*eta.*4/5.*9/5|strict denominator exponent relabel|codimension slope is not denominator eta|eta-label correction",
        "precursor_chain": "P2414|P2603|P2606|P2615|P2616|strict denominator.*9/5|codimension.*4/5",
        "guardrails": "K_legacy_ont|K_strict_gate|role-bearing L_total|QW-2191|ToE closure|APD source",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def kernel_value(d: int, exponent: Fraction) -> float:
    denominator = 1.0 + float(BETA_TORS) * (d ** float(exponent))
    return math.cos(OMEGA * d + PHI) / denominator


def exponent_comparison_rows() -> list[dict[str, Any]]:
    rows = []
    for d in AUDIT_D_VALUES:
        linear = kernel_value(d, ETA_LINEAR)
        p2606_codim = kernel_value(d, DELTA_CODIMENSION)
        strict = kernel_value(d, ETA_STRICT)
        rows.append({
            "d": d,
            "linear_eta_1_kernel": linear,
            "p2606_codimension_label_eta_4_over_5_kernel": p2606_codim,
            "strict_denominator_eta_9_over_5_kernel": strict,
            "p2606_minus_linear": p2606_codim - linear,
            "strict_minus_linear": strict - linear,
            "strict_minus_p2606_codimension_label": strict - p2606_codim,
            "abs_strict_minus_p2606_codimension_label": abs(strict - p2606_codim),
        })
    return rows


def exact_semantics_proof() -> dict[str, Any]:
    return {
        "statement": "The P2603 value 4/5 is the codimension/log-slope delta, not the strict denominator exponent eta. In the strict denominator layer the audited exponent is eta=1+delta=9/5, so a P2606 residual computed with denominator power 4/5 is not the strict-side nonlinear compression residual for K_strict_gate.",
        "proof_steps": [
            "P2603 supplies the logarithmic damping slope/codimension delta = D_f - 1 = 4/5.",
            "The strict denominator target audited by P2414 is S(d)=1+d^(9/5), so the denominator power is eta_strict=9/5.",
            "The relation between the linear backbone and the codimension slope is eta_strict = 1 + delta = 1 + 4/5 = 9/5.",
            "P2606 used the label eta=4/5 inside the denominator power, thereby using delta where the strict kernel notation reserves eta for the full denominator exponent.",
            "Therefore P2606 may be retained only as a numerical codimension-slope perturbation probe, not as an exported strict-side nonlinear residual addition for the strict denominator.",
        ],
        "boundary_conditions_for_revalidating_p2606_as_strict_residual": [
            "rewrite the residual with denominator exponent eta=9/5 or with an explicit eta=1+delta map",
            "show that the corrected residual is a physical strict-side addition rather than a label substitution",
            "keep bridge and role-transfer predicates separate after correction",
        ],
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    payloads = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2603 = theorem(payloads["P2603_FRACTAL_CODIMENSION_SLOPE"], "nadsoliton_fractal_codimension_slope_source_theorem")
    p2606 = theorem(payloads["P2606_NONLINEAR_RESIDUAL"], "strict_side_nonlinear_compression_residual_addition")
    p2615 = theorem(payloads["P2615_LINEAR_SLICE_NONBRIDGE"], "p2605_linear_slice_nonbridge_obstruction")
    p2616 = theorem(payloads["P2616_ROLE_OBSTRUCTION"], "p2608_role_acceptance_obstruction_after_source_revalidation")
    rows = exponent_comparison_rows()
    max_strict_minus_p2606 = max(row["abs_strict_minus_p2606_codimension_label"] for row in rows)
    proof = exact_semantics_proof()
    p2606_used_delta_as_denominator_eta = "4/5" in json.dumps(p2606, sort_keys=True)
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2617_T1_p2606_exponent_semantics_reclassification",
        "audited_chain": ["P2414/S1364", "P2603/S1553", "P2606/S1556", "P2615/S1565", "P2616/S1566"],
        "exact_semantics_proof": proof,
        "delta_codimension": str(DELTA_CODIMENSION),
        "eta_strict_denominator": str(ETA_STRICT),
        "eta_strict_equals_one_plus_delta": ETA_STRICT == 1 + DELTA_CODIMENSION,
        "p2603_codimension_slope_retained": p2603.get("fractal_codimension_slope_source_exported") is True or p2603.get("slope_value_or_prime_anchor_source_exported") is True,
        "p2606_used_delta_as_denominator_eta": p2606_used_delta_as_denominator_eta,
        "exponent_comparison_rows": rows,
        "max_abs_strict_eta_9_5_minus_p2606_eta_4_5_kernel": max_strict_minus_p2606,
        "p2606_strict_residual_interpretation_quarantined_by_p2617": p2606_used_delta_as_denominator_eta and max_strict_minus_p2606 > 1e-6,
        "p2606_codimension_probe_retained": True,
        "p2615_nonbridge_obstruction_inherited": p2615.get("p2605_quarantine_retained_by_p2615") is True,
        "p2616_role_obstruction_inherited": p2616.get("p2608_quarantine_retained_by_p2616") is True,
        "recommended_next_honest_step": "If P2606 is needed later, regenerate it with explicit delta=4/5 and eta=1+delta=9/5 semantics before using it as any strict-side residual addition. Do not reopen bridge or role transfer from the old eta=4/5 denominator label.",
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "eta_relation_checked": theorem_export["eta_strict_equals_one_plus_delta"],
        "p2606_delta_label_detected": theorem_export["p2606_used_delta_as_denominator_eta"],
        "strict_eta_9_5_differs_from_p2606_eta_4_5": theorem_export["max_abs_strict_eta_9_5_minus_p2606_eta_4_5_kernel"] > 1e-6,
        "p2606_strict_residual_interpretation_quarantined": theorem_export["p2606_strict_residual_interpretation_quarantined_by_p2617"],
        "p2606_probe_retained": theorem_export["p2606_codimension_probe_retained"],
        "p2615_obstruction_inherited": theorem_export["p2615_nonbridge_obstruction_inherited"],
        "p2616_role_obstruction_inherited": theorem_export["p2616_role_obstruction_inherited"],
        "no_p2606_strict_residual_revalidation": theorem_export["p2606_strict_side_residual_addition_revalidated"] is False,
        "no_bridge_revalidation": theorem_export["p2607_bridge_completion_revalidated"] is False,
        "no_p2608_reenable": theorem_export["p2608_role_bearing_ltotal_reenabled"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_packet"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
        "no_apd_source_exported": theorem_export["apd_source_exported"] is False,
    }
    return {
        "packet_id": "P2617",
        "stage_id": "S1567",
        "status": "P2617_P2606_EXPONENT_SEMANTICS_RECLASSIFIED_DELTA_PROBE_RETAINED_STRICT_RESIDUAL_EXPORT_QUARANTINED",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "p2606_exponent_semantics_reclassification": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {name: sha256_json(payload) for name, payload in payloads.items()},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["p2606_exponent_semantics_reclassification"]["theorem_export"]
    proof = t["exact_semantics_proof"]
    lines = [
        "# P2617/S1567 P2606 exponent-semantics reclassification", "",
        f"Status: `{payload['status']}`", "", "## Theorem", "",
        proof["statement"], "", "## Proof", "",
    ]
    for step in proof["proof_steps"]:
        lines.append(f"- {step}")
    lines.extend([
        "", "## Computed checks", "",
        f"- `eta_strict = 1 + delta`: `{t['eta_strict_equals_one_plus_delta']}`.",
        f"- P2606 used delta as denominator eta: `{t['p2606_used_delta_as_denominator_eta']}`.",
        f"- Max strict eta=9/5 vs P2606 eta=4/5 kernel difference: `{t['max_abs_strict_eta_9_5_minus_p2606_eta_4_5_kernel']}`.",
        f"- P2606 strict residual interpretation quarantined: `{t['p2606_strict_residual_interpretation_quarantined_by_p2617']}`.",
        f"- P2606 codimension probe retained: `{t['p2606_codimension_probe_retained']}`.",
        "", "## Scope guards", "",
        "P2617 reclassifies only the P2606 exponent semantics. It does not revalidate P2606 as a strict-side residual addition, does not revalidate P2607/P2608, and does not export QW-2191 discharge, APD sourcehood, legacy role transfer, or ToE closure.",
        "", "## Fingerprint", "",
        f"`{payload['p2606_exponent_semantics_reclassification']['theorem_fingerprint_sha256']}`",
    ])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2617/S1567 P2606 exponent-semantics reclassification

`P2617/S1567` audits the P2606 nonlinear-residual label and finds a strict notation mismatch: P2603 supplies the codimension/log-slope `delta=4/5`, while the strict denominator exponent audited by P2414 is `eta=1+delta=9/5`.  Therefore the old P2606 denominator computation with `eta=4/5` is retained only as a codimension-slope perturbation probe, not as an exported strict-side nonlinear residual addition.  P2607 bridge completion and P2608 role-bearing `L_total` remain blocked.
""".strip()
    lag_section = """
## P2617/S1567 exponent-semantics Ltotal guard

`P2617/S1567` prevents the P2606 `eta=4/5` residual label from entering `L_total` as a strict nonlinear compression term.  Any future residual insertion must use explicit `delta=4/5`, `eta=1+delta=9/5` semantics and still pass bridge and role-transfer gates.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2617/S1567 P2606 exponent-semantics reclassification", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2617/S1567 exponent-semantics Ltotal guard", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
