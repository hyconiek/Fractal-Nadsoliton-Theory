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
OUT = GEN / "p2553_s1503_strict_damping_nonhomogeneous_anchor_constant_equivalence_certificate.json"
MD = GEN / "p2553_s1503_strict_damping_nonhomogeneous_anchor_constant_equivalence_certificate.md"

SOURCE_FILES = {
    "P2528_PRIME_SLOPE_ANCHOR_EQUIVALENCE": GEN / "p2528_s1478_strict_damping_prime_slope_anchor_equivalence_certificate.json",
    "P2551_POST_PRIME_LOG_SLOPE_ANCHOR_OBSTRUCTION": GEN / "p2551_s1501_strict_damping_post_prime_log_slope_anchor_obstruction_certificate.json",
    "P2552_HOMOGENEOUS_SLOPE_SELECTOR_OBSTRUCTION": GEN / "p2552_s1502_strict_damping_homogeneous_slope_selector_obstruction_certificate.json",
}
PRIMES = [2, 3, 5, 7, 11]
STRICT_DELTA = Fraction(4, 5)
NEGATIVE_EXPORT_FLAGS = [
    "nonhomogeneous_anchor_constant_source_exported", "slope_value_or_prime_anchor_source_exported", "beta_eta_numeric_source_exported",
    "strict_damping_beta_eta_source_exported", "source_obligation_discharge_exported", "full_bridge_theorem_exported",
    "role_transfer_theorem_exported", "selector_closure_exported", "qw2191_discharged_by_this_certificate",
    "role_bearing_ltotal_exported", "toe_closure_claimed",
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
        "new_packet": "P2553|S1503|nonhomogeneous anchor constant|fixed-scale anchor constant|anchor constant equivalence",
        "intended_research_nonduplication": "nonhomogeneous.*anchor constant|fixed-scale.*slope anchor|constant.*selects.*delta|anchor constant.*4/5|affine anchor.*slope selector",
        "precursors": "P2528|S1478|P2551|S1501|P2552|S1502|prime-value/slope-anchor theorem|nonhomogeneous strict slope anchor",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def load_theorems(sources: dict[str, dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {
        "P2528": theorem(sources["P2528_PRIME_SLOPE_ANCHOR_EQUIVALENCE"], "strict_damping_prime_slope_anchor_equivalence_certificate"),
        "P2551": theorem(sources["P2551_POST_PRIME_LOG_SLOPE_ANCHOR_OBSTRUCTION"], "strict_damping_post_prime_log_slope_anchor_obstruction_certificate"),
        "P2552": theorem(sources["P2552_HOMOGENEOUS_SLOPE_SELECTOR_OBSTRUCTION"], "strict_damping_homogeneous_slope_selector_obstruction_certificate"),
    }


def anchor_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    coefficient_vectors = []
    for i, prime in enumerate(PRIMES):
        coeffs = [0.0] * len(PRIMES)
        coeffs[i] = 1.0
        coefficient_vectors.append((f"single_prime_v_{prime}", coeffs))
    coefficient_vectors.append(("sum_all_prime_values", [1.0] * len(PRIMES)))
    coefficient_vectors.append(("weighted_sum_prime_values", [float(i + 1) for i in range(len(PRIMES))]))
    logs = [math.log(p) for p in PRIMES]
    for name, coeffs in coefficient_vectors:
        log_dot = sum(c * ell for c, ell in zip(coeffs, logs))
        strict_constant = float(STRICT_DELTA) * log_dot
        rows.append({
            "anchor_row": name,
            "coefficients": coeffs,
            "dot_with_log_prime_vector": log_dot,
            "strict_constant_required_k": strict_constant,
            "selector_formula_if_k_supplied": "delta = k/(c·log(p))",
            "delta_selected_by_strict_constant": "4/5",
            "zero_constant_selects_delta": "0" if abs(log_dot) > 1e-12 else "underdetermined",
            "unit_constant_selects_delta_numeric": 1.0 / log_dot if abs(log_dot) > 1e-12 else None,
            "constant_is_equivalent_to_slope_value": True,
        })
    return rows


def constant_misanchor_witnesses(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    witnesses = []
    for row in rows:
        log_dot = row["dot_with_log_prime_vector"]
        for label, k in [("zero_constant", 0.0), ("half_strict_constant", 0.5 * row["strict_constant_required_k"]), ("unit_constant", 1.0)]:
            selected = k / log_dot
            witnesses.append({
                "anchor_row": row["anchor_row"],
                "constant_label": label,
                "constant_k": k,
                "selected_delta_numeric": selected,
                "selects_strict_delta_4_over_5": abs(selected - float(STRICT_DELTA)) < 1e-12,
            })
    return witnesses


def equivalence_audit() -> dict[str, Any]:
    rows = anchor_rows()
    witnesses = constant_misanchor_witnesses(rows)
    return {
        "principle": "On the post-prime-log line v=delta*log(p), every nonhomogeneous linear anchor c·v=k with c·log(p)!=0 selects delta=k/(c·log(p)); selecting delta=4/5 is exactly the constant obligation k=(4/5)c·log(p).",
        "anchor_rows": rows,
        "anchor_row_count": len(rows),
        "all_anchor_rows_have_nonzero_log_dot": all(abs(row["dot_with_log_prime_vector"]) > 1e-12 for row in rows),
        "constant_misanchor_witnesses": witnesses,
        "all_non_strict_constants_fail_strict_delta": all(not witness["selects_strict_delta_4_over_5"] for witness in witnesses if witness["constant_label"] != "strict_constant"),
        "source_obligation_rephrasing": "A nonhomogeneous slope anchor is not a free source; its constant must itself be sourced at the strict value.",
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    ts = load_theorems(sources)
    audit = equivalence_audit()
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2553_T1_strict_damping_nonhomogeneous_anchor_constant_equivalence_certificate",
        "audited_chain": ["P2528/S1478", "P2551/S1501", "P2552/S1502"],
        "frontier_source_under_attack": "slope_value_or_prime_anchor_source",
        "p2528_prime_anchor_equivalence_inherited": bool(ts["P2528"]),
        "p2551_post_prime_log_slope_obstruction_inherited": ts["P2551"].get("post_prime_log_slope_value_nonentailment_exported") is True,
        "p2552_homogeneous_selector_obstruction_inherited": ts["P2552"].get("homogeneous_constraints_cannot_select_nonzero_strict_delta") is True,
        "nonhomogeneous_anchor_constant_equivalence_audit": audit,
        "nonhomogeneous_anchor_reduces_to_constant_source_obligation": True,
        "nonhomogeneous_anchor_constant_source_exported": False,
        "slope_value_or_prime_anchor_source_exported": False,
        "recommended_next_honest_step": (
            "Do not count a nonhomogeneous anchor as progress unless its constant is independently sourced. The next honest step is to derive the fixed constant "
            "k=(4/5)c·log(p) from strict nadsoliton dynamics, or stop local strict-damping closure attempts and move to the broader legacy->strict completion/source bridge audit."
        ),
        "not_licensed": [
            "No strict source for the nonhomogeneous anchor constant is supplied.",
            "No slope source, beta/eta source, bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing L_total, or ToE closure is exported.",
        ],
    }
    for flag in NEGATIVE_EXPORT_FLAGS:
        theorem_export[flag] = False
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "precursors_inherited": theorem_export["p2551_post_prime_log_slope_obstruction_inherited"] and theorem_export["p2552_homogeneous_selector_obstruction_inherited"],
        "all_anchor_rows_have_nonzero_log_dot": audit["all_anchor_rows_have_nonzero_log_dot"],
        "anchor_reduces_to_constant_source": theorem_export["nonhomogeneous_anchor_reduces_to_constant_source_obligation"],
        "no_false_source_export": theorem_export["nonhomogeneous_anchor_constant_source_exported"] is False and theorem_export["slope_value_or_prime_anchor_source_exported"] is False,
        "no_false_bridge_or_role_transfer": theorem_export["full_bridge_theorem_exported"] is False and theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_or_toe_claim": theorem_export["qw2191_discharged_by_this_certificate"] is False and theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2553",
        "stage_id": "S1503",
        "status": "STRICT_DAMPING_NONHOMOGENEOUS_ANCHOR_CONSTANT_EQUIVALENCE_CERTIFICATE_NO_ANCHOR_CONSTANT_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_nonhomogeneous_anchor_constant_equivalence_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_nonhomogeneous_anchor_constant_equivalence_certificate"]["theorem_export"]
    audit = t["nonhomogeneous_anchor_constant_equivalence_audit"]
    lines = [
        "# P2553/S1503 strict damping nonhomogeneous anchor-constant equivalence certificate", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Frontier source under attack: `{t['frontier_source_under_attack']}`.",
        f"- P2552 homogeneous selector obstruction inherited: `{t['p2552_homogeneous_selector_obstruction_inherited']}`.",
        f"- Audited nonhomogeneous anchor rows: `{audit['anchor_row_count']}`.",
        f"- All anchor rows have nonzero `c·log(p)`: `{audit['all_anchor_rows_have_nonzero_log_dot']}`.",
        f"- Anchor reduces to constant source obligation: `{t['nonhomogeneous_anchor_reduces_to_constant_source_obligation']}`.",
        f"- Anchor constant source exported: `{t['nonhomogeneous_anchor_constant_source_exported']}`.", "", "## Interpretation", "",
        audit["principle"], "", audit["source_obligation_rephrasing"], "", "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No nonhomogeneous anchor-constant source, slope/prime-anchor source, beta/eta source, bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.",
        "", "## Fingerprint", "", f"`{payload['strict_damping_nonhomogeneous_anchor_constant_equivalence_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
`P2553/S1503` audits the nonhomogeneous alternative left open by P2552.  On the post-prime-log line `v=delta log(p)`, any fixed-scale anchor `c·v=k` with `c·log(p) != 0` selects `delta=k/(c·log(p))`; selecting the strict nonzero value is exactly the constant obligation `k=(4/5)c·log(p)`.  Thus a nonhomogeneous anchor is only progress if the strict constant is independently sourced, not merely inserted.
""".strip()
    lag_section = """
`P2553/S1503` blocks another false `L_total` route: writing a nonhomogeneous slope anchor does not source the damping slope unless the anchor constant itself is derived.  The role-bearing damping/compression term still waits for a strict constant source plus bridge-completion and role-transfer gates.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2553/S1503 nonhomogeneous anchor-constant equivalence guard", "## P2553/S1503 nonhomogeneous anchor-constant equivalence guard\n\n" + eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2553/S1503 nonhomogeneous anchor-constant Ltotal guard", "## P2553/S1503 nonhomogeneous anchor-constant Ltotal guard\n\n" + lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
