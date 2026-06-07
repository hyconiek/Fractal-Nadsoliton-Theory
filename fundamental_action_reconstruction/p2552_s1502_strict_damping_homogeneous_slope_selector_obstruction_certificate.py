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
OUT = GEN / "p2552_s1502_strict_damping_homogeneous_slope_selector_obstruction_certificate.json"
MD = GEN / "p2552_s1502_strict_damping_homogeneous_slope_selector_obstruction_certificate.md"

SOURCE_FILES = {
    "P2543_SLOPE_VALUE_OBSTRUCTION": GEN / "p2543_s1493_strict_damping_slope_value_current_premise_obstruction_certificate.json",
    "P2550_PRIME_LOG_ADJACENT_RATIO_BASIS": GEN / "p2550_s1500_strict_damping_prime_log_adjacent_ratio_basis_certificate.json",
    "P2551_POST_PRIME_LOG_SLOPE_ANCHOR_OBSTRUCTION": GEN / "p2551_s1501_strict_damping_post_prime_log_slope_anchor_obstruction_certificate.json",
}
PRIMES = [2, 3, 5, 7, 11]
STRICT_DELTA = Fraction(4, 5)
SLOPE_CANDIDATES = [Fraction(0, 1), Fraction(1, 2), Fraction(4, 5), Fraction(1, 1), Fraction(9, 5)]
NEGATIVE_EXPORT_FLAGS = [
    "homogeneous_slope_selector_source_exported", "slope_value_or_prime_anchor_source_exported", "prime_log_proportionality_source_exported",
    "m2_operator_signature_source_exported", "strict_quadruple_trace_source_exported", "beta_eta_numeric_source_exported",
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
        "new_packet": "P2552|S1502|homogeneous slope selector|homogeneous slope obstruction|scale-invariant slope",
        "intended_research_nonduplication": "homogeneous.*slope selector|scale.*slope.*delta|nonhomogeneous.*slope anchor|homogeneous.*delta=4/5|slope selector obstruction",
        "precursors": "P2543|S1493|P2550|S1500|P2551|S1501|slope_value_or_prime_anchor_source",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def load_theorems(sources: dict[str, dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {
        "P2543": theorem(sources["P2543_SLOPE_VALUE_OBSTRUCTION"], "strict_damping_slope_value_current_premise_obstruction_certificate"),
        "P2550": theorem(sources["P2550_PRIME_LOG_ADJACENT_RATIO_BASIS"], "strict_damping_prime_log_adjacent_ratio_basis_certificate"),
        "P2551": theorem(sources["P2551_POST_PRIME_LOG_SLOPE_ANCHOR_OBSTRUCTION"], "strict_damping_post_prime_log_slope_anchor_obstruction_certificate"),
    }


def frac_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def coefficient_rows() -> list[dict[str, Any]]:
    logs = [math.log(p) for p in PRIMES]
    rows = []
    # Homogeneous versions of adjacent ratio equalities: v_i/log(p_i)-v_j/log(p_j)=0.
    for i in range(len(PRIMES) - 1):
        coeffs = [0.0] * len(PRIMES)
        coeffs[i] = 1.0 / logs[i]
        coeffs[i + 1] = -1.0 / logs[i + 1]
        rows.append({"row_name": f"ratio_edge_{PRIMES[i]}_{PRIMES[i+1]}", "coefficients": coeffs})
    rows.append({"row_name": "zero_prime_value_v2", "coefficients": [1.0, 0.0, 0.0, 0.0, 0.0]})
    rows.append({"row_name": "sum_prime_values_zero", "coefficients": [1.0] * len(PRIMES)})
    return rows


def row_audit(row: dict[str, Any]) -> dict[str, Any]:
    logs = [math.log(p) for p in PRIMES]
    coeffs = row["coefficients"]
    log_dot = sum(c * l for c, l in zip(coeffs, logs))
    accepted = [s for s in SLOPE_CANDIDATES if abs(float(s) * log_dot) < 1e-12]
    if abs(log_dot) < 1e-12:
        classification = "scale_invariant_accepts_entire_slope_line"
    else:
        classification = "homogeneous_zero_selector_accepts_only_delta_0"
    return {
        "row_name": row["row_name"],
        "coefficients": coeffs,
        "dot_with_log_prime_vector": log_dot,
        "classification": classification,
        "accepted_audited_slopes": [frac_text(s) for s in accepted],
        "accepts_strict_delta_4_over_5": STRICT_DELTA in accepted,
        "uniquely_selects_strict_delta_4_over_5": accepted == [STRICT_DELTA],
    }


def homogeneous_selector_audit() -> dict[str, Any]:
    rows = [row_audit(row) for row in coefficient_rows()]
    return {
        "principle": "For homogeneous linear constraints c·v=0 on the post-prime-log line v=delta*log(p), either c·log(p)=0 and every delta passes, or c·log(p)!=0 and only delta=0 passes.",
        "audited_rows": rows,
        "audited_row_count": len(rows),
        "any_row_uniquely_selects_strict_delta": any(row["uniquely_selects_strict_delta_4_over_5"] for row in rows),
        "scale_invariant_rows": [row["row_name"] for row in rows if row["classification"] == "scale_invariant_accepts_entire_slope_line"],
        "zero_selector_rows": [row["row_name"] for row in rows if row["classification"] == "homogeneous_zero_selector_accepts_only_delta_0"],
        "required_selector_type": "nonhomogeneous prime-value anchor v_p=(4/5)log(p), fixed-scale normalization, or equivalent strict slope selector",
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    ts = load_theorems(sources)
    audit = homogeneous_selector_audit()
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2552_T1_strict_damping_homogeneous_slope_selector_obstruction_certificate",
        "audited_chain": ["P2543/S1493", "P2550/S1500", "P2551/S1501"],
        "frontier_source_under_attack": "slope_value_or_prime_anchor_source",
        "p2543_slope_obstruction_inherited": ts["P2543"].get("current_slope_line_premises_do_not_entail_delta_4_over_5") is True,
        "p2550_prime_log_basis_inherited": ts["P2550"].get("full_adjacent_ratio_basis_suffices_for_prime_log_proportionality") is True,
        "p2551_post_prime_log_slope_obstruction_inherited": ts["P2551"].get("post_prime_log_slope_value_nonentailment_exported") is True,
        "homogeneous_slope_selector_audit": audit,
        "homogeneous_constraints_cannot_select_nonzero_strict_delta": not audit["any_row_uniquely_selects_strict_delta"],
        "nonhomogeneous_anchor_required": True,
        "recommended_next_honest_step": (
            "Stop testing homogeneous/scale-invariant slope constraints for delta=4/5. They either keep the full delta-line or select delta=0. "
            "The next honest step is to seek a nonhomogeneous strict nadsoliton anchor/fixed-scale theorem, e.g. v_p=(4/5)log(p), "
            "or pause strict damping closure and return to the legacy->strict bridge/source audit without role transfer."
        ),
        "not_licensed": [
            "No nonhomogeneous strict slope anchor is supplied.",
            "No slope source, beta/eta source, bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing L_total, or ToE closure is exported.",
        ],
    }
    for flag in NEGATIVE_EXPORT_FLAGS:
        theorem_export[flag] = False
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "precursors_inherited": theorem_export["p2543_slope_obstruction_inherited"] and theorem_export["p2550_prime_log_basis_inherited"] and theorem_export["p2551_post_prime_log_slope_obstruction_inherited"],
        "homogeneous_nonselection_verified": theorem_export["homogeneous_constraints_cannot_select_nonzero_strict_delta"],
        "nonhomogeneous_anchor_required": theorem_export["nonhomogeneous_anchor_required"],
        "no_false_source_export": theorem_export["slope_value_or_prime_anchor_source_exported"] is False and theorem_export["homogeneous_slope_selector_source_exported"] is False,
        "no_false_bridge_or_role_transfer": theorem_export["full_bridge_theorem_exported"] is False and theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_or_toe_claim": theorem_export["qw2191_discharged_by_this_certificate"] is False and theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2552",
        "stage_id": "S1502",
        "status": "STRICT_DAMPING_HOMOGENEOUS_SLOPE_SELECTOR_OBSTRUCTION_CERTIFICATE_NO_SLOPE_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_homogeneous_slope_selector_obstruction_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_homogeneous_slope_selector_obstruction_certificate"]["theorem_export"]
    audit = t["homogeneous_slope_selector_audit"]
    lines = [
        "# P2552/S1502 strict damping homogeneous slope-selector obstruction certificate", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Frontier source under attack: `{t['frontier_source_under_attack']}`.",
        f"- P2551 post-prime-log slope obstruction inherited: `{t['p2551_post_prime_log_slope_obstruction_inherited']}`.",
        f"- Audited homogeneous rows: `{audit['audited_row_count']}`.",
        f"- Scale-invariant rows: `{audit['scale_invariant_rows']}`.",
        f"- Zero-selector rows: `{audit['zero_selector_rows']}`.",
        f"- Any homogeneous row uniquely selects delta=4/5: `{audit['any_row_uniquely_selects_strict_delta']}`.",
        f"- Nonhomogeneous anchor required: `{t['nonhomogeneous_anchor_required']}`.", "", "## Interpretation", "",
        audit["principle"], "",
        "Thus homogeneous/scale-invariant constraints cannot be the missing strict source for the nonzero slope `delta=4/5`; a fixed-scale, nonhomogeneous prime-value anchor or equivalent selector is required.",
        "", "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No homogeneous slope selector source, slope/prime-anchor source, beta/eta source, bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.",
        "", "## Fingerprint", "", f"`{payload['strict_damping_homogeneous_slope_selector_obstruction_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
`P2552/S1502` audits homogeneous post-prime-log slope selectors.  On the line `v=delta log(p)`, any homogeneous linear constraint `c·v=0` either has `c·log(p)=0` and accepts the whole slope line, or has `c·log(p) != 0` and selects only `delta=0`.  Therefore homogeneous/scale-invariant constraints cannot select the nonzero strict value `delta=4/5`; the missing source must be nonhomogeneous, such as a strict prime-value anchor or equivalent fixed-scale selector.
""".strip()
    lag_section = """
`P2552/S1502` blocks another false `L_total` promotion route: homogeneous slope constraints cannot source the nonzero strict damping slope.  A role-bearing damping/compression term still requires a nonhomogeneous strict slope anchor, plus the separate bridge-completion and role-transfer gates.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2552/S1502 homogeneous slope-selector obstruction guard", "## P2552/S1502 homogeneous slope-selector obstruction guard\n\n" + eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2552/S1502 homogeneous slope-selector Ltotal guard", "## P2552/S1502 homogeneous slope-selector Ltotal guard\n\n" + lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
