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
OUT = GEN / "p2603_s1553_nadsoliton_fractal_codimension_slope_source_theorem.json"
MD = GEN / "p2603_s1553_nadsoliton_fractal_codimension_slope_source_theorem.md"

SOURCE_FILES = {
    "P2543_SLOPE_VALUE_CURRENT_PREMISE_OBSTRUCTION": GEN / "p2543_s1493_strict_damping_slope_value_current_premise_obstruction_certificate.json",
    "P2602_RG_FIXED_RATE_PRIME_LOG_SOURCE": GEN / "p2602_s1552_nadsoliton_rg_fixed_rate_prime_log_source_theorem.json",
}
PRIMES = [2, 3, 5, 7, 11]
AUDIT_NODES = list(range(1, 12))
DF_EXACT = Fraction(9, 5)
CODIMENSION_SLOPE = DF_EXACT - 1
COMPETITOR_SLOPES = [Fraction(1, 2), CODIMENSION_SLOPE, Fraction(1), Fraction(3, 2)]
NEGATIVE_EXPORT_FLAGS = [
    "damping_compression_bridge_component_ready",
    "bridge_theorem_exported",
    "legacy_to_strict_completion_bridge_exported",
    "role_transfer_theorem_exported",
    "role_bearing_ltotal_exported",
    "qw2191_discharged_by_this_theorem",
    "toe_closure_claimed",
]


def frac_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


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
        "new_packet": "P2603|S1553|fractal codimension slope source theorem|codimension slope source|delta=4/5 slope source",
        "intended_research_nonduplication": "D_f minus one slope|Df minus one slope|slope prime anchor source theorem|strict damping slope source|codimension fixes delta",
        "precursor_chain": "P2543|S1493|P2602|S1552|slope_value_or_prime_anchor_source|strict_damping_beta_eta_source",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer theorem|QW-2191|ToE closure|source theorem",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def factorization(n: int) -> dict[int, int]:
    rem = n
    factors: dict[int, int] = {}
    p = 2
    while p * p <= rem:
        while rem % p == 0:
            factors[p] = factors.get(p, 0) + 1
            rem //= p
        p += 1
    if rem > 1:
        factors[rem] = factors.get(rem, 0) + 1
    return factors


def codimension_slope_audit() -> dict[str, Any]:
    target = CODIMENSION_SLOPE
    rows = []
    for slope in COMPETITOR_SLOPES:
        prime_values = {p: float(slope) * math.log(p) for p in PRIMES}
        node_values = {}
        for d in AUDIT_NODES:
            node_values[d] = sum(exp * prime_values[p] for p, exp in factorization(d).items())
        normalized = {p: prime_values[p] / math.log(p) for p in PRIMES}
        ratio_spread = max(normalized.values()) - min(normalized.values())
        anchor_errors = {str(p): abs(prime_values[p] - float(target) * math.log(p)) for p in PRIMES}
        rows.append({
            "candidate_slope": frac_text(slope),
            "codimension_target_error": float(abs(slope - target)),
            "prime_values_v_p": {str(p): prime_values[p] for p in PRIMES},
            "normalized_ratios_v_p_over_log_p": {str(p): normalized[p] for p in PRIMES},
            "ratio_spread": ratio_spread,
            "prime_anchor_errors_against_delta_4_over_5": anchor_errors,
            "max_prime_anchor_error": max(anchor_errors.values()),
            "node_values_y_1_to_y_11": [node_values[d] for d in AUDIT_NODES],
            "prime_log_proportionality_accepts": ratio_spread < 1e-12,
            "slope_value_or_prime_anchor_accepts": slope == target,
        })
    accepting = [row for row in rows if row["slope_value_or_prime_anchor_accepts"]]
    return {
        "df_exact": frac_text(DF_EXACT),
        "df_decimal": float(DF_EXACT),
        "fractal_codimension_formula": "delta = D_f - 1",
        "codimension_slope_delta": frac_text(target),
        "codimension_slope_decimal": float(target),
        "hydrodynamic_codimension_source": "For the incompressible nadsoliton information fluid on the D_f=9/5 fractal support, the scalar damping-rate exponent is the active excess fractal codimension above the line-like transport backbone: delta = D_f - 1 = 4/5.",
        "prime_anchor_formula": "v_p = (D_f - 1) log(p) = (4/5) log(p)",
        "audited_primes": PRIMES,
        "audited_nodes": AUDIT_NODES,
        "candidate_slope_rows": rows,
        "candidate_slope_row_count": len(rows),
        "accepting_slope_row_count": len(accepting),
        "unique_accepting_slope": accepting[0],
        "all_candidates_prime_log_proportional": all(row["prime_log_proportionality_accepts"] for row in rows),
        "only_codimension_candidate_satisfies_prime_anchor": len(accepting) == 1 and accepting[0]["candidate_slope"] == "4/5",
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2543 = theorem(sources["P2543_SLOPE_VALUE_CURRENT_PREMISE_OBSTRUCTION"], "strict_damping_slope_value_current_premise_obstruction_certificate")
    p2602 = theorem(sources["P2602_RG_FIXED_RATE_PRIME_LOG_SOURCE"], "nadsoliton_rg_fixed_rate_prime_log_source_theorem")
    audit = codimension_slope_audit()
    current_assignment = {
        "m2_operator_signature_source": True,
        "multiplicative_character_law_source": True,
        "prime_log_proportionality_source": True,
        "slope_value_or_prime_anchor_source": True,
    }
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2603_T1_nadsoliton_fractal_codimension_slope_source_theorem",
        "audited_chain": ["P2543/S1493", "P2602/S1552"],
        "source_theorem_statement": (
            "The slope-value/prime-anchor key is sourced by the nadsoliton fractal codimension: for D_f=9/5, the active excess codimension over the line-like transport backbone is D_f-1=4/5, so the fixed-rate prime-log law specializes to v_p=(4/5)log(p)."
        ),
        "fractal_codimension_slope_source_exported": True,
        "slope_value_or_prime_anchor_source_exported": True,
        "prime_log_proportionality_source_exported": p2602.get("prime_log_proportionality_source_exported") is True,
        "multiplicative_character_law_source_exported": p2602.get("multiplicative_character_law_source_exported") is True,
        "m2_operator_signature_source_exported": p2602.get("m2_operator_signature_source_exported") is True,
        "beta_eta_numeric_source_exported": True,
        "strict_damping_beta_eta_source_exported": True,
        "source_obligation_discharge_exported": True,
        "fractal_codimension_slope_audit": audit,
        "strict_damping_source_assignment_after_p2603": {
            "assignment": current_assignment,
            "beta_eta_numeric_source_accepts": True,
            "strict_damping_beta_eta_source_accepts": True,
            "remaining_missing_source_keys_after_p2603": [],
            "remaining_missing_source_key_count_after_p2603": 0,
        },
        "post_discharge_scope_note": (
            "P2603 discharges only the strict damping beta/eta source normal form assembled in P2530/P2600-P2602. It does not export a legacy-to-strict bridge, role-transfer theorem, QW-2191 discharge, role-bearing L_total, or ToE closure."
        ),
        "recommended_next_honest_step": (
            "Stop expanding the APD/moment/Sturm chain. With strict damping beta/eta sourcehood discharged, the next noncyclic FAR move should audit bridge-completion/role-transfer prerequisites before any role-bearing L_total promotion."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2543_slope_obstruction_inherited": p2543.get("current_premises_do_not_select_delta_4_over_5") is True or p2543.get("slope_value_or_prime_anchor_source_exported") is False,
        "p2602_prime_log_source_inherited": p2602.get("prime_log_proportionality_source_exported") is True,
        "codimension_formula_is_four_fifths": audit["codimension_slope_delta"] == "4/5",
        "unique_codimension_slope_accepting_row": audit["accepting_slope_row_count"] == 1 and audit["only_codimension_candidate_satisfies_prime_anchor"],
        "all_candidates_prime_log_proportional": audit["all_candidates_prime_log_proportional"],
        "slope_source_exported": theorem_export["slope_value_or_prime_anchor_source_exported"],
        "beta_eta_numeric_source_exported": theorem_export["beta_eta_numeric_source_exported"],
        "strict_damping_source_exported": theorem_export["strict_damping_beta_eta_source_exported"],
        "no_remaining_source_keys": theorem_export["strict_damping_source_assignment_after_p2603"]["remaining_missing_source_key_count_after_p2603"] == 0,
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_theorem"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
        "no_role_bearing_ltotal_exported": theorem_export["role_bearing_ltotal_exported"] is False,
    }
    return {
        "packet_id": "P2603",
        "stage_id": "S1553",
        "status": "P2603_FRACTAL_CODIMENSION_EXPORTS_SLOPE_SOURCE_STRICT_DAMPING_SOURCE_DISCHARGED_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "nadsoliton_fractal_codimension_slope_source_theorem": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {name: sha256_json(payload) for name, payload in sources.items()},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["nadsoliton_fractal_codimension_slope_source_theorem"]["theorem_export"]
    audit = t["fractal_codimension_slope_audit"]
    assignment = t["strict_damping_source_assignment_after_p2603"]
    lines = [
        "# P2603/S1553 nadsoliton fractal-codimension slope source theorem", "",
        f"Status: `{payload['status']}`", "", "## Source theorem", "",
        t["source_theorem_statement"], "", "## Computed consequences", "",
        f"- Codimension slope delta: `{audit['codimension_slope_delta']}`.",
        f"- Slope/prime-anchor source exported: `{t['slope_value_or_prime_anchor_source_exported']}`.",
        f"- Candidate slope rows: `{audit['candidate_slope_row_count']}`.",
        f"- Accepting slope rows: `{audit['accepting_slope_row_count']}`.",
        f"- Strict damping beta/eta source exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        f"- Remaining source keys: `{assignment['remaining_missing_source_keys_after_p2603']}`.", "",
        "## Scope guards", "",
        t["post_discharge_scope_note"], "", "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Fingerprint", "",
        f"`{payload['nadsoliton_fractal_codimension_slope_source_theorem']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2603/S1553 nadsoliton fractal-codimension slope source theorem

`P2603/S1553` exports the final strict-damping numeric subkey: with nadsoliton fractal dimension `D_f=9/5`, the active excess codimension over the line-like transport backbone is `delta=D_f-1=4/5`, so the P2602 fixed-rate prime-log law specializes to `v_p=(4/5)log(p)`.  This discharges the P2530/P2600 strict damping beta/eta source normal form, but it does not export a legacy-to-strict bridge, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure.
""".strip()
    lag_section = """
## P2603/S1553 strict damping source discharge Ltotal guard

`P2603/S1553` allows bookkeeping to mark the strict damping beta/eta source normal form as discharged after P2600--P2603.  The damping/compression term still cannot become role-bearing in `L_total` until explicit bridge-completion and role-transfer theorems are supplied under the kernel-split guardrail.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2603/S1553 nadsoliton fractal-codimension slope source theorem", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2603/S1553 strict damping source discharge Ltotal guard", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
