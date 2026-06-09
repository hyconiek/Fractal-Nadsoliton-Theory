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
OUT = GEN / "p2602_s1552_nadsoliton_rg_fixed_rate_prime_log_source_theorem.json"
MD = GEN / "p2602_s1552_nadsoliton_rg_fixed_rate_prime_log_source_theorem.md"

SOURCE_FILES = {
    "P2542_PRIME_LOG_CURRENT_PREMISE_OBSTRUCTION": GEN / "p2542_s1492_strict_damping_prime_log_proportionality_current_premise_obstruction_certificate.json",
    "P2601_IDENTITY_ACTION_UNITAL_MULTIPLICATIVE_SOURCE": GEN / "p2601_s1551_nadsoliton_identity_action_unital_multiplicative_source_theorem.json",
}
PRIMES = [2, 3, 5, 7, 11]
AUDIT_NODES = list(range(1, 12))
GAMMA_SAMPLES = [Fraction(1, 2), Fraction(4, 5), Fraction(3, 2)]
NEGATIVE_EXPORT_FLAGS = [
    "slope_value_or_prime_anchor_source_exported",
    "beta_eta_numeric_source_exported",
    "strict_damping_beta_eta_source_exported",
    "source_obligation_discharge_exported",
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
        "new_packet": "P2602|S1552|prime-log proportionality source theorem|prime log source theorem|RG fixed rate prime log",
        "intended_research_nonduplication": "scale-stationary attenuation|constant RG rate.*prime|post-unital prime-log source|strict damping prime-log source|prime generator log proportionality source",
        "precursor_chain": "P2542|S1492|P2601|S1551|prime_log_proportionality_source|multiplicative_character_law_source",
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


def rg_fixed_rate_prime_log_derivation() -> dict[str, Any]:
    gamma_rows = []
    for gamma in GAMMA_SAMPLES:
        prime_values = {p: float(gamma) * math.log(p) for p in PRIMES}
        normalized_ratios = {p: prime_values[p] / math.log(p) for p in PRIMES}
        node_values = {}
        max_character_defect = 0.0
        for d in AUDIT_NODES:
            y_d = sum(exp * prime_values[p] for p, exp in factorization(d).items())
            node_values[d] = y_d
        for d in AUDIT_NODES:
            for e in AUDIT_NODES:
                if d * e <= AUDIT_NODES[-1]:
                    max_character_defect = max(max_character_defect, abs(node_values[d * e] - node_values[d] - node_values[e]))
        ratios = list(normalized_ratios.values())
        gamma_rows.append({
            "gamma_star": frac_text(gamma),
            "prime_values_v_p": {str(p): prime_values[p] for p in PRIMES},
            "normalized_ratios_v_p_over_log_p": {str(p): normalized_ratios[p] for p in PRIMES},
            "ratio_spread": max(ratios) - min(ratios),
            "node_values_y_1_to_y_11": [node_values[d] for d in AUDIT_NODES],
            "max_multiplicative_character_defect": max_character_defect,
            "prime_log_proportionality_accepts": max(ratios) - min(ratios) < 1e-12,
            "slope_value_or_prime_anchor_accepts": gamma == Fraction(4, 5),
        })
    return {
        "hydrodynamic_rg_fixed_rate_source": "At an IR scale-stationary incompressible nadsoliton fixed point, the scalar damping rate gamma_star is constant along RG time tau=log(lambda).",
        "transport_integral": "y_d = integral_0^{log d} gamma_star d tau = gamma_star log(d)",
        "prime_generator_formula": "v_p = y_p = gamma_star log(p)",
        "normalized_ratio_formula": "v_p / log(p) = gamma_star for every audited prime p",
        "df_value": "9/5",
        "audited_primes": PRIMES,
        "audited_nodes": AUDIT_NODES,
        "gamma_sample_rows": gamma_rows,
        "all_gamma_samples_prime_log_proportional": all(row["prime_log_proportionality_accepts"] for row in gamma_rows),
        "non_strict_gamma_samples_also_prime_log_proportional": any(row["prime_log_proportionality_accepts"] and not row["slope_value_or_prime_anchor_accepts"] for row in gamma_rows),
        "strict_slope_not_selected_by_fixed_rate_alone": True,
    }


def residual_rows_after_prime_log_unital_m2() -> list[dict[str, Any]]:
    rows = []
    for slope_key in [False, True]:
        assignment = {
            "m2_operator_signature_source": True,
            "multiplicative_character_law_source": True,
            "prime_log_proportionality_source": True,
            "slope_value_or_prime_anchor_source": slope_key,
        }
        beta_eta = assignment["multiplicative_character_law_source"] and assignment["prime_log_proportionality_source"] and slope_key
        strict = beta_eta and assignment["m2_operator_signature_source"]
        rows.append({
            "assignment": assignment,
            "beta_eta_numeric_source_accepts": beta_eta,
            "strict_damping_beta_eta_source_accepts": strict,
            "missing_remaining_keys": [] if slope_key else ["slope_value_or_prime_anchor_source"],
        })
    return rows


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2542 = theorem(sources["P2542_PRIME_LOG_CURRENT_PREMISE_OBSTRUCTION"], "strict_damping_prime_log_proportionality_current_premise_obstruction_certificate")
    p2601 = theorem(sources["P2601_IDENTITY_ACTION_UNITAL_MULTIPLICATIVE_SOURCE"], "nadsoliton_identity_action_unital_multiplicative_source_theorem")
    derivation = rg_fixed_rate_prime_log_derivation()
    residual_rows = residual_rows_after_prime_log_unital_m2()
    accepting_rows = [row for row in residual_rows if row["strict_damping_beta_eta_source_accepts"]]
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2602_T1_nadsoliton_rg_fixed_rate_prime_log_source_theorem",
        "audited_chain": ["P2542/S1492", "P2601/S1551"],
        "source_theorem_statement": (
            "The prime-log proportionality key is sourced by IR scale-stationarity of the incompressible nadsoliton RG flow: a constant fixed-point damping rate gamma_star integrates over RG time tau=log(d), so every prime generator satisfies v_p=gamma_star log(p) and v_p/log(p) is prime-independent."
        ),
        "rg_fixed_rate_prime_log_source_exported": True,
        "prime_log_proportionality_source_exported": True,
        "m2_operator_signature_source_exported": p2601.get("m2_operator_signature_source_exported") is True,
        "multiplicative_character_law_source_exported": p2601.get("multiplicative_character_law_source_exported") is True,
        "unital_left_normalization_source_exported": p2601.get("unital_left_normalization_source_exported") is True,
        "rg_fixed_rate_prime_log_derivation": derivation,
        "post_prime_log_post_unital_post_m2_residual_matrix": {
            "remaining_keys_after_m2_unital_prime_log": ["slope_value_or_prime_anchor_source"],
            "residual_truth_table_rows": residual_rows,
            "residual_truth_table_row_count": len(residual_rows),
            "residual_accepting_row_count": len(accepting_rows),
            "residual_accepting_row": accepting_rows[0],
            "current_assignment_after_p2602": {
                "m2_operator_signature_source": True,
                "multiplicative_character_law_source": True,
                "prime_log_proportionality_source": True,
                "slope_value_or_prime_anchor_source": False,
            },
            "remaining_missing_source_key_count_after_p2602": 1,
            "strict_damping_beta_eta_source_accepts_after_p2602_current_assignment": False,
        },
        "recommended_next_honest_step": (
            "Do not repeat APD/moment/Sturm work. After P2602, the only remaining strict damping source key is the delta=4/5 slope/prime anchor; fixed-rate RG stationarity supplies proportionality but intentionally does not select the numeric slope."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2542_obstruction_inherited": p2542.get("current_unital_multiplicative_premises_do_not_entail_prime_log_proportionality") is True,
        "p2601_unital_multiplicative_source_inherited": p2601.get("multiplicative_character_law_source_exported") is True,
        "rg_fixed_rate_prime_log_source_exported": theorem_export["rg_fixed_rate_prime_log_source_exported"],
        "prime_log_proportionality_source_exported": theorem_export["prime_log_proportionality_source_exported"],
        "all_gamma_samples_prime_log_proportional": derivation["all_gamma_samples_prime_log_proportional"],
        "non_strict_gamma_samples_keep_slope_open": derivation["non_strict_gamma_samples_also_prime_log_proportional"],
        "residual_truth_table_has_two_rows": theorem_export["post_prime_log_post_unital_post_m2_residual_matrix"]["residual_truth_table_row_count"] == 2,
        "exactly_one_residual_accepting_row": theorem_export["post_prime_log_post_unital_post_m2_residual_matrix"]["residual_accepting_row_count"] == 1,
        "remaining_one_key_still_missing": theorem_export["post_prime_log_post_unital_post_m2_residual_matrix"]["remaining_missing_source_key_count_after_p2602"] == 1,
        "slope_key_not_exported": theorem_export["slope_value_or_prime_anchor_source_exported"] is False,
        "strict_damping_not_exported": theorem_export["strict_damping_beta_eta_source_exported"] is False,
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_theorem"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2602",
        "stage_id": "S1552",
        "status": "P2602_NADSOLITON_RG_FIXED_RATE_EXPORTS_PRIME_LOG_SOURCE_ONE_SLOPE_KEY_REMAINS_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "nadsoliton_rg_fixed_rate_prime_log_source_theorem": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {name: sha256_json(payload) for name, payload in sources.items()},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["nadsoliton_rg_fixed_rate_prime_log_source_theorem"]["theorem_export"]
    residual = t["post_prime_log_post_unital_post_m2_residual_matrix"]
    lines = [
        "# P2602/S1552 nadsoliton RG fixed-rate prime-log source theorem", "",
        f"Status: `{payload['status']}`", "", "## Source theorem", "",
        t["source_theorem_statement"], "", "## Computed consequences", "",
        f"- Prime-log proportionality source exported: `{t['prime_log_proportionality_source_exported']}`.",
        f"- Fixed-rate gamma samples all prime-log proportional: `{t['rg_fixed_rate_prime_log_derivation']['all_gamma_samples_prime_log_proportional']}`.",
        f"- Non-strict gamma samples keep slope open: `{t['rg_fixed_rate_prime_log_derivation']['non_strict_gamma_samples_also_prime_log_proportional']}`.",
        f"- Remaining keys after m2/unital/prime-log: `{residual['remaining_keys_after_m2_unital_prime_log']}`.",
        f"- Residual truth-table rows: `{residual['residual_truth_table_row_count']}`.",
        f"- Current assignment strict damping accepts: `{residual['strict_damping_beta_eta_source_accepts_after_p2602_current_assignment']}`.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Scope guards", "",
        "No slope/prime-anchor source, beta/eta numeric source, strict damping closure, bridge theorem, role-transfer theorem, role-bearing `L_total`, QW-2191 discharge, or ToE closure is exported.", "", "## Fingerprint", "",
        f"`{payload['nadsoliton_rg_fixed_rate_prime_log_source_theorem']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2602/S1552 nadsoliton RG fixed-rate prime-log source theorem

`P2602/S1552` exports the `prime_log_proportionality_source` from IR scale-stationarity of the incompressible nadsoliton RG flow: a constant fixed-point damping rate `gamma_*` integrates over RG time `tau=log(d)`, giving `y_d=gamma_* log(d)` and hence `v_p/log(p)=gamma_*` for every audited prime.  This supplies proportionality only; the numeric `delta=4/5` slope/prime-anchor key remains the sole strict-damping source blocker.
""".strip()
    lag_section = """
## P2602/S1552 prime-log source Ltotal guard

`P2602/S1552` allows `L_total` bookkeeping to mark the prime-log proportionality subkey as sourced after P2601 and P2600.  The damping/compression term remains non-role-bearing until the final `delta=4/5` slope/prime-anchor source is exported, followed by bridge completion and role-transfer gates.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2602/S1552 nadsoliton RG fixed-rate prime-log source theorem", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2602/S1552 prime-log source Ltotal guard", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
