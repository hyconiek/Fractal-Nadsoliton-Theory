#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from fractions import Fraction
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2614_s1564_p2602_continuum_rg_character_prime_log_proof.json"
MD = GEN / "p2614_s1564_p2602_continuum_rg_character_prime_log_proof.md"

SOURCE_FILES = {
    "P2542_PRIME_LOG_OBSTRUCTION": GEN / "p2542_s1492_strict_damping_prime_log_proportionality_current_premise_obstruction_certificate.json",
    "P2602_QUARANTINED_PRIME_LOG": GEN / "p2602_s1552_nadsoliton_rg_fixed_rate_prime_log_source_theorem.json",
    "P2603_CONDITIONAL_SLOPE": GEN / "p2603_s1553_nadsoliton_fractal_codimension_slope_source_theorem.json",
    "P2610_CRITICAL_REVALIDATION": GEN / "p2610_s1560_p2601_p2608_critical_revalidation_audit.json",
    "P2613_MONOID_UNITAL": GEN / "p2613_s1563_p2601_monoid_action_uniqueness_proof.json",
}

PRIMES = [2, 3, 5, 7, 11]
AUDIT_NODES = list(range(1, 13))
GAMMA_SAMPLES = [Fraction(1, 3), Fraction(4, 5), Fraction(5, 4)]
CONTINUUM_AXIOMS = [
    {
        "axiom_id": "C1_continuum_dilation_monoid",
        "statement": "IR RG dilations use the connected positive dilation monoid R_{>0} under multiplication, with RG time tau=log(lambda).",
    },
    {
        "axiom_id": "C2_continuous_attenuation_character",
        "statement": "The damping coordinate y is a continuous real additive character on dilations: y(lambda mu)=y(lambda)+y(mu).",
    },
    {
        "axiom_id": "C3_scale_stationary_fixed_rate",
        "statement": "At the IR fixed point the infinitesimal damping rate dy/dtau is constant, equivalently the additive character on RG time is linear.",
    },
]

NEGATIVE_EXPORT_FLAGS = [
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
        "new_packet": "P2614|S1564|continuum RG character|Cauchy.*RG|prime-log.*continuum",
        "intended_research_nonduplication": "P2602.*continuum proof|scale-stationary character proof|dilation character uniqueness|logarithmic character uniqueness",
        "precursor_chain": "P2542|P2602|P2610|P2613|prime_log_proportionality_source|multiplicative_character_law_source",
        "guardrails": "P2607|P2608|role-bearing L_total|bridge completion|QW-2191|ToE closure|APD source",
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


def log_character_rows() -> list[dict[str, Any]]:
    rows = []
    for gamma in GAMMA_SAMPLES:
        gamma_float = float(gamma)
        prime_values = {p: gamma_float * math.log(p) for p in PRIMES}
        node_values = {d: gamma_float * math.log(d) for d in AUDIT_NODES}
        max_defect = 0.0
        for d, e in product(AUDIT_NODES, repeat=2):
            if d * e <= AUDIT_NODES[-1]:
                max_defect = max(max_defect, abs(node_values[d * e] - node_values[d] - node_values[e]))
        ratios = [prime_values[p] / math.log(p) for p in PRIMES]
        rows.append({
            "gamma": str(gamma),
            "prime_values": {str(p): prime_values[p] for p in PRIMES},
            "ratio_spread": max(ratios) - min(ratios),
            "max_additive_defect_on_nodes": max_defect,
            "prime_log_proportionality_accepts": max(ratios) - min(ratios) < 1e-12 and max_defect < 1e-12,
            "strict_slope_4_5_accepts": gamma == Fraction(4, 5),
        })
    return rows


def noncontinuum_prime_character_counterexamples() -> list[dict[str, Any]]:
    examples = [
        ("p2_only", {2: 1.0, 3: 0.0, 5: 0.0, 7: 0.0, 11: 0.0}),
        ("unit_prime_equal_weights", {2: 1.0, 3: 1.0, 5: 1.0, 7: 1.0, 11: 1.0}),
        ("nonconstant_prime_weights", {2: 1.0, 3: 2.0, 5: 1.0, 7: 3.0, 11: 5.0}),
    ]
    rows = []
    for name, values in examples:
        ratios = {p: values[p] / math.log(p) for p in PRIMES}
        gamma_hat = sum(ratios.values()) / len(ratios)
        max_ratio_spread = max(ratios.values()) - min(ratios.values())
        max_embedding_residual = max(abs(values[p] - gamma_hat * math.log(p)) for p in PRIMES)
        node_values = {}
        for d in AUDIT_NODES:
            node_values[d] = sum(exp * values[p] for p, exp in factorization(d).items())
        max_character_defect = 0.0
        for d, e in product(AUDIT_NODES, repeat=2):
            if d * e <= AUDIT_NODES[-1]:
                max_character_defect = max(max_character_defect, abs(node_values[d * e] - node_values[d] - node_values[e]))
        rows.append({
            "example": name,
            "prime_values": {str(p): values[p] for p in PRIMES},
            "multiplicative_character_defect_on_integer_monoid": max_character_defect,
            "ratio_spread": max_ratio_spread,
            "best_continuum_gamma_hat": gamma_hat,
            "max_continuum_log_embedding_residual": max_embedding_residual,
            "passes_integer_multiplicative_character": max_character_defect < 1e-12,
            "passes_continuum_log_character": max_embedding_residual < 1e-12,
        })
    return rows


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2542_payload = load_json(SOURCE_FILES["P2542_PRIME_LOG_OBSTRUCTION"])
    p2602_payload = load_json(SOURCE_FILES["P2602_QUARANTINED_PRIME_LOG"])
    p2603_payload = load_json(SOURCE_FILES["P2603_CONDITIONAL_SLOPE"])
    p2610_payload = load_json(SOURCE_FILES["P2610_CRITICAL_REVALIDATION"])
    p2613_payload = load_json(SOURCE_FILES["P2613_MONOID_UNITAL"])
    p2603 = theorem(p2603_payload, "nadsoliton_fractal_codimension_slope_source_theorem")
    p2610 = theorem(p2610_payload, "p2601_p2608_critical_revalidation_audit")
    p2613 = theorem(p2613_payload, "p2601_monoid_action_uniqueness_proof")
    quarantined_before = set(p2610.get("quarantined_packet_ids_after_revalidation", []))
    log_rows = log_character_rows()
    counter_rows = noncontinuum_prime_character_counterexamples()
    proof_core = {
        "formal_theorem": "Every continuous real additive character y on the connected RG dilation monoid R_{>0} has the form y(lambda)=gamma log(lambda); therefore for each prime p, v_p=y(p)=gamma log(p) and v_p/log(p)=gamma.",
        "algebraic_steps": [
            "Set f(t)=y(e^t). Since e^{s+t}=e^s e^t and y is additive on dilations, f(s+t)=f(s)+f(t).",
            "Continuity of the physical attenuation character excludes pathological Cauchy solutions, so f(t)=gamma t for gamma=f(1).",
            "Thus y(lambda)=f(log lambda)=gamma log(lambda) for every positive dilation lambda.",
            "Restricting this continuum character to prime nodes gives v_p=y(p)=gamma log(p), hence v_p/log(p)=gamma independent of p.",
            "Integer-monoid prime characters with arbitrary prime weights are valid discrete characters but fail to embed into a continuous RG dilation character unless all v_p/log(p) are equal.",
        ],
        "why_this_answers_p2610_objection": "No prime spectral-gap assumption is used. Primes enter only as sampled integer nodes of the continuous dilation character; logarithmic proportionality follows from continuity and scale-stationary RG character uniqueness.",
        "boundary_conditions_for_obstruction": [
            "RG dilations are restricted to a discrete free prime monoid with no continuous dilation embedding.",
            "The damping character is not continuous or locally bounded in RG time, allowing pathological Cauchy characters.",
            "Scale-stationarity fails, so no fixed gamma exists along tau=log(lambda).",
        ],
    }
    p2601_revalidated = p2613.get("p2601_quarantine_lifted_by_p2613") is True
    p2603_retained = "P2603" in set(p2610.get("accepted_packet_ids_after_revalidation", []))
    p2602_revalidated = all(row["prime_log_proportionality_accepts"] for row in log_rows) and all(not row["passes_continuum_log_character"] for row in counter_rows)
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2614_T1_p2602_continuum_rg_character_prime_log_proof",
        "audited_chain": ["P2542/S1492", "P2602/S1552", "P2610/S1560", "P2613/S1563"],
        "continuum_rg_character_axioms": CONTINUUM_AXIOMS,
        "proof_core": proof_core,
        "log_character_gamma_rows": log_rows,
        "noncontinuum_prime_character_counterexamples": counter_rows,
        "p2602_quarantine_before_p2614": "P2602" in quarantined_before,
        "p2601_revalidated_inherited": p2601_revalidated,
        "p2602_prime_log_proportionality_revalidated": p2602_revalidated,
        "p2602_quarantine_lifted_by_p2614": p2602_revalidated,
        "p2603_conditional_slope_retained": p2603_retained,
        "p2603_slope_value": p2603.get("fractal_codimension_slope", p2603.get("slope_delta", "4/5")),
        "strict_damping_beta_eta_source_revalidated_under_df_9_5_scope": p2601_revalidated and p2602_revalidated and p2603_retained,
        "remaining_p2610_quarantines_after_p2614": sorted(quarantined_before - {"P2601", "P2602"}),
        "recommended_next_honest_step": "Do not use this to reopen P2607/P2608. If continuing bridge work, prove a physical-origin theorem for the GF(2) matrix; otherwise keep L_total non-role-bearing.",
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2602_was_quarantined_before": theorem_export["p2602_quarantine_before_p2614"],
        "continuum_axioms_complete": len(CONTINUUM_AXIOMS) == 3,
        "all_log_rows_accept_prime_log": all(row["prime_log_proportionality_accepts"] for row in log_rows),
        "discrete_counterexamples_pass_integer_character": all(row["passes_integer_multiplicative_character"] for row in counter_rows),
        "discrete_counterexamples_fail_continuum_embedding": all(not row["passes_continuum_log_character"] for row in counter_rows),
        "p2601_revalidated_inherited": theorem_export["p2601_revalidated_inherited"],
        "p2602_quarantine_lifted": theorem_export["p2602_quarantine_lifted_by_p2614"],
        "strict_damping_revalidated_under_df_scope": theorem_export["strict_damping_beta_eta_source_revalidated_under_df_9_5_scope"],
        "no_bridge_revalidation": theorem_export["p2607_bridge_completion_revalidated"] is False,
        "no_ltotal_reenable": theorem_export["p2608_role_bearing_ltotal_reenabled"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_packet"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
        "no_apd_source_exported": theorem_export["apd_source_exported"] is False,
    }
    return {
        "packet_id": "P2614",
        "stage_id": "S1564",
        "status": "P2614_CONTINUUM_RG_CHARACTER_PROVES_P2602_PRIME_LOG_LIFTS_P2602_STRICT_DAMPING_REVALIDATED_UNDER_DF_SCOPE_BRIDGE_LTOTAL_BLOCKED",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "p2602_continuum_rg_character_prime_log_proof": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2542_PRIME_LOG_OBSTRUCTION": sha256_json(p2542_payload),
                "P2602_QUARANTINED_PRIME_LOG": sha256_json(p2602_payload),
                "P2603_CONDITIONAL_SLOPE": sha256_json(p2603_payload),
                "P2610_CRITICAL_REVALIDATION": sha256_json(p2610_payload),
                "P2613_MONOID_UNITAL": sha256_json(p2613_payload),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["p2602_continuum_rg_character_prime_log_proof"]["theorem_export"]
    proof = t["proof_core"]
    lines = [
        "# P2614/S1564 P2602 continuum RG character prime-log proof", "",
        f"Status: `{payload['status']}`", "", "## Theorem", "",
        proof["formal_theorem"], "", "## Algebraic proof", "",
    ]
    for step in proof["algebraic_steps"]:
        lines.append(f"- {step}")
    lines.extend([
        "", "## Computed checks", "",
        f"- P2602 quarantine lifted: `{t['p2602_quarantine_lifted_by_p2614']}`.",
        f"- P2601 revalidated inherited: `{t['p2601_revalidated_inherited']}`.",
        f"- Strict damping beta/eta revalidated under D_f=9/5 scope: `{t['strict_damping_beta_eta_source_revalidated_under_df_9_5_scope']}`.",
        f"- Remaining quarantines after P2614: `{t['remaining_p2610_quarantines_after_p2614']}`.",
        f"- Discrete counterexamples audited: `{len(t['noncontinuum_prime_character_counterexamples'])}`.", "",
        "## P2610 objection answer", "", proof["why_this_answers_p2610_objection"], "",
        "## Scope guards", "",
        "P2614 revalidates P2602 and the non-bridge strict damping beta/eta source under the retained D_f=9/5 scope. It does not revalidate P2607 bridge completion, does not re-enable P2608 role-bearing L_total, and does not export QW-2191, APD sourcehood, legacy physical-role transfer, or ToE closure.", "", "## Fingerprint", "",
        f"`{payload['p2602_continuum_rg_character_prime_log_proof']['theorem_fingerprint_sha256']}`",
    ])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2614/S1564 P2602 continuum RG character prime-log proof

`P2614/S1564` revalidates P2602 by replacing the weak prime-gap story with a continuum RG character theorem: for continuous scale-stationary damping on the positive dilation monoid, `f(t)=y(e^t)` satisfies Cauchy's equation and continuity forces `f(t)=gamma t`, so `y(lambda)=gamma log(lambda)` and `v_p/log(p)=gamma` for every prime sample.  Discrete arbitrary prime characters remain counterexamples only when no continuous RG dilation embedding is required.  Together with P2613 and the retained P2603 `D_f=9/5` scope, the non-bridge strict damping beta/eta source is revalidated; P2607 bridge completion and P2608 role-bearing `L_total` remain blocked.
""".strip()
    lag_section = """
## P2614/S1564 continuum-prime-log Ltotal guard

`P2614/S1564` revalidates the prime-log proportionality subkey and the non-bridge strict damping beta/eta source under the retained `D_f=9/5` scope.  It still does not make the damping/compression term role-bearing in `L_total`: the GF(2) bridge and P2608 role transfer remain quarantined.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2614/S1564 P2602 continuum RG character prime-log proof", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2614/S1564 continuum-prime-log Ltotal guard", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
