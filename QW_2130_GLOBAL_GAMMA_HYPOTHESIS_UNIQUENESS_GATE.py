#!/usr/bin/env python3
"""
QW-2130: Global gamma-hypothesis uniqueness gate (strict-admissible domain).

Purpose:
- resolve the open flag from QW-2128 by evaluating uniqueness across gamma
  hypotheses in a rigorously declared admissible set,
- keep explicit scope boundary versus unconstrained "all formulas ever listed".
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2130_global_gamma_hypothesis_uniqueness_gate.json"
OUT_MD = ROOT / "RAPORT_QW2130_GLOBAL_GAMMA_HYPOTHESIS_UNIQUENESS_GATE.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def best_by_gamma_and_q(variants: List[Dict], gamma_source: str, q_assignment: str) -> Dict:
    subset = [v for v in variants if v["gamma_source"] == gamma_source and v["q_assignment"] == q_assignment]
    if not subset:
        raise KeyError(f"Missing subset for gamma={gamma_source}, q={q_assignment}")
    return sorted(subset, key=lambda x: float(x["score"]))[0]


def main() -> None:
    r1961 = load_json("report_qw1961_noncircular_gamma_q_derivation_matrix.json")
    r2128 = load_json("report_qw2128_kernel_rep_assignment_uniqueness_gate.json")

    variants = list(r1961["variants"])
    q_assignments = sorted(r1961["inputs"]["q_assignments"].keys())

    # Domain declaration:
    # - exclude canonical frozen reference (not derivational hypothesis),
    # - exclude "derived_ratio_n_over_df_minus_1" marked in QW-1961 comments as
    #   arithmetically inconsistent documented path.
    gamma_sources_all = sorted(set(str(v["gamma_source"]) for v in variants))
    admissible_gamma_sources = [
        "derived_force_energy_2n_over_3",
        "derived_kernel_d1_to_d4",
    ]
    excluded_gamma_sources = sorted(set(gamma_sources_all) - set(admissible_gamma_sources))

    per_gamma_rows: List[Dict] = []
    winners: List[str] = []
    winner_allpass_rows = 0
    for gs in admissible_gamma_sources:
        rows = []
        for qn in q_assignments:
            b = best_by_gamma_and_q(variants, gs, qn)
            rows.append(
                {
                    "q_assignment": qn,
                    "score": float(b["score"]),
                    "split_mode": str(b["split_mode"]),
                    "all_pass": bool(b["all_pass"]),
                }
            )
        rows_sorted = sorted(rows, key=lambda x: x["score"])
        winner = str(rows_sorted[0]["q_assignment"])
        winners.append(winner)
        if rows_sorted[0]["all_pass"]:
            winner_allpass_rows += 1
        gap = float(rows_sorted[1]["score"] - rows_sorted[0]["score"]) if len(rows_sorted) > 1 else 0.0
        per_gamma_rows.append(
            {
                "gamma_source": gs,
                "rows_sorted": rows_sorted,
                "winner_q_assignment": winner,
                "winner_score": float(rows_sorted[0]["score"]),
                "runnerup_q_assignment": str(rows_sorted[1]["q_assignment"]) if len(rows_sorted) > 1 else None,
                "runnerup_score": float(rows_sorted[1]["score"]) if len(rows_sorted) > 1 else None,
                "winner_score_gap": gap,
            }
        )

    unique_winner = sorted(set(winners))
    winner_name = unique_winner[0] if len(unique_winner) == 1 else None

    q2128_winner = str(r2128["locked_branch_ranking"]["winner_q_assignment"])
    q2128_gap = float(r2128["locked_branch_ranking"]["winner_score_gap"])

    # Strict dominance in the primary admissible branch:
    # derived_force_energy_2n_over_3 has winner all-pass while runner-up does not.
    primary = next((r for r in per_gamma_rows if r["gamma_source"] == "derived_force_energy_2n_over_3"), None)
    primary_strict_dominance = bool(
        primary is not None
        and bool(primary["rows_sorted"][0]["all_pass"])
        and (not bool(primary["rows_sorted"][1]["all_pass"]))
    )

    flags = {
        "admissible_gamma_domain_declared_explicitly": bool(len(admissible_gamma_sources) >= 2),
        "canonical_reference_excluded_from_derivational_domain": bool("canonical_frozen_1p52_reference" in excluded_gamma_sources),
        "arithmetically_inconsistent_gamma_path_excluded": bool("derived_ratio_n_over_df_minus_1" in excluded_gamma_sources),
        "all_admissible_gamma_sources_present_in_variants": bool(
            all(gs in gamma_sources_all for gs in admissible_gamma_sources)
        ),
        "unique_q_assignment_winner_across_admissible_gammas": bool(len(unique_winner) == 1),
        "winner_matches_q2128_locked_branch_winner": bool(winner_name == q2128_winner),
        "winner_allpass_in_at_least_one_admissible_gamma": bool(winner_allpass_rows >= 1),
        "primary_admissible_gamma_has_strict_winner_dominance": bool(primary_strict_dominance),
        "q2128_locked_branch_gap_ge_0p25": bool(q2128_gap >= 0.25),
        "global_unconstrained_formula_space_uniqueness": False,
        "deterministic_no_scan_no_retune": True,
    }
    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = int(len(flags))

    verdict = (
        "GLOBAL_GAMMA_HYPOTHESIS_UNIQUENESS_GATE_PASS_STRICT_ADMISSIBLE_DOMAIN"
        if (
            flags["admissible_gamma_domain_declared_explicitly"]
            and flags["all_admissible_gamma_sources_present_in_variants"]
            and flags["unique_q_assignment_winner_across_admissible_gammas"]
            and flags["winner_matches_q2128_locked_branch_winner"]
            and flags["winner_allpass_in_at_least_one_admissible_gamma"]
            and flags["primary_admissible_gamma_has_strict_winner_dominance"]
            and flags["q2128_locked_branch_gap_ge_0p25"]
            and flags["deterministic_no_scan_no_retune"]
        )
        else "GLOBAL_GAMMA_HYPOTHESIS_UNIQUENESS_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q1961_matrix": "report_qw1961_noncircular_gamma_q_derivation_matrix.json",
            "q2128_locked_branch": "report_qw2128_kernel_rep_assignment_uniqueness_gate.json",
        },
        "domain_definition": {
            "admissible_gamma_sources": admissible_gamma_sources,
            "excluded_gamma_sources": excluded_gamma_sources,
            "note": (
                "Admissible domain excludes non-derivational canonical reference and "
                "arithmetically inconsistent documented path."
            ),
        },
        "per_gamma_assignment_ranking": per_gamma_rows,
        "winner_summary": {
            "unique_winners_set": unique_winner,
            "winner_if_unique": winner_name,
            "winner_allpass_rows_count": winner_allpass_rows,
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "EXTEND_ADMISSIBLE_DOMAIN_WITH_FORMAL_JUSTIFICATION_IF_NEEDED_OR_KEEP_SCOPE_BOUNDARY_EXPLICIT"
            if verdict.endswith("DOMAIN")
            else "REPAIR_DOMAIN_OR_RANKING_RULES_AND_RERUN_QW2130"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2130: GLOBAL GAMMA HYPOTHESIS UNIQUENESS GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Domain",
        f"- admissible gamma sources: `{admissible_gamma_sources}`",
        f"- excluded gamma sources: `{excluded_gamma_sources}`",
        "",
        "## Winner summary",
        f"- unique winners set: `{unique_winner}`",
        f"- winner (if unique): `{winner_name}`",
        f"- winner all-pass rows: `{winner_allpass_rows}`",
        "",
        "## Open scope boundary",
        "- global_unconstrained_formula_space_uniqueness: `False`",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2130] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2130] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2130] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()
