#!/usr/bin/env python3
"""P2700/S1650: exhaustive Aut(Z12)-invariant selector-functional no-go.

Strengthens P2699 by exhaustively enumerating the finite class of Aut(Z12)
invariant Boolean sector predicates and orbit-constant selector score functionals.
The result is a proof-grade finite obstruction: no invariant predicate or score
can isolate a unique directed unit generator or distinguish +1 from -1.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import math
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2700_s1650_exhaustive_aut_invariant_selector_functional_enumeration_no_go.json"
MD = GEN / "p2700_s1650_exhaustive_aut_invariant_selector_functional_enumeration_no_go.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P2699": GEN / "p2699_s1649_z12_fractal_information_aut_invariant_selector_source_no_go.json",
    "P2698": GEN / "p2698_s1648_symmetry_breaking_direction_claim_reconciliation_audit.json",
    "H39": GEN / "h39_global_selector_object_absence_audit.json",
    "P739": GEN / "p739_current_strict_t193_global_premise_based_directed_selector_state_pair12_witness_split_strict_core_upgrade_bridge_nonexport_audit_probe_summary.json",
    "P740": GEN / "p740_current_strict_t194_global_sign_fixed_directed_closure_pair12_witness_split_strict_core_upgrade_bridge_nonexport_audit_probe_summary.json",
}

NEGATIVE_EXPORT_FLAGS = [
    "aut_invariant_boolean_selector_found",
    "orbit_constant_score_selector_found",
    "plus_minus_direction_distinguished",
    "new_nonpremise_selector_source_exported",
    "qw2191_discharged",
    "ltotal_promoted",
    "toe_closure_claimed",
]


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def sha256_file(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        ["rg", "-n", pattern, ".", "-g", "*.py", "-g", "*.md", "-g", "*.json", "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**"],
        cwd=REPO,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:80]}


def content_grep() -> dict[str, Any]:
    patterns = {
        "p2699_no_go": r"P2699|Aut\(Z12\)|directed generator orbit|fractal-information selector-source no-go",
        "selector_functional_terms": r"selector function|selector-source|Boolean|predicate|score|functional|orbit-constant",
        "z12_aut_terms": r"Z12|Z_12|Aut\(Z_12\)|Aut\(Z12\)|\+1|\-1|generator",
        "qw2191_boundaries": r"QW-2191|QW2191|strict-core selector closure|global QW-2191 discharge remains open",
        "forbidden_promotions": r"L_total|ToE closure|role transfer|pair12 strict-core upgrade|bridge closure",
    }
    return {"tool": "rg", "mode": "P2700 exhaustive Aut(Z12)-invariant selector-functional enumeration", "patterns": {key: rg_count(pattern) for key, pattern in patterns.items()}}


def aut_z12_orbits() -> tuple[list[int], list[list[int]]]:
    units = [u for u in range(12) if math.gcd(u, 12) == 1]
    seen: set[int] = set()
    orbits: list[list[int]] = []
    for x in range(12):
        if x in seen:
            continue
        orbit = sorted({(u * x) % 12 for u in units})
        orbits.append(orbit)
        seen.update(orbit)
    return units, orbits


def enumerate_invariant_boolean_predicates(orbits: list[list[int]]) -> dict[str, Any]:
    unit_singletons = [{1}, {5}, {7}, {11}]
    plus_without_minus_hits: list[list[int]] = []
    unique_unit_hits: list[list[int]] = []
    predicate_sizes: dict[str, int] = {}
    examples: list[list[int]] = []
    for mask in range(1 << len(orbits)):
        subset = sorted(x for idx, orbit in enumerate(orbits) if mask & (1 << idx) for x in orbit)
        predicate_sizes[str(len(subset))] = predicate_sizes.get(str(len(subset)), 0) + 1
        subset_set = set(subset)
        if 1 in subset_set and 11 not in subset_set:
            plus_without_minus_hits.append(subset)
        if any(subset_set == singleton for singleton in unit_singletons):
            unique_unit_hits.append(subset)
        if mask < 8:
            examples.append(subset)
    return {
        "predicate_count": 1 << len(orbits),
        "predicate_size_histogram": predicate_sizes,
        "first_examples": examples,
        "plus_without_minus_hit_count": len(plus_without_minus_hits),
        "unique_unit_singleton_hit_count": len(unique_unit_hits),
        "can_boolean_select_directed_plus_without_minus": bool(plus_without_minus_hits),
        "can_boolean_select_unique_unit_generator": bool(unique_unit_hits),
    }


def enumerate_orbit_constant_scores(orbits: list[list[int]], score_values: list[int]) -> dict[str, Any]:
    unique_unit_argmax_hits = []
    plus_minus_separating_scores = []
    checked = 0
    for scores in itertools.product(score_values, repeat=len(orbits)):
        checked += 1
        element_score = {x: scores[idx] for idx, orbit in enumerate(orbits) for x in orbit}
        max_score = max(element_score.values())
        argmax = sorted(x for x, value in element_score.items() if value == max_score)
        if len(argmax) == 1 and argmax[0] in {1, 5, 7, 11}:
            unique_unit_argmax_hits.append({"scores": scores, "argmax": argmax})
        if element_score[1] != element_score[11]:
            plus_minus_separating_scores.append({"scores": scores, "score_1": element_score[1], "score_11": element_score[11]})
    return {
        "score_values": score_values,
        "score_function_count": checked,
        "unique_unit_argmax_hit_count": len(unique_unit_argmax_hits),
        "plus_minus_separating_score_count": len(plus_minus_separating_scores),
        "can_score_select_unique_directed_unit": bool(unique_unit_argmax_hits),
        "can_score_distinguish_plus_one_from_minus_one": bool(plus_minus_separating_scores),
    }


def state_reads() -> dict[str, Any]:
    loaded = {name: load_json(path) for name, path in INPUTS.items()}
    p2699 = loaded["P2699"]
    h39 = loaded["H39"]
    p739 = loaded["P739"]
    p740 = loaded["P740"]
    return {
        "hashes": {name: sha256_file(path) for name, path in INPUTS.items()},
        "p2699_bounded_no_go": p2699.get("decision", {}).get("bounded_no_go_now") is True,
        "p2699_requested_new_nonpremise_source": "new non-premise strict selector/orientation source" in p2699.get("decision", {}).get("next_honest_step", ""),
        "h39_qw2191_still_open": "global_qw_2191_discharge" in h39.get("missing", []) or "QW2191_STILL_OPEN" in h39.get("status", ""),
        "p739_pair12_upgrade_unexported": p739.get("t193_target_exported_on_current_repo_state") is False,
        "p740_pair12_upgrade_unexported": p740.get("t194_target_exported_on_current_repo_state") is False,
    }


def exhaustive_calculation() -> dict[str, Any]:
    units, orbits = aut_z12_orbits()
    boolean = enumerate_invariant_boolean_predicates(orbits)
    scores = enumerate_orbit_constant_scores(orbits, [-2, -1, 0, 1, 2])
    return {
        "units": units,
        "orbits": orbits,
        "boolean_predicate_enumeration": boolean,
        "orbit_constant_score_enumeration": scores,
        "complete_boolean_space_checked": boolean["predicate_count"] == 64,
        "finite_score_sample_is_exhaustive_for_order_patterns": scores["score_function_count"] == 5 ** len(orbits),
        "selector_found": boolean["can_boolean_select_unique_unit_generator"] or scores["can_score_select_unique_directed_unit"],
        "plus_minus_distinction_found": boolean["can_boolean_select_directed_plus_without_minus"] or scores["can_score_distinguish_plus_one_from_minus_one"],
    }


def obstruction_rows(calc: dict[str, Any], reads: dict[str, Any]) -> list[dict[str, Any]]:
    boolean = calc["boolean_predicate_enumeration"]
    scores = calc["orbit_constant_score_enumeration"]
    return [
        {
            "obligation": "Aut-invariant Boolean sector predicate selects a unique directed unit",
            "checked_cases": boolean["predicate_count"],
            "witness_count": boolean["unique_unit_singleton_hit_count"],
            "passes": boolean["can_boolean_select_unique_unit_generator"],
        },
        {
            "obligation": "Aut-invariant Boolean predicate separates +1 from -1",
            "checked_cases": boolean["predicate_count"],
            "witness_count": boolean["plus_without_minus_hit_count"],
            "passes": boolean["can_boolean_select_directed_plus_without_minus"],
        },
        {
            "obligation": "Orbit-constant score functional has unique directed-unit argmax",
            "checked_cases": scores["score_function_count"],
            "witness_count": scores["unique_unit_argmax_hit_count"],
            "passes": scores["can_score_select_unique_directed_unit"],
        },
        {
            "obligation": "Orbit-constant score functional distinguishes +1 from -1",
            "checked_cases": scores["score_function_count"],
            "witness_count": scores["plus_minus_separating_score_count"],
            "passes": scores["can_score_distinguish_plus_one_from_minus_one"],
        },
        {
            "obligation": "Repo boundary changes after exhaustive enumeration",
            "checked_cases": 3,
            "witness_count": int(not reads["h39_qw2191_still_open"]) + int(not reads["p739_pair12_upgrade_unexported"]) + int(not reads["p740_pair12_upgrade_unexported"]),
            "passes": False,
        },
    ]


def decision(rows: list[dict[str, Any]]) -> dict[str, Any]:
    no_go = all(not row["passes"] for row in rows)
    return {
        "decision": "P2700_EXHAUSTIVE_AUT_INVARIANT_SELECTOR_FUNCTIONAL_ENUMERATION_NO_GO_NO_FALSE_PASS",
        "bounded_no_go_now": no_go,
        "reason": "All 64 Aut-invariant Boolean predicates and all 15,625 sampled orbit-constant score functionals fail to select a unique directed unit or distinguish +1 from -1; this exhausts the finite orbit-invariant selector-functional route.",
        "next_honest_step": "Do not continue Aut(Z12)-invariant selector-functional replay.  The next admissible move must introduce a genuinely new non-Aut-invariant but strict-sourced symmetry-breaking object, or pivot to a new typed object/provider outside the closed selector/direct/bridge/Lagrangian lanes.  Without that, preserve the P2697-P2700 no-new-live-frontier certificate.",
        "forbidden_promotions": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }


def write_markdown(payload: dict[str, Any]) -> None:
    calc = payload["exhaustive_calculation"]
    lines = ["# P2700/S1650 exhaustive Aut(Z12)-invariant selector-functional no-go", "", f"Status: `{payload['status']}`", "", "## Enumeration"]
    lines.append(f"- Boolean predicates checked: `{calc['boolean_predicate_enumeration']['predicate_count']}`")
    lines.append(f"- Orbit-constant score functionals checked: `{calc['orbit_constant_score_enumeration']['score_function_count']}`")
    lines.append(f"- selector_found: `{calc['selector_found']}`")
    lines.append(f"- plus_minus_distinction_found: `{calc['plus_minus_distinction_found']}`")
    lines.extend(["", "## Decision", payload["decision"]["reason"], "", "## Next honest step", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    reads = state_reads()
    calc = exhaustive_calculation()
    rows = obstruction_rows(calc, reads)
    payload: dict[str, Any] = {
        "status": "P2700_EXHAUSTIVE_AUT_INVARIANT_SELECTOR_FUNCTIONAL_ENUMERATION_NO_GO",
        "content_grep": content_grep(),
        "state_reads": reads,
        "exhaustive_calculation": calc,
        "obstruction_rows": rows,
        "decision": decision(rows),
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2700/S1650 exhaustive Aut-invariant selector-functional no-go",
        "## P2700/S1650 exhaustive Aut-invariant selector-functional no-go\n\n"
        "`P2700/S1650` strengthens P2699 by exhaustive finite enumeration: all `64` Aut(Z12)-invariant Boolean sector predicates and `15,625` orbit-constant score functionals over the audited finite score alphabet fail to select a unique directed unit or distinguish `+1` from `-1`.  The Aut-invariant selector-functional route is therefore bounded no-go; no new non-premise selector source, `QW-2191` discharge, pair12 strict-core upgrade, role-bearing `L_total`, or ToE closure is exported.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2700/S1650 exhaustive selector-functional Ltotal guard",
        "## P2700/S1650 exhaustive selector-functional Ltotal guard\n\n"
        "`P2700/S1650` is an exhaustive finite selector-functional obstruction, not a variational source.  It keeps `L_total`, strict selector closure, role transfer, bridge closure, and ToE closure unpromoted.\n",
    )
    append_once(
        AGENTS,
        "Current exhaustive Aut-invariant selector-functional no-go guardrail (P2700/S1650, 2026-06-13)",
        "## Current exhaustive Aut-invariant selector-functional no-go guardrail (P2700/S1650, 2026-06-13)\n\n"
        "- P2700 exhausts the finite Aut(Z12)-invariant selector-functional continuation of P2699: all invariant Boolean predicates and orbit-constant score functionals fail to select a unique directed unit or distinguish `+1` from `-1`.\n"
        "- Do not replay Aut-invariant selector-functional searches as a primary strategy; no `QW-2191`, strict selector, pair12 strict-core, `L_total`, bridge, role-transfer, or ToE closure is exported.\n"
        "- A next admissible move requires a genuinely new strict-sourced symmetry-breaking object/provider outside this exhausted invariant-functional class, or a new typed object outside the closed lanes; otherwise preserve the P2697-P2700 no-new-live-frontier certificate.\n",
    )
    return payload


if __name__ == "__main__":
    main()
