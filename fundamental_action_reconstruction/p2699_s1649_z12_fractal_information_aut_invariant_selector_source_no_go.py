#!/usr/bin/env python3
"""P2699/S1649: Z12 fractal-information Aut-invariant selector-source no-go.

Executes the P2698 recommendation in the only honest way available without a
new external premise: test a candidate *internal* selector source built from the
nadsoliton-as-fractal-information ontology and require it to be Aut(Z12)-safe.
The finite calculation proves that Aut-invariant information on Z12 cannot pick
a directed +1 generator, distinguish +1 from -1, or export QW-2191 closure.
"""
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2699_s1649_z12_fractal_information_aut_invariant_selector_source_no_go.json"
MD = GEN / "p2699_s1649_z12_fractal_information_aut_invariant_selector_source_no_go.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P2698": GEN / "p2698_s1648_symmetry_breaking_direction_claim_reconciliation_audit.json",
    "H36": GEN / "h36_directed_axis_orientation_audit.json",
    "H37": GEN / "h37_sign_distinction_state_audit.json",
    "H39": GEN / "h39_global_selector_object_absence_audit.json",
    "P739": GEN / "p739_current_strict_t193_global_premise_based_directed_selector_state_pair12_witness_split_strict_core_upgrade_bridge_nonexport_audit_probe_summary.json",
    "P740": GEN / "p740_current_strict_t194_global_sign_fixed_directed_closure_pair12_witness_split_strict_core_upgrade_bridge_nonexport_audit_probe_summary.json",
    "KAPPA": GEN / "kappa_z12_generator_orientation_canonical_fixing_datum_strict_provenance_v1.json",
}

NEGATIVE_EXPORT_FLAGS = [
    "new_nonpremise_selector_source_exported",
    "directed_generator_selected",
    "plus_minus_generator_distinguished_aut_invariantly",
    "qw2191_discharged",
    "strict_core_selector_closure_exported",
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
        "p2698_next_step": r"P2699|one-object strict-internal selector-source|non-premise-based strict selector/orientation source",
        "fractal_information_ontology": r"nadsoliton.*fractal information|primordial fractal information|fractal|information in a solitonic state",
        "z12_aut_invariant_boundary": r"Aut\(Z_12\)|Aut\(Z12\)|generator/orientation|Aut-invariant canonicity|Z_12/Aut",
        "premise_based_direction_support": r"H36|H37|T164|premise-based|directed/sign-sensitive|strict provenance fixing datum",
        "qw2191_still_open": r"QW-2191.*open|QW2191_STILL_OPEN|global QW-2191 discharge remains open|qw2191_closed.*false",
        "forbidden_promotions": r"strict-core selector closure|L_total|ToE closure|role transfer|pair12 strict-core upgrade",
    }
    return {"tool": "rg", "mode": "P2699 finite Aut(Z12) selector-source candidate no-go", "patterns": {key: rg_count(pattern) for key, pattern in patterns.items()}}


def aut_z12() -> list[int]:
    return [u for u in range(12) if math.gcd(u, 12) == 1]


def orbit(seed: int, units: list[int]) -> list[int]:
    return sorted({(u * seed) % 12 for u in units})


def z12_aut_calculation() -> dict[str, Any]:
    units = aut_z12()
    elements = list(range(12))
    seen: set[int] = set()
    orbits: list[list[int]] = []
    for element in elements:
        if element not in seen:
            row = orbit(element, units)
            orbits.append(row)
            seen.update(row)
    fixed_points = [x for x in elements if all((u * x) % 12 == x for u in units)]
    directed_generators = units[:]
    generator_orbit = orbit(1, units)
    plus_minus_same_orbit = 11 in generator_orbit
    invariant_translation_strides = fixed_points[:]
    invariant_scalar_basis_dimension = len(orbits)
    orbit_by_element = {x: row for row in orbits for x in row}
    valuation_signature = {
        str(x): {
            "gcd_with_12": math.gcd(x, 12),
            "orbit": orbit_by_element[x],
            "unit_generator": x in directed_generators,
        }
        for x in elements
    }
    return {
        "group": "Aut(Z12)=U(12)",
        "units": units,
        "elements": elements,
        "orbits": orbits,
        "fixed_points": fixed_points,
        "directed_generators": directed_generators,
        "directed_generator_orbit": generator_orbit,
        "plus_one_and_minus_one_same_orbit": plus_minus_same_orbit,
        "invariant_translation_strides": invariant_translation_strides,
        "invariant_scalar_basis_dimension": invariant_scalar_basis_dimension,
        "valuation_signature_by_element": valuation_signature,
        "can_select_unique_directed_generator_aut_invariantly": len(generator_orbit) == 1,
        "can_distinguish_plus_one_from_minus_one_aut_invariantly": not plus_minus_same_orbit,
        "can_select_nonzero_stride_aut_invariantly_without_premise": invariant_translation_strides == [6],
    }


def state_reads() -> dict[str, Any]:
    loaded = {name: load_json(path) for name, path in INPUTS.items()}
    p2698 = loaded["P2698"]
    h36 = loaded["H36"]
    h37 = loaded["H37"]
    h39 = loaded["H39"]
    p739 = loaded["P739"]
    p740 = loaded["P740"]
    kappa = loaded["KAPPA"]
    return {
        "hashes": {name: sha256_file(path) for name, path in INPUTS.items()},
        "p2698_recommended_p2699": "P2699" in p2698.get("decision", {}).get("next_honest_step", ""),
        "ontology_candidate_is_internal_only": True,
        "ontology_candidate_description": "Use only nadsoliton-as-primordial-fractal-information symmetry data on the typed Z12 carrier; do not import T164/KAPPA as a source theorem.",
        "existing_direction_support_real_but_premise_based": "PREMISE_BASED_T164" in h36.get("status", "") and "PREMISE_BASED_T164" in h37.get("status", ""),
        "h39_qw2191_still_open": "global_qw_2191_discharge" in h39.get("missing", []) or "QW2191_STILL_OPEN" in h39.get("status", ""),
        "p739_pair12_upgrade_unexported": p739.get("t193_target_exported_on_current_repo_state") is False,
        "p740_pair12_upgrade_unexported": p740.get("t194_target_exported_on_current_repo_state") is False,
        "kappa_is_premise_based_not_used_as_nonpremise_source": "premise-based" in json.dumps(kappa).lower() and "does not claim Aut-invariant canonicity" in json.dumps(kappa),
    }


def no_go_matrix(calc: dict[str, Any], reads: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "candidate_obligation": "unique_directed_generator_from_aut_invariant_fractal_information",
            "finite_test": "Aut(Z12) orbit of +1",
            "computed_witness": calc["directed_generator_orbit"],
            "passes": calc["can_select_unique_directed_generator_aut_invariantly"],
            "no_go_reason": "+1 has a four-element Aut orbit, so invariant information cannot select it as a unique directed generator.",
        },
        {
            "candidate_obligation": "distinguish_plus_one_from_minus_one_without_premise",
            "finite_test": "membership of -1=11 in the +1 Aut orbit",
            "computed_witness": {"plus_one_orbit": calc["directed_generator_orbit"], "minus_one_in_orbit": calc["plus_one_and_minus_one_same_orbit"]},
            "passes": calc["can_distinguish_plus_one_from_minus_one_aut_invariantly"],
            "no_go_reason": "+1 and -1 lie in the same Aut orbit, so a scalar Aut-invariant fractal-information source cannot provide directed sign choice.",
        },
        {
            "candidate_obligation": "derive_existing_premise_based_direction_as_nonpremise_source",
            "finite_test": "P2698/H36/H37/P739/P740 scope read plus Aut calculation",
            "computed_witness": {
                "existing_direction_support_real_but_premise_based": reads["existing_direction_support_real_but_premise_based"],
                "kappa_is_premise_based_not_used_as_nonpremise_source": reads["kappa_is_premise_based_not_used_as_nonpremise_source"],
            },
            "passes": False,
            "no_go_reason": "The existing directed lane is real but premise-based; P2699 does not add a new non-premise source theorem.",
        },
        {
            "candidate_obligation": "strict_core_selector_or_qw2191_closure",
            "finite_test": "QW-2191/P739/P740 boundary read",
            "computed_witness": {
                "h39_qw2191_still_open": reads["h39_qw2191_still_open"],
                "p739_pair12_upgrade_unexported": reads["p739_pair12_upgrade_unexported"],
                "p740_pair12_upgrade_unexported": reads["p740_pair12_upgrade_unexported"],
            },
            "passes": False,
            "no_go_reason": "The finite Aut calculation supplies no missing strict-core selector source and does not change P739/P740 nonexport boundaries.",
        },
    ]


def decision(rows: list[dict[str, Any]]) -> dict[str, Any]:
    bounded_no_go = all(not row["passes"] for row in rows)
    return {
        "decision": "P2699_Z12_FRACTAL_INFORMATION_AUT_INVARIANT_SELECTOR_SOURCE_NO_GO_NO_FALSE_PASS",
        "bounded_no_go_now": bounded_no_go,
        "reason": "A pure internal fractal-information candidate constrained to Aut(Z12)-invariant data cannot choose a unique directed generator or distinguish +1 from -1; the real directed lane remains premise-based rather than a new strict source theorem.",
        "next_honest_step": "No proof-grade closure move is unlocked by P2699.  The next admissible move must introduce a genuinely new non-premise strict selector/orientation source beyond Aut-invariant Z12/fractal-information data, or pivot to a new typed object outside the closed selector/direct/bridge/Lagrangian lanes.  Otherwise preserve the P2697-P2699 no-new-live-frontier certificate.",
        "forbidden_promotions": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = ["# P2699/S1649 Z12 fractal-information Aut-invariant selector-source no-go", "", f"Status: `{payload['status']}`", "", "## Finite Aut(Z12) calculation"]
    calc = payload["z12_aut_calculation"]
    lines.append(f"- units: `{calc['units']}`")
    lines.append(f"- orbits: `{calc['orbits']}`")
    lines.append(f"- fixed points / invariant translation strides: `{calc['fixed_points']}`")
    lines.append(f"- directed generator orbit of +1: `{calc['directed_generator_orbit']}`")
    lines.extend(["", "## No-go matrix"])
    for row in payload["no_go_matrix"]:
        lines.append(f"- `{row['candidate_obligation']}`: passes=`{row['passes']}`; {row['no_go_reason']}")
    lines.extend(["", "## Decision", payload["decision"]["reason"], "", "## Next honest step", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    reads = state_reads()
    calc = z12_aut_calculation()
    rows = no_go_matrix(calc, reads)
    payload: dict[str, Any] = {
        "status": "P2699_Z12_FRACTAL_INFORMATION_AUT_INVARIANT_SELECTOR_SOURCE_NO_GO",
        "content_grep": content_grep(),
        "state_reads": reads,
        "z12_aut_calculation": calc,
        "no_go_matrix": rows,
        "decision": decision(rows),
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2699/S1649 Z12 fractal-information selector-source no-go",
        "## P2699/S1649 Z12 fractal-information selector-source no-go\n\n"
        "`P2699/S1649` executes the one-object P2698 follow-up as a finite Aut(Z12) calculation.  Treating the nadsoliton ontology as pure primordial fractal information does not by itself break Aut(Z12): the directed generator orbit is `{1,5,7,11}`, `+1` and `-1` are in the same orbit, and the only Aut-fixed Z12 elements/translation strides are `0` and `6`.  Thus no non-premise strict selector/orientation source, `QW-2191` discharge, pair12 strict-core upgrade, role-bearing `L_total`, or ToE closure is exported.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2699/S1649 Aut-invariant selector-source Ltotal guard",
        "## P2699/S1649 Aut-invariant selector-source Ltotal guard\n\n"
        "`P2699/S1649` is a finite group-theoretic no-go, not a variational source.  Aut-invariant fractal-information data on Z12 cannot supply the missing non-premise directed selector source, so `L_total`, strict selector closure, role transfer, and ToE closure remain unpromoted.\n",
    )
    append_once(
        AGENTS,
        "Current Z12 fractal-information selector-source no-go guardrail (P2699/S1649, 2026-06-13)",
        "## Current Z12 fractal-information selector-source no-go guardrail (P2699/S1649, 2026-06-13)\n\n"
        "- P2699 tests the honest one-object follow-up to P2698: an internal nadsoliton-as-fractal-information selector source constrained by Aut(Z12)-invariance.  The finite calculation finds no unique directed generator and no Aut-invariant distinction of `+1` from `-1`; existing directed support remains premise-based.\n"
        "- Do not claim `QW-2191` discharge, strict selector closure, pair12 strict-core upgrade, role-bearing `L_total`, role transfer, bridge closure, or ToE closure from this candidate.\n"
        "- No proof-grade closure move is unlocked unless a genuinely new non-premise strict selector/orientation source or a new typed object outside the closed lanes is exported; otherwise preserve the P2697-P2699 no-new-live-frontier certificate.\n",
    )
    return payload


if __name__ == "__main__":
    main()
