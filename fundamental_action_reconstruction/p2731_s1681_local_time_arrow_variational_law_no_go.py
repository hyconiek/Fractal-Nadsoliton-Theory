#!/usr/bin/env python3
"""P2731/S1681: local time-arrow variational source-law no-go.

P2730 closed bare Z12 tau fields.  P2731 tests a more dynamical candidate:
translation-invariant nearest-neighbour variational laws for tau: Z12->{-1,+1}.
It exhausts all integer local pair-energy tables e(tau_i,tau_{i+1}) in {-1,0,+1}
and asks whether an admissible time-reversal-even local Hamiltonian can choose
one global time-arrow sign rather than the paired +/- constants.
"""
from __future__ import annotations

import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2731_s1681_local_time_arrow_variational_law_no_go.json"
MD = GEN / "p2731_s1681_local_time_arrow_variational_law_no_go.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

Z12 = 12
SPINS = (-1, 1)
PAIR_ORDER = ((-1, -1), (-1, 1), (1, -1), (1, 1))
FIELDS = tuple(itertools.product(SPINS, repeat=Z12))
PLUS_FIELD = (1,) * Z12
MINUS_FIELD = (-1,) * Z12
INPUTS = {
    "P2730_TIME_ARROW_SOURCE_FIELD_NO_GO": GEN / "p2730_s1680_time_arrow_source_field_equivariance_no_go.json",
    "P2729_TIME_ARROW_AUDIT": GEN / "p2729_s1679_time_arrow_orientation_coupling_law_audit.json",
    "P2721_SIGN_TORSOR_COUPLING": GEN / "p2721_s1671_chiral_bispectrum_sign_to_orientation_torsor_coupling_audit.json",
}
NEGATIVE_EXPORT_FLAGS = [
    "strict_time_arrow_variational_source_exported",
    "unique_nonpremise_tau_sign_selected",
    "time_reversal_even_local_law_selects_single_tau",
    "time_arrow_p2721_polarity_coupling_exported",
    "strict_mechanism_fixing_lambda_exported",
    "qw2191_discharged",
    "pair12_strict_core_upgrade_exported",
    "bridge_closure_exported",
    "role_transfer_started",
    "ltotal_promoted",
    "toe_closure_exported",
]


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def neg(field: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(-x for x in field)


def energy(field: tuple[int, ...], table: dict[tuple[int, int], int]) -> int:
    return sum(table[(field[i], field[(i + 1) % Z12])] for i in range(Z12))


def time_reversal_even(table: dict[tuple[int, int], int]) -> bool:
    return all(table[pair] == table[(-pair[0], -pair[1])] for pair in PAIR_ORDER)


def table_label(table: dict[tuple[int, int], int]) -> str:
    return ",".join(f"e{a}{b}={table[(a, b)]}" for a, b in PAIR_ORDER)


def ground_state_summary(table: dict[tuple[int, int], int]) -> dict[str, Any]:
    energies = [energy(field, table) for field in FIELDS]
    min_energy = min(energies)
    ground_states = [field for field, value in zip(FIELDS, energies) if value == min_energy]
    plus_ground = PLUS_FIELD in ground_states
    minus_ground = MINUS_FIELD in ground_states
    unpaired_ground_states = [field for field in ground_states if neg(field) not in ground_states]
    return {
        "min_energy": min_energy,
        "ground_state_count": len(ground_states),
        "plus_constant_ground": plus_ground,
        "minus_constant_ground": minus_ground,
        "selects_plus_constant_only": plus_ground and not minus_ground,
        "selects_minus_constant_only": minus_ground and not plus_ground,
        "unpaired_ground_state_count": len(unpaired_ground_states),
        "ground_magnetizations": sorted(set(sum(field) for field in ground_states)),
    }


def audit_laws() -> dict[str, Any]:
    laws: list[dict[str, Any]] = []
    for values in itertools.product((-1, 0, 1), repeat=len(PAIR_ORDER)):
        table = {pair: value for pair, value in zip(PAIR_ORDER, values)}
        summary = ground_state_summary(table)
        even = time_reversal_even(table)
        selector = summary["selects_plus_constant_only"] or summary["selects_minus_constant_only"] or summary["unpaired_ground_state_count"] > 0
        laws.append({
            "label": table_label(table),
            "values": {f"{a},{b}": table[(a, b)] for a, b in PAIR_ORDER},
            "time_reversal_even": even,
            "selects_unpaired_tau_ground_state": selector,
            **summary,
        })

    even_laws = [law for law in laws if law["time_reversal_even"]]
    odd_or_mixed_selectors = [law for law in laws if not law["time_reversal_even"] and law["selects_unpaired_tau_ground_state"]]
    even_selectors = [law for law in even_laws if law["selects_unpaired_tau_ground_state"]]
    ferro_even_laws = [law for law in even_laws if law["plus_constant_ground"] and law["minus_constant_ground"]]
    ground_count_hist = Counter(law["ground_state_count"] for law in even_laws)
    return {
        "field_count": len(FIELDS),
        "local_pair_energy_table_count": len(laws),
        "time_reversal_even_law_count": len(even_laws),
        "time_reversal_even_unpaired_selector_count": len(even_selectors),
        "non_time_reversal_even_unpaired_selector_count": len(odd_or_mixed_selectors),
        "time_reversal_even_laws_with_both_constant_ground_states": len(ferro_even_laws),
        "time_reversal_even_ground_state_count_histogram": dict(sorted(ground_count_hist.items())),
        "sample_non_time_reversal_even_selectors": odd_or_mixed_selectors[:5],
        "sample_time_reversal_even_laws": even_laws[:5],
        "obstruction": "Every time-reversal-even translation-invariant nearest-neighbour tau Hamiltonian has sign-paired ground states.  Laws that select +tau over -tau exist only when the local energy table itself contains a time-reversal-odd sign bias, which is exactly the missing strict source rather than a derived non-premise selector.",
    }


def acceptance_matrix(audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "all_local_pair_energy_tables_enumerated": audit["local_pair_energy_table_count"] == 81,
        "all_tau_fields_scanned_per_law": audit["field_count"] == 4096,
        "time_reversal_even_law_class_nonempty": audit["time_reversal_even_law_count"] == 9,
        "time_reversal_even_unpaired_selector_exists": audit["time_reversal_even_unpaired_selector_count"] > 0,
        "strict_nonpremise_time_arrow_variational_source_exported": False,
        "p2721_polarity_coupling_theorem_exported": False,
    }
    missing = [k for k, v in facts.items() if not v]
    return {
        "facts": facts,
        "missing_criteria": missing,
        "accepted_as_strict_time_arrow_source_law": not missing,
        "blocker": "The admissible even local variational class cannot choose one tau sign.  Odd tables can choose a sign only by inserting the missing signed source into the energy law.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    a = payload["local_variational_time_arrow_audit"]
    lines = [
        "# P2731/S1681 local time-arrow variational source-law no-go",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Exhaustive local variational audit",
        f"- field_count={a['field_count']}",
        f"- local_pair_energy_table_count={a['local_pair_energy_table_count']}",
        f"- time_reversal_even_law_count={a['time_reversal_even_law_count']}",
        f"- time_reversal_even_unpaired_selector_count={a['time_reversal_even_unpaired_selector_count']}",
        f"- non_time_reversal_even_unpaired_selector_count={a['non_time_reversal_even_unpaired_selector_count']}",
        f"- time_reversal_even_ground_state_count_histogram={a['time_reversal_even_ground_state_count_histogram']}",
        a["obstruction"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    inputs = {key: read_json(path) for key, path in INPUTS.items()}
    audit = audit_laws()
    acceptance = acceptance_matrix(audit)
    no_go = not acceptance["accepted_as_strict_time_arrow_source_law"]
    payload = {
        "status": "P2731_LOCAL_TIME_ARROW_VARIATIONAL_LAW_NO_GO" if no_go else "P2731_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "input_statuses": {key: value.get("status") for key, value in inputs.items()},
        "audited_candidate_class": "translation-invariant nearest-neighbour tau Hamiltonians H=sum_i e(tau_i,tau_{i+1}) with e-values in {-1,0,+1}",
        "local_variational_time_arrow_audit": audit,
        "acceptance_matrix": acceptance,
        "decision": {
            "strict_time_arrow_variational_source_exported": False,
            "unique_nonpremise_tau_sign_selected": False,
            "time_arrow_p2721_polarity_coupling_exported": False,
            "strict_mechanism_fixing_lambda_exported": False,
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": "P2731 exhausts finite nearest-neighbour local tau Hamiltonians.  The time-reversal-even admissible subclass has no unpaired tau-sign selector; sign-selecting tables exist only by inserting a time-reversal-odd bias into the law, which is the missing strict source premise.",
            "next_honest_step": "Do not repeat local nearest-neighbour tau variational laws.  A next admissible proof-grade move must provide a concrete strict time-reversal-odd source term with a nonpremise coefficient sign and an exported P2721 polarity-coupling theorem; otherwise preserve the P2697-P2731 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2731/S1681 local time-arrow variational source-law no-go", "## P2731/S1681 local time-arrow variational source-law no-go\n\n`P2731/S1681` exhausts the finite class of translation-invariant nearest-neighbour time-arrow Hamiltonians `H=sum_i e(tau_i,tau_{i+1})` with `e in {-1,0,+1}`.  All 81 local tables are checked against all `2^12=4096` fields; the 9 time-reversal-even laws have no unpaired `tau`-sign ground selector.  Sign-selecting tables require a time-reversal-odd bias in the law itself, so no strict non-premise time-arrow source value, P2721 coupling, `lambda` fixing, `QW-2191` discharge, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2731/S1681 local time-arrow variational Ltotal guard", "## P2731/S1681 local time-arrow variational Ltotal guard\n\n`P2731/S1681` tests a genuine local variational source-law class for `tau`, but every admissible time-reversal-even nearest-neighbour Hamiltonian keeps the ground sector sign-paired.  Odd sign-selecting tables merely import the missing signed source coefficient, so this audit adds no strict variational source term and cannot promote `L_total`.\n")
    append_once(AGENTS, "Current local time-arrow variational source-law no-go guardrail (P2731/S1681, 2026-06-14)", "## Current local time-arrow variational source-law no-go guardrail (P2731/S1681, 2026-06-14)\n\n- P2731 exhausts all 81 translation-invariant nearest-neighbour local tau energy tables `e(tau_i,tau_{i+1}) in {-1,0,+1}` across all `2^12=4096` Z12 tau fields.\n- The 9 time-reversal-even laws have no unpaired tau-sign ground selector; sign-selecting tables exist only by inserting a time-reversal-odd bias into the law, which is the missing source premise.\n- Do not replay local nearest-neighbour tau variational laws as `QW-2191` discharge, selector closure, pair12 strict-core upgrade, role transfer, bridge closure, `L_total`, or ToE.  A next admissible move must export a concrete strict time-reversal-odd source term with a nonpremise coefficient sign and a P2721 coupling theorem, or preserve the P2697-P2731 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
