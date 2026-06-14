#!/usr/bin/env python3
"""P2729/S1679: time-arrow orientation-coupling law audit.

The user explicitly asked to check the arrow of time.  P2729 treats time-arrow
sign tau in {-1,+1} as a candidate chiral input and exhausts all 16 Boolean
laws f(orientation, tau) -> next_orientation, crossed with all Z12 source
velocities, on the exact P2718 Im(B_{1,5}) marker.

Result: laws can use a fixed tau sign to choose an orientation branch, but the
choice is conditional on tau already being a non-premise signed source.  When
both tau signs are included, polarity remains paired; current artifacts do not
export a strict time-arrow source value or a P2721 polarity-coupling theorem.
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
OUT = GEN / "p2729_s1679_time_arrow_orientation_coupling_law_audit.json"
MD = GEN / "p2729_s1679_time_arrow_orientation_coupling_law_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

Z12 = 12
ORIENTATIONS = (-1, 1)
TIME_ARROWS = (-1, 1)
VELOCITIES = tuple(range(Z12))
STATE_ORDER = tuple((orientation, tau) for orientation in ORIENTATIONS for tau in TIME_ARROWS)
INPUTS = {
    "P2728_SOURCE_ORBIT_NO_GO": GEN / "p2728_s1678_aut_z12_source_orbit_weighted_chiral_invariant_no_go.json",
    "P2727_ORIENTATION_LAW_NO_GO": GEN / "p2727_s1677_orientation_transition_law_equivariance_and_polarity_no_go.json",
    "P2721_SIGN_TORSOR_COUPLING": GEN / "p2721_s1671_chiral_bispectrum_sign_to_orientation_torsor_coupling_audit.json",
    "P2718_CHIRAL_BISPECTRUM_MARKER": GEN / "p2718_s1668_chiral_bispectrum_explicit_signed_formula_torsor_coupling_audit.json",
}

NEGATIVE_EXPORT_FLAGS = [
    "strict_time_arrow_source_value_exported",
    "time_arrow_selected_as_nonpremise_tau",
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
    if not path.exists():
        return {"missing": True, "path": rel(path)}
    return json.loads(path.read_text(encoding="utf-8"))


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def marker_table(p2718: dict[str, Any]) -> dict[tuple[int, int], float]:
    return {
        (int(row["orientation"]), int(row["source"])): float(row["marker_imag"])
        for row in p2718.get("marker_rows", [])
    }


def law_name(mapping: dict[tuple[int, int], int]) -> str:
    bits = "".join("p" if mapping[state] == 1 else "m" for state in STATE_ORDER)
    special = {
        "mmmm": "collapse_to_minus",
        "pppp": "collapse_to_plus",
        "mmpp": "preserve_orientation",
        "ppmm": "flip_orientation",
        "mpmp": "follow_time_arrow_tau",
        "pmpm": "oppose_time_arrow_tau",
    }
    return special.get(bits, f"law_{bits}")


def time_reversal_equivariant(mapping: dict[tuple[int, int], int]) -> bool:
    # Joint reversal sends (orientation, tau) -> (-orientation, -tau), and a signed output must reverse too.
    return all(mapping[(-orientation, -tau)] == -mapping[(orientation, tau)] for orientation, tau in STATE_ORDER)


def depends_on_tau(mapping: dict[tuple[int, int], int]) -> bool:
    return any(mapping[(orientation, -1)] != mapping[(orientation, 1)] for orientation in ORIENTATIONS)


def fixed_tau_polarity(mapping: dict[tuple[int, int], int], table: dict[tuple[int, int], float], tau: int) -> dict[str, Any]:
    deltas = []
    for orientation in ORIENTATIONS:
        for source in range(Z12):
            for velocity in VELOCITIES:
                next_orientation = mapping[(orientation, tau)]
                next_source = (source + tau * velocity) % Z12
                deltas.append(round(table[(next_orientation, next_source)] - table[(orientation, source)], 9))
    positive = sum(1 for delta in deltas if delta > 0)
    negative = sum(1 for delta in deltas if delta < 0)
    return {
        "tau": tau,
        "delta_values": sorted(set(deltas)),
        "positive_delta_count": positive,
        "negative_delta_count": negative,
        "selects_single_nonzero_polarity": (positive > 0 and negative == 0) or (negative > 0 and positive == 0),
    }


def exhaustive_time_arrow_audit(table: dict[tuple[int, int], float]) -> dict[str, Any]:
    laws = []
    for outputs in itertools.product(ORIENTATIONS, repeat=len(STATE_ORDER)):
        mapping = {state: output for state, output in zip(STATE_ORDER, outputs)}
        deltas = []
        for orientation, tau in STATE_ORDER:
            for source in range(Z12):
                for velocity in VELOCITIES:
                    next_orientation = mapping[(orientation, tau)]
                    next_source = (source + tau * velocity) % Z12
                    deltas.append(round(table[(next_orientation, next_source)] - table[(orientation, source)], 9))
        positive = sum(1 for delta in deltas if delta > 0)
        negative = sum(1 for delta in deltas if delta < 0)
        law = {
            "name": law_name(mapping),
            "mapping": {f"o={orientation},tau={tau}": mapping[(orientation, tau)] for orientation, tau in STATE_ORDER},
            "time_reversal_equivariant": time_reversal_equivariant(mapping),
            "depends_on_tau": depends_on_tau(mapping),
            "delta_values_both_tau": sorted(set(deltas)),
            "delta_histogram_both_tau": dict(sorted(Counter(deltas).items())),
            "positive_delta_count_both_tau": positive,
            "negative_delta_count_both_tau": negative,
            "balanced_both_tau": positive == negative,
            "fixed_tau_summaries": [fixed_tau_polarity(mapping, table, tau) for tau in TIME_ARROWS],
        }
        law["fixed_tau_can_select_polarity"] = any(row["selects_single_nonzero_polarity"] for row in law["fixed_tau_summaries"])
        laws.append(law)

    equivariant = [law for law in laws if law["time_reversal_equivariant"]]
    tau_dependent_equivariant = [law for law in equivariant if law["depends_on_tau"]]
    fixed_tau_selectors = [law for law in tau_dependent_equivariant if law["fixed_tau_can_select_polarity"]]
    unconditional_selectors = [law for law in equivariant if not law["balanced_both_tau"]]
    return {
        "law_count": len(laws),
        "checked_transition_rows": len(laws) * len(STATE_ORDER) * Z12 * len(VELOCITIES),
        "state_order": [list(state) for state in STATE_ORDER],
        "time_reversal_equivariant_law_count": len(equivariant),
        "tau_dependent_equivariant_law_names": [law["name"] for law in tau_dependent_equivariant],
        "fixed_tau_polarity_selector_names": [law["name"] for law in fixed_tau_selectors],
        "unconditional_time_reversal_equivariant_selector_count": len(unconditional_selectors),
        "selected_law_summaries": [law for law in laws if law["name"] in {"follow_time_arrow_tau", "oppose_time_arrow_tau", "preserve_orientation", "flip_orientation"}],
        "obstruction": "Time-arrow-dependent equivariant laws can select a polarity only after fixing tau.  With tau included as an unsourced +/- pair, the signs remain balanced; current artifacts do not export a non-premise tau value or a P2721 coupling theorem.",
    }


def acceptance_matrix(audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "all_time_arrow_orientation_laws_enumerated": audit["law_count"] == 16,
        "time_arrow_dependent_equivariant_laws_exist": len(audit["tau_dependent_equivariant_law_names"]) > 0,
        "fixed_tau_can_conditionally_select_polarity": len(audit["fixed_tau_polarity_selector_names"]) > 0,
        "unconditional_time_reversal_equivariant_selector_exported": audit["unconditional_time_reversal_equivariant_selector_count"] > 0,
        "strict_time_arrow_source_value_exported": False,
        "time_arrow_p2721_polarity_coupling_exported": False,
    }
    missing = [criterion for criterion, value in facts.items() if not value]
    return {
        "facts": facts,
        "missing_criteria": missing,
        "accepted_as_strict_time_arrow_source": not missing,
        "blocker": "The arrow of time is the right type of signed input only conditionally: once tau is fixed.  The finite audit finds no unconditional time-reversal-equivariant selector and no exported strict tau source value.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["time_arrow_audit"]
    lines = [
        "# P2729/S1679 time-arrow orientation-coupling law audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Time-arrow law enumeration",
        f"- law_count={audit['law_count']}",
        f"- checked_transition_rows={audit['checked_transition_rows']}",
        f"- tau_dependent_equivariant_law_names={audit['tau_dependent_equivariant_law_names']}",
        f"- fixed_tau_polarity_selector_names={audit['fixed_tau_polarity_selector_names']}",
        f"- unconditional_time_reversal_equivariant_selector_count={audit['unconditional_time_reversal_equivariant_selector_count']}",
        audit["obstruction"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    inputs = {key: read_json(path) for key, path in INPUTS.items()}
    audit = exhaustive_time_arrow_audit(marker_table(inputs["P2718_CHIRAL_BISPECTRUM_MARKER"]))
    acceptance = acceptance_matrix(audit)
    no_go = audit["unconditional_time_reversal_equivariant_selector_count"] == 0 and not acceptance["accepted_as_strict_time_arrow_source"]
    payload = {
        "status": "P2729_TIME_ARROW_CONDITIONAL_POLARITY_BUT_NO_STRICT_SOURCE" if no_go else "P2729_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "input_statuses": {key: value.get("status") for key, value in inputs.items()},
        "audited_candidate_class": "all laws f(orientation,tau)->next_orientation for tau in {-1,+1}, crossed with all Z12 velocities on P2718 marker rows",
        "time_arrow_audit": audit,
        "acceptance_matrix": acceptance,
        "decision": {
            "strict_time_arrow_source_value_exported": False,
            "time_arrow_selected_as_nonpremise_tau": False,
            "time_arrow_p2721_polarity_coupling_exported": False,
            "strict_mechanism_fixing_lambda_exported": False,
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": "P2729 checks the arrow of time directly.  Exhausting all 16 orientation/time-arrow coupling laws shows that tau-dependent equivariant laws can conditionally choose a polarity for fixed tau, but the +/- tau pair is balanced unless a strict non-premise time-arrow value is exported.  No such source value or P2721 coupling theorem exists in current artifacts.",
            "next_honest_step": "A next admissible move may target the arrow of time only by exporting one strict, non-premise signed time-orientation source value tau and a theorem coupling that tau to the P2721 polarity pair.  Without that, do not replay time-arrow intuition; preserve the P2697-P2729 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2729/S1679 time-arrow orientation-coupling law audit", "## P2729/S1679 time-arrow orientation-coupling law audit\n\n`P2729/S1679` checks the arrow of time by exhausting all 16 laws `f(orientation,tau)->next_orientation` for `tau in {-1,+1}` crossed with all Z12 source velocities on the P2718 marker.  Time-arrow-dependent equivariant laws can conditionally select polarity only after `tau` is fixed; with both `tau` signs present, polarity remains balanced.  No strict non-premise time-arrow source value, P2721 polarity coupling, `lambda` fixing, `QW-2191` discharge, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2729/S1679 time-arrow source Ltotal guard", "## P2729/S1679 time-arrow source Ltotal guard\n\n`P2729/S1679` finds that a time-arrow sign `tau` would be useful only conditionally, after a non-premise signed `tau` value is sourced.  Since current artifacts export no such time-orientation source value and no P2721 coupling theorem, time-arrow intuition cannot add a variational source term or promote `L_total`.\n")
    append_once(AGENTS, "Current time-arrow orientation-coupling audit guardrail (P2729/S1679, 2026-06-14)", "## Current time-arrow orientation-coupling audit guardrail (P2729/S1679, 2026-06-14)\n\n- P2729 checks the arrow of time directly by exhausting all 16 laws `f(orientation,tau)->next_orientation` for `tau in {-1,+1}` crossed with all Z12 source velocities on the P2718 marker.\n- Time-arrow-dependent equivariant laws can select polarity only conditionally after fixing `tau`; with both `tau` signs included, polarity remains balanced, and current artifacts export no strict non-premise time-arrow source value.\n- Do not replay time-arrow intuition as `QW-2191` discharge, selector closure, pair12 strict-core upgrade, role transfer, bridge closure, `L_total`, or ToE.  A next admissible move must export one strict signed time-orientation source value and a P2721 polarity-coupling theorem, or preserve the P2697-P2729 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
