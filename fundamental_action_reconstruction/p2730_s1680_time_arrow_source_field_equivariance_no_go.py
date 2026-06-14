#!/usr/bin/env python3
"""P2730/S1680: time-arrow source-field equivariance no-go.

After P2729 showed that a fixed time-arrow sign tau would be useful only if it
were strictly sourced, P2730 tests the simplest concrete source class: all
Z12-internal signed time-arrow fields tau: Z12 -> {-1,+1}.  The audit exhausts
all 2^12 fields and checks whether translation/Aut-invariant data can export a
non-premise signed value rather than a paired convention.
"""
from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2730_s1680_time_arrow_source_field_equivariance_no_go.json"
MD = GEN / "p2730_s1680_time_arrow_source_field_equivariance_no_go.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

Z12 = 12
UNITS = (1, 5, 7, 11)
INPUTS = {
    "P2729_TIME_ARROW_AUDIT": GEN / "p2729_s1679_time_arrow_orientation_coupling_law_audit.json",
    "P2721_SIGN_TORSOR_COUPLING": GEN / "p2721_s1671_chiral_bispectrum_sign_to_orientation_torsor_coupling_audit.json",
}
NEGATIVE_EXPORT_FLAGS = [
    "strict_time_arrow_source_value_exported",
    "nonpremise_tau_sign_selected",
    "translation_gauge_safe_tau_exported",
    "aut_invariant_unpaired_tau_exported",
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


def shift(field: tuple[int, ...], k: int) -> tuple[int, ...]:
    return tuple(field[(i - k) % Z12] for i in range(Z12))


def aut_action(field: tuple[int, ...], unit: int) -> tuple[int, ...]:
    return tuple(field[(pow(unit, -1, Z12) * i) % Z12] for i in range(Z12))


def neg(field: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(-x for x in field)


def source_orbits() -> list[list[int]]:
    unseen = set(range(Z12))
    orbits = []
    while unseen:
        start = min(unseen)
        orbit = sorted({(u * start) % Z12 for u in UNITS})
        orbits.append(orbit)
        unseen -= set(orbit)
    return orbits


def audit_fields() -> dict[str, Any]:
    fields = list(itertools.product((-1, 1), repeat=Z12))
    translation_invariant = [f for f in fields if all(shift(f, k) == f for k in range(Z12))]
    aut_invariant = [f for f in fields if all(aut_action(f, u) == f for u in UNITS)]
    both_invariant = [f for f in translation_invariant if f in aut_invariant]
    sign_unpaired_aut = [f for f in aut_invariant if neg(f) not in aut_invariant]
    sign_unpaired_translation = [f for f in translation_invariant if neg(f) not in translation_invariant]
    nonconstant_aut = [f for f in aut_invariant if len(set(f)) > 1]
    signed_totals = sorted(set(sum(f) for f in fields))
    aut_signed_total_histogram = {str(total): sum(1 for f in aut_invariant if sum(f) == total) for total in signed_totals}
    translation_signed_total_histogram = {str(total): sum(1 for f in translation_invariant if sum(f) == total) for total in signed_totals}
    canonical_positive_attempts = [f for f in both_invariant if sum(f) > 0]
    canonical_negative_attempts = [f for f in both_invariant if sum(f) < 0]
    return {
        "field_count": len(fields),
        "source_orbits_under_aut_z12": source_orbits(),
        "aut_invariant_field_count": len(aut_invariant),
        "translation_invariant_field_count": len(translation_invariant),
        "translation_and_aut_invariant_field_count": len(both_invariant),
        "nonconstant_aut_invariant_field_count": len(nonconstant_aut),
        "sign_unpaired_aut_invariant_field_count": len(sign_unpaired_aut),
        "sign_unpaired_translation_invariant_field_count": len(sign_unpaired_translation),
        "aut_signed_total_histogram": aut_signed_total_histogram,
        "translation_signed_total_histogram": translation_signed_total_histogram,
        "translation_and_aut_invariant_fields": [list(f) for f in both_invariant],
        "canonical_positive_attempt_count": len(canonical_positive_attempts),
        "canonical_negative_attempt_count": len(canonical_negative_attempts),
        "obstruction": "The only translation-gauge-safe tau fields are the two constant fields +1 and -1.  They are exchanged by global time reversal, and no current strict artifact selects one sign or couples it to one P2721 polarity.  Aut-invariant nonconstant fields exist, but they retain source labels/orbits and remain paired by tau -> -tau.",
    }


def acceptance_matrix(audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "all_z12_tau_fields_enumerated": audit["field_count"] == 4096,
        "translation_gauge_safe_tau_fields_exist": audit["translation_invariant_field_count"] == 2,
        "translation_gauge_safe_tau_sign_unpaired": audit["sign_unpaired_translation_invariant_field_count"] > 0,
        "aut_invariant_tau_sign_unpaired": audit["sign_unpaired_aut_invariant_field_count"] > 0,
        "strict_nonpremise_tau_sign_selected": False,
        "p2721_polarity_coupling_theorem_exported": False,
    }
    missing = [k for k, v in facts.items() if not v]
    return {
        "facts": facts,
        "missing_criteria": missing,
        "accepted_as_strict_time_arrow_source": not missing,
        "blocker": "Finite source-field enumeration leaves only sign-paired tau candidates; selecting +tau over -tau would be an extra premise unless a new strict signed source law is exported.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    a = payload["time_arrow_source_field_audit"]
    lines = [
        "# P2730/S1680 time-arrow source-field equivariance no-go",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Exhaustive Z12 tau-field audit",
        f"- field_count={a['field_count']}",
        f"- source_orbits_under_aut_z12={a['source_orbits_under_aut_z12']}",
        f"- aut_invariant_field_count={a['aut_invariant_field_count']}",
        f"- translation_invariant_field_count={a['translation_invariant_field_count']}",
        f"- translation_and_aut_invariant_field_count={a['translation_and_aut_invariant_field_count']}",
        f"- sign_unpaired_aut_invariant_field_count={a['sign_unpaired_aut_invariant_field_count']}",
        f"- sign_unpaired_translation_invariant_field_count={a['sign_unpaired_translation_invariant_field_count']}",
        a["obstruction"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    inputs = {key: read_json(path) for key, path in INPUTS.items()}
    audit = audit_fields()
    acceptance = acceptance_matrix(audit)
    no_go = not acceptance["accepted_as_strict_time_arrow_source"]
    payload = {
        "status": "P2730_TIME_ARROW_SOURCE_FIELD_EQUIVARIANCE_NO_GO" if no_go else "P2730_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "input_statuses": {key: value.get("status") for key, value in inputs.items()},
        "audited_candidate_class": "all Z12-internal signed time-arrow fields tau: Z12 -> {-1,+1}, tested under translation gauge and Aut(Z12) actions",
        "time_arrow_source_field_audit": audit,
        "acceptance_matrix": acceptance,
        "decision": {
            "strict_time_arrow_source_value_exported": False,
            "nonpremise_tau_sign_selected": False,
            "time_arrow_p2721_polarity_coupling_exported": False,
            "strict_mechanism_fixing_lambda_exported": False,
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": "P2730 exhausts the finite class of Z12-internal signed tau fields.  Translation-safe candidates reduce to the paired constants +tau and -tau; Aut-invariant source-dependent candidates also remain paired by global sign.  No strict law selects one sign or couples it to the P2721 polarity pair.",
            "next_honest_step": "Do not repeat bare time-arrow field searches.  The next admissible proof-grade move must provide one exported strict law that fixes the sign of a time-orientation source and proves a polarity coupling to P2721; otherwise preserve the P2697-P2730 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2730/S1680 time-arrow source-field equivariance no-go", "## P2730/S1680 time-arrow source-field equivariance no-go\n\n`P2730/S1680` exhausts all `2^12=4096` Z12-internal signed time-arrow fields `tau: Z12 -> {-1,+1}`.  Translation-gauge-safe fields are only the paired constants `+tau` and `-tau`; Aut-invariant source-dependent fields exist but remain paired by global sign and do not select a P2721 polarity.  No strict non-premise time-arrow source value, `lambda` fixing, `QW-2191` discharge, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2730/S1680 time-arrow source-field Ltotal guard", "## P2730/S1680 time-arrow source-field Ltotal guard\n\n`P2730/S1680` finds no unpaired strict signed `tau` source among finite Z12 time-arrow fields.  Because selecting `+tau` over `-tau` remains an unsourced premise and no P2721 coupling theorem is exported, this audit adds no variational source term and cannot promote `L_total`.\n")
    append_once(AGENTS, "Current time-arrow source-field equivariance no-go guardrail (P2730/S1680, 2026-06-14)", "## Current time-arrow source-field equivariance no-go guardrail (P2730/S1680, 2026-06-14)\n\n- P2730 exhausts all `2^12=4096` Z12-internal signed time-arrow fields `tau: Z12 -> {-1,+1}` as the finite source class requested after P2729.\n- Translation-gauge-safe candidates reduce to the paired constants `+tau` and `-tau`; Aut-invariant source-dependent candidates also remain paired by global sign, so no non-premise strict time-arrow sign or P2721 polarity coupling is exported.\n- Do not replay bare time-arrow field searches as `QW-2191` discharge, selector closure, pair12 strict-core upgrade, role transfer, bridge closure, `L_total`, or ToE.  A next admissible move must export a strict law fixing the time-orientation sign and coupling it to P2721, or preserve the P2697-P2730 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
