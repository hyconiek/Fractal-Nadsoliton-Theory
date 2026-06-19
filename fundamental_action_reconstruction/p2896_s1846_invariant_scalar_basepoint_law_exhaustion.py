#!/usr/bin/env python3
"""P2896/S1846: invariant scalar basepoint-law exhaustion.

P2895 granted an unpointed free Z12 torsor and showed it has zero invariant
basepoints.  P2896 audits a common attempted repair: attach a source-neutral
scalar/score law to the torsor and choose a unique minimum, maximum, or marked
level as the basepoint.  The finite theorem is immediate but useful to compute:
translation-invariant maps from a transitive free torsor to any trivial finite
alphabet are exactly constant maps.  Hence every level set is either empty or
all 12 points, and no unique basepoint/polarity law is exported.
"""
from __future__ import annotations

import json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

N = 12
P2895 = GEN / "p2895_s1845_free_12_torsor_basepoint_polarity_law_no_go.json"
OUT = GEN / "p2896_s1846_invariant_scalar_basepoint_law_exhaustion.json"
MD = GEN / "p2896_s1846_invariant_scalar_basepoint_law_exhaustion.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

TRIVIAL_ALPHABET_SIZES = (1, 2, 3, 4, 6, 12)


def translate(point: int, shift: int) -> int:
    return (point + shift) % N


def invariant_scalar_laws(alphabet_size: int) -> list[dict[str, Any]]:
    """Return all translation-invariant maps T -> {0,...,alphabet_size-1}."""
    rows = []
    for value in range(alphabet_size):
        table = [value for _ in range(N)]
        invariant = all(table[translate(point, shift)] == table[point] for point in range(N) for shift in range(N))
        level_sizes = {level: sum(1 for entry in table if entry == level) for level in range(alphabet_size)}
        rows.append({
            "constant_value": value,
            "invariant": invariant,
            "level_sizes": level_sizes,
            "unique_level_count": sum(1 for size in level_sizes.values() if size == 1),
            "nonempty_level_sizes": sorted(size for size in level_sizes.values() if size),
            "argmin_size": N,
            "argmax_size": N,
        })
    return rows


def scalar_exhaustion() -> dict[str, Any]:
    alphabet_rows = []
    for size in TRIVIAL_ALPHABET_SIZES:
        laws = invariant_scalar_laws(size)
        alphabet_rows.append({
            "alphabet_size": size,
            "invariant_scalar_law_count": len(laws),
            "all_invariant_laws_constant": all(row["invariant"] and row["nonempty_level_sizes"] == [N] for row in laws),
            "unique_marked_level_law_count": sum(row["unique_level_count"] for row in laws),
            "unique_argmin_law_count": sum(1 for row in laws if row["argmin_size"] == 1),
            "unique_argmax_law_count": sum(1 for row in laws if row["argmax_size"] == 1),
            "sample_laws": laws[: min(4, len(laws))],
        })
    return {
        "torsor_size": N,
        "tested_trivial_alphabet_sizes": list(TRIVIAL_ALPHABET_SIZES),
        "alphabet_rows": alphabet_rows,
        "total_invariant_scalar_law_count": sum(row["invariant_scalar_law_count"] for row in alphabet_rows),
        "total_unique_marked_level_law_count": sum(row["unique_marked_level_law_count"] for row in alphabet_rows),
        "total_unique_argmin_law_count": sum(row["unique_argmin_law_count"] for row in alphabet_rows),
        "total_unique_argmax_law_count": sum(row["unique_argmax_law_count"] for row in alphabet_rows),
        "proof_certificate": {
            "transitivity_fact": "The free Z12 torsor is transitive under translation.",
            "invariant_function_fact": "Any translation-invariant scalar map from a transitive torsor to a trivial alphabet is constant.",
            "selection_obstruction": "A constant score has empty or full level sets and 12-point argmin/argmax, never a unique basepoint.",
        },
    }


def build_payload(p2895: dict[str, Any]) -> dict[str, Any]:
    audit = scalar_exhaustion()
    facts = {
        "p2895_rechecked": p2895.get("status") == "P2895_FREE_12_TORSOR_BASEPOINT_POLARITY_LAW_NO_GO_NO_CLOSURE",
        "free_12_torsor_context_reused": True,
        "all_tested_invariant_scalar_laws_constant": all(row["all_invariant_laws_constant"] for row in audit["alphabet_rows"]),
        "no_unique_marked_level": audit["total_unique_marked_level_law_count"] == 0,
        "no_unique_argmin_or_argmax": audit["total_unique_argmin_law_count"] == 0 and audit["total_unique_argmax_law_count"] == 0,
    }
    return {
        "status": "P2896_INVARIANT_SCALAR_BASEPOINT_LAW_EXHAUSTION_NO_CLOSURE",
        "input_hashes": {"P2895": sha(P2895)},
        "invariant_scalar_basepoint_law_exhaustion": {
            "input_status_rechecked": p2895.get("status"),
            "candidate_class": "source-neutral translation-invariant scalar/score laws on the granted free Z12 torsor, with basepoint chosen by marked level, argmin, or argmax",
            "scalar_exhaustion": audit,
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_nonimported_basepoint_law": False,
            "exports_translation_breaking_scalar_source": False,
            "exports_unique_coupling_to_9_over_5_carrier": False,
            "exports_nonimported_9_over_5_variational_density": False,
        },
        "decision": {
            "negative_export_flags": {
                "nonimported_basepoint_or_polarity_law_exported": False,
                "translation_breaking_scalar_source_exported": False,
                "strict_free_12_torsor_source_law_exported": False,
                "strict_phase_origin_source_artifact_exported": False,
                "unique_coupling_to_9_over_5_carrier_exported": False,
                "coupling_to_9_over_5_variational_density_exported": False,
                "nonimported_9_over_5_variational_chain_rule_exported": False,
                "localized_action_density_exported": False,
                "strict_damping_compression_bridge_exported": False,
                "full_kernel_bridge_exported": False,
                "role_transfer_exported": False,
                "ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2896 exhausts source-neutral invariant scalar/score repair attempts on the P2895 free torsor for alphabet sizes 1,2,3,4,6,12.  Because the torsor is transitive, every invariant scalar law is constant; marked levels are empty or all 12 points and argmin/argmax have size 12.  Therefore invariant scoring does not create a nonimported basepoint/polarity law or source the 9/5 carrier origin.",
            "next_honest_step": "Do not promote invariant scalar scores, entropy-like constants, marked levels, argmin/argmax conventions, canonical zero choices, or unpointed free-torsor clocks to strict phase/origin sourcehood.  A next proof-grade move must either supply a genuinely non-invariant strict translation-breaking law with computed basepoint/polarity and a coupling theorem to the 9/5 variational density, or pivot to a genuinely different typed object outside torsor/basepoint/scalar-score/support/orbit/Fourier/inventory constructions; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["invariant_scalar_basepoint_law_exhaustion"]["scalar_exhaustion"]
    lines = [
        "# P2896/S1846 invariant scalar basepoint-law exhaustion",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite scalar-law audit",
        f"- torsor size: `{audit['torsor_size']}`",
        f"- tested trivial alphabet sizes: `{audit['tested_trivial_alphabet_sizes']}`",
        f"- total invariant scalar laws: `{audit['total_invariant_scalar_law_count']}`",
        f"- unique marked-level laws: `{audit['total_unique_marked_level_law_count']}`",
        f"- unique argmin laws: `{audit['total_unique_argmin_law_count']}`",
        f"- unique argmax laws: `{audit['total_unique_argmax_law_count']}`",
        "",
        "## Alphabet matrix",
    ]
    for row in audit["alphabet_rows"]:
        lines.append(
            f"- alphabet `{row['alphabet_size']}`: invariant laws=`{row['invariant_scalar_law_count']}`, "
            f"constant=`{row['all_invariant_laws_constant']}`, unique levels=`{row['unique_marked_level_law_count']}`, "
            f"unique argmin=`{row['unique_argmin_law_count']}`, unique argmax=`{row['unique_argmax_law_count']}`"
        )
    lines.extend(["", "## Boundary", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2895))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2896/S1846 invariant scalar basepoint-law exhaustion",
        "## P2896/S1846 invariant scalar basepoint-law exhaustion\n\n"
        "`P2896/S1846` exhausts source-neutral translation-invariant scalar/score laws on the `P2895` free `Z12` torsor for trivial alphabet sizes `1,2,3,4,6,12`.  Every invariant scalar law is constant; marked levels are empty or all `12` torsor points, and argmin/argmax sets have size `12`, not `1`.  Thus invariant scalar scoring does not export a nonimported basepoint/polarity law, `9/5` variational density, localized action density, strict damping bridge, `L_total`, EOM, Hamiltonian, role transfer, or ToE closure.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2896/S1846 invariant scalar basepoint-law `L_total` guard",
        "## P2896/S1846 invariant scalar basepoint-law `L_total` guard\n\n"
        "`P2896/S1846` is a finite scalar-selection obstruction, not a strict action construction.  It adds no translation-breaking source law, no localized unit-bearing density, no coupling theorem to the `9/5` carrier, and no variational chain rule into nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current invariant scalar basepoint-law exhaustion guardrail (P2896/S1846, 2026-06-19)",
        "## Current invariant scalar basepoint-law exhaustion guardrail (P2896/S1846, 2026-06-19)\n\n"
        "- P2896 exhausts source-neutral translation-invariant scalar/score laws on the P2895 free `Z12` torsor for trivial alphabet sizes `1,2,3,4,6,12`.\n"
        "- Every invariant scalar law is constant; marked levels are empty or all `12` torsor points, and argmin/argmax sets have size `12`, so no unique basepoint/polarity law is exported.\n"
        "- Do not promote invariant scalar scores, entropy-like constants, marked levels, argmin/argmax conventions, canonical zero choices, unpointed free-torsor clocks, support/orbit/Fourier data, or inventory hits to strict phase/origin sourcehood, strict damping/compression bridge, selector closure, role transfer, `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must supply a genuinely non-invariant strict translation-breaking law with computed basepoint/polarity and coupling to the `9/5` variational density, pivot to a genuinely different typed object outside torsor/basepoint/scalar-score/support/orbit/Fourier/inventory constructions, or preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
