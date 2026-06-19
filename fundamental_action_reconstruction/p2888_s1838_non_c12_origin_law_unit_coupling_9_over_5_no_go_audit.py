#!/usr/bin/env python3
"""P2888/S1838: non-C12 origin-law/unit-coupling 9/5 no-go audit.

P2887 showed that C12-neutral unit measures and ternary localized densities are
only global and cannot produce 9/5.  This packet tests the next honest
candidate: allow non-C12-neutral finite supports, but demand that any origin is
selected by a source-neutral intrinsic rule rather than by an endpoint label.

The finite rule audited here is deliberately explicit and bounded.  For every
nonempty binary support S in Z12, compute each point's cyclic distance profile
inside S.  A point is intrinsically selectable only if its profile is unique in
S.  This catches the strongest source-neutral information available from the
unlabeled support geometry without importing a chart origin.  The packet then
enumerates bounded integer localized densities on the selected support with
values {-2,-1,0,1,2}; unit coupling is the uniform support average.

Result: 9/5 is arithmetically representable in this non-C12 family, but every
record still depends on a non-C12 support/origin and a bounded coefficient
assignment.  The intrinsic distance-profile rule is equivariant, not a strict
source law for which support/origin/density the nadsoliton must choose.  Thus it
is a representation audit, not a strict unit-bearing 9/5 source.
"""
from __future__ import annotations

import itertools
import json
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2887 = GEN / "p2887_s1837_c12_unit_measure_localized_density_9_over_5_no_go_audit.json"
OUT = GEN / "p2888_s1838_non_c12_origin_law_unit_coupling_9_over_5_no_go_audit.json"
MD = GEN / "p2888_s1838_non_c12_origin_law_unit_coupling_9_over_5_no_go_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

N = 12
TARGET = Fraction(9, 5)
DENSITY_VALUES = (-2, -1, 0, 1, 2)


def cyclic_distance(a: int, b: int) -> int:
    d = abs(a - b) % N
    return min(d, N - d)


def support_points(mask: int) -> tuple[int, ...]:
    return tuple(i for i in range(N) if mask & (1 << i))


def distance_profile(point: int, points: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted(cyclic_distance(point, q) for q in points if q != point))


def intrinsic_unique_origins(points: tuple[int, ...]) -> list[int]:
    profiles = {p: distance_profile(p, points) for p in points}
    counts: dict[tuple[int, ...], int] = {}
    for profile in profiles.values():
        counts[profile] = counts.get(profile, 0) + 1
    return [p for p, profile in profiles.items() if counts[profile] == 1]


def enumerate_support_origin_records() -> dict[str, Any]:
    nonempty_count = 0
    supports_with_unique_origin = 0
    unique_origin_records: list[dict[str, Any]] = []
    support_size_histogram: dict[str, int] = {}
    unique_origin_size_histogram: dict[str, int] = {}
    for mask in range(1, 1 << N):
        points = support_points(mask)
        nonempty_count += 1
        support_size_histogram[str(len(points))] = support_size_histogram.get(str(len(points)), 0) + 1
        origins = intrinsic_unique_origins(points)
        if len(origins) == 1:
            supports_with_unique_origin += 1
            unique_origin_size_histogram[str(len(points))] = unique_origin_size_histogram.get(str(len(points)), 0) + 1
            if len(unique_origin_records) < 12:
                unique_origin_records.append({"support": list(points), "selected_origin": origins[0], "profile": list(distance_profile(origins[0], points))})
    return {
        "nonempty_binary_support_count": nonempty_count,
        "supports_with_unique_intrinsic_distance_profile_origin": supports_with_unique_origin,
        "support_size_histogram": support_size_histogram,
        "unique_origin_size_histogram": unique_origin_size_histogram,
        "sample_unique_origin_records": unique_origin_records,
    }


def count_assignments_with_sum(k: int, target_sum: int) -> int:
    counts: dict[int, int] = {0: 1}
    for _ in range(k):
        nxt: dict[int, int] = {}
        for partial, count in counts.items():
            for value in DENSITY_VALUES:
                nxt[partial + value] = nxt.get(partial + value, 0) + count
        counts = nxt
    return counts.get(target_sum, 0)


def sample_assignment_with_sum(k: int, target_sum: int) -> tuple[int, ...] | None:
    values: list[int] = []
    remaining = target_sum
    for slot in range(k):
        rest = k - slot - 1
        for value in DENSITY_VALUES:
            low = min(DENSITY_VALUES) * rest
            high = max(DENSITY_VALUES) * rest
            if low <= remaining - value <= high:
                values.append(value)
                remaining -= value
                break
        else:
            return None
    return tuple(values) if remaining == 0 else None


def enumerate_unit_couplings() -> dict[str, Any]:
    candidate_density_assignments = 0
    target_records = 0
    target_supports: set[tuple[int, ...]] = set()
    target_samples: list[dict[str, Any]] = []
    target_by_size: dict[str, int] = {}
    for mask in range(1, 1 << N):
        points = support_points(mask)
        origins = intrinsic_unique_origins(points)
        if len(origins) != 1:
            continue
        k = len(points)
        candidate_density_assignments += len(DENSITY_VALUES) ** k
        if (TARGET.numerator * k) % TARGET.denominator != 0:
            continue
        target_sum = TARGET.numerator * k // TARGET.denominator
        assignment_count = count_assignments_with_sum(k, target_sum)
        if assignment_count == 0:
            continue
        target_records += assignment_count
        target_supports.add(points)
        target_by_size[str(k)] = target_by_size.get(str(k), 0) + assignment_count
        if len(target_samples) < 12:
            sample = sample_assignment_with_sum(k, target_sum)
            target_samples.append({
                "support": list(points),
                "selected_origin": origins[0],
                "density_values_on_support_order": list(sample or ()),
                "unit_average_coupling": str(TARGET),
            })
    return {
        "bounded_density_value_set": list(DENSITY_VALUES),
        "candidate_density_assignments_on_unique_origin_supports": candidate_density_assignments,
        "target_9_over_5_record_count": target_records,
        "target_9_over_5_support_count": len(target_supports),
        "target_9_over_5_record_count_by_support_size": target_by_size,
        "sample_target_records": target_samples,
    }

def build_payload(p2887: dict[str, Any]) -> dict[str, Any]:
    support_audit = enumerate_support_origin_records()
    coupling_audit = enumerate_unit_couplings()
    facts = {
        "p2887_rechecked": p2887.get("status") == "P2887_C12_UNIT_MEASURE_LOCALIZED_DENSITY_9_OVER_5_NO_GO_AUDIT_NO_CLOSURE",
        "all_nonempty_binary_supports_checked": support_audit["nonempty_binary_support_count"] == 2**N - 1,
        "non_c12_unique_origins_exist": support_audit["supports_with_unique_intrinsic_distance_profile_origin"] > 0,
        "target_9_over_5_representable_in_bounded_non_c12_family": coupling_audit["target_9_over_5_record_count"] > 0,
        "target_remains_nonunique": coupling_audit["target_9_over_5_record_count"] > 1 and coupling_audit["target_9_over_5_support_count"] > 1,
        "strict_source_law_missing": True,
    }
    return {
        "status": "P2888_NON_C12_ORIGIN_LAW_UNIT_COUPLING_9_OVER_5_NO_GO_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2887": sha(P2887)},
        "non_c12_origin_law_unit_coupling_9_over_5_no_go_audit": {
            "input_status_rechecked": p2887.get("status"),
            "candidate_class": "non-C12 binary supports with intrinsic distance-profile origin selector and bounded integer localized densities {-2,-1,0,1,2}",
            "support_origin_audit": support_audit,
            "unit_coupling_audit": coupling_audit,
            "proof_certificate": {
                "positive_representation_result": "Some non-C12 supports have a unique intrinsic distance-profile point, and bounded integer densities on such supports can average to 9/5.",
                "negative_source_result": "The finite family supplies many support/density carriers and no strict law selecting one support, one origin, one unit measure, or one density assignment from nadsoliton data.",
                "import_boundary": "Non-C12 support choice is already an origin/source input; the 9/5 target records additionally depend on coefficient assignments.  Representation is therefore not sourcehood.",
            },
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_non_c12_strict_origin_source_law": False,
            "accepted_as_unit_bearing_9_over_5_coupling_source": False,
            "exports_selected_support_origin_density_triple": False,
            "exports_variational_chain_rule_to_ltotal": False,
            "exports_strict_damping_compression_bridge": False,
        },
        "decision": {
            "negative_export_flags": {
                "non_c12_strict_origin_source_law_exported": False,
                "unit_bearing_9_over_5_coupling_source_exported": False,
                "selected_support_origin_density_triple_exported": False,
                "localized_action_density_exported": False,
                "variational_chain_rule_exported": False,
                "strict_damping_compression_bridge_exported": False,
                "full_kernel_bridge_exported": False,
                "role_transfer_exported": False,
                "ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2888 allows non-C12 supports and audits an intrinsic distance-profile origin selector plus bounded integer unit couplings.  The target 9/5 is representable, but nonunique and still depends on an unsourced non-C12 support/origin/density triple rather than an exported strict source law.",
            "next_honest_step": "Do not replay non-C12 support choice, intrinsic distance-profile origin selection, bounded density coefficient assignment, C12-neutral unit measures, ratio algebra, or scalar Euler transmission as sourcehood.  A next proof-grade move must either supply one explicit strict source law selecting a unique support-origin-density triple with a nonimported 9/5 coupling and variational chain rule, or pivot to a genuinely different typed object; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["non_c12_origin_law_unit_coupling_9_over_5_no_go_audit"]
    support = audit["support_origin_audit"]
    coupling = audit["unit_coupling_audit"]
    lines = [
        "# P2888/S1838 non-C12 origin-law/unit-coupling 9/5 no-go audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite origin/coupling audit",
        f"- nonempty binary supports checked: `{support['nonempty_binary_support_count']}`",
        f"- supports with unique intrinsic distance-profile origin: `{support['supports_with_unique_intrinsic_distance_profile_origin']}`",
        f"- bounded density assignments on unique-origin supports: `{coupling['candidate_density_assignments_on_unique_origin_supports']}`",
        f"- target 9/5 records: `{coupling['target_9_over_5_record_count']}`",
        f"- target 9/5 supports: `{coupling['target_9_over_5_support_count']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload(read_json(P2887))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2888/S1838 non-C12 origin-law/unit-coupling 9/5 no-go audit",
        "## P2888/S1838 non-C12 origin-law/unit-coupling 9/5 no-go audit\n\n"
        "`P2888/S1838` tests the post-`P2887` non-`C12` route by enumerating all `4095` nonempty binary supports, selecting any unique intrinsic distance-profile origin, and then enumerating bounded integer localized densities `{-2,-1,0,1,2}` on those selected supports.  The target `9/5` is representable, but it is nonunique and depends on an unsourced non-`C12` support/origin/density triple plus coefficient assignment.  Therefore the packet exports no strict origin/source law, no localized unit-bearing action density, no variational chain rule, no strict damping bridge, no `L_total`, EOM, Hamiltonian, role transfer, or ToE closure.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2888/S1838 non-C12 origin-law/unit-coupling `L_total` guard",
        "## P2888/S1838 non-C12 origin-law/unit-coupling `L_total` guard\n\n"
        "`P2888/S1838` adds no strict action term.  Non-`C12` distance-profile origins and bounded densities can carry `9/5`, but the audited carriers are not selected by a strict source law and do not provide a unit-bearing localized action density or variational chain rule into nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current non-C12 origin-law/unit-coupling 9/5 no-go guardrail (P2888/S1838, 2026-06-18)",
        "## Current non-C12 origin-law/unit-coupling 9/5 no-go guardrail (P2888/S1838, 2026-06-18)\n\n"
        "- P2888 tests the post-P2887 non-`C12` route by enumerating all nonempty binary supports, intrinsic distance-profile origin selectors, and bounded integer density assignments on selected supports.\n"
        "- The `9/5` unit coupling is representable, but nonunique and still depends on an unsourced non-`C12` support/origin/density triple plus coefficient assignment; representation is not strict sourcehood.\n"
        "- Do not promote non-`C12` support choice, distance-profile origin selection, bounded density coefficients, `C12`-neutral unit measures, ratio algebra, or scalar Euler transmission to strict damping/compression bridge, selector closure, role transfer, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must supply one explicit strict source law selecting a unique support-origin-density triple with nonimported `9/5` coupling and variational chain rule, pivot to a genuinely different typed object, or preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
