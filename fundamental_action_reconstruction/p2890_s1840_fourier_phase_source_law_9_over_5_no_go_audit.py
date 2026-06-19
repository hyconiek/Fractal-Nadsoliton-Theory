#!/usr/bin/env python3
"""P2890/S1840: Fourier phase-source law 9/5 no-go audit.

P2889 showed that the P2888 9/5 carriers form free Z12 translation orbits, so a
translation-neutral law can select at most an orbit.  This packet tests the next
honest source candidate: a finite Fourier/phase law on the carrier density.

The audit regenerates the P2888 target triples, quotients by translation as in
P2889, and computes exact cyclic autocorrelation signatures.  These signatures
are the phase-blind Fourier power data: they are translation invariant and do
not choose an embedded representative.  The packet also records the phaseful
character boundary: a nonzero Fourier phase would rotate under translation and
therefore requires an external phase/origin pin before it can select a unique
representative.
"""
from __future__ import annotations

import json
from collections import Counter
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2889_s1839_translation_orbit_source_law_9_over_5_no_go_audit import N, P2888, orbit, regenerate_p2888_target_triples

P2889 = GEN / "p2889_s1839_translation_orbit_source_law_9_over_5_no_go_audit.json"
OUT = GEN / "p2890_s1840_fourier_phase_source_law_9_over_5_no_go_audit.json"
MD = GEN / "p2890_s1840_fourier_phase_source_law_9_over_5_no_go_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def orbit_representatives() -> list[Any]:
    unseen = set(regenerate_p2888_target_triples())
    reps = []
    while unseen:
        current = min(unseen)
        orb = orbit(current)
        reps.append(min(orb))
        unseen.difference_update(orb)
    return reps


def density_vector(triple: Any) -> tuple[int, ...]:
    values = [0] * N
    for node, value in triple.density_by_node:
        values[node] = value
    return tuple(values)


def cyclic_autocorrelation(values: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sum(values[j] * values[(j + shift) % N] for j in range(N)) for shift in range(N))


def phase_blind_power_audit(reps: list[Any]) -> dict[str, Any]:
    signatures: dict[tuple[int, ...], list[Any]] = {}
    for rep in reps:
        signatures.setdefault(cyclic_autocorrelation(density_vector(rep)), []).append(rep)
    class_sizes = sorted(len(members) for members in signatures.values())
    multiplicity_histogram = Counter(class_sizes)
    sample_collision_classes = []
    for signature, members in signatures.items():
        if len(members) > 1 and len(sample_collision_classes) < 6:
            sample_collision_classes.append(
                {
                    "autocorrelation_signature": list(signature),
                    "orbit_class_count": len(members),
                    "sample_representatives": [member.primitive_record() for member in members[:3]],
                }
            )
    return {
        "translation_orbit_representative_count": len(reps),
        "phase_blind_autocorrelation_signature_count": len(signatures),
        "phase_blind_signature_multiplicity_histogram": {str(size): count for size, count in sorted(multiplicity_histogram.items())},
        "phase_blind_signatures_unique_for_all_orbits": len(signatures) == len(reps),
        "sample_collision_classes": sample_collision_classes,
    }


def character_phase_boundary() -> dict[str, Any]:
    rows = []
    for k in range(1, N):
        orbit_size = N // gcd(k, N)
        phase_stabilizer_size = gcd(k, N)
        rows.append(
            {
                "character_mode_k": k,
                "phase_orbit_size_under_translation": orbit_size,
                "phase_stabilizer_size": phase_stabilizer_size,
                "needs_external_phase_pin_for_representative": orbit_size > 1,
            }
        )
    return {
        "nontrivial_character_mode_count": N - 1,
        "character_rows": rows,
        "all_nontrivial_characters_need_phase_pin_or_have_nontrivial_translation_phase_orbit": all(row["needs_external_phase_pin_for_representative"] for row in rows),
    }


def gcd(a: int, b: int) -> int:
    while b:
        a, b = b, a % b
    return abs(a)


def build_payload(p2889: dict[str, Any]) -> dict[str, Any]:
    reps = orbit_representatives()
    power = phase_blind_power_audit(reps)
    phase = character_phase_boundary()
    facts = {
        "p2889_rechecked": p2889.get("status") == "P2889_TRANSLATION_ORBIT_SOURCE_LAW_9_OVER_5_NO_GO_AUDIT_NO_CLOSURE",
        "p2889_orbit_count_reproduced": power["translation_orbit_representative_count"] == p2889["translation_orbit_source_law_9_over_5_no_go_audit"]["orbit_audit"]["translation_orbit_count"] == 50,
        "phase_blind_power_not_unique_for_all_orbits": not power["phase_blind_signatures_unique_for_all_orbits"],
        "phaseful_characters_need_external_phase_pin": phase["all_nontrivial_characters_need_phase_pin_or_have_nontrivial_translation_phase_orbit"],
        "strict_phase_source_law_missing": True,
    }
    return {
        "status": "P2890_FOURIER_PHASE_SOURCE_LAW_9_OVER_5_NO_GO_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2888": sha(P2888), "P2889": sha(P2889)},
        "fourier_phase_source_law_9_over_5_no_go_audit": {
            "input_status_rechecked": p2889.get("status"),
            "candidate_class": "phase-blind Fourier power/autocorrelation invariants and phaseful Z12 characters on P2888/P2889 9/5 carrier orbits",
            "phase_blind_power_audit": power,
            "phaseful_character_boundary": phase,
            "proof_certificate": {
                "phase_blind_result": "Autocorrelation/Fourier-power data compresses 50 translation orbits to 29 signatures, so it neither chooses embedded representatives nor uniquely classifies all carrier orbits.",
                "phaseful_result": "Every nontrivial character phase transforms under translation; using it to select a representative requires a phase/origin pin not exported by the current artifacts.",
                "sourcehood_obstruction": "Fourier power is too invariant, while Fourier phase is too origin-dependent unless a new strict translation-breaking phase source is supplied.",
            },
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_strict_translation_breaking_phase_source": False,
            "exports_unique_9_over_5_carrier_orbit_or_representative": False,
            "exports_nonimported_9_over_5_variational_chain_rule": False,
            "exports_strict_damping_compression_bridge": False,
        },
        "decision": {
            "negative_export_flags": {
                "strict_translation_breaking_phase_source_exported": False,
                "unique_9_over_5_carrier_orbit_exported": False,
                "unique_embedded_representative_exported": False,
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
            "reason": "P2890 tests the natural Fourier/phase follow-up to P2889.  Phase-blind autocorrelation/Fourier-power invariants are source-neutral but collapse 50 carrier orbits to 29 signatures; phaseful characters can distinguish representatives only after an external phase/origin pin is supplied.  No strict translation-breaking phase source or variational 9/5 density is exported.",
            "next_honest_step": "Do not replay Fourier-power signatures, Fourier-character phase pins, P2888/P2889 carrier representatives, bounded coefficients, or C12-neutral unit measures as strict sourcehood.  A next proof-grade move must supply one explicit strict phase/origin source artifact with a nonconventional sign/phase and a coupling theorem to the 9/5 variational density, or pivot to a genuinely different typed object; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["fourier_phase_source_law_9_over_5_no_go_audit"]
    power = audit["phase_blind_power_audit"]
    phase = audit["phaseful_character_boundary"]
    lines = [
        "# P2890/S1840 Fourier phase-source law 9/5 no-go audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite Fourier/power audit",
        f"- translation orbit representatives: `{power['translation_orbit_representative_count']}`",
        f"- phase-blind autocorrelation signatures: `{power['phase_blind_autocorrelation_signature_count']}`",
        f"- phase-blind multiplicity histogram: `{power['phase_blind_signature_multiplicity_histogram']}`",
        f"- nontrivial character modes checked: `{phase['nontrivial_character_mode_count']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2889))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2890/S1840 Fourier phase-source law 9/5 no-go audit",
        "## P2890/S1840 Fourier phase-source law 9/5 no-go audit\n\n"
        "`P2890/S1840` tests the natural Fourier/phase follow-up to `P2889`.  Exact cyclic autocorrelation/Fourier-power data is translation-neutral but compresses the `50` carrier orbits to `29` signatures, so it does not uniquely classify all orbits and cannot choose embedded representatives.  Phaseful `Z12` characters can distinguish representatives only after an external phase/origin pin is supplied.  No strict translation-breaking phase source, nonimported `9/5` variational density, localized action density, strict damping bridge, `L_total`, EOM, Hamiltonian, role transfer, or ToE closure follows.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2890/S1840 Fourier phase-source law `L_total` guard",
        "## P2890/S1840 Fourier phase-source law `L_total` guard\n\n"
        "`P2890/S1840` adds no strict action term.  Fourier-power/autocorrelation data is too invariant, while phaseful character data requires an external phase pin; neither supplies a unit-bearing localized action density or variational chain rule into nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current Fourier phase-source law 9/5 no-go guardrail (P2890/S1840, 2026-06-19)",
        "## Current Fourier phase-source law 9/5 no-go guardrail (P2890/S1840, 2026-06-19)\n\n"
        "- P2890 tests Fourier-power/autocorrelation invariants and phaseful `Z12` character data on the P2888/P2889 `9/5` carrier orbits.\n"
        "- Phase-blind power data compresses `50` carrier orbits to `29` signatures and cannot select embedded representatives; phaseful characters need an external phase/origin pin before they can select a representative.\n"
        "- Do not promote Fourier-power signatures, Fourier-character phase pins, P2888/P2889 carrier representatives, bounded coefficients, or `C12`-neutral unit measures to strict damping/compression bridge, selector closure, role transfer, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must supply one explicit strict phase/origin source artifact with nonconventional sign/phase and coupling theorem to the `9/5` variational density, pivot to a genuinely different typed object, or preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
