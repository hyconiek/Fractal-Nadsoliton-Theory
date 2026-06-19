#!/usr/bin/env python3
"""P2892/S1842: post phase/origin inventory state-map no-new-live-frontier certificate.

P2891 closed the current generated-artifact intake gate for the exact missing
post-P2890 object: a strict phase/origin source artifact with nonconventional
sign/phase plus a coupling theorem to the 9/5 variational density.  P2892 is the
honest next step when no such new artifact has been supplied: reconcile the
recent finite lane and check whether any current JSON artifact exports an unlock
that would justify replaying carrier/orbit/Fourier/inventory moves.

This is a current-artifact certificate, not a future impossibility theorem.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2888 = GEN / "p2888_s1838_non_c12_origin_law_unit_coupling_9_over_5_no_go_audit.json"
P2889 = GEN / "p2889_s1839_translation_orbit_source_law_9_over_5_no_go_audit.json"
P2890 = GEN / "p2890_s1840_fourier_phase_source_law_9_over_5_no_go_audit.json"
P2891 = GEN / "p2891_s1841_strict_phase_origin_source_artifact_inventory_no_go_audit.json"
OUT = GEN / "p2892_s1842_post_phase_origin_inventory_state_map_no_new_live_frontier_certificate.json"
MD = GEN / "p2892_s1842_post_phase_origin_inventory_state_map_no_new_live_frontier_certificate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {"P2888": P2888, "P2889": P2889, "P2890": P2890, "P2891": P2891}

UNLOCK_TERMS = (
    "strict_phase_origin_source_artifact_exported",
    "nonconventional_phase_or_sign_source_exported",
    "coupling_to_9_over_5_variational_density_exported",
    "nonimported_9_over_5_variational_chain_rule_exported",
    "localized_action_density_exported",
    "translation_neutral_strict_source_law_exported",
    "strict_translation_breaking_phase_source_exported",
    "strict_damping_compression_bridge_exported",
    "full_kernel_bridge_exported",
    "role_transfer_exported",
    "ltotal_exported",
    "eom_closure_exported",
    "hamiltonian_closure_exported",
    "toe_closure_exported",
)

LANES = [
    ("non_c12_support_origin_density", "P2888", "9/5 representable but nonunique; support/origin/density triple remains unsourced."),
    ("translation_orbit_representatives", "P2889", "600 embedded carriers quotient to 50 free Z12 orbits; no embedded representative is sourced."),
    ("fourier_power_and_phase", "P2890", "Power/autocorrelation is too invariant; phaseful characters need an external phase/origin pin."),
    ("phase_origin_artifact_inventory", "P2891", "Generated-artifact inventory has no coupled positive strict phase/origin + 9/5 variational-density export."),
]


def json_files() -> list[Path]:
    return sorted(path for path in GEN.glob("*.json") if path.is_file() and path != OUT)


def walk(value: Any, prefix: str = "") -> list[tuple[str, Any]]:
    rows: list[tuple[str, Any]] = []
    if isinstance(value, dict):
        for key, child in value.items():
            rows.extend(walk(child, f"{prefix}.{key}" if prefix else str(key)))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            rows.extend(walk(child, f"{prefix}[{index}]"))
    else:
        rows.append((prefix, value))
    return rows


def negative_flags(payload: dict[str, Any]) -> dict[str, Any]:
    decision = payload.get("decision") if isinstance(payload.get("decision"), dict) else {}
    flags = decision.get("negative_export_flags") if isinstance(decision.get("negative_export_flags"), dict) else {}
    return flags


def lane_matrix(loaded: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for lane, key, reason in LANES:
        payload = loaded[key]
        flags = negative_flags(payload)
        checked = {term: flags.get(term) for term in UNLOCK_TERMS if term in flags}
        unlock_true = any(value is True for value in checked.values())
        rows.append({
            "lane": lane,
            "artifact": key,
            "status": payload.get("status"),
            "checked_negative_export_flags": checked,
            "unlock_flag_true": unlock_true,
            "blocked_now": not unlock_true,
            "reason": reason,
        })
    return rows


def broad_unlock_scan() -> dict[str, Any]:
    hits = []
    true_positive_like_hits = []
    unquarantined_true_positive_like_hits = []
    for file_path in json_files():
        try:
            payload = read_json(file_path)
        except Exception:
            continue
        for path, value in walk(payload):
            lower = path.lower()
            if any(term in lower for term in UNLOCK_TERMS):
                hits.append({"file": str(file_path.relative_to(ROOT)), "path": path, "value": value})
                if value is True and "negative_export_flags" not in lower and "no_" not in lower and "not_" not in lower:
                    hit = {"file": str(file_path.relative_to(ROOT)), "path": path}
                    true_positive_like_hits.append(hit)
                    if not ("p2608" in hit["file"] or "p2611" in hit["file"]):
                        unquarantined_true_positive_like_hits.append(hit)
    return {
        "generated_json_file_count": len(json_files()),
        "unlock_term_occurrence_count": len(hits),
        "stale_or_raw_true_positive_like_unlock_count": len(true_positive_like_hits),
        "unquarantined_true_positive_like_unlock_count": len(unquarantined_true_positive_like_hits),
        "sample_unlock_term_occurrences": hits[:20],
        "stale_or_raw_true_positive_like_unlock_hits": true_positive_like_hits[:20],
        "unquarantined_true_positive_like_unlock_hits": unquarantined_true_positive_like_hits[:20],
    }


def build_payload() -> dict[str, Any]:
    loaded = {key: read_json(path) for key, path in INPUTS.items()}
    lanes = lane_matrix(loaded)
    scan = broad_unlock_scan()
    no_frontier = all(row["blocked_now"] for row in lanes) and scan["unquarantined_true_positive_like_unlock_count"] == 0
    return {
        "status": "P2892_POST_PHASE_ORIGIN_INVENTORY_STATE_MAP_NO_NEW_LIVE_FRONTIER_CERTIFICATE" if no_frontier else "P2892_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "post_p2891_state_map": {
            "candidate_class": "current post-P2891 9/5 phase/origin sourcehood frontier across support, translation, Fourier, and generated-artifact inventory lanes",
            "lane_count": len(lanes),
            "blocked_lane_count": sum(1 for row in lanes if row["blocked_now"]),
            "broad_unlock_scan": scan,
        },
        "lane_matrix": lanes,
        "acceptance_matrix": {
            "p2891_rechecked": loaded["P2891"].get("status") == "P2891_STRICT_PHASE_ORIGIN_SOURCE_ARTIFACT_INVENTORY_NO_GO_AUDIT_NO_CLOSURE",
            "all_recent_lanes_blocked": all(row["blocked_now"] for row in lanes),
            "no_positive_unlock_found_in_broad_scan": scan["unquarantined_true_positive_like_unlock_count"] == 0,
            "accepted_as_no_new_live_frontier_certificate": no_frontier,
        },
        "decision": {
            "no_new_live_frontier_certificate": no_frontier,
            "negative_export_flags": {term: False for term in UNLOCK_TERMS},
            "reason": "P2892 reconciles the post-P2888/P2891 9/5 phase/origin sourcehood lane.  Non-C12 carriers, translation orbits, Fourier power/phase, and generated-artifact inventory are all blocked on current artifacts; the broad generated-JSON scan finds no unquarantined positive unlock for strict phase/origin sourcehood, 9/5 variational density, localized action density, strict damping bridge, L_total, EOM, Hamiltonian, role transfer, or ToE closure.",
            "next_honest_step": "Do not replay non-C12 support choices, translation representatives, Fourier phase/power, bounded density coefficients, C12-neutral unit measures, ratio algebra, scalar Euler transmission, or bare inventory hits as sourcehood.  A next proof-grade move must supply one genuinely new explicit strict source formula/artifact with computed nonconventional sign/phase and a coupling theorem to the 9/5 variational density, or pivot to a genuinely different typed object outside the support/origin/density/orbit/Fourier/inventory family.  If neither is supplied, preserve this no-new-live-frontier certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2892/S1842 post phase/origin inventory state-map no-new-live-frontier certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Computed state-map facts",
        f"- recent lane count: `{payload['post_p2891_state_map']['lane_count']}`",
        f"- blocked recent lanes: `{payload['post_p2891_state_map']['blocked_lane_count']}`",
        f"- generated JSON files scanned: `{payload['post_p2891_state_map']['broad_unlock_scan']['generated_json_file_count']}`",
        f"- unquarantined positive-like unlock hits: `{payload['post_p2891_state_map']['broad_unlock_scan']['unquarantined_true_positive_like_unlock_count']}`",
        "",
        "## Lane matrix",
    ]
    for row in payload["lane_matrix"]:
        lines.append(f"- `{row['lane']}`: blocked_now=`{row['blocked_now']}`. {row['reason']}")
    lines.extend(["", "## Boundary", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2892/S1842 post phase/origin inventory state-map no-new-live-frontier certificate",
        "## P2892/S1842 post phase/origin inventory state-map no-new-live-frontier certificate\n\n"
        "`P2892/S1842` reconciles the post-`P2888`/`P2891` `9/5` phase/origin sourcehood frontier.  Non-`C12` carriers, translation orbits, Fourier power/phase, and generated-artifact inventory remain blocked on current artifacts, and the broad generated-JSON scan finds no unquarantined positive unlock for strict phase/origin sourcehood, nonimported `9/5` variational density, localized action density, strict damping bridge, `L_total`, EOM, Hamiltonian, role transfer, or ToE closure.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2892/S1842 post phase/origin inventory state-map `L_total` guard",
        "## P2892/S1842 post phase/origin inventory state-map `L_total` guard\n\n"
        "`P2892/S1842` is a no-new-live-frontier certificate, not a new action construction.  It does not add a strict action term, localized unit-bearing density, variational chain rule into nonproxy `L_total`, EOM, Hamiltonian, role transfer, bridge closure, or ToE closure.\n",
    )
    append_once(
        AGENTS,
        "Current post-P2891 phase/origin state-map no-new-live-frontier guardrail (P2892/S1842, 2026-06-19)",
        "## Current post-P2891 phase/origin state-map no-new-live-frontier guardrail (P2892/S1842, 2026-06-19)\n\n"
        "- P2892 reconciles the post-P2888/P2891 `9/5` phase/origin sourcehood lane across non-`C12` carriers, translation orbits, Fourier power/phase, and generated-artifact inventory.\n"
        "- All four recent lanes remain blocked on current artifacts, and no unquarantined positive generated-JSON unlock is found for strict phase/origin sourcehood, `9/5` variational density, localized action density, strict damping bridge, `L_total`, EOM, Hamiltonian, role transfer, or ToE closure.\n"
        "- Do not promote non-`C12` support choices, translation representatives, Fourier phase/power, bounded density coefficients, `C12`-neutral unit measures, ratio algebra, scalar Euler transmission, or bare inventory hits to sourcehood or closure.\n"
        "- A next admissible proof-grade move must supply one genuinely new explicit strict source formula/artifact with computed nonconventional sign/phase and a coupling theorem to the `9/5` variational density, pivot to a genuinely different typed object outside the support/origin/density/orbit/Fourier/inventory family, or preserve this no-new-live-frontier certificate.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
