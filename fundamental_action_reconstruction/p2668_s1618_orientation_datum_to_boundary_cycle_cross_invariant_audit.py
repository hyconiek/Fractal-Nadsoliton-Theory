#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2668_s1618_orientation_datum_to_boundary_cycle_cross_invariant_audit.json"
MD = GEN / "p2668_s1618_orientation_datum_to_boundary_cycle_cross_invariant_audit.md"
P2667 = GEN / "p2667_s1617_pair12_boundary_orientation_reversal_no_go_audit.json"
N546 = ROOT / "N546_CURRENT_EXPORTED_S_SEL_INT_STRICT_CORE_SOURCE_OBJECT_ADMISSIBLE_ORIENTATION_EXPORT_THEOREM.md"
N500 = ROOT / "N500_CURRENT_FIRST_ACTUAL_STRICT_SHANNON_ELEMENT_ORDER_REFERENCE_INTERNAL_ORIENTATION_DATUM_EXPORT_THEOREM.md"
N493 = ROOT / "N493_CURRENT_FIRST_STRICT_QW2191_RESIDUAL_Z2_SIGN_FLIP_GAUGE_EQUIVALENCE_THEOREM.md"
F467 = ROOT / "F467_CURRENT_STRICT_PAIR12345_SELECTOR_ATLAS_ORIENTED_TRANSPORT_LIFT_PACKET.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"

NEGATIVE_EXPORT_FLAGS = [
    "orientation_datum_to_boundary_cycle_cross_invariant_exported",
    "canonical_pair12_boundary_orientation_map_exported",
    "orientation_reversal_forbidden_internally",
    "pair12_to_boundary_phase_sector_descent_exported",
    "boundary_phase_bit_target_exported_unconditionally",
    "intrinsic_entropy_level_exported",
    "bit_to_action_map_sourced_unconditionally",
    "bit_to_length_map_sourced_unconditionally",
    "target_independent_beta_source_exported",
    "q_w_2191_discharged",
    "role_bearing_ltotal_reenabled",
    "toe_closure_claimed",
]


def sha256_file(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def file_contains(path: Path, needle: str) -> bool:
    return path.exists() and needle in path.read_text(encoding="utf-8")


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run([
        "rg", "-n", pattern, ".",
        "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "*.json",
        "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**",
    ], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:120]}


def semantic_rg_audit() -> dict[str, Any]:
    patterns = {
        "orientation_datum_content": "orientation datum|internal orientation|admissible orientation|physical orientation datum",
        "boundary_cycle_cross_content": "boundary-cycle|boundary cycle|cross-invariant|sector swap|boundary sector|square holonomy",
        "lane_scope_limits_content": "lane-scoped|axis-only|residual Z2|sign-gauge|projector level|convention layer",
        "selector_pair12_content": "pair1|pair2|pair12|w_break|pair2.*positive|selector witness",
        "nonclosure_guard_content": "QW-2191|role-bearing L_total|ToE closure|beta source|nonexport|future-route",
    }
    return {"tool": "rg", "mode": "content-first orientation-datum to boundary-cycle cross-invariant audit", "patterns": {key: rg_count(pattern) for key, pattern in patterns.items()}}


def upstream_consistency() -> dict[str, Any]:
    p2667 = load_json(P2667)
    decision = p2667.get("closure_decision", {})
    return {
        "p2667_sector_one_mapping_exists": decision.get("sector_one_mapping_exists") is True,
        "p2667_reversal_unforbidden": decision.get("all_reversal_pairs_unforbidden") is True,
        "p2667_no_canonical_orientation_map": decision.get("canonical_orientation_map_exported_now") is False,
        "p2667_no_boundary_phase_bit_target": decision.get("boundary_phase_bit_target_exported_now") is False,
    }


def existing_orientation_sources() -> list[dict[str, Any]]:
    return [
        {
            "source": "N546_S_sel_int_pair1_orientation",
            "path": rel(N546),
            "present": N546.exists(),
            "exports_orientation": file_contains(N546, "exports one admissible orientation datum"),
            "scope_limit": "pair1 frame / not admissible S_sel_int or downstream completion",
            "ties_pair2_positive_to_boundary_sector1": False,
            "forbids_sector_swap": False,
        },
        {
            "source": "N500_Shannon_axis_only_orientation",
            "path": rel(N500),
            "present": N500.exists(),
            "exports_orientation": file_contains(N500, "internal orientation datum exists"),
            "scope_limit": "axis-only; residual Z2 sign remains",
            "ties_pair2_positive_to_boundary_sector1": False,
            "forbids_sector_swap": False,
        },
        {
            "source": "N493_residual_Z2_sign_flip_gauge_equivalence",
            "path": rel(N493),
            "present": N493.exists(),
            "exports_orientation": False,
            "scope_limit": "sign flips are gauge equivalences; no global physical uniqueness",
            "ties_pair2_positive_to_boundary_sector1": False,
            "forbids_sector_swap": False,
        },
        {
            "source": "F467_oriented_transport_convention_lift",
            "path": rel(F467),
            "present": F467.exists(),
            "exports_orientation": file_contains(F467, "oriented transport"),
            "scope_limit": "sign-tracked convention layer, not physical orientation datum",
            "ties_pair2_positive_to_boundary_sector1": False,
            "forbids_sector_swap": False,
        },
    ]


def cross_invariant_witness(sources: list[dict[str, Any]]) -> dict[str, Any]:
    rows = []
    for source in sources:
        changes_under_sector_swap = source["forbids_sector_swap"]
        tied_to_pair2 = source["ties_pair2_positive_to_boundary_sector1"]
        strict_current = source["present"] and source["exports_orientation"] and tied_to_pair2 and changes_under_sector_swap
        rows.append({
            **source,
            "changes_under_sector_swap": changes_under_sector_swap,
            "is_boundary_cycle_functor": False,
            "passes_cross_invariant_acceptance_now": strict_current,
        })
    return {
        "statement": "P2668 audits whether existing orientation-datum exports already provide the P2667 cross-invariant.  They do not: current orientation exports are pair1-frame, axis-only, gauge-equivalence, or convention-layer objects.  None is a boundary-cycle functor, none changes under the P2667 sector swap, and none ties pair2-positive to boundary sector 1.",
        "rows": rows,
        "orientation_sources_present": [row["source"] for row in rows if row["present"]],
        "passing_cross_invariant_sources": [row["source"] for row in rows if row["passes_cross_invariant_acceptance_now"]],
        "any_existing_orientation_export_present": any(row["present"] and row["exports_orientation"] for row in rows),
        "any_source_forbids_sector_swap": any(row["forbids_sector_swap"] for row in rows),
        "any_source_ties_pair2_to_sector1": any(row["ties_pair2_positive_to_boundary_sector1"] for row in rows),
    }


def closure_decision(witness: dict[str, Any]) -> dict[str, Any]:
    return {
        "decision": "P2668_ORIENTATION_DATUM_TO_BOUNDARY_CYCLE_CROSS_INVARIANT_AUDIT__NO_TRANSFER",
        "professorial_verdict": "P2668 checks the obvious possible rescue after P2667: maybe an existing orientation datum already supplies the missing cross-invariant.  Content grep and file-level audit find real orientation-datum material, but it is scoped as pair1-frame, axis-only with residual Z2, gauge-equivalence, or convention-layer transport.  None is a boundary-cycle functor; none forbids the P2667 sector swap; none ties pair2-positive to boundary sector 1.  Therefore no boundary-phase bit target, UV unit, beta source, QW-2191 discharge, L_total reopening, or ToE closure follows.",
        "next_honest_step": "Construct a genuinely new boundary-cycle cross-invariant or stop this route.  The invariant must be a current strict theorem, not a convention layer; must evaluate on the boundary cycle; must change under sector swap; and must tie pair2-positive to sector 1.  If that cannot be done, promote P2667/P2668 into a no-go for using existing orientation-datum exports as the entropy-bit source.",
        "passing_cross_invariant_sources": witness["passing_cross_invariant_sources"],
        "existing_orientation_export_present": witness["any_existing_orientation_export_present"],
        "source_forbids_sector_swap": witness["any_source_forbids_sector_swap"],
        "source_ties_pair2_to_sector1": witness["any_source_ties_pair2_to_sector1"],
        "boundary_phase_bit_target_exported_now": False,
        "unconditional_uv_unit_selected_now": False,
        "beta_source_exported_now": False,
        "qw2191_discharged_now": False,
        "role_bearing_ltotal_now": False,
        "toe_closure_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    witness = payload["cross_invariant_witness"]
    decision = payload["closure_decision"]
    lines = [
        "# P2668/S1618 orientation-datum to boundary-cycle cross-invariant audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first repo grep",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines += ["", "## Upstream consistency"]
    for key, value in payload["upstream_consistency"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines += [
        "",
        "## Cross-invariant witness",
        witness["statement"],
        f"Existing orientation export present? `{witness['any_existing_orientation_export_present']}`.",
        f"Any source forbids sector swap? `{witness['any_source_forbids_sector_swap']}`.",
        f"Any source ties pair2-positive to sector 1? `{witness['any_source_ties_pair2_to_sector1']}`.",
        f"Passing cross-invariant sources: `{witness['passing_cross_invariant_sources']}`.",
        "",
        "| source | present | exports orientation? | scope limit | boundary-cycle functor? | forbids swap? | ties pair2->sector1? | passes? |",
        "| --- | ---: | ---: | --- | ---: | ---: | ---: | ---: |",
    ]
    for row in witness["rows"]:
        lines.append(f"| `{row['source']}` | `{row['present']}` | `{row['exports_orientation']}` | {row['scope_limit']} | `{row['is_boundary_cycle_functor']}` | `{row['forbids_sector_swap']}` | `{row['ties_pair2_positive_to_boundary_sector1']}` | `{row['passes_cross_invariant_acceptance_now']}` |")
    lines += [
        "",
        "## Verdict",
        decision["professorial_verdict"],
        f"Decision: `{decision['decision']}`.",
        f"Boundary-phase bit target exported now? `{decision['boundary_phase_bit_target_exported_now']}`.",
        f"Beta source exported now? `{decision['beta_source_exported_now']}`.",
        f"QW-2191 discharged now? `{decision['qw2191_discharged_now']}`.",
        f"Role-bearing L_total now? `{decision['role_bearing_ltotal_now']}`.",
        f"ToE closure now? `{decision['toe_closure_now']}`.",
        "",
        "## Next honest step",
        decision["next_honest_step"],
        "",
        "## Negative exports",
    ]
    for key, value in payload["negative_export_flags"].items():
        lines.append(f"- `{key}`: `{value}`")
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    audit = semantic_rg_audit()
    upstream = upstream_consistency()
    sources = existing_orientation_sources()
    witness = cross_invariant_witness(sources)
    decision = closure_decision(witness)
    payload: dict[str, Any] = {
        "status": "P2668_ORIENTATION_DATUM_TO_BOUNDARY_CYCLE_CROSS_INVARIANT_AUDIT_NO_FALSE_PASS",
        "semantic_rg_antiduplication_audit": audit,
        "source_hashes": {"P2667": sha256_file(P2667), "N546": sha256_file(N546), "N500": sha256_file(N500), "N493": sha256_file(N493), "F467": sha256_file(F467), "STRICT_EQUATION_SHEET": sha256_file(STRICT_EQUATION_SHEET), "STRICT_LAGRANGIAN_DRAFT": sha256_file(STRICT_LAGRANGIAN_DRAFT)},
        "upstream_consistency": upstream,
        "existing_orientation_sources": sources,
        "cross_invariant_witness": witness,
        "closure_decision": decision,
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2668/S1618 orientation-datum to boundary-cycle cross-invariant guard", "## P2668/S1618 orientation-datum to boundary-cycle cross-invariant guard\n\n`P2668/S1618` audits whether existing orientation-datum exports can supply the P2667 cross-invariant.  They cannot: the current orientation material is pair1-frame, axis-only with residual `Z2`, sign-flip gauge-equivalence, or convention-layer oriented transport.  None is a boundary-cycle functor, none forbids sector swap, and none ties `pair2_positive` to boundary sector `1`.  Therefore no boundary-phase bit target, intrinsic UV unit, beta source, `QW-2191` discharge, bridge completion, role transfer, role-bearing `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2668/S1618 orientation-datum cross-invariant Ltotal guard", "## P2668/S1618 orientation-datum cross-invariant Ltotal guard\n\n`P2668/S1618` keeps `L_total` closed to orientation-datum-derived boundary entropy terms.  Existing orientation exports may guide proof search, but a variational coefficient still requires a current strict boundary-cycle cross-invariant that changes under sector swap and ties `pair2_positive` to sector `1`.\n")
    return payload


if __name__ == "__main__":
    result = main()
    print(rel(OUT))
    print(rel(MD))
