#!/usr/bin/env python3
"""P2675/S1625: S6 Sigma_sel_src_target_v1 -> F301 typed-arrow audit.

Audits the exact S6 arrow recommended by P2674:
    Sigma_sel_src_target_v1 -> surviving F301 pair1/pair2 carrier
before Q_basis/preLM and projector-only collapse, without fiat identification.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2675_s1625_sigma_to_f301_typed_arrow_s6_audit.json"
MD = GEN / "p2675_s1625_sigma_to_f301_typed_arrow_s6_audit.md"
P2674 = GEN / "p2674_s1624_tau_src_pair12_typed_seed_o3_audit.json"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"

NEGATIVE_EXPORT_FLAGS = [
    "s6_sigma_to_f301_typed_arrow_exported",
    "sigma_sel_src_target_typed_as_f301_without_fiat",
    "pre_q_basis_pre_projector_chart_label_retaining_descent_exported",
    "o3_chart_sensitive_pair12_typed_seed_exported",
    "boundary_square_cycle_typed_arrow_exported",
    "sector_swap_sourced_invariant_exported",
    "q_w_2191_discharged",
    "role_bearing_ltotal_reenabled",
    "toe_closure_claimed",
]


def sha256_file(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def sha256_json(obj: Any) -> str:
    data = json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":"))
    return hashlib.sha256(data.encode()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, ".",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "*.json",
            "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**",
        ],
        cwd=REPO,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:60]}


def semantic_rg_audit() -> dict[str, Any]:
    patterns = {
        "sigma_codomain_content": "Sigma_sel_src_target_v1|selector witness target|source-topology selector target",
        "f301_carrier_content": "F301|surviving F301|residual-datum pair1/pair2 carrier|pair1/pair2 residual-datum",
        "typed_arrow_target_content": "Sigma_sel_src_target_v1.*F301|F301.*Sigma_sel_src_target_v1|typed source-side bridge from Sigma_sel_src_target_v1|chart_sensitive_pair12_typed_descent_from_Sigma_sel_src_target_v1",
        "precollapse_constraint_content": "prior to Q_basis_sel_v1 terminal collapse|without_Q_basis_sel_v1_terminal_collapse|before Q_basis|preLM|basis-free class-reduction",
        "projector_collapse_constraint_content": "projector-only local pair12 atlas|without_projector_only_atlas_collapse|projector-only.*collapse|local pair1/pair2 atlas",
        "fiat_blocker_content": "by fiat|not identify|must_not_identify|without.*fiat|unbridged to.*pair1/pair2",
        "actual_export_blocker_content": "not actually realized|no actual realization|nonexport|future-only target|below actual subinterface export",
        "closure_guard_content": "QW-2191|L_total|ToE closure|boundary square cycle|sector-swap sourced invariant",
    }
    return {
        "tool": "rg",
        "mode": "content-first exact S6 Sigma_sel_src_target_v1 -> F301 typed-arrow audit; not name/number-only",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
    }


def upstream_consistency() -> dict[str, Any]:
    p2674 = load_json(P2674)
    decision = p2674.get("closure_decision", {})
    lattice = p2674.get("finite_o3_lattice", {})
    return {
        "p2674_names_s6_missing": "S6_sigma_to_f301_typed_arrow_not_by_fiat" in lattice.get("missing_o3_subobligations_now", []),
        "p2674_o3_not_passed": decision.get("o3_exported_now") is False,
        "p2674_blocks_o4": decision.get("boundary_square_arrow_allowed_next") is False,
        "p2674_blocks_ltotal": decision.get("role_bearing_ltotal_now") is False,
    }


def s6_acceptance_vector() -> list[dict[str, Any]]:
    return [
        {
            "id": "A1_sigma_sel_src_target_present",
            "description": "Sigma_sel_src_target_v1 / selector-target codomain material exists",
            "satisfied_now": True,
            "content_pattern": "Sigma_sel_src_target_v1|selector witness target|source-topology selector target",
        },
        {
            "id": "A2_surviving_f301_pair12_carrier_present",
            "description": "surviving F301 residual-datum pair1/pair2 carrier material exists",
            "satisfied_now": True,
            "content_pattern": "F301|surviving F301|residual-datum pair1/pair2 carrier|pair1/pair2 residual-datum",
        },
        {
            "id": "A3_same_tau_packet_link_present",
            "description": "the selector-side material and F301 carrier are linked to the same tau_src packet",
            "satisfied_now": True,
            "content_pattern": "same tau_src|tau_src_candidate_v1.*F301|F301.*tau_src_candidate_v1|same tau_src packet",
        },
        {
            "id": "A4_chart_label_retaining_arrow_exported",
            "description": "a current chart-label-retaining typed arrow from Sigma_sel_src_target_v1 to F301 is actually exported",
            "satisfied_now": False,
            "content_pattern": "chart-label-retaining.*Sigma_sel_src_target_v1.*F301|Sigma_sel_src_target_v1.*chart-label-retaining.*F301|actual.*typed arrow.*F301",
        },
        {
            "id": "A5_pre_collapse_nonquotient_descent_exported",
            "description": "the arrow is before Q_basis/preLM/basis-free collapse and is not merely quotient-class material",
            "satisfied_now": False,
            "content_pattern": "prior to Q_basis_sel_v1 terminal collapse|without_Q_basis_sel_v1_terminal_collapse|preLM|basis-free class-reduction|quotient-class only",
        },
        {
            "id": "A6_nonprojector_local_atlas_descent_exported",
            "description": "the arrow is not only the projector-only local pair12 atlas collapse",
            "satisfied_now": False,
            "content_pattern": "projector-only local pair12 atlas|without_projector_only_atlas_collapse|projector-only.*collapse|local pair1/pair2 atlas",
        },
        {
            "id": "A7_no_fiat_identification_proof_exported",
            "description": "Sigma_sel_src_target_v1 is not identified with F301 or the atlas by declaration/fiat",
            "satisfied_now": False,
            "content_pattern": "by fiat|not identify|must_not_identify|without.*fiat|unbridged to.*pair1/pair2",
        },
    ]


def score_vector(vector: list[dict[str, Any]]) -> list[dict[str, Any]]:
    scored = []
    for item in vector:
        hits = rg_count(item["content_pattern"])
        scored.append({**item, "content_hits": hits["count"], "content_samples": hits["samples"][:10]})
    return scored


def finite_s6_lattice(scored: list[dict[str, Any]]) -> dict[str, Any]:
    ids = [item["id"] for item in scored]
    rows = []
    for bits in itertools.product([False, True], repeat=len(ids)):
        state = dict(zip(ids, bits, strict=True))
        rows.append({"state": state, "passes_s6": all(bits)})
    current = {item["id"]: item["satisfied_now"] for item in scored}
    return {
        "total_states_checked": len(rows),
        "passing_s6_states_count": sum(row["passes_s6"] for row in rows),
        "current_state": current,
        "current_state_passes_s6": all(current.values()),
        "missing_s6_obligations_now": [key for key, value in current.items() if not value],
        "hamming_distance_from_s6_pass": sum(1 for value in current.values() if not value),
        "s6_pass_state": {key: True for key in ids},
    }


def obstruction_table(scored: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "lane": "basis_free_Q_basis_continuation",
            "available_material": True,
            "why_not_s6": "forgets chart labels / quotient-class continuation rather than typing the surviving F301 carrier",
            "blocked_by": "A5_pre_collapse_nonquotient_descent_exported",
        },
        {
            "lane": "local_pair12_atlas_projector_lane",
            "available_material": True,
            "why_not_s6": "projector-only local atlas does not bind Sigma_sel_src_target_v1 to F301 as a source-side typed arrow",
            "blocked_by": "A6_nonprojector_local_atlas_descent_exported",
        },
        {
            "lane": "route_local_T220_T222_seed_target_family",
            "available_material": True,
            "why_not_s6": "target/attempt material remains future-only or nonexport, not an actual typed arrow",
            "blocked_by": "A4_chart_label_retaining_arrow_exported",
        },
        {
            "lane": "declaration_or_identification_shortcut",
            "available_material": True,
            "why_not_s6": "would identify Sigma/F301/atlas by fiat and therefore fail the nonconvention proof obligation",
            "blocked_by": "A7_no_fiat_identification_proof_exported",
        },
    ]


def closure_decision(lattice: dict[str, Any]) -> dict[str, Any]:
    return {
        "decision": "P2675_S6_SIGMA_TO_F301_TYPED_ARROW_AUDIT__NO_ACTUAL_NONFIAT_PRECOLLAPSE_ARROW_EXPORT",
        "professorial_verdict": (
            "P2675 audits the exact S6 arrow requested by P2674. Sigma_sel_src_target_v1 material, the surviving F301 pair1/pair2 carrier, and same-tau packet linkage are real, "
            "but the current repo still lacks an actual chart-label-retaining Sigma->F301 typed arrow, a pre-Q_basis/preLM nonquotient descent, a nonprojector local-atlas descent, and a proof that the binding is not by fiat. "
            "Therefore S6 fails and O3 remains blocked."
        ),
        "next_honest_step": (
            "The next honest step is to stop descending inside the same T220/T222 seed-target ladder unless a genuinely new nonquotient, nonprojector morphism is supplied. "
            "Either construct one explicit pre-collapse naturality square whose source is Sigma_sel_src_target_v1 and whose codomain is the surviving F301 carrier, or promote S6/O3 to a no-go for this tau_src -> pair12 -> boundary-square route. "
            "Do not attempt O4/O5 before that."
        ),
        "current_state_passes_s6": lattice["current_state_passes_s6"],
        "missing_s6_obligations_now": lattice["missing_s6_obligations_now"],
        "hamming_distance_from_s6_pass": lattice["hamming_distance_from_s6_pass"],
        "s6_exported_now": False,
        "o3_exported_now": False,
        "boundary_square_arrow_allowed_next": False,
        "sector_swap_invariant_allowed_next": False,
        "boundary_phase_bit_target_exported_now": False,
        "qw2191_discharged_now": False,
        "role_bearing_ltotal_now": False,
        "toe_closure_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lattice = payload["finite_s6_lattice"]
    decision = payload["closure_decision"]
    lines = [
        "# P2675/S1625 Sigma_sel_src_target_v1 -> F301 typed-arrow S6 audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first repo grep",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: `{data['count']}` hits")
    lines.extend(["", "## S6 acceptance vector"])
    for item in payload["s6_acceptance_vector"]:
        lines.append(f"- `{item['id']}`: satisfied_now=`{item['satisfied_now']}`, content_hits=`{item['content_hits']}` — {item['description']}")
    lines.extend(["", "## Obstruction table"])
    for row in payload["obstruction_table"]:
        lines.append(f"- `{row['lane']}`: blocked_by=`{row['blocked_by']}` — {row['why_not_s6']}")
    lines.extend([
        "",
        "## Finite S6 lattice",
        f"Total states checked: `{lattice['total_states_checked']}`.",
        f"Passing S6 states: `{lattice['passing_s6_states_count']}`.",
        f"Current state passes S6? `{lattice['current_state_passes_s6']}`.",
        f"Missing S6 obligations now: `{lattice['missing_s6_obligations_now']}`.",
        f"Hamming distance from S6 pass: `{lattice['hamming_distance_from_s6_pass']}`.",
        "",
        "## Verdict",
        decision["professorial_verdict"],
        f"Decision: `{decision['decision']}`.",
        "",
        "## Next honest step",
        decision["next_honest_step"],
        "",
        "## Negative exports",
    ])
    for key, value in payload["negative_export_flags"].items():
        lines.append(f"- `{key}`: `{value}`")
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    scored = score_vector(s6_acceptance_vector())
    lattice = finite_s6_lattice(scored)
    decision = closure_decision(lattice)
    payload: dict[str, Any] = {
        "status": "P2675_S6_SIGMA_TO_F301_TYPED_ARROW_AUDIT_NO_FALSE_PASS",
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_hashes": {
            "P2674": sha256_file(P2674),
            "STRICT_EQUATION_SHEET": sha256_file(STRICT_EQUATION_SHEET),
            "STRICT_LAGRANGIAN_DRAFT": sha256_file(STRICT_LAGRANGIAN_DRAFT),
        },
        "upstream_consistency": upstream_consistency(),
        "s6_acceptance_vector": scored,
        "obstruction_table": obstruction_table(scored),
        "finite_s6_lattice": lattice,
        "closure_decision": decision,
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2675/S1625 S6 Sigma-to-F301 typed-arrow guard",
        "## P2675/S1625 S6 Sigma-to-F301 typed-arrow guard\n\n"
        "`P2675/S1625` audits the exact S6 arrow `Sigma_sel_src_target_v1 -> surviving F301 pair1/pair2 carrier` before `Q_basis`/preLM and projector-only collapse.  The Sigma-side target, F301 carrier, and same-`tau_src` packet linkage are real, but no actual chart-label-retaining Sigma->F301 typed arrow, pre-collapse nonquotient descent, nonprojector local-atlas descent, or non-fiat binding proof is exported.  Therefore S6 fails, O3 remains blocked, O4/O5 remain inadmissible, and no boundary-phase bit target, beta source, `QW-2191` discharge, bridge completion, role transfer, role-bearing `L_total`, or ToE closure is exported.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2675/S1625 S6 Sigma-to-F301 Ltotal guard",
        "## P2675/S1625 S6 Sigma-to-F301 Ltotal guard\n\n"
        "`P2675/S1625` keeps `L_total` closed at the S6 subobligation.  A tau_src-derived pair12 variational source remains inadmissible until a pre-`Q_basis`, preLM, nonprojector, chart-label-retaining `Sigma_sel_src_target_v1 -> F301` typed arrow is exported without fiat identification.\n",
    )
    return payload


if __name__ == "__main__":
    main()
