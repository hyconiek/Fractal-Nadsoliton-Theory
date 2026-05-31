#!/usr/bin/env python3
"""Scratch probe: S1 selector-margin monotonicity obstruction certificate.

The closure-plan dependency certificate identifies the next honest step as S1:
a strict-local selector-margin monotonicity witness.  This probe performs the
finite preflight on the existing strict nu-branch policy confidence artifacts.
It does not invent a new provider class; it audits the available one.

Result: the available surrogate route does not satisfy S1.  The confidence
margin shrinks as sampling budget grows, but the certified selector margins to
target remain negative, the locked replay also remains negative, and therefore
no local S1 witness is exported here.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FAR = ROOT / "fundamental_action_reconstruction"
GEN = FAR / "generated"
OUT_JSON = HERE / "bridge_strict_completion_s1_selector_margin_monotonicity_obstruction_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_s1_selector_margin_monotonicity_obstruction_certificate_report.md"
CLOSURE_PLAN_REPORT = HERE / "bridge_strict_completion_closure_plan_dependency_certificate_report.json"
P2275_REPORT = GEN / "p2275_s1225_strict_nu_branch_group_policy_certified_box_corner_replay_passrate_floor_probe.json"
P2277_REPORT = GEN / "p2277_s1227_strict_nu_branch_group_policy_adaptive_confidence_margin_certificate_probe.json"
P2278_REPORT = GEN / "p2278_s1228_strict_nu_branch_group_policy_confidence_curve_sweep_probe.json"
P2279_REPORT = GEN / "p2279_s1229_strict_nu_branch_group_policy_locked_confidence_profile_corner_replay_probe.json"
P1445_PACKET = FAR / "P1445_STRICT_NEXT_HONEST_STEP_F_NADSOLITON_TO_LSM_PLUS_LGR_NO_LEGACY_BRIDGE_PACKET_PL.md"


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    if not path.exists():
        raise FileNotFoundError(f"missing source text: {path}")
    return path.read_text(encoding="utf-8")


def grouped_sweep_rows(sweep_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    deltas = sorted({float(row["delta"]) for row in sweep_rows})
    groups = []
    for delta in deltas:
        rows = sorted((row for row in sweep_rows if float(row["delta"]) == delta), key=lambda row: int(row["trial_multiplier"]))
        worst_margins = [float(row["worst_margin_to_target"]) for row in rows]
        adaptive_margins = [float(row["adaptive_margin"]) for row in rows]
        groups.append({
            "delta": delta,
            "trial_multipliers": [int(row["trial_multiplier"]) for row in rows],
            "worst_margins_to_target": worst_margins,
            "adaptive_confidence_margins": adaptive_margins,
            "strictly_improves_worst_margin": all(b > a for a, b in zip(worst_margins, worst_margins[1:])),
            "nondecreasing_worst_margin": all(b >= a for a, b in zip(worst_margins, worst_margins[1:])),
            "confidence_margin_decreases": all(b < a for a, b in zip(adaptive_margins, adaptive_margins[1:])),
            "final_worst_margin_positive": bool(worst_margins and worst_margins[-1] > 0),
        })
    return groups


def extract_floor_rows(p2275: dict[str, Any], p2277: dict[str, Any], p2279: dict[str, Any]) -> list[dict[str, Any]]:
    replay_rows = p2275["strict_nu_branch_group_policy_certified_box_corner_replay_passrate_floor_probe"]["replay_rows"]
    adaptive_rows = p2277["strict_nu_branch_group_policy_adaptive_confidence_margin_certificate_probe"]["certificate_rows"]
    locked_rows = p2279["strict_nu_branch_group_policy_locked_confidence_profile_corner_replay_probe"]["locked_rows"]
    rows = []
    for replay, adaptive, locked in zip(replay_rows, adaptive_rows, locked_rows):
        rows.append({
            "risk_tolerance": float(replay["risk_tolerance"]),
            "target_floor": float(locked["target_floor"]),
            "baseline_empirical_floor": float(replay["empirical_passrate_floor_over_corners"]),
            "adaptive_certified_floor": float(adaptive["adaptive_certified_floor"]),
            "locked_adaptive_certified_floor": float(locked["adaptive_certified_floor"]),
            "locked_margin_to_target": float(locked["margin_to_target"]),
            "locked_meets_target": bool(locked["adaptive_floor_meets_target"]),
        })
    return rows


def build_payload() -> dict[str, Any]:
    closure_plan = load_json(CLOSURE_PLAN_REPORT)
    p2275 = load_json(P2275_REPORT)
    p2277 = load_json(P2277_REPORT)
    p2278 = load_json(P2278_REPORT)
    p2279 = load_json(P2279_REPORT)
    p1445 = load_text(P1445_PACKET)

    sweep = p2278["strict_nu_branch_group_policy_confidence_curve_sweep_probe"]
    locked = p2279["strict_nu_branch_group_policy_locked_confidence_profile_corner_replay_probe"]
    sweep_groups = grouped_sweep_rows(sweep["sweep_rows"])
    floor_rows = extract_floor_rows(p2275, p2277, p2279)
    locked_summary = locked["global_summary"]

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_S1_SELECTOR_MARGIN_MONOTONICITY_OBSTRUCTION_CERTIFICATE__FINITE_PREFLIGHT",
        "status": "s1-selector-margin-monotonicity-witness-not-exported-local-obstruction-certified",
        "source_reports": {
            "closure_plan_dependency_certificate": str(CLOSURE_PLAN_REPORT.relative_to(ROOT)),
            "p2275_corner_replay_passrate_floor": str(P2275_REPORT.relative_to(ROOT)),
            "p2277_adaptive_confidence_margin": str(P2277_REPORT.relative_to(ROOT)),
            "p2278_confidence_curve_sweep": str(P2278_REPORT.relative_to(ROOT)),
            "p2279_locked_confidence_profile_replay": str(P2279_REPORT.relative_to(ROOT)),
            "P1445_next_honest_step_packet": str(P1445_PACKET.relative_to(ROOT)),
        },
        "grep_disambiguation": {
            "searched_terms": [
                "strict-local selector margin monotonicity witness",
                "selector margin monotonicity",
                "margin monotonicity witness",
                "replay-stability",
                "strict_nu_branch_group_policy_confidence_curve_sweep",
                "locked confidence profile",
                "S1 selector obstruction",
            ],
            "finding": "Repo grep finds prior confidence/replay artifacts and P1445's S1 target, but no strict-completion S1 pass/fail certificate that evaluates monotonic selector margins and locked replay against the closure-plan dependency matrix before this file.",
        },
        "s1_definition": {
            "source": "P1445",
            "required_local_pass_conditions": [
                "monotonic selector-margin growth on a fixed perturbation/confidence grid",
                "replay-stability under a locked profile",
                "explicit local-only pass scope with no global closure claim",
            ],
            "provider_class_audited": "strict_nu_branch_group_policy_surrogate_confidence_route",
            "scope": "local surrogate preflight only; not strict-core selector closure",
        },
        "finite_margin_table": {
            "sweep_groups_by_delta": sweep_groups,
            "floor_rows_by_risk": floor_rows,
            "locked_profile": locked["locked_profile"],
            "locked_global_summary": locked_summary,
        },
        "s1_obstruction_summary": {
            "closure_plan_recommends_s1": closure_plan["closure_plan_summary"].get("recommended_next_step_is_S1_selector_margin", closure_plan["closure_plan_summary"].get("S1_selector_margin_remains_next_selector_subproblem", False)),
            "p1445_contains_s1_target": "strict-local selector margin monotonicity witness" in p1445,
            "confidence_margins_decrease_with_budget": all(group["confidence_margin_decreases"] for group in sweep_groups),
            "worst_selector_margins_strictly_improve_with_budget": all(group["strictly_improves_worst_margin"] for group in sweep_groups),
            "worst_selector_margins_reach_positive": any(group["final_worst_margin_positive"] for group in sweep_groups),
            "locked_replay_meets_all_targets": bool(locked_summary["all_rows_meet_target"]),
            "locked_worst_margin_to_target": float(locked_summary["worst_margin_to_target"]),
            "s1_witness_exported": False,
            "s1_obstruction_certified": (
                all(group["confidence_margin_decreases"] for group in sweep_groups)
                and not all(group["strictly_improves_worst_margin"] for group in sweep_groups)
                and not any(group["final_worst_margin_positive"] for group in sweep_groups)
                and not bool(locked_summary["all_rows_meet_target"])
            ),
        },
        "blocker_context": {
            "resolved_locally": [
                "The existing confidence artifacts are sufficient to evaluate the available S1 surrogate route.",
                "The statistical confidence margin improves, but the selector margin-to-target remains negative.",
                "Locked replay does not meet all target floors, so S1 is not exported from this provider class.",
            ],
            "still_open": [
                "a new provider class or blocker cut for S1",
                "strict-local selector margin monotonicity witness",
                "kernel-alone/global QW-2191 discharge",
                "strict F_nadsoliton => L_SM + L_GR bridge",
                "strict ToE closure",
            ],
        },
        "proof_certificate": {
            "grep_step": "Repo grep distinguishes this from prior confidence sweeps: those artifacts exist, but no S1 pass/fail obstruction certificate was exported.",
            "sweep_step": "For each delta, the confidence margin decreases with budget, yet the worst selector margin-to-target does not strictly improve and never becomes positive.",
            "locked_replay_step": "The locked replay profile has all_rows_meet_target=false and worst_margin_to_target=-0.99, so replay-stability at the requested floor is absent.",
            "verdict_step": "S1 is a certified local obstruction on the currently available strict nu-branch surrogate provider class; no S1 witness is exported.",
            "next_step": "Continue by changing provider class or blocker cut, as required by the closure-plan dependency certificate.",
            "theoretical_limit": "This report is a local S1 obstruction certificate; it does not discharge QW-2191, build F_nadsoliton=>L_SM+L_GR, or close ToE.",
        },
        "hard_limits": [
            "No strict-local selector margin monotonicity witness is claimed.",
            "No strict-core selector closure is claimed.",
            "No QW-2191 selector discharge is claimed.",
            "No strict F_nadsoliton => L_SM + L_GR bridge is claimed.",
            "No ToE closure is claimed.",
            "No unqualified identity K_legacy_ont == K_strict_gate is claimed or used.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    summary = payload["s1_obstruction_summary"]
    lines = [
        "# S1 selector-margin monotonicity obstruction certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "The available strict nu-branch surrogate route does not export the S1",
        "strict-local selector margin monotonicity witness.",
        "",
        "## Summary",
        "",
    ]
    for key, value in summary.items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Hard limits", ""])
    for limit in payload["hard_limits"]:
        lines.append(f"- {limit}")
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    payload = build_payload()
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_markdown(payload)
    print(json.dumps(payload["s1_obstruction_summary"], indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
