#!/usr/bin/env python3
"""P2681/S1631: broad research-state map and AGENTS reorientation audit.

This packet answers the methodological correction that AGENTS.md had become a
stale steering source.  It performs broad content-first repo scans, records
which lanes have already been exhausted/repeated, updates AGENTS.md to require
state-map-first prioritisation, and then gives the next honest proof-grade move.
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
OUT = GEN / "p2681_s1631_professorial_research_state_map_and_agents_reorientation_audit.json"
MD = GEN / "p2681_s1631_professorial_research_state_map_and_agents_reorientation_audit.md"
AGENTS = REPO / "AGENTS.md"
P2680 = GEN / "p2680_s1630_legacy_strict_bridge_source_inventory_no_selector_replay_audit.json"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"

NEGATIVE_EXPORT_FLAGS = [
    "toe_closure_claimed",
    "q_w_2191_discharged",
    "legacy_strict_bridge_completed",
    "role_transfer_started",
    "selector_reopened_without_new_object",
    "tau_pair12_reopened_without_new_object",
    "beta_tors_chi11_reopened_without_new_object",
    "local_direct_mass_declared_main_bottleneck",
]


def sha256_file(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


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
    return {"count": len(lines), "samples": lines[:80]}


def broad_state_grep() -> dict[str, Any]:
    patterns = {
        "legacy_strict_bridge_and_role_transfer": r"legacy -> strict|legacy-to-strict|completion bridge|K_legacy_ont|K_strict_gate|role transfer|role-transfer",
        "selector_orientation_qw2191": r"QW-2191|selector source|orientation source|orientation-odd|XOR|XNOR|Sigma_sel_src_target_v1|F301|symmetry-breaking|spin/Pin|torsor",
        "tau_pair12_boundary_square": r"tau_src|pair12|pair1/pair2|typed seed|typed-seed|boundary-square|sector-swap",
        "beta_tors_chi11_and_damping": r"beta_tors.*chi|chi_11|chi11|torsion damping|d\^eta|eta = 9/5|eta=9/5|Z_beta|positive beta|canonical length|UV unit",
        "strict_lagrangian_eom_qg": r"L_total|Lagrangian|Euler|EOM|stress-energy|Bianchi|ADM|Riemann|QG closure|variational",
        "direct_route_mass_defects_p46": r"P46|N49|c1s1|m2 psi4|m2_psi4|psi4\*\*2/2|defect polynomial|zero witness|pairwise matching",
        "alpha_s_and_couplings": r"alpha_s|alpha_EM|sin\^2\(theta_W\)|Weinberg|coupling|lawful refined|SIR_DA|ERFS",
        "noncyclic_new_provider_constraints": r"noncyclic|new provider|provider class|new object|blocker-cut|reopen|repetition gate|no-go",
    }
    return {
        "tool": "rg",
        "mode": "broad content-first state-map grep before choosing a new research priority",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
    }


def prior_lane_ledger() -> list[dict[str, Any]]:
    return [
        {
            "lane": "selector/orientation/QW-2191",
            "already_done": ["P2619-P2623", "P2676-P2678", "P2679 repetition gate"],
            "current_status": "blocked unless a genuinely new strict internal source object is exported",
            "repeat_now_admissible": False,
            "scientific_leverage_now": "low without new object; high only if a new source theorem appears",
        },
        {
            "lane": "tau_src -> pair12 -> boundary-square / S6-O3 typed seed",
            "already_done": ["P2673-P2677", "P2679 repetition gate"],
            "current_status": "bounded no-go for current typed-seed route",
            "repeat_now_admissible": False,
            "scientific_leverage_now": "low; repetition would not advance ToE closure",
        },
        {
            "lane": "beta_tors -> chi11",
            "already_done": ["P2390", "P2435", "P2619", "P2677-P2679"],
            "current_status": "beta_tors remains scalar/damping input, not an orientation-odd sign source",
            "repeat_now_admissible": False,
            "scientific_leverage_now": "low unless a new bridge source atom changes the type",
        },
        {
            "lane": "legacy -> strict bridge-source atoms",
            "already_done": ["K1/K2/F2/F3/S2", "P2679", "P2680"],
            "current_status": "formal inventory exists; three missing source atoms remain after P2680",
            "repeat_now_admissible": False,
            "scientific_leverage_now": "medium only for one specific missing atom; not as a generic loop",
        },
        {
            "lane": "strict Lagrangian/EOM/QG nonproxy closure",
            "already_done": ["P1563-P1747 family", "P1833-P1907 family", "P1964-P1988 family", "P2314-P2329 family"],
            "current_status": "large executable scaffold exists; remaining blockers are concrete tensor/variation and reverse-closure obligations",
            "repeat_now_admissible": True,
            "scientific_leverage_now": "high if targeted at a finite obstruction matrix rather than broad restatement",
        },
        {
            "lane": "kernel-split-robust direct route P46/N49 mass-defect witnesses",
            "already_done": ["F3 classification", "R33-R83 direct defect packets", "multiple boundary theorem packets"],
            "current_status": "live frontier under F3, but not automatically the global bottleneck",
            "repeat_now_admissible": True,
            "scientific_leverage_now": "high for finite zero-witness/no-go computation on one named defect",
        },
        {
            "lane": "coupling-law/alpha_s/alpha_EM role-bearing transfer",
            "already_done": ["many alpha_s lawful-refined packets", "K1/K2/S2 role-transfer guard"],
            "current_status": "role-bearing physical claims remain downstream of bridge or independent strict source theorem",
            "repeat_now_admissible": False,
            "scientific_leverage_now": "premature unless converted into non-role-bearing finite obstruction audit",
        },
    ]


def p2680_consistency() -> dict[str, Any]:
    p2680 = load_json(P2680)
    decision = p2680.get("closure_decision", {})
    lattice = p2680.get("bridge_completion_lattice", {})
    return {
        "p2680_exists": not p2680.get("missing", False),
        "p2680_no_selector_replay": decision.get("selector_replay_used") is False,
        "p2680_bridge_not_completed": decision.get("full_bridge_completed_now") is False,
        "p2680_missing_atoms": lattice.get("missing_current_obligations", []),
        "p2680_hash": sha256_file(P2680),
    }


def opportunity_matrix() -> list[dict[str, Any]]:
    rows = [
        ("finite P46/N49 zero-witness obstruction matrix", 5, 5, 4, 2, "Compute one live direct-route defect as pass/no-go without touching stale selector/bridge loops."),
        ("strict Lagrangian/EOM reverse-closure obstruction matrix", 5, 4, 5, 3, "Use existing executable Lagrangian/EOM scaffolds to isolate a concrete nonproxy reverse-closure gap."),
        ("one P2680 missing source atom", 4, 3, 3, 4, "Try target-independent positive beta/Z_beta or canonical UV unit, but only as one atom, not generic bridge repetition."),
        ("generic legacy->strict bridge loop", 2, 2, 2, 5, "Stale AGENTS priority made this too attractive; P2680 already inventories the current gap."),
        ("selector/tau_pair12/beta_tors_chi11 replay", 1, 1, 1, 5, "Already bounded or repetition-gated; no new proof object exists."),
        ("role-transfer audit now", 1, 1, 2, 5, "Premature before bridge/source closure."),
    ]
    return [
        {
            "candidate": name,
            "toe_closure_leverage_1_to_5": leverage,
            "proof_computability_1_to_5": computability,
            "evidence_base_1_to_5": evidence,
            "loop_risk_1_to_5": loop_risk,
            "net_score": leverage + computability + evidence - loop_risk,
            "rationale": rationale,
        }
        for name, leverage, computability, evidence, loop_risk, rationale in rows
    ]


def finite_priority_lattice() -> dict[str, Any]:
    obligations = [
        "not_a_repetition_lane",
        "finite_executable_or_grep_auditable",
        "touches_live_toe_blocker",
        "does_not_need_role_transfer_first",
        "does_not_discharge_qw2191_by_fiat",
        "can_end_in_pass_or_no_go",
        "updates_stale_steering_before_continuing",
    ]
    candidates = {
        "finite_p46_n49_zero_witness_matrix": [True, True, True, True, True, True, True],
        "strict_lagrangian_eom_reverse_closure_matrix": [True, True, True, True, True, True, True],
        "single_p2680_missing_source_atom": [True, True, True, True, True, True, True],
        "generic_bridge_loop": [False, True, True, True, True, False, True],
        "selector_tau_beta_replay": [False, True, True, True, False, False, True],
        "role_transfer_now": [False, True, True, False, True, False, True],
    }
    pass_count = 0
    for bits in itertools.product([False, True], repeat=len(obligations)):
        pass_count += int(all(bits))
    return {
        "obligations": obligations,
        "total_states": 2 ** len(obligations),
        "passing_states": pass_count,
        "candidate_states": {name: dict(zip(obligations, values)) for name, values in candidates.items()},
        "passing_candidates": [name for name, values in candidates.items() if all(values)],
        "selected_next_lane": "finite_p46_n49_zero_witness_matrix",
        "selected_reason": "It is F3-live, finite, computational, avoids stale bridge/selector loops, and can produce either a real witness or a no-go certificate.",
    }


def closure_decision(matrix: list[dict[str, Any]], lattice: dict[str, Any]) -> dict[str, Any]:
    ranked = sorted(matrix, key=lambda row: row["net_score"], reverse=True)
    return {
        "decision": "P2681_STATE_MAP_REORIENTS_AGENTS_AND_SELECTS_FINITE_P46_N49_ZERO_WITNESS_FRONTIER_NO_FALSE_PASS",
        "professorial_verdict": (
            "Broad state reading shows that the stale AGENTS priority caused an over-focus on generic legacy->strict bridge repetition.  P2679-P2680 already made that lane precise: only specific missing source atoms remain, while selector, tau_src->pair12, and beta_tors->chi11 are not reopened.  The highest near-term scientific chance is therefore not another generic bridge audit, but a finite, executable live-frontier obstruction/witness computation."
        ),
        "top_ranked_candidates": ranked[:3],
        "next_honest_step": (
            "Run P2682 as a finite P46/N49 zero-witness/no-go matrix, preferably the direct m2 psi4 target action coefficient defect on common psi4**2/2 support named by F3, with explicit guards excluding selector replay, tau_src->pair12, beta_tors->chi11, role transfer, and ToE closure.  If that route is blocked by unavailable symbolic carriers, fall back to the strict Lagrangian/EOM reverse-closure obstruction matrix, not to generic bridge looping."
        ),
        "selected_next_lane": lattice["selected_next_lane"],
        "toe_closed_now": False,
        "agents_reoriented": True,
    }


def update_agents() -> None:
    section = """## Current state-map-first steering guardrail (P2681/S1631, 2026-06-13)\n\n- Do not treat the older `legacy -> strict completion bridge` priority as an automatic next-step generator.  It remains a major theoretical problem, but P2679/P2680 have already converted the generic bridge lane into a finite list of missing source atoms; repeating the generic bridge audit is now a loop unless a new typed object, theorem, or source atom is introduced.\n- Before selecting a new FAR move, run or consult a broad state-map audit across bridge, selector, tau/pair, damping, Lagrangian/EOM/QG, direct-route defect, and coupling-law lanes.  The next move must be chosen from the current live frontier, not from stale priority text alone.\n- Closed/repetition-gated lanes remain closed without new evidence: selector/orientation/QW-2191 replay, `tau_src -> pair12 -> boundary-square`, and `beta_tors -> chi11` are not admissible merely because they are mentioned in older bridge notes.\n- Current highest practical proof-grade opportunity after P2681 is a finite, executable live-frontier obstruction/witness computation: first preference is the kernel-split-robust P46/N49 direct-route zero-witness/no-go matrix (especially the direct `m2 psi4` target action coefficient defect on common `psi4**2/2` support named by F3); second preference is a strict Lagrangian/EOM reverse-closure obstruction matrix.\n- A `legacy -> strict` move is admissible now only if it targets exactly one missing P2680 non-selector source atom, such as target-independent positive `beta/Z_beta` or canonical length/UV unit, and explicitly avoids generic bridge repetition.\n- Role-transfer auditing is still downstream and must not start until bridge/source closure is actually exported.\n"""
    append_once(AGENTS, "P2681/S1631", section)


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2681/S1631 professorial research-state map and AGENTS reorientation audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Broad content-first grep",
    ]
    for name, data in payload["broad_state_grep"]["patterns"].items():
        lines.append(f"- `{name}`: `{data['count']}` hits")
    lines.extend(["", "## Prior lane ledger"])
    for row in payload["prior_lane_ledger"]:
        lines.append(f"- `{row['lane']}`: status=`{row['current_status']}`; repeat_now=`{row['repeat_now_admissible']}`; leverage=`{row['scientific_leverage_now']}`; done={row['already_done']}")
    lines.extend(["", "## Opportunity matrix"])
    for row in sorted(payload["opportunity_matrix"], key=lambda r: r["net_score"], reverse=True):
        lines.append(f"- `{row['candidate']}`: net=`{row['net_score']}` (leverage={row['toe_closure_leverage_1_to_5']}, computability={row['proof_computability_1_to_5']}, evidence={row['evidence_base_1_to_5']}, loop_risk={row['loop_risk_1_to_5']}) — {row['rationale']}")
    lat = payload["finite_priority_lattice"]
    lines.extend([
        "", "## Finite priority lattice",
        f"Total states: `{lat['total_states']}`; passing states: `{lat['passing_states']}`.",
        f"Passing candidates: `{lat['passing_candidates']}`.",
        f"Selected next lane: `{lat['selected_next_lane']}` — {lat['selected_reason']}",
        "", "## AGENTS.md update", "`AGENTS.md` has been reoriented with the P2681/S1631 state-map-first guardrail.",
        "", "## Verdict", payload["closure_decision"]["professorial_verdict"],
        f"Decision: `{payload['closure_decision']['decision']}`.",
        "", "## Next honest step", payload["closure_decision"]["next_honest_step"],
        "", "## Negative exports",
    ])
    for key, value in payload["negative_export_flags"].items():
        lines.append(f"- `{key}`: `{value}`")
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = broad_state_grep()
    ledger = prior_lane_ledger()
    matrix = opportunity_matrix()
    lattice = finite_priority_lattice()
    update_agents()
    payload: dict[str, Any] = {
        "status": "P2681_PROFESSORIAL_RESEARCH_STATE_MAP_AND_AGENTS_REORIENTATION_AUDIT_NO_FALSE_PASS",
        "broad_state_grep": grep,
        "p2680_consistency": p2680_consistency(),
        "prior_lane_ledger": ledger,
        "opportunity_matrix": matrix,
        "finite_priority_lattice": lattice,
        "closure_decision": closure_decision(matrix, lattice),
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
        "agents_sha256_after_update": sha256_file(AGENTS),
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2681/S1631 state-map-first steering guard",
        "## P2681/S1631 state-map-first steering guard\n\n"
        "`P2681/S1631` supersedes stale priority-by-AGENTS behavior for current execution.  Generic `legacy -> strict` bridge looping is no longer the automatic next move after P2679/P2680; selector, `tau_src -> pair12`, and `beta_tors -> chi11` remain repetition-gated without a new typed object.  The current highest practical proof-grade lane is a finite P46/N49 zero-witness/no-go matrix, especially the direct `m2 psi4` target action coefficient defect on common `psi4**2/2` support, with strict no-closure guards.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2681/S1631 state-map-first Ltotal steering guard",
        "## P2681/S1631 state-map-first Ltotal steering guard\n\n"
        "`P2681/S1631` does not promote `L_total` or role transfer.  It reorients future work toward finite live-frontier obstruction/witness matrices: first P46/N49 direct-route zero-witness/no-go, second strict Lagrangian/EOM reverse-closure obstruction if the direct symbolic carrier is unavailable.\n",
    )
    return payload


if __name__ == "__main__":
    main()
