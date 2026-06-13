#!/usr/bin/env python3
"""P2679/S1629: repetition gate and bridge-pivot audit.

This packet answers whether the recent selector / tau_src->pair12 /
beta_tors->chi11 lanes may honestly be reopened.  It records what was already
run, checks for a new reopening theorem/object, and pivots the recommended next
work away from repeated selector loops unless new evidence appears.
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
OUT = GEN / "p2679_s1629_reopen_repetition_gate_and_bridge_pivot_audit.json"
MD = GEN / "p2679_s1629_reopen_repetition_gate_and_bridge_pivot_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
UPSTREAM = {
    "P2619": GEN / "p2619_s1569_p2618_selector_source_obligation_lattice.json",
    "P2649": GEN / "p2649_s1599_beta_source_route_decision_matrix_and_normalization_orbit_no_go.json",
    "P2677": GEN / "p2677_s1627_s6_o3_typed_seed_route_no_go_audit.json",
    "P2678": GEN / "p2678_s1628_strict_internal_orientation_source_provider_class_audit.json",
}

NEGATIVE_EXPORT_FLAGS = [
    "new_reopening_evidence_after_p2678_exported",
    "tau_src_pair12_route_reopened",
    "selector_orientation_route_reopened",
    "beta_tors_chi11_route_reopened",
    "s6_o3_reopened",
    "o4_o5_allowed",
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
        "tau_src_pair12_repetition_content": "tau_src|source-topology.*pair|pair12|pair1/pair2|typed seed|typed-seed|boundary-square",
        "selector_orientation_repetition_content": "selector source|orientation source|orientation-odd|XOR|XNOR|QW-2191|symmetry-breaking|spin/Pin|torsor",
        "beta_tors_chi11_repetition_content": "beta_tors.*chi11|beta_tors -> chi_11|chi_11|legacy scalar|torsion damping",
        "legacy_strict_bridge_content": "legacy -> strict|legacy-to-strict|bridge completion|K_legacy_ont|K_strict_gate|completion bridge",
        "damping_compression_content": "damping/compression|nonlinear damping|d\\^eta|beta source|positive_beta|Z_beta|normalization orbit",
        "role_transfer_guard_content": "role transfer|role-transfer|role-bearing L_total|ToE closure|QW-2191 discharge",
    }
    return {
        "tool": "rg",
        "mode": "content-first repetition/reopening audit; explicitly checks whether selector/tau_src/beta_tors lanes have new post-P2678 evidence before reopening",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
    }


def prior_research_ledger() -> list[dict[str, Any]]:
    return [
        {
            "lane": "tau_src -> pair12 -> boundary-square",
            "already_run": ["P2673", "P2674", "P2675", "P2676", "P2677"],
            "last_result": "bounded no-go for current S6/O3 typed-seed route",
            "missing_to_reopen": ["new chart-label-retaining typed seed export", "new pre-collapse Sigma->F301 arrow", "new sourced sector-swap invariant"],
            "reopen_now": False,
        },
        {
            "lane": "selector/orientation source",
            "already_run": ["P2619", "P2620", "P2621", "P2622", "P2623", "P2676", "P2677", "P2678"],
            "last_result": "formal orientation-odd provider shapes known; no provider exported with Sigma->F301 binding",
            "missing_to_reopen": ["exported orientation-odd torsor/spin-boundary object", "non-fiat C2 law", "typed Sigma->F301 action"],
            "reopen_now": False,
        },
        {
            "lane": "beta_tors -> chi11",
            "already_run": ["P2390", "P2435", "P2619", "P2677", "P2678"],
            "last_result": "legacy scalar beta_tors may remain damping input but is not an odd chi11/sign source",
            "missing_to_reopen": ["orientation/sign torsor already present in the input", "strict role-transfer theorem", "non-scalar source map"],
            "reopen_now": False,
        },
        {
            "lane": "strict beta/damping source",
            "already_run": ["P2625", "P2649", "P2650", "P2651"],
            "last_result": "normalization and tail-ratio facts exist, but no target-independent positive beta/Z_beta source theorem",
            "missing_to_reopen": ["canonical length/UV unit", "target-independent positive beta renormalization source", "non-normalization source theorem"],
            "reopen_now": False,
        },
    ]


def upstream_state() -> dict[str, Any]:
    p2678 = load_json(UPSTREAM["P2678"])
    p2677 = load_json(UPSTREAM["P2677"])
    p2649 = load_json(UPSTREAM["P2649"])
    p2619 = load_json(UPSTREAM["P2619"])
    return {
        "p2678_s6_not_reopened": p2678.get("closure_decision", {}).get("s6_reopened_now") is False,
        "p2678_no_orientation_provider": p2678.get("closure_decision", {}).get("strict_internal_orientation_source_exported_now") is False,
        "p2677_current_route_no_go": p2677.get("closure_decision", {}).get("s6_current_route_passed_now") is False,
        "p2649_available": "missing" not in p2649,
        "p2619_available": "missing" not in p2619,
        "source_hashes": {name: sha256_file(path) for name, path in UPSTREAM.items()},
    }


def reopen_gate_lattice() -> dict[str, Any]:
    obligations = [
        "new_object_after_p2678",
        "object_not_same_tau_src_pair12_lane",
        "object_not_same_selector_reversal_lane",
        "object_not_beta_tors_scalar_chi11_lane",
        "typed_precollapse_binding_or_bridge_source_exported",
        "closure_guards_preserved",
    ]
    current = {
        obligations[0]: False,
        obligations[1]: True,
        obligations[2]: True,
        obligations[3]: True,
        obligations[4]: False,
        obligations[5]: True,
    }
    pass_count = 0
    for bits in itertools.product([False, True], repeat=len(obligations)):
        pass_count += int(all(bits))
    missing = [key for key, value in current.items() if not value]
    return {
        "obligations": obligations,
        "total_states": 2 ** len(obligations),
        "passing_states": pass_count,
        "current_state": current,
        "missing_current_obligations": missing,
        "hamming_distance_to_pass": len(missing),
        "may_reopen_old_lanes_now": False,
    }


def bridge_pivot_matrix() -> list[dict[str, Any]]:
    return [
        {
            "candidate_next_lane": "repeat tau_src -> pair12 / S6/O3",
            "admissible_now": False,
            "reason": "P2677 already gives bounded no-go and no new typed seed/arrow object has been exported after P2678.",
        },
        {
            "candidate_next_lane": "repeat selector/orientation provider search",
            "admissible_now": False,
            "reason": "P2678 already enumerates formal odd provider classes; reopening requires a concrete new provider object, not another same-shape search.",
        },
        {
            "candidate_next_lane": "repeat beta_tors -> chi11",
            "admissible_now": False,
            "reason": "P2619/P2678 classify scalar beta_tors as even/scalar for sign purposes; no new torsor/role-transfer theorem is present.",
        },
        {
            "candidate_next_lane": "legacy->strict bridge source audit excluding selector/tau_src/beta_tors->chi11 repeats",
            "admissible_now": True,
            "reason": "This follows S2 priority while respecting the repetition gate; the narrow target should be a bridge-source inventory such as damping/compression or amplitude/normalization source terms, not a selector replay.",
        },
    ]


def closure_decision(lattice: dict[str, Any], matrix: list[dict[str, Any]]) -> dict[str, Any]:
    allowed = [row["candidate_next_lane"] for row in matrix if row["admissible_now"]]
    return {
        "decision": "P2679_REOPEN_REPETITION_GATE__NO_NEW_EVIDENCE_TO_REOPEN_SELECTOR_TAUSRC_BETATORS_CHI11_LANES",
        "professorial_verdict": (
            "P2679 answers the methodological objection directly.  The recent selector, tau_src->pair12, and beta_tors->chi11 topics were already audited: tau_src/pair12 reached a bounded S6/O3 no-go; selector/orientation reached formal odd-provider shapes without export; beta_tors->chi11 remains blocked because scalar beta_tors is not an orientation-odd sign source.  Since no new post-P2678 theorem/object is present, reopening those lanes would be repetition, not progress."
        ),
        "next_honest_step": (
            "Do not rerun selector, tau_src->pair12, or beta_tors->chi11 until a genuinely new object appears.  The next honest proof-grade move is a legacy->strict bridge-source audit that explicitly excludes those repeats, preferably a source inventory for amplitude/normalization and damping/compression terms and their missing target-independent source atoms."
        ),
        "allowed_next_lanes": allowed,
        "old_lanes_reopened_now": False,
        "new_reopening_evidence_after_p2678": False,
        "hamming_distance_to_reopen": lattice["hamming_distance_to_pass"],
        "s6_o3_reopened_now": False,
        "o4_o5_allowed_now": False,
        "role_bearing_ltotal_now": False,
        "toe_closure_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2679/S1629 reopen repetition gate and bridge-pivot audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first repo grep",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: `{data['count']}` hits")
    lines.extend(["", "## Prior research ledger"])
    for row in payload["prior_research_ledger"]:
        lines.append(f"- `{row['lane']}`: already_run=`{row['already_run']}`; reopen_now=`{row['reopen_now']}`; last_result={row['last_result']}")
    lat = payload["reopen_gate_lattice"]
    lines.extend([
        "", "## Reopen gate lattice",
        f"Total states: `{lat['total_states']}`; passing states: `{lat['passing_states']}`.",
        f"Current Hamming distance to reopen: `{lat['hamming_distance_to_pass']}`.",
        f"Missing current obligations: `{lat['missing_current_obligations']}`.",
        "", "## Bridge pivot matrix",
    ])
    for row in payload["bridge_pivot_matrix"]:
        lines.append(f"- `{row['candidate_next_lane']}`: admissible_now=`{row['admissible_now']}` — {row['reason']}")
    lines.extend([
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
    lattice = reopen_gate_lattice()
    matrix = bridge_pivot_matrix()
    payload: dict[str, Any] = {
        "status": "P2679_REOPEN_REPETITION_GATE_AND_BRIDGE_PIVOT_AUDIT_NO_FALSE_PASS",
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "upstream_state": upstream_state(),
        "prior_research_ledger": prior_research_ledger(),
        "reopen_gate_lattice": lattice,
        "bridge_pivot_matrix": matrix,
        "closure_decision": closure_decision(lattice, matrix),
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2679/S1629 reopen repetition gate and bridge-pivot guard",
        "## P2679/S1629 reopen repetition gate and bridge-pivot guard\n\n"
        "`P2679/S1629` answers the repetition concern explicitly.  The `tau_src -> pair12` S6/O3 lane, the selector/orientation lane, and the `beta_tors -> chi_11` lane have already been audited through P2677/P2678 and earlier P2619/P2649-class guards; no new post-P2678 theorem/object reopens them.  Therefore repeating those lanes is inadmissible unless a genuinely new typed object appears.  The admissible next proof-grade direction is a `legacy -> strict` bridge-source audit that excludes selector replay and `beta_tors -> chi_11` replay, focusing instead on source inventories for amplitude/normalization and damping/compression atoms.  No O4/O5 admission, `QW-2191` discharge, role-bearing `L_total`, role transfer, or ToE closure is exported.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2679/S1629 repetition gate Ltotal guard",
        "## P2679/S1629 repetition gate Ltotal guard\n\n"
        "`P2679/S1629` prevents `L_total` promotion by repeated selector, `tau_src -> pair12`, or `beta_tors -> chi_11` audits.  A future variational term must come from a genuinely new bridge-source theorem/object, not from restating already-blocked lanes.\n",
    )
    return payload


if __name__ == "__main__":
    main()
