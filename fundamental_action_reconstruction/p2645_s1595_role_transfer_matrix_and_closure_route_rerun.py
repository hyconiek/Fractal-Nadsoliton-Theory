#!/usr/bin/env python3
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
OUT = GEN / "p2645_s1595_role_transfer_matrix_and_closure_route_rerun.json"
MD = GEN / "p2645_s1595_role_transfer_matrix_and_closure_route_rerun.md"

SOURCES = {
    "P2629_ZBETA_GAUGE": GEN / "p2629_s1579_zbeta_normalization_gauge_obstruction.json",
    "P2642_NODE_DEMOTION": GEN / "p2642_s1592_universal_affine_zero_lattice_source_no_go_and_node_role_demotion.json",
    "P2643_BETA_THRESHOLD": GEN / "p2643_s1593_inverse_hierarchy_beta_threshold_role_rejection_certificate.json",
    "P2644_COMPRESSED_SUCCESSOR": GEN / "p2644_s1594_modified_compressed_inverse_hierarchy_successor_theorem.json",
    "STRICT_EQUATION_SHEET": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "STRICT_LAGRANGIAN_DRAFT": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

NEGATIVE_EXPORT_FLAGS = [
    "full_legacy_role_transfer_revalidated",
    "legacy_integer_node_gauge_role_reopened",
    "unchanged_inverse_hierarchy_role_reopened",
    "beta_source_exported",
    "alpha_geo_role_source_exported",
    "ew_angle_role_transferred",
    "alpha_em_role_transferred",
    "strict_kernel_full_kernel_claimed",
    "role_bearing_ltotal_reenabled",
    "qw2191_discharged",
    "toe_closure_claimed",
]

ROLE_CANDIDATES = [
    {
        "role": "legacy_integer_node_gauge_role",
        "legacy_claim": "integer node/gauge coordinate repair from zero-lattice / affine node list",
        "required_atoms": {"origin_source", "selector_free_node_source", "strict_topology_safe_coordinate"},
    },
    {
        "role": "unchanged_inverse_hierarchy_role",
        "legacy_claim": "distant-octave / beta^N hierarchy amplification transfers unchanged",
        "required_atoms": {"beta_below_inverse_threshold", "role_preserving_damping", "node_role_not_demoted"},
    },
    {
        "role": "modified_compressed_successor_role",
        "legacy_claim": "strict denominator is a compression/locality-bias successor, not unchanged amplification",
        "required_atoms": {"monotone_compression_theorem", "explicit_successor_semantics"},
    },
    {
        "role": "strict_beta_source_role",
        "legacy_claim": "beta=1 / Z_beta=100 is sourced target-independently rather than normalized by convention",
        "required_atoms": {"target_independent_beta_identity", "normalization_gauge_fixed", "micro_strict_mismatch_removed"},
    },
    {
        "role": "legacy_ew_angle_alpha_geo_role",
        "legacy_claim": "sin^2(theta_W)=alpha_geo/12 transfers to strict layer",
        "required_atoms": {"alpha_geo_survives_completion", "ew_operator_source", "full_bridge_completion", "role_transfer_theorem"},
    },
    {
        "role": "legacy_alpha_em_beta_tors_role",
        "legacy_claim": "alpha_EM^-1 = alpha_geo/(2*beta_tors)*(1-beta_tors) transfers to strict layer",
        "required_atoms": {"alpha_geo_survives_completion", "beta_tors_maps_to_strict_beta", "em_operator_source", "full_bridge_completion", "role_transfer_theorem"},
    },
]

CURRENT_ATOMS = {
    "origin_source": False,
    "selector_free_node_source": False,
    "strict_topology_safe_coordinate": False,
    "beta_below_inverse_threshold": False,
    "role_preserving_damping": False,
    "node_role_not_demoted": False,
    "monotone_compression_theorem": True,
    "explicit_successor_semantics": True,
    "target_independent_beta_identity": False,
    "normalization_gauge_fixed": False,
    "micro_strict_mismatch_removed": False,
    "alpha_geo_survives_completion": False,
    "ew_operator_source": False,
    "full_bridge_completion": False,
    "role_transfer_theorem": False,
    "beta_tors_maps_to_strict_beta": False,
    "em_operator_source": False,
    "selector_source_qw2191": False,
    "blind_frozen_empirical_confirmation": False,
}


def sha256_file(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, ".",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "*.json", "-g", "*.lean",
            "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**",
        ],
        cwd=REPO,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:40]}


def semantic_rg_audit() -> dict[str, Any]:
    patterns = {
        "role_transfer_matrix_content": (
            "role-transfer matrix|role transfer audit|role-transfer theorem|legacy physical role|"
            "survives unchanged|modified/compressed|successor statement|rejected"
        ),
        "node_demotion_origin_source_content": (
            "zero-lattice origin|origin-source|hidden selector|node/gauge|affine.*repair|"
            "coordinate repair|demoted"
        ),
        "inverse_hierarchy_compression_content": (
            "inverse hierarchy|distant-octave|modified/compressed successor|compression/locality-bias|"
            "unchanged inverse-hierarchy|beta_crit"
        ),
        "beta_alpha_role_content": (
            r"beta source|Z_beta|normalization gauge|alpha_geo|sin\^2\(theta_W\)|alpha_EM|beta_tors"
        ),
        "closure_guard_content": (
            "full bridge completion|role-bearing L_total|QW-2191|selector source|ToE closure|"
            "blind frozen-kernel empirical"
        ),
    }
    return {"tool": "rg", "mode": "content-first semantic audit for post-P2644 role-transfer matrix rerun", "patterns": {k: rg_count(v) for k, v in patterns.items()}}


def upstream_consistency() -> dict[str, Any]:
    p2629 = load_json(SOURCES["P2629_ZBETA_GAUGE"])
    p2642 = load_json(SOURCES["P2642_NODE_DEMOTION"])
    p2643 = load_json(SOURCES["P2643_BETA_THRESHOLD"])
    p2644 = load_json(SOURCES["P2644_COMPRESSED_SUCCESSOR"])
    return {
        "p2629_beta_source_exported": p2629.get("exact_source_gate", {}).get("all_exact_source_gates_passed"),
        "p2629_invariant_ratio": p2629.get("normalization_orbit_certificate", {}).get("invariant_ratio"),
        "p2642_node_role_status": p2642.get("demotion_decision", {}).get("legacy_integer_node_gauge_role_status"),
        "p2643_inverse_role_status": p2643.get("role_verdict", {}).get("unchanged_inverse_hierarchy_role_status"),
        "p2644_successor_verdict": p2644.get("role_successor_decision", {}).get("legacy_role_transfer_verdict"),
        "p2644_suppression_ratio_d7_over_d1": p2644.get("compression_successor_theorem", {}).get("d7_over_d1_suppression_ratio"),
    }


def evaluate_roles() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for role in ROLE_CANDIDATES:
        missing = sorted(atom for atom in role["required_atoms"] if not CURRENT_ATOMS[atom])
        passed = not missing
        if role["role"] == "modified_compressed_successor_role" and passed:
            verdict = "SURVIVES_ONLY_AS_MODIFIED_COMPRESSED_SUCCESSOR_DESCRIPTIVE_NOT_FULL_TRANSFER"
        elif role["role"] in {"legacy_integer_node_gauge_role", "unchanged_inverse_hierarchy_role"}:
            verdict = "REJECTED_OR_DEMOTED_UNDER_CURRENT_STRICT_COMPLETION_AUDIT"
        else:
            verdict = "BLOCKED_PENDING_EXPLICIT_SOURCE_AND_ROLE_TRANSFER_THEOREM"
        rows.append({
            "role": role["role"],
            "legacy_claim": role["legacy_claim"],
            "required_atoms": sorted(role["required_atoms"]),
            "missing_atoms": missing,
            "passes_current_gate": passed,
            "verdict": verdict,
        })
    return rows


def minimal_unblock_sets(rows: list[dict[str, Any]]) -> dict[str, Any]:
    blocked_roles = [row for row in rows if not row["passes_current_gate"]]
    universe = sorted({atom for row in blocked_roles for atom in row["missing_atoms"]})
    # Compute exact minimal atom additions that would make at least one blocked role pass.
    per_role_minimal = {row["role"]: row["missing_atoms"] for row in blocked_roles}
    # For full legacy transfer all blocked legacy rows would need every missing atom; this is not a recommendation.
    full_transfer_union = sorted(set(itertools.chain.from_iterable(per_role_minimal.values())))
    priority_sets = [
        {
            "route": "beta_source_route",
            "new_atoms_needed": ["target_independent_beta_identity", "normalization_gauge_fixed", "micro_strict_mismatch_removed"],
            "unblocks": ["strict_beta_source_role"],
        },
        {
            "route": "empirical_compression_route",
            "new_atoms_needed": ["blind_frozen_empirical_confirmation"],
            "unblocks": ["empirical support for modified_compressed_successor_role only; not role transfer"],
        },
        {
            "route": "alpha_role_transfer_route",
            "new_atoms_needed": ["alpha_geo_survives_completion", "full_bridge_completion", "role_transfer_theorem", "ew_operator_source", "em_operator_source", "beta_tors_maps_to_strict_beta"],
            "unblocks": ["legacy_ew_angle_alpha_geo_role", "legacy_alpha_em_beta_tors_role"],
        },
    ]
    return {
        "blocked_role_count": len(blocked_roles),
        "passing_role_count": len(rows) - len(blocked_roles),
        "all_missing_atoms_union": universe,
        "per_role_minimal_atom_sets": per_role_minimal,
        "full_legacy_transfer_union_not_recommended": full_transfer_union,
        "professorial_priority_routes": priority_sets,
    }


def closure_decision(rows: list[dict[str, Any]], unblock: dict[str, Any]) -> dict[str, Any]:
    full_transfer = all(row["passes_current_gate"] for row in rows)
    beta_ready = CURRENT_ATOMS["target_independent_beta_identity"] and CURRENT_ATOMS["normalization_gauge_fixed"] and CURRENT_ATOMS["micro_strict_mismatch_removed"]
    return {
        "full_legacy_role_transfer_now": full_transfer,
        "role_bearing_ltotal_now": False,
        "toe_closure_now": False,
        "q_w_2191_discharged_now": False,
        "accepted_current_strict_role_entry": "modified_compressed_successor_role only, as descriptive compression/locality-bias semantics",
        "rejected_current_entries": [row["role"] for row in rows if row["verdict"].startswith("REJECTED")],
        "blocked_current_entries": [row["role"] for row in rows if row["verdict"].startswith("BLOCKED")],
        "beta_source_route_ready": beta_ready,
        "professorial_verdict": (
            "P2645 reruns the post-P2644 role-transfer matrix as gates, not prose.  Exactly one audited entry passes now: "
            "modified/compressed successor semantics for the strict denominator.  The old node/gauge role and unchanged inverse-hierarchy role are not recoverable under current evidence.  "
            "The alpha_geo/EW and alpha_EM/beta_tors roles remain blocked because the completion map, operator sources, and beta_tors->strict_beta map are absent.  "
            "Therefore strict remains a robust working kernel, not a full ToE kernel or role-bearing L_total source."
        ),
        "next_honest_step": (
            "Do not reopen node-lifts or unchanged inverse-hierarchy.  The next proof-grade move is a target-independent beta-source identity if one exists; "
            "if not, preregister a blind frozen-kernel empirical compression test while marking beta/alpha legacy roles as blocked and the compression role as descriptive only."
        ),
        "unblock_summary": unblock,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2645/S1595 role-transfer matrix and closure-route rerun",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first anti-duplication audit",
        "",
        "This audit greps role-transfer, origin-source/node demotion, inverse-hierarchy compression, beta/alpha role, and closure-guard content before rerunning the matrix.",
        "",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines.extend(["", "## Gate matrix", "", "| role | gate | missing atoms | verdict |", "| --- | --- | --- | --- |"])
    for row in payload["role_transfer_matrix"]:
        gate = "PASS" if row["passes_current_gate"] else "FAIL"
        missing = ", ".join(row["missing_atoms"]) if row["missing_atoms"] else "none"
        lines.append(f"| `{row['role']}` | `{gate}` | `{missing}` | `{row['verdict']}` |")
    dec = payload["closure_decision"]
    lines.extend([
        "",
        "## Closure decision",
        "",
        dec["professorial_verdict"],
        "",
        f"Full legacy role transfer now: `{dec['full_legacy_role_transfer_now']}`.",
        f"Role-bearing L_total now: `{dec['role_bearing_ltotal_now']}`.",
        f"QW-2191 discharged now: `{dec['q_w_2191_discharged_now']}`.",
        f"ToE closure now: `{dec['toe_closure_now']}`.",
        "",
        "## Professorial route map",
        "",
    ])
    for route in dec["unblock_summary"]["professorial_priority_routes"]:
        lines.append(f"- `{route['route']}` needs `{', '.join(route['new_atoms_needed'])}` and unblocks `{', '.join(route['unblocks'])}`.")
    lines.extend(["", "## Next honest step", "", dec["next_honest_step"], "", "## Negative exports", ""])
    for key, value in payload["negative_export_flags"].items():
        lines.append(f"- `{key}`: `{value}`")
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    rows = evaluate_roles()
    unblock = minimal_unblock_sets(rows)
    decision = closure_decision(rows, unblock)
    payload: dict[str, Any] = {
        "status": "P2645_ROLE_TRANSFER_MATRIX_RERUN_ONE_MODIFIED_SUCCESSOR_PASS_NO_FULL_TRANSFER_NO_LTOTAL_NO_QW2191_NO_TOE",
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_hashes": {name: sha256_file(path) for name, path in SOURCES.items()},
        "upstream_consistency": upstream_consistency(),
        "current_atoms": CURRENT_ATOMS,
        "role_transfer_matrix": rows,
        "closure_decision": decision,
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({k: v for k, v in payload.items() if k != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        SOURCES["STRICT_EQUATION_SHEET"],
        "P2645/S1595 role-transfer matrix rerun guard",
        "## P2645/S1595 role-transfer matrix rerun guard\n\n"
        "`P2645/S1595` reruns the legacy-role transfer table after P2642-P2644.  The only current pass is the modified/compressed successor reading of the strict denominator as monotone compression/locality-bias.  The legacy integer node/gauge role and unchanged inverse-hierarchy transfer are rejected/demoted, while `alpha_geo` EW and `alpha_EM/beta_tors` roles remain blocked pending explicit completion-map, operator-source, and beta-source theorems.  Thus bridge completion, full role transfer, role-bearing `L_total`, `QW-2191`, and ToE closure remain closed.\n",
    )
    append_once(
        SOURCES["STRICT_LAGRANGIAN_DRAFT"],
        "P2645/S1595 role-transfer matrix Ltotal guard",
        "## P2645/S1595 role-transfer matrix Ltotal guard\n\n"
        "`P2645/S1595` prevents promoting the post-P2644 compression semantics into a full variational source.  The matrix has only a descriptive modified-successor pass and still lacks target-independent beta sourcehood, alpha-role transfer, completion-map closure, and a `QW-2191` selector source; `L_total` remains non-role-bearing.\n",
    )
    return payload


if __name__ == "__main__":
    result = main()
    print(json.dumps({"status": result["status"], "out": rel(OUT), "md": rel(MD), "next": result["closure_decision"]["next_honest_step"]}, indent=2, ensure_ascii=False, sort_keys=True))
