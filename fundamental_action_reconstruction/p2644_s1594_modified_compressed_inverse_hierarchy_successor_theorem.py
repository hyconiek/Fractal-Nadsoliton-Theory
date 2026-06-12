#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import re
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2644_s1594_modified_compressed_inverse_hierarchy_successor_theorem.json"
MD = GEN / "p2644_s1594_modified_compressed_inverse_hierarchy_successor_theorem.md"

BETA_TORS = 0.01
STRICT_BETA = 1.0
ETA = 9.0 / 5.0
OMEGA = math.pi / 4.0
PHI = math.pi / 6.0
GRID = [1, 2, 3, 5, 7, 8, 10, 11, 12]

SOURCE_FILES = {
    "P2633_DIAGRAM_RETENTION": GEN / "p2633_s1583_diagram_grounded_strict_kernel_characteristic_preservation_audit.json",
    "P2642_NODE_DEMOTION": GEN / "p2642_s1592_universal_affine_zero_lattice_source_no_go_and_node_role_demotion.json",
    "P2643_BETA_THRESHOLD": GEN / "p2643_s1593_inverse_hierarchy_beta_threshold_role_rejection_certificate.json",
    "STRICT_EQUATION_SHEET": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "STRICT_LAGRANGIAN_DRAFT": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

NEGATIVE_EXPORT_FLAGS = [
    "unchanged_inverse_hierarchy_role_restored",
    "modified_successor_promoted_to_full_legacy_transfer",
    "beta_source_exported",
    "legacy_integer_node_gauge_role_reopened",
    "legacy_to_strict_bridge_completion_exported",
    "legacy_role_transfer_revalidated",
    "strict_kernel_full_kernel_claimed",
    "toe_closure_claimed",
    "selector_source_exported",
    "q_w_2191_discharged",
    "role_bearing_ltotal_reenabled",
    "blind_empirical_confirmation_claimed",
]


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
    return {"count": len(lines), "samples": lines[:35]}


def semantic_rg_audit() -> dict[str, Any]:
    patterns = {
        "modified_compressed_successor_content": (
            "modified/compressed successor|compressed successor|modified.*successor|fractal compression|"
            "strict compression|heavy-tailed compression|UV/local compression"
        ),
        "inverse_hierarchy_rejection_content": (
            "inverse hierarchy|distant octave|Wilson-loop|Wilson loop|unchanged.*transfer|rejected unchanged|"
            "K\\(7\\).*K\\(1\\)"
        ),
        "attention_suppression_content": (
            "attention denominator|strict attention|legacy attention|suppression factor|compression factor|"
            "1\\+beta_tors|d\\^eta|d\\^\\(9/5\\)"
        ),
        "source_and_ltotal_guard_content": (
            "beta source|target-independent|role-bearing L_total|QW-2191|selector source|ToE closure|"
            "bridge completion"
        ),
    }
    return {"tool": "rg", "mode": "content-first semantic audit for modified/compressed inverse-hierarchy successor theorem", "patterns": {k: rg_count(v) for k, v in patterns.items()}}


def latest_commit_audit() -> list[dict[str, Any]]:
    proc = subprocess.run(["git", "log", "-n", "8", "--oneline", "--name-only"], cwd=REPO, text=True, stdout=subprocess.PIPE, check=True)
    rows: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None
    for line in proc.stdout.splitlines():
        if not line.strip():
            continue
        if re.match(r"^[0-9a-f]{7,12} ", line):
            if current:
                rows.append(current)
            sha, subject = line.split(" ", 1)
            current = {"sha": sha, "subject": subject, "files": []}
        elif current is not None:
            current["files"].append(line)
    if current:
        rows.append(current)
    return rows


def phase_abs(d: float) -> float:
    return abs(math.cos(OMEGA * d + PHI))


def legacy_attention(d: float) -> float:
    return 1.0 / (1.0 + BETA_TORS * d)


def strict_attention(d: float) -> float:
    return 1.0 / (1.0 + STRICT_BETA * d**ETA)


def suppression_factor(d: float) -> float:
    return strict_attention(d) / legacy_attention(d)


def suppression_derivative_numerator(d: float) -> float:
    # derivative of (1+b_t d)/(1+d^eta), multiplied by positive denominator^2
    return BETA_TORS * (1.0 + d**ETA) - ETA * d ** (ETA - 1.0) * (1.0 + BETA_TORS * d)


def compression_theorem() -> dict[str, Any]:
    rows = []
    for d in GRID:
        rows.append({
            "d": d,
            "phase_abs": phase_abs(d),
            "legacy_attention": legacy_attention(d),
            "strict_attention": strict_attention(d),
            "strict_over_legacy_attention_suppression": suppression_factor(d),
            "derivative_numerator": suppression_derivative_numerator(d),
        })
    return {
        "suppression_factor": "S(d)=(1+0.01*d)/(1+d^(9/5)) = A_strict(d)/A_legacy(d)",
        "analytic_derivative_numerator": "S'(d) has sign 0.01*(1+d^(9/5))-(9/5)*d^(4/5)*(1+0.01*d)",
        "monotonicity_proof_on_d_ge_1": "For d>=1, numerator <= 0.01 - 9/5 + 0.01*(1-9/5)*d^(9/5) < 0, so S is strictly decreasing.",
        "grid_rows": rows,
        "all_grid_suppression_below_one": all(row["strict_over_legacy_attention_suppression"] < 1.0 for row in rows),
        "all_grid_derivative_negative": all(row["derivative_numerator"] < 0.0 for row in rows),
        "d1_suppression": suppression_factor(1.0),
        "d7_suppression": suppression_factor(7.0),
        "d7_over_d1_suppression_ratio": suppression_factor(7.0) / suppression_factor(1.0),
    }


def upstream_consistency() -> dict[str, Any]:
    p2633 = load_json(SOURCE_FILES["P2633_DIAGRAM_RETENTION"])
    p2642 = load_json(SOURCE_FILES["P2642_NODE_DEMOTION"])
    p2643 = load_json(SOURCE_FILES["P2643_BETA_THRESHOLD"])
    return {
        "p2633_strict_ratio": p2633.get("finite_witness", {}).get("inverse_hierarchy_ratio_abs_k7_over_abs_k1", {}).get("strict"),
        "p2633_legacy_ratio": p2633.get("finite_witness", {}).get("inverse_hierarchy_ratio_abs_k7_over_abs_k1", {}).get("legacy_amplitude_normalized"),
        "p2642_node_role_status": p2642.get("demotion_decision", {}).get("legacy_integer_node_gauge_role_status"),
        "p2643_inverse_role_status": p2643.get("role_verdict", {}).get("unchanged_inverse_hierarchy_role_status"),
        "p2643_allowed_successor_status": p2643.get("role_verdict", {}).get("allowed_successor_status"),
        "p2643_beta_crit": p2643.get("beta_threshold_theorem", {}).get("beta_critical_exact_role_boundary"),
    }


def role_successor_decision(theorem: dict[str, Any]) -> dict[str, Any]:
    gates = {
        "suppression_factor_proven_monotone_decreasing_on_d_ge_1": theorem["all_grid_derivative_negative"],
        "strict_attention_below_legacy_attention_on_grid": theorem["all_grid_suppression_below_one"],
        "distant_octave_suppressed_relative_to_near": theorem["d7_over_d1_suppression_ratio"] < 1.0,
        "unchanged_inverse_hierarchy_reopened": False,
        "beta_source_available": False,
    }
    return {
        "gates": gates,
        "modified_successor_statement": "Strict beta=1, eta=9/5 converts the legacy distant-octave amplification role into a monotone fractal-compression / locality-bias role on d>=1.",
        "legacy_role_transfer_verdict": "UNCHANGED_INVERSE_HIERARCHY_REJECTED__MODIFIED_COMPRESSED_SUCCESSOR_ACCEPTED_AS_DESCRIPTIVE_NOT_FULL_ROLE_TRANSFER",
        "professorial_verdict": (
            "P2644 proves the honest successor semantics left open by P2643: the strict denominator is not an inverse-hierarchy carrier; "
            "it is a monotone compression operator relative to the legacy hyperbolic denominator.  This supports a modified/compressed successor reading, "
            "but it still does not export beta sourcehood, role-bearing L_total, QW-2191 discharge, or ToE closure."
        ),
        "next_honest_step": (
            "Use this as the role-transfer table entry: inverse hierarchy is rejected unchanged and replaced by monotone compression/locality-bias. "
            "The next proof-grade move is a target-independent beta-source theorem or, if beta=1 is retained, a blind frozen-kernel empirical test of the compression signature."
        ),
        "full_kernel_now": False,
        "toe_closure_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    theorem = payload["compression_successor_theorem"]
    decision = payload["role_successor_decision"]
    lines = [
        "# P2644/S1594 modified compressed inverse-hierarchy successor theorem",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first anti-duplication audit",
        "",
        "This audit greps modified/compressed successor, inverse-hierarchy rejection, attention suppression, beta/source, and L_total/ToE guard content before adding the theorem.",
        "",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines.extend([
        "",
        "## Compression theorem",
        "",
        f"Suppression factor: `{theorem['suppression_factor']}`.",
        f"Derivative numerator: `{theorem['analytic_derivative_numerator']}`.",
        theorem["monotonicity_proof_on_d_ge_1"],
        f"S(1)=`{theorem['d1_suppression']:.12f}`, S(7)=`{theorem['d7_suppression']:.12f}`, S(7)/S(1)=`{theorem['d7_over_d1_suppression_ratio']:.12f}`.",
        "",
        "## Grid witness",
        "",
        "| d | strict/legacy attention | derivative numerator |",
        "| ---: | ---: | ---: |",
    ])
    for row in theorem["grid_rows"]:
        lines.append(f"| {row['d']} | `{row['strict_over_legacy_attention_suppression']:.12f}` | `{row['derivative_numerator']:.12f}` |")
    lines.extend([
        "",
        "## Verdict",
        "",
        decision["professorial_verdict"],
        "",
        f"Role-transfer verdict: `{decision['legacy_role_transfer_verdict']}`.",
        f"Full kernel now? `{decision['full_kernel_now']}`.",
        f"ToE closure now? `{decision['toe_closure_now']}`.",
        "",
        "## Next honest step",
        "",
        decision["next_honest_step"],
        "",
        "## Negative exports",
        "",
    ])
    for key, value in payload["negative_export_flags"].items():
        lines.append(f"- `{key}`: `{value}`")
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    theorem = compression_theorem()
    decision = role_successor_decision(theorem)
    payload: dict[str, Any] = {
        "status": "P2644_MODIFIED_COMPRESSED_INVERSE_HIERARCHY_SUCCESSOR_THEOREM_NO_FULL_TRANSFER_NO_FALSE_PASS",
        "latest_commit_audit": latest_commit_audit(),
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_hashes": {name: sha256_file(path) for name, path in SOURCE_FILES.items()},
        "upstream_consistency": upstream_consistency(),
        "compression_successor_theorem": theorem,
        "role_successor_decision": decision,
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({k: v for k, v in payload.items() if k != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        SOURCE_FILES["STRICT_EQUATION_SHEET"],
        "P2644/S1594 modified compressed inverse-hierarchy successor guard",
        "## P2644/S1594 modified compressed inverse-hierarchy successor guard\n\n"
        "`P2644/S1594` supplies the modified successor theorem left open by P2643: relative to the legacy hyperbolic denominator, the strict `beta=1, eta=9/5` denominator has suppression factor `S(d)=(1+0.01d)/(1+d^(9/5))`, and `S` is strictly decreasing for `d>=1`.  Thus the unchanged distant-octave/inverse-hierarchy role is rejected, while a modified strict compression/locality-bias successor is admissible as descriptive semantics only.  This still does not export beta sourcehood, bridge completion, role-bearing `L_total`, `QW-2191` discharge, or ToE closure.\n",
    )
    append_once(
        SOURCE_FILES["STRICT_LAGRANGIAN_DRAFT"],
        "P2644/S1594 modified compressed inverse-hierarchy Ltotal guard",
        "## P2644/S1594 modified compressed inverse-hierarchy Ltotal guard\n\n"
        "`P2644/S1594` keeps `L_total` closed: the strict denominator can be read as a monotone compression/locality-bias successor, not as unchanged legacy inverse-hierarchy transfer.  A role-bearing variational term still requires a target-independent beta source and typed source semantics for the compression operator.\n",
    )
    print(rel(OUT))
    print(rel(MD))


if __name__ == "__main__":
    main()
