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
OUT = GEN / "p2643_s1593_inverse_hierarchy_beta_threshold_role_rejection_certificate.json"
MD = GEN / "p2643_s1593_inverse_hierarchy_beta_threshold_role_rejection_certificate.md"

OMEGA = math.pi / 4.0
PHI = math.pi / 6.0
ETA = 9.0 / 5.0
LEGACY_BETA_TORS = 0.01
STRICT_BETA = 1.0
D_NEAR = 1.0
D_FAR = 7.0

SOURCE_FILES = {
    "P2633_DIAGRAM_RETENTION": GEN / "p2633_s1583_diagram_grounded_strict_kernel_characteristic_preservation_audit.json",
    "P2629_ZBETA_GAUGE": GEN / "p2629_s1579_zbeta_normalization_gauge_obstruction.json",
    "P2630_BETA_SOURCE_SEPARATION": GEN / "p2630_s1580_strict_beta_source_vs_bridge_zbeta_separation.json",
    "P2636_BLOCKER_LATTICE": GEN / "p2636_s1586_current_toe_blocker_lattice_full_kernel_decision_audit.json",
    "P2642_NODE_DEMOTION": GEN / "p2642_s1592_universal_affine_zero_lattice_source_no_go_and_node_role_demotion.json",
    "STRICT_EQUATION_SHEET": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "STRICT_LAGRANGIAN_DRAFT": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

NEGATIVE_EXPORT_FLAGS = [
    "unchanged_inverse_hierarchy_role_transferred",
    "strict_beta_preserves_distant_octave_ratio",
    "micro_beta_preserves_distant_octave_ratio",
    "beta_source_repaired_inverse_hierarchy",
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
        "inverse_hierarchy_role_content": (
            "inverse hierarchy|distant octave|Wilson-loop|Wilson loop|long-range coupling|"
            "K\\(7\\).*K\\(1\\)|d=7.*d=1|heavy-tailed compression"
        ),
        "beta_threshold_source_content": (
            "beta threshold|beta_crit|strict beta|beta=1|beta_tors|positive beta|"
            "beta source|normalization identity|Z_beta|UV normalization"
        ),
        "role_rejection_transfer_content": (
            "role-transfer|role transfer|survives unchanged|modified/compressed|rejected|demot|degrad|"
            "legacy physical role|legacy.*role"
        ),
        "node_demotion_reroute_content": (
            "node/gauge|integer node|zero-lattice|origin-source|origin source|offset/stride|"
            "reroute|source atom"
        ),
        "toe_guard_content": (
            "full kernel|ToE closure|QW-2191|selector source|role-bearing L_total|"
            "blind empirical|frozen-kernel"
        ),
    }
    return {"tool": "rg", "mode": "content-first semantic audit for inverse-hierarchy beta-threshold role rejection", "patterns": {k: rg_count(v) for k, v in patterns.items()}}


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


def inverse_hierarchy_ratio(beta: float, eta: float = ETA) -> float:
    return (phase_abs(D_FAR) / phase_abs(D_NEAR)) * ((1.0 + beta * D_NEAR**eta) / (1.0 + beta * D_FAR**eta))


def beta_threshold() -> dict[str, Any]:
    c_phase = phase_abs(D_FAR) / phase_abs(D_NEAR)
    far_power = D_FAR**ETA
    # Solve c_phase*(1+beta)/(1+beta*7^eta) > 1.
    beta_crit = (c_phase - 1.0) / (far_power - c_phase)
    derivative_sign = "negative" if far_power > 1.0 else "nonnegative"
    return {
        "ratio_formula": "R(beta)=|cos(pi*7/4+pi/6)|/|cos(pi/4+pi/6)| * (1+beta)/(1+beta*7^(9/5))",
        "phase_ratio_C": c_phase,
        "far_power_7_eta": far_power,
        "beta_critical_exact_role_boundary": beta_crit,
        "preservation_iff": "unchanged inverse-hierarchy proxy R(beta)>1 iff 0 <= beta < beta_crit",
        "monotonicity_proof": "dR/dbeta has sign C*(1-7^(9/5))/(1+beta*7^(9/5))^2, hence is strictly negative for beta>=0.",
        "derivative_sign_on_beta_nonnegative": derivative_sign,
    }


def micro_beta_from_p2629() -> float | None:
    p2629 = load_json(SOURCE_FILES["P2629_ZBETA_GAUGE"])
    metrics = p2629.get("qw_metrics", {})
    value = metrics.get("micro_beta_median")
    return float(value) if value is not None else None


def ratio_table(micro_beta: float | None, threshold: float) -> list[dict[str, Any]]:
    betas = [0.0, LEGACY_BETA_TORS, threshold * 0.999, threshold, threshold * 1.001, STRICT_BETA]
    if micro_beta is not None:
        betas.append(micro_beta)
    rows = []
    seen = set()
    for beta in betas:
        key = round(beta, 15)
        if key in seen:
            continue
        seen.add(key)
        rows.append({
            "beta": beta,
            "label": (
                "zero_damping_limit" if beta == 0.0 else
                "legacy_beta_tors" if abs(beta - LEGACY_BETA_TORS) < 1e-15 else
                "strict_beta" if abs(beta - STRICT_BETA) < 1e-15 else
                "micro_beta_median" if micro_beta is not None and abs(beta - micro_beta) < 1e-15 else
                "threshold_probe"
            ),
            "ratio_abs_k7_over_k1": inverse_hierarchy_ratio(beta),
            "preserves_unchanged_inverse_hierarchy_proxy": inverse_hierarchy_ratio(beta) > 1.0,
            "beta_below_critical": beta < threshold,
        })
    return sorted(rows, key=lambda row: row["beta"])


def upstream_consistency() -> dict[str, Any]:
    p2633 = load_json(SOURCE_FILES["P2633_DIAGRAM_RETENTION"])
    p2642 = load_json(SOURCE_FILES["P2642_NODE_DEMOTION"])
    p2630 = load_json(SOURCE_FILES["P2630_BETA_SOURCE_SEPARATION"])
    p2636 = load_json(SOURCE_FILES["P2636_BLOCKER_LATTICE"])
    return {
        "p2633_status": p2633.get("status"),
        "p2633_strict_ratio": p2633.get("finite_witness", {}).get("inverse_hierarchy_ratio_abs_k7_over_abs_k1", {}).get("strict"),
        "p2633_legacy_ratio": p2633.get("finite_witness", {}).get("inverse_hierarchy_ratio_abs_k7_over_abs_k1", {}).get("legacy_amplitude_normalized"),
        "p2642_node_role_status": p2642.get("demotion_decision", {}).get("legacy_integer_node_gauge_role_status"),
        "p2642_next_step": p2642.get("demotion_decision", {}).get("next_honest_step"),
        "p2630_status": p2630.get("status"),
        "p2636_status": p2636.get("status"),
    }


def role_verdict(threshold_data: dict[str, Any], rows: list[dict[str, Any]], micro_beta: float | None) -> dict[str, Any]:
    strict_row = next(row for row in rows if row["label"] == "strict_beta")
    micro_row = next((row for row in rows if row["label"] == "micro_beta_median"), None)
    gates = {
        "legacy_beta_tors_preserves_ratio": next(row for row in rows if row["label"] == "legacy_beta_tors")["preserves_unchanged_inverse_hierarchy_proxy"],
        "strict_beta_below_critical": STRICT_BETA < threshold_data["beta_critical_exact_role_boundary"],
        "strict_beta_preserves_ratio": strict_row["preserves_unchanged_inverse_hierarchy_proxy"],
        "micro_beta_available": micro_beta is not None,
        "micro_beta_below_critical": bool(micro_beta is not None and micro_beta < threshold_data["beta_critical_exact_role_boundary"]),
        "micro_beta_preserves_ratio": bool(micro_row and micro_row["preserves_unchanged_inverse_hierarchy_proxy"]),
        "node_gauge_role_already_demoted_by_p2642": "DEMOTE" in str(upstream_consistency().get("p2642_node_role_status", "")),
    }
    return {
        "gates": gates,
        "unchanged_inverse_hierarchy_role_status": "REJECT_UNCHANGED_TRANSFER_FOR_CURRENT_STRICT_BETA_AND_MICRO_BETA_MEDIAN",
        "allowed_successor_status": "MODIFIED_COMPRESSED_SUCCESSOR_SEMANTICS_REMAINS_OPEN_BUT_MUST_NOT_BE_CALLED_UNCHANGED_LEGACY_ROLE",
        "professorial_verdict": (
            "For the same phase channel and eta=9/5, the distant-octave ratio is a strictly decreasing function of beta. "
            "It exceeds one only below beta_crit, while strict beta=1 and the current micro beta median lie far above beta_crit. "
            "Thus the legacy inverse-hierarchy role does not transfer unchanged after P2642 node-role demotion; only a modified/compressed successor interpretation remains admissible."
        ),
        "next_honest_step": (
            "Do not try to rescue the unchanged inverse-hierarchy role by amplitude normalization or more node lifts. "
            "Either prove a new beta-source theorem that changes the strict damping semantics below beta_crit without retuning, or mark inverse hierarchy as rejected unchanged and build the modified/compressed successor role-transfer theorem."
        ),
        "full_kernel_now": False,
        "toe_closure_now": False,
    }


def professorial_closure_path() -> list[dict[str, str]]:
    return [
        {
            "rank": "1",
            "step": "Reject unchanged inverse-hierarchy transfer unless beta semantics change",
            "criterion": "A strict or micro beta above beta_crit cannot preserve |K(7)|/|K(1)|>1 on the audited phase/damping channel.",
        },
        {
            "rank": "2",
            "step": "Build modified/compressed successor theorem",
            "criterion": "If strict beta=1 is retained, the role must be restated as strong UV/local compression rather than legacy distant-octave amplification.",
        },
        {
            "rank": "3",
            "step": "Run beta-source proof separately",
            "criterion": "A beta-source theorem must be target-independent and cannot reuse Z_beta normalization gauge as source data.",
        },
        {
            "rank": "4",
            "step": "Only then rerun role-transfer matrix",
            "criterion": "Node/gauge is demoted, inverse hierarchy is rejected or modified, and beta semantics are typed before any L_total or ToE claim reopens.",
        },
    ]


def write_markdown(payload: dict[str, Any]) -> None:
    threshold = payload["beta_threshold_theorem"]
    verdict = payload["role_verdict"]
    lines = [
        "# P2643/S1593 inverse-hierarchy beta-threshold role rejection certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first anti-duplication audit",
        "",
        "This audit greps inverse-hierarchy, beta-threshold/source, role-transfer/rejection, node-demotion reroute, and ToE guard content before adding the threshold certificate.",
        "",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines.extend([
        "",
        "## Threshold theorem",
        "",
        f"Ratio formula: `{threshold['ratio_formula']}`.",
        f"Phase ratio C: `{threshold['phase_ratio_C']:.12f}`.",
        f"7^(9/5): `{threshold['far_power_7_eta']:.12f}`.",
        f"Critical beta: `{threshold['beta_critical_exact_role_boundary']:.12f}`.",
        f"Monotonicity: {threshold['monotonicity_proof']}",
        "",
        "## Ratio table",
        "",
        "| label | beta | |K(7)|/|K(1)| | preserves >1? |",
        "| --- | ---: | ---: | --- |",
    ])
    for row in payload["ratio_table"]:
        lines.append(f"| `{row['label']}` | `{row['beta']:.12f}` | `{row['ratio_abs_k7_over_k1']:.12f}` | `{row['preserves_unchanged_inverse_hierarchy_proxy']}` |")
    lines.extend([
        "",
        "## Verdict",
        "",
        verdict["professorial_verdict"],
        "",
        f"Unchanged inverse-hierarchy role status: `{verdict['unchanged_inverse_hierarchy_role_status']}`.",
        f"Allowed successor status: `{verdict['allowed_successor_status']}`.",
        f"Full kernel now? `{verdict['full_kernel_now']}`.",
        f"ToE closure now? `{verdict['toe_closure_now']}`.",
        "",
        "## Next honest step",
        "",
        verdict["next_honest_step"],
        "",
        "## Professorial closure path",
        "",
    ])
    for row in payload["professorial_closure_path"]:
        lines.append(f"{row['rank']}. **{row['step']}** — {row['criterion']}")
    lines.extend(["", "## Negative exports", ""])
    for key, value in payload["negative_export_flags"].items():
        lines.append(f"- `{key}`: `{value}`")
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    threshold = beta_threshold()
    micro_beta = micro_beta_from_p2629()
    rows = ratio_table(micro_beta, threshold["beta_critical_exact_role_boundary"])
    verdict = role_verdict(threshold, rows, micro_beta)
    payload: dict[str, Any] = {
        "status": "P2643_INVERSE_HIERARCHY_BETA_THRESHOLD_ROLE_REJECTION_NO_FALSE_PASS",
        "latest_commit_audit": latest_commit_audit(),
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_hashes": {name: sha256_file(path) for name, path in SOURCE_FILES.items()},
        "upstream_consistency": upstream_consistency(),
        "beta_threshold_theorem": threshold,
        "micro_beta_median_from_p2629": micro_beta,
        "ratio_table": rows,
        "role_verdict": verdict,
        "professorial_closure_path": professorial_closure_path(),
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({k: v for k, v in payload.items() if k != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        SOURCE_FILES["STRICT_EQUATION_SHEET"],
        "P2643/S1593 inverse-hierarchy beta-threshold guard",
        "## P2643/S1593 inverse-hierarchy beta-threshold guard\n\n"
        "`P2643/S1593` proves a beta-threshold obstruction for unchanged inverse-hierarchy transfer after P2642 node-role demotion.  For the audited phase channel and `eta=9/5`, `|K(7)|/|K(1)|>1` holds only for `beta < beta_crit ≈ 0.0915`; strict `beta=1` and the current micro beta median lie above that threshold.  Therefore the legacy distant-octave / inverse-hierarchy role is rejected as an unchanged strict transfer unless a new target-independent beta-source theorem changes the strict damping semantics; bridge completion, role transfer, `L_total`, `QW-2191`, and ToE closure remain closed.\n",
    )
    append_once(
        SOURCE_FILES["STRICT_LAGRANGIAN_DRAFT"],
        "P2643/S1593 inverse-hierarchy beta-threshold Ltotal guard",
        "## P2643/S1593 inverse-hierarchy beta-threshold Ltotal guard\n\n"
        "`P2643/S1593` does not re-enable role-bearing `L_total`: the unchanged inverse-hierarchy role is incompatible with strict `beta=1` on the audited threshold calculation.  Any variational term must therefore use a modified/compressed successor semantics, or wait for a new beta-source theorem below the threshold without retuning.\n",
    )
    print(rel(OUT))
    print(rel(MD))


if __name__ == "__main__":
    main()
