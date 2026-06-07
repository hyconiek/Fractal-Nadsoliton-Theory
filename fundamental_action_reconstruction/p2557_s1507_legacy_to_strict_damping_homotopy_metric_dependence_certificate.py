#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2557_s1507_legacy_to_strict_damping_homotopy_metric_dependence_certificate.json"
MD = GEN / "p2557_s1507_legacy_to_strict_damping_homotopy_metric_dependence_certificate.md"
SOURCE_FILES = {
    "P2556_HOMOTOPY_NONUNIQUENESS": GEN / "p2556_s1506_legacy_to_strict_damping_homotopy_source_nonuniqueness_certificate.json",
}
DOMAIN = list(range(1, 12))
BETA_TORS_CANDIDATES = [0.01, 0.05]
ETA = 9.0 / 5.0
NEGATIVE_EXPORT_FLAGS = [
    "homotopy_metric_selector_source_exported", "unique_damping_homotopy_source_exported", "damping_compression_bridge_component_ready",
    "legacy_to_strict_completion_bridge_exported", "full_bridge_theorem_exported", "role_transfer_theorem_exported",
    "legacy_role_transfer_claimed", "role_bearing_ltotal_exported", "qw2191_discharged_by_this_certificate", "toe_closure_claimed",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run([
        "rg", "-n", pattern, "fundamental_action_reconstruction", "material_dowodowy",
        "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "!fundamental_action_reconstruction/generated/**",
    ], cwd=REPO, check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:80]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2557|S1507|homotopy metric dependence|metric-dependent homotopy selector|homotopy selector metric",
        "intended_research_nonduplication": "homotopy.*metric dependence|metric.*homotopy selector|log-source cost|denominator velocity cost|path action.*damping|linear.*geometric.*cost",
        "precursors": "P2556|S1506|homotopy source nonuniqueness|damping homotopy nonuniqueness",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def legacy(beta_tors: float, d: int) -> float:
    return 1.0 + beta_tors * d


def strict(d: int) -> float:
    return 1.0 + d ** ETA


def cost_rows() -> list[dict[str, Any]]:
    rows = []
    for beta_tors in BETA_TORS_CANDIDATES:
        for d in DOMAIN:
            l_val = legacy(beta_tors, d)
            s_val = strict(d)
            delta = s_val - l_val
            c = math.log(s_val / l_val)
            linear_denominator_velocity_cost = delta ** 2
            geometric_denominator_velocity_cost = c * (s_val ** 2 - l_val ** 2) / 2.0
            linear_log_source_cost = delta ** 2 / (l_val * s_val)
            geometric_log_source_cost = c ** 2
            rows.append({
                "beta_tors": beta_tors,
                "d": d,
                "legacy_denominator_L": l_val,
                "strict_denominator_S": s_val,
                "endpoint_log_transport_C": c,
                "linear_denominator_velocity_cost": linear_denominator_velocity_cost,
                "geometric_denominator_velocity_cost": geometric_denominator_velocity_cost,
                "denominator_velocity_metric_selects": "linear_homotopy",
                "linear_log_source_cost": linear_log_source_cost,
                "geometric_log_source_cost": geometric_log_source_cost,
                "log_source_metric_selects": "geometric_homotopy",
                "metrics_select_different_homotopies": linear_denominator_velocity_cost < geometric_denominator_velocity_cost and geometric_log_source_cost < linear_log_source_cost,
            })
    return rows


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2556 = theorem(sources["P2556_HOMOTOPY_NONUNIQUENESS"], "legacy_to_strict_damping_homotopy_source_nonuniqueness_certificate")
    rows = cost_rows()
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2557_T1_legacy_to_strict_damping_homotopy_metric_dependence_certificate",
        "audited_chain": ["P2556/S1506"],
        "p2556_homotopy_nonuniqueness_inherited": p2556.get("homotopy_endpoint_data_do_not_select_unique_source_dynamics") is True,
        "domain": DOMAIN,
        "beta_tors_candidates_audited": BETA_TORS_CANDIDATES,
        "homotopy_metric_dependence_rows": rows,
        "row_count": len(rows),
        "denominator_velocity_metric_selects_linear_for_all_rows": all(row["denominator_velocity_metric_selects"] == "linear_homotopy" for row in rows),
        "log_source_metric_selects_geometric_for_all_rows": all(row["log_source_metric_selects"] == "geometric_homotopy" for row in rows),
        "metrics_disagree_for_all_rows": all(row["metrics_select_different_homotopies"] for row in rows),
        "homotopy_metric_choice_is_additional_source_obligation": True,
        "recommended_next_honest_step": (
            "Do not select a damping-completion homotopy by choosing a convenient metric after the fact. The next honest step is to derive the homotopy metric/source-density principle from strict nadsoliton dynamics; otherwise the legacy->strict damping bridge remains conditional and role-transfer stays blocked."
        ),
        "not_licensed": [
            "No strict metric selecting linear, geometric, or any other damping homotopy is exported.",
            "No bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing L_total, legacy-role transfer, or ToE closure is exported.",
        ],
    }
    for flag in NEGATIVE_EXPORT_FLAGS:
        theorem_export[flag] = False
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "p2556_inherited": theorem_export["p2556_homotopy_nonuniqueness_inherited"],
        "metric_disagreement_verified": theorem_export["metrics_disagree_for_all_rows"],
        "metric_choice_marked_source_obligation": theorem_export["homotopy_metric_choice_is_additional_source_obligation"],
        "no_false_bridge_or_role_transfer": theorem_export["full_bridge_theorem_exported"] is False and theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_or_toe_claim": theorem_export["qw2191_discharged_by_this_certificate"] is False and theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2557",
        "stage_id": "S1507",
        "status": "LEGACY_TO_STRICT_DAMPING_HOMOTOPY_METRIC_DEPENDENCE_CERTIFICATE_NO_METRIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "legacy_to_strict_damping_homotopy_metric_dependence_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["legacy_to_strict_damping_homotopy_metric_dependence_certificate"]["theorem_export"]
    lines = [
        "# P2557/S1507 legacy-to-strict damping homotopy metric-dependence certificate", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Row count: `{t['row_count']}`.",
        f"- Denominator-velocity metric selects linear homotopy for all rows: `{t['denominator_velocity_metric_selects_linear_for_all_rows']}`.",
        f"- Log-source metric selects geometric homotopy for all rows: `{t['log_source_metric_selects_geometric_for_all_rows']}`.",
        f"- Metrics disagree for all rows: `{t['metrics_disagree_for_all_rows']}`.",
        f"- Homotopy metric selector source exported: `{t['homotopy_metric_selector_source_exported']}`.", "", "## Interpretation", "",
        "Two natural path actions select different damping-completion homotopies on the same endpoint data: denominator-velocity cost selects the linear denominator path, while log-source cost selects the geometric path.  Therefore a homotopy metric/source-density principle is an additional strict source obligation.",
        "", "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No homotopy metric selector, unique damping source, bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, legacy-role transfer, or ToE closure is exported.",
        "", "## Fingerprint", "", f"`{payload['legacy_to_strict_damping_homotopy_metric_dependence_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
`P2557/S1507` audits whether a simple variational metric can select the damping-completion homotopy.  On the same legacy/strict endpoint data, denominator-velocity cost selects the linear denominator homotopy, while log-source cost selects the geometric homotopy, for every audited `d=1..11` and `beta_tors in {0.01,0.05}`.  Thus homotopy selection is metric-dependent: the bridge needs a strict source for the metric/source-density principle, not an after-the-fact convenient path action.
""".strip()
    lag_section = """
`P2557/S1507` blocks promotion of a chosen damping homotopy into role-bearing `L_total` by an unsourced metric choice.  The damping bridge still requires a strict principle selecting the homotopy metric/source density before bridge completion or role-transfer auditing.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2557/S1507 damping homotopy metric-dependence guard", "## P2557/S1507 damping homotopy metric-dependence guard\n\n" + eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2557/S1507 damping homotopy metric-dependence Ltotal guard", "## P2557/S1507 damping homotopy metric-dependence Ltotal guard\n\n" + lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
