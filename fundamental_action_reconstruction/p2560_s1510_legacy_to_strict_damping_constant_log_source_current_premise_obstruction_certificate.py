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
OUT = GEN / "p2560_s1510_legacy_to_strict_damping_constant_log_source_current_premise_obstruction_certificate.json"
MD = GEN / "p2560_s1510_legacy_to_strict_damping_constant_log_source_current_premise_obstruction_certificate.md"
SOURCE_FILES = {
    "P2559_CONSTANT_LOG_SOURCE_CONDITIONAL": GEN / "p2559_s1509_legacy_to_strict_damping_constant_log_source_conditional_selector_certificate.json",
}
DOMAIN = list(range(1, 12))
BETA_TORS_CANDIDATES = [0.01, 0.05]
ETA = 9.0 / 5.0
AUDITED_Q_VALUES = [-2.0, -1.0, 0.0, 1.0, 2.0]
NONCONSTANT_Q_VALUES = [-2.0, -1.0, 1.0, 2.0]
SAMPLES = [0.0, 0.5, 1.0]
TOL = 1.0e-12
NEGATIVE_EXPORT_FLAGS = [
    "constant_log_source_law_exported", "geometric_homotopy_source_exported", "unique_damping_homotopy_source_exported",
    "damping_compression_bridge_component_ready", "legacy_to_strict_completion_bridge_exported", "full_bridge_theorem_exported",
    "role_transfer_theorem_exported", "legacy_role_transfer_claimed", "role_bearing_ltotal_exported",
    "qw2191_discharged_by_this_certificate", "toe_closure_claimed",
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
        "new_packet": "P2560|S1510|constant log-source obstruction|constant-log-source obstruction|constant log source nonentailment",
        "intended_research_nonduplication": "geometric selector obstruction|current premise.*constant.*log|current-premise.*constant.*log|constant-log-source nonentailment",
        "precursors": "P2559|S1509|constant-log-source conditional selector|constant log-source conditional",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def legacy(beta_tors: float, d: int) -> float:
    return 1.0 + beta_tors * d


def strict(d: int) -> float:
    return 1.0 + d ** ETA


def power_mean_denominator(l_val: float, s_val: float, q_value: float, s_param: float) -> float:
    if abs(q_value) <= TOL:
        return math.exp((1.0 - s_param) * math.log(l_val) + s_param * math.log(s_val))
    return ((1.0 - s_param) * (l_val ** q_value) + s_param * (s_val ** q_value)) ** (1.0 / q_value)


def log_source_density(l_val: float, s_val: float, q_value: float, s_param: float) -> float:
    if abs(q_value) <= TOL:
        return math.log(s_val / l_val)
    mixture = (1.0 - s_param) * (l_val ** q_value) + s_param * (s_val ** q_value)
    return (s_val ** q_value - l_val ** q_value) / (q_value * mixture)


def countermodel_rows() -> list[dict[str, Any]]:
    rows = []
    for beta_tors in BETA_TORS_CANDIDATES:
        for d in DOMAIN:
            l_val = legacy(beta_tors, d)
            s_val = strict(d)
            for q_value in NONCONSTANT_Q_VALUES:
                values = [log_source_density(l_val, s_val, q_value, sample) for sample in SAMPLES]
                rows.append({
                    "beta_tors": beta_tors,
                    "d": d,
                    "q": q_value,
                    "legacy_denominator_L": l_val,
                    "strict_denominator_S": s_val,
                    "u_0_equals_L": abs(power_mean_denominator(l_val, s_val, q_value, 0.0) - l_val) <= TOL,
                    "u_1_equals_S": abs(power_mean_denominator(l_val, s_val, q_value, 1.0) - s_val) <= TOL,
                    "u_midpoint_positive": power_mean_denominator(l_val, s_val, q_value, 0.5) > 0.0,
                    "sampled_log_source_values": values,
                    "sampled_log_source_spread": max(values) - min(values),
                    "violates_constant_log_source": max(values) - min(values) > TOL,
                    "passes_current_endpoint_premises": True,
                })
    return rows


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2559 = theorem(sources["P2559_CONSTANT_LOG_SOURCE_CONDITIONAL"], "legacy_to_strict_damping_constant_log_source_conditional_selector_certificate")
    rows = countermodel_rows()
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2560_T1_constant_log_source_current_premise_obstruction_certificate",
        "audited_chain": ["P2559/S1509"],
        "p2559_conditional_selector_inherited": p2559.get("conditional_geometric_selector_exported") is True,
        "domain": DOMAIN,
        "beta_tors_candidates_audited": BETA_TORS_CANDIDATES,
        "audited_q_values": AUDITED_Q_VALUES,
        "countermodel_q_values": NONCONSTANT_Q_VALUES,
        "source_density_samples": SAMPLES,
        "countermodel_count": len(rows),
        "countermodel_rows": rows,
        "all_countermodels_pass_endpoint_premises": all(row["passes_current_endpoint_premises"] and row["u_0_equals_L"] and row["u_1_equals_S"] for row in rows),
        "all_countermodels_have_positive_midpoint_denominator": all(row["u_midpoint_positive"] for row in rows),
        "all_countermodels_violate_constant_log_source": all(row["violates_constant_log_source"] for row in rows),
        "current_endpoint_premises_do_not_entail_constant_log_source": True,
        "p2559_conditional_selector_cannot_be_promoted_to_source_by_current_premises": True,
        "recommended_next_honest_step": (
            "Do not rerun endpoint or power-mean scans. The next honest step is to derive a strict dynamic source for constant log-denominator source density, or else abandon q=0/geometric selection as a bridge theorem and move to another bridge component such as phase/topological-bit passage."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    theorem_export["constant_log_source_current_premise_obstruction_exported"] = True
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2559_inherited": theorem_export["p2559_conditional_selector_inherited"],
        "countermodels_present": theorem_export["countermodel_count"] == len(DOMAIN) * len(BETA_TORS_CANDIDATES) * len(NONCONSTANT_Q_VALUES),
        "no_constant_log_source_exported": theorem_export["constant_log_source_law_exported"] is False,
        "no_bridge_completion_exported": theorem_export["full_bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2560",
        "stage_id": "S1510",
        "status": "P2560_CURRENT_PREMISE_CONSTANT_LOG_SOURCE_OBSTRUCTION_NO_FALSE_PASS",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "legacy_to_strict_damping_constant_log_source_current_premise_obstruction_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["legacy_to_strict_damping_constant_log_source_current_premise_obstruction_certificate"]["theorem_export"]
    lines = [
        "# P2560/S1510 legacy-to-strict damping constant-log-source current-premise obstruction certificate", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Countermodel count: `{t['countermodel_count']}`.",
        f"- Countermodel q-values: `{t['countermodel_q_values']}`.",
        f"- All countermodels pass endpoint premises: `{t['all_countermodels_pass_endpoint_premises']}`.",
        f"- All countermodels violate constant log-source: `{t['all_countermodels_violate_constant_log_source']}`.",
        f"- Current endpoint premises entail constant log-source: `{not t['current_endpoint_premises_do_not_entail_constant_log_source']}`.", "", "## Interpretation", "",
        "P2559 gives only a conditional selector.  P2560 shows that the current endpoint/power-mean premises do not supply its missing premise: every audited nonzero q path is an explicit countermodel to constant log-source while preserving the legacy/strict endpoints.",
        "", "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No constant-log-source law, geometric homotopy source, unique damping source, bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, legacy-role transfer, or ToE closure is exported.",
        "", "## Fingerprint", "", f"`{payload['legacy_to_strict_damping_constant_log_source_current_premise_obstruction_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
`P2560/S1510` audits whether the P2559 constant-log-source premise follows from the current endpoint/power-mean bridge premises.  It does not: the audited nonzero power-mean parameters `q in {-2,-1,1,2}` give `88` countermodels that preserve the legacy/strict endpoints and positive denominators while violating constant log-source density.  Thus the P2559 `q=0` selector remains conditional and cannot be promoted into a damping bridge source by current premises.
""".strip()
    lag_section = """
`P2560/S1510` blocks role-bearing `L_total` promotion of the geometric damping path from current endpoint premises alone.  The missing constant-log-source law is still a strict source obligation, not a consequence of the audited bridge bookkeeping.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2560/S1510 constant-log-source current-premise obstruction guard", "## P2560/S1510 constant-log-source current-premise obstruction guard\n\n" + eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2560/S1510 constant-log-source current-premise obstruction Ltotal guard", "## P2560/S1510 constant-log-source current-premise obstruction Ltotal guard\n\n" + lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
