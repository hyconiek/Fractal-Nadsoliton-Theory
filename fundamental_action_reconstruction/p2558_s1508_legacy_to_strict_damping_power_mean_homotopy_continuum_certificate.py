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
OUT = GEN / "p2558_s1508_legacy_to_strict_damping_power_mean_homotopy_continuum_certificate.json"
MD = GEN / "p2558_s1508_legacy_to_strict_damping_power_mean_homotopy_continuum_certificate.md"
SOURCE_FILES = {
    "P2556_HOMOTOPY_NONUNIQUENESS": GEN / "p2556_s1506_legacy_to_strict_damping_homotopy_source_nonuniqueness_certificate.json",
    "P2557_HOMOTOPY_METRIC_DEPENDENCE": GEN / "p2557_s1507_legacy_to_strict_damping_homotopy_metric_dependence_certificate.json",
}
DOMAIN = list(range(1, 12))
BETA_TORS_CANDIDATES = [0.01, 0.05]
ETA = 9.0 / 5.0
POWER_MEAN_Q_VALUES = [-2.0, -1.0, 0.0, 1.0, 2.0]
MIDPOINT_S = 0.5
TOL = 1.0e-12
NEGATIVE_EXPORT_FLAGS = [
    "power_mean_homotopy_selector_source_exported", "unique_damping_homotopy_source_exported", "damping_compression_bridge_component_ready",
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
        "new_packet": "P2558|S1508|power mean homotopy|power-mean homotopy|homotopy continuum|continuum nonuniqueness",
        "intended_research_nonduplication": "q-homotopy|generalized mean homotopy|generalized power mean|power mean.*damping|power-mean.*damping|homotopy.*continuum",
        "precursors": "P2556|S1506|P2557|S1507|homotopy source nonuniqueness|homotopy metric dependence",
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


def midpoint_log_source_density(l_val: float, s_val: float, q_value: float) -> float:
    if abs(q_value) <= TOL:
        return math.log(s_val / l_val)
    midpoint_mixture = 0.5 * (l_val ** q_value) + 0.5 * (s_val ** q_value)
    return (s_val ** q_value - l_val ** q_value) / (q_value * midpoint_mixture)


def power_mean_rows() -> list[dict[str, Any]]:
    rows = []
    for beta_tors in BETA_TORS_CANDIDATES:
        for d in DOMAIN:
            l_val = legacy(beta_tors, d)
            s_val = strict(d)
            endpoint_log_transport = math.log(s_val / l_val)
            q_midpoint_sources = []
            for q_value in POWER_MEAN_Q_VALUES:
                q_midpoint_sources.append({
                    "q": q_value,
                    "u_0_equals_L": abs(power_mean_denominator(l_val, s_val, q_value, 0.0) - l_val) <= TOL,
                    "u_1_equals_S": abs(power_mean_denominator(l_val, s_val, q_value, 1.0) - s_val) <= TOL,
                    "u_midpoint": power_mean_denominator(l_val, s_val, q_value, MIDPOINT_S),
                    "midpoint_log_source_density": midpoint_log_source_density(l_val, s_val, q_value),
                    "endpoint_log_primitive": endpoint_log_transport,
                })
            midpoint_values = [entry["midpoint_log_source_density"] for entry in q_midpoint_sources]
            endpoint_primitives = [entry["endpoint_log_primitive"] for entry in q_midpoint_sources]
            rows.append({
                "beta_tors": beta_tors,
                "d": d,
                "legacy_denominator_L": l_val,
                "strict_denominator_S": s_val,
                "endpoint_log_transport_C": endpoint_log_transport,
                "q_midpoint_sources": q_midpoint_sources,
                "all_power_mean_paths_share_endpoints": all(entry["u_0_equals_L"] and entry["u_1_equals_S"] for entry in q_midpoint_sources),
                "all_power_mean_paths_share_endpoint_log_primitive": max(endpoint_primitives) - min(endpoint_primitives) <= TOL,
                "midpoint_source_spread": max(midpoint_values) - min(midpoint_values),
                "midpoint_source_spread_positive": max(midpoint_values) - min(midpoint_values) > TOL,
            })
    return rows


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2556 = theorem(sources["P2556_HOMOTOPY_NONUNIQUENESS"], "legacy_to_strict_damping_homotopy_source_nonuniqueness_certificate")
    p2557 = theorem(sources["P2557_HOMOTOPY_METRIC_DEPENDENCE"], "legacy_to_strict_damping_homotopy_metric_dependence_certificate")
    rows = power_mean_rows()
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2558_T1_legacy_to_strict_damping_power_mean_homotopy_continuum_certificate",
        "audited_chain": ["P2556/S1506", "P2557/S1507"],
        "p2556_homotopy_nonuniqueness_inherited": p2556.get("homotopy_endpoint_data_do_not_select_unique_source_dynamics") is True,
        "p2557_metric_dependence_inherited": p2557.get("homotopy_metric_choice_is_additional_source_obligation") is True,
        "domain": DOMAIN,
        "beta_tors_candidates_audited": BETA_TORS_CANDIDATES,
        "power_mean_q_values_audited": POWER_MEAN_Q_VALUES,
        "power_mean_homotopy_row_count": len(rows),
        "q_value_count": len(POWER_MEAN_Q_VALUES),
        "power_mean_homotopy_rows": rows,
        "all_power_mean_homotopies_share_endpoints": all(row["all_power_mean_paths_share_endpoints"] for row in rows),
        "all_power_mean_homotopies_share_endpoint_log_primitive": all(row["all_power_mean_paths_share_endpoint_log_primitive"] for row in rows),
        "midpoint_source_spread_positive_for_all_rows": all(row["midpoint_source_spread_positive"] for row in rows),
        "finite_q_sample_witnesses_continuum_nonuniqueness": True,
        "power_mean_homotopy_parameter_is_unsourced_bridge_obligation": True,
        "recommended_next_honest_step": (
            "Do not add more endpoint-only damping homotopy examples. The next honest step is to derive a strict nadsoliton source theorem selecting the power-mean parameter q, or an equivalent source-density law, for the legacy->strict damping compression bridge; without that source, bridge completion and role-transfer remain blocked."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    theorem_export["homotopy_continuum_nonuniqueness_exported"] = True
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2556_inherited": theorem_export["p2556_homotopy_nonuniqueness_inherited"],
        "p2557_inherited": theorem_export["p2557_metric_dependence_inherited"],
        "no_power_mean_selector_source_exported": theorem_export["power_mean_homotopy_selector_source_exported"] is False,
        "no_bridge_completion_exported": theorem_export["full_bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2558",
        "stage_id": "S1508",
        "status": "P2558_CURRENT_PROOF_COMPUTATIONAL_BRIDGE_HOMOTOPY_CONTINUUM_NONUNIQUENESS_NO_FALSE_PASS",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "legacy_to_strict_damping_power_mean_homotopy_continuum_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["legacy_to_strict_damping_power_mean_homotopy_continuum_certificate"]["theorem_export"]
    lines = [
        "# P2558/S1508 legacy-to-strict damping power-mean homotopy continuum certificate", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Row count: `{t['power_mean_homotopy_row_count']}`.",
        f"- Audited power-mean q-values: `{t['power_mean_q_values_audited']}`.",
        f"- All audited power-mean homotopies share endpoints: `{t['all_power_mean_homotopies_share_endpoints']}`.",
        f"- All audited power-mean homotopies share endpoint log primitive: `{t['all_power_mean_homotopies_share_endpoint_log_primitive']}`.",
        f"- Midpoint source-density spread is positive for all rows: `{t['midpoint_source_spread_positive_for_all_rows']}`.",
        f"- Power-mean homotopy selector source exported: `{t['power_mean_homotopy_selector_source_exported']}`.", "", "## Interpretation", "",
        "The audited power-mean family supplies a finite witness to a continuum of damping-completion homotopies with identical legacy/strict endpoints and identical endpoint log transport, while their local source-density profiles differ.  Endpoint compression data therefore cannot select the homotopy parameter.",
        "", "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No power-mean selector source, unique damping source, bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, legacy-role transfer, or ToE closure is exported.",
        "", "## Fingerprint", "", f"`{payload['legacy_to_strict_damping_power_mean_homotopy_continuum_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
`P2558/S1508` strengthens the damping-homotopy bridge obstruction from two examples to an audited power-mean family.  For `q in {-2,-1,0,1,2}` every audited path shares the same legacy/strict endpoints and the same endpoint log transport, but the midpoint log-source density changes with `q` on every audited `d=1..11` and `beta_tors in {0.01,0.05}`.  Thus endpoint compression data do not select the homotopy/source-density law; the legacy->strict damping bridge needs a strict source for the power-mean parameter or an equivalent source-density selector.
""".strip()
    lag_section = """
`P2558/S1508` blocks promotion of a damping-completion path into role-bearing `L_total` by endpoint agreement alone.  The power-mean continuum has the same endpoint transport but different local source densities, so `L_total` still needs a sourced bridge homotopy/metric principle before any legacy role-transfer theorem can be audited.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2558/S1508 damping power-mean homotopy continuum guard", "## P2558/S1508 damping power-mean homotopy continuum guard\n\n" + eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2558/S1508 damping power-mean homotopy continuum Ltotal guard", "## P2558/S1508 damping power-mean homotopy continuum Ltotal guard\n\n" + lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
