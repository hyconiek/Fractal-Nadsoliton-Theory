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
OUT = GEN / "p2559_s1509_legacy_to_strict_damping_constant_log_source_conditional_selector_certificate.json"
MD = GEN / "p2559_s1509_legacy_to_strict_damping_constant_log_source_conditional_selector_certificate.md"
SOURCE_FILES = {
    "P2558_POWER_MEAN_CONTINUUM": GEN / "p2558_s1508_legacy_to_strict_damping_power_mean_homotopy_continuum_certificate.json",
}
DOMAIN = list(range(1, 12))
BETA_TORS_CANDIDATES = [0.01, 0.05]
ETA = 9.0 / 5.0
POWER_MEAN_Q_VALUES = [-2.0, -1.0, 0.0, 1.0, 2.0]
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
        "new_packet": "P2559|S1509|constant log-source|constant log source|log-source constant|constant source density",
        "intended_research_nonduplication": "source-density selector|geometric homotopy selector|constant.*log.*homotopy|power-mean.*selector",
        "precursors": "P2558|S1508|power mean homotopy|power-mean homotopy|homotopy continuum",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def legacy(beta_tors: float, d: int) -> float:
    return 1.0 + beta_tors * d


def strict(d: int) -> float:
    return 1.0 + d ** ETA


def log_source_density(l_val: float, s_val: float, q_value: float, s_param: float) -> float:
    if abs(q_value) <= TOL:
        return math.log(s_val / l_val)
    mixture = (1.0 - s_param) * (l_val ** q_value) + s_param * (s_val ** q_value)
    return (s_val ** q_value - l_val ** q_value) / (q_value * mixture)


def q_source_row(l_val: float, s_val: float, q_value: float) -> dict[str, Any]:
    samples = [{"s": sample, "rho_q_s": log_source_density(l_val, s_val, q_value, sample)} for sample in SAMPLES]
    values = [sample["rho_q_s"] for sample in samples]
    spread = max(values) - min(values)
    return {
        "q": q_value,
        "sampled_log_source_density": samples,
        "sampled_log_source_spread": spread,
        "sampled_constant_log_source": spread <= TOL,
        "analytic_constant_log_source_condition": "q=0 in the audited power-mean family, since rho_q(s)=(S^q-L^q)/(q*((1-s)L^q+sS^q)) varies in s for q!=0 and S!=L; q=0 gives rho_0=log(S/L).",
    }


def selector_rows() -> list[dict[str, Any]]:
    rows = []
    for beta_tors in BETA_TORS_CANDIDATES:
        for d in DOMAIN:
            l_val = legacy(beta_tors, d)
            s_val = strict(d)
            q_rows = [q_source_row(l_val, s_val, q_value) for q_value in POWER_MEAN_Q_VALUES]
            constant_q_values = [entry["q"] for entry in q_rows if entry["sampled_constant_log_source"]]
            rows.append({
                "beta_tors": beta_tors,
                "d": d,
                "legacy_denominator_L": l_val,
                "strict_denominator_S": s_val,
                "endpoint_log_transport_C": math.log(s_val / l_val),
                "q_source_rows": q_rows,
                "constant_log_source_q_values": constant_q_values,
                "constant_log_source_selects_geometric_q0": constant_q_values == [0.0],
            })
    return rows


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2558 = theorem(sources["P2558_POWER_MEAN_CONTINUUM"], "legacy_to_strict_damping_power_mean_homotopy_continuum_certificate")
    rows = selector_rows()
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2559_T1_constant_log_source_conditional_geometric_selector_certificate",
        "audited_chain": ["P2558/S1508"],
        "p2558_power_mean_continuum_inherited": p2558.get("finite_q_sample_witnesses_continuum_nonuniqueness") is True,
        "domain": DOMAIN,
        "beta_tors_candidates_audited": BETA_TORS_CANDIDATES,
        "power_mean_q_values_audited": POWER_MEAN_Q_VALUES,
        "source_density_samples": SAMPLES,
        "row_count": len(rows),
        "selector_rows": rows,
        "constant_log_source_selects_geometric_q0_for_all_rows": all(row["constant_log_source_selects_geometric_q0"] for row in rows),
        "conditional_selector_theorem": "Within the audited power-mean denominator homotopy family from legacy L to strict S, adding the premise d/ds log u_s = constant selects the geometric q=0 path.",
        "constant_log_source_premise_is_unsourced_bridge_obligation": True,
        "recommended_next_honest_step": (
            "Do not treat the conditional q=0/geometric selection as a bridge theorem. The next honest step is to derive, from strict nadsoliton dynamics, why the damping compression source density must be constant in log-denominator time; without that source law, the legacy->strict damping bridge and any role-transfer audit remain conditional."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    theorem_export["conditional_geometric_selector_exported"] = True
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2558_inherited": theorem_export["p2558_power_mean_continuum_inherited"],
        "conditional_only_no_constant_source_exported": theorem_export["constant_log_source_law_exported"] is False,
        "no_geometric_homotopy_source_exported": theorem_export["geometric_homotopy_source_exported"] is False,
        "no_bridge_completion_exported": theorem_export["full_bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2559",
        "stage_id": "S1509",
        "status": "P2559_CONDITIONAL_GEOMETRIC_HOMOTOPY_SELECTOR_NO_CONSTANT_LOG_SOURCE_EXPORT_NO_FALSE_PASS",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "legacy_to_strict_damping_constant_log_source_conditional_selector_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["legacy_to_strict_damping_constant_log_source_conditional_selector_certificate"]["theorem_export"]
    lines = [
        "# P2559/S1509 legacy-to-strict damping constant-log-source conditional selector certificate", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Row count: `{t['row_count']}`.",
        f"- Audited power-mean q-values: `{t['power_mean_q_values_audited']}`.",
        f"- Constant log-source selects geometric `q=0` for all rows: `{t['constant_log_source_selects_geometric_q0_for_all_rows']}`.",
        f"- Constant log-source law exported: `{t['constant_log_source_law_exported']}`.",
        f"- Geometric homotopy source exported: `{t['geometric_homotopy_source_exported']}`.", "", "## Interpretation", "",
        "This is a conditional theorem, not a bridge completion: if strict dynamics supplies a constant log-denominator source-density law, then the power-mean continuum collapses to the geometric homotopy `q=0`.  The source law itself remains unsourced.",
        "", "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No constant-log-source law, geometric homotopy source, unique damping source, bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, legacy-role transfer, or ToE closure is exported.",
        "", "## Fingerprint", "", f"`{payload['legacy_to_strict_damping_constant_log_source_conditional_selector_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
`P2559/S1509` turns the P2558 power-mean continuum obstruction into a conditional selector audit.  In the audited family, the additional premise `d/ds log u_s = constant` selects exactly the geometric path `q=0`; every audited `q != 0` row has nonconstant sampled log-source density.  This is only a conditional selector: the constant-log-source law is not exported by current strict nadsoliton dynamics, so the legacy->strict damping bridge remains unclosed.
""".strip()
    lag_section = """
`P2559/S1509` records that a constant log-denominator source law would select the geometric damping homotopy, but the law itself is still an unsourced bridge premise.  Therefore the selected path still cannot be promoted into role-bearing `L_total` or role-transfer evidence.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2559/S1509 constant-log-source conditional selector guard", "## P2559/S1509 constant-log-source conditional selector guard\n\n" + eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2559/S1509 constant-log-source conditional selector Ltotal guard", "## P2559/S1509 constant-log-source conditional selector Ltotal guard\n\n" + lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
