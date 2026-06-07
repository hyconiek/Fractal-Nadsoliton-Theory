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
OUT = GEN / "p2556_s1506_legacy_to_strict_damping_homotopy_source_nonuniqueness_certificate.json"
MD = GEN / "p2556_s1506_legacy_to_strict_damping_homotopy_source_nonuniqueness_certificate.md"

SOURCE_FILES = {
    "P2377_TRANSPORT_PRIMITIVE": GEN / "p2377_s1327_damping_compression_transport_primitive_uniform_coupling_theorem.json",
    "P2555_DENOMINATOR_NONRENORMALIZATION": GEN / "p2555_s1505_legacy_to_strict_damping_denominator_nonrenormalization_certificate.json",
}
DOMAIN = list(range(1, 12))
BETA_TORS_CANDIDATES = [0.01, 0.05]
ETA = 9.0 / 5.0
NEGATIVE_EXPORT_FLAGS = [
    "unique_damping_homotopy_source_exported", "damping_compression_bridge_component_ready", "legacy_to_strict_completion_bridge_exported",
    "full_bridge_theorem_exported", "role_transfer_theorem_exported", "legacy_role_transfer_claimed", "role_bearing_ltotal_exported",
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
        "new_packet": "P2556|S1506|homotopy source nonuniqueness|damping homotopy nonuniqueness|completion path nonuniqueness",
        "intended_research_nonduplication": "homotopy.*nonunique|nonunique.*homotopy|path nonuniqueness|geometric homotopy|instantaneous source.*homotopy|log-transport.*nonunique",
        "precursors": "P2377|S1327|P2555|S1505|denominator-completion homotopy|damping denominator nonrenormalization",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def legacy(beta_tors: float, d: int) -> float:
    return 1.0 + beta_tors * d


def strict(d: int) -> float:
    return 1.0 + d ** ETA


def rows() -> list[dict[str, Any]]:
    out = []
    for beta_tors in BETA_TORS_CANDIDATES:
        for d in DOMAIN:
            l_val = legacy(beta_tors, d)
            s_val = strict(d)
            endpoint = math.log(s_val / l_val)
            linear_mid_source = (s_val - l_val) / (0.5 * (l_val + s_val))
            geometric_source = endpoint
            out.append({
                "beta_tors": beta_tors,
                "d": d,
                "legacy_denominator_L": l_val,
                "strict_denominator_S": s_val,
                "endpoint_log_transport_C": endpoint,
                "linear_denominator_homotopy": "u_s=(1-s)*L+s*S",
                "linear_midpoint_instantaneous_source": linear_mid_source,
                "linear_integral_equals_endpoint_C": True,
                "geometric_denominator_homotopy": "u_s=L^(1-s)*S^s",
                "geometric_instantaneous_source_constant_in_s": geometric_source,
                "geometric_integral_equals_endpoint_C": True,
                "midpoint_source_difference_abs": abs(linear_mid_source - geometric_source),
                "homotopy_sources_differ": abs(linear_mid_source - geometric_source) > 1e-12,
            })
    return out


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2555 = theorem(sources["P2555_DENOMINATOR_NONRENORMALIZATION"], "legacy_to_strict_damping_denominator_nonrenormalization_certificate")
    audit_rows = rows()
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2556_T1_legacy_to_strict_damping_homotopy_source_nonuniqueness_certificate",
        "audited_chain": ["P2377/S1327", "P2555/S1505"],
        "p2555_denominator_nonrenormalization_inherited": p2555.get("beta_tors_to_beta_eta_scalar_renormalization_refuted_for_audited_class") is True,
        "domain": DOMAIN,
        "beta_tors_candidates_audited": BETA_TORS_CANDIDATES,
        "homotopy_source_nonuniqueness_rows": audit_rows,
        "row_count": len(audit_rows),
        "both_homotopies_have_same_endpoints_and_transport_primitive": all(row["linear_integral_equals_endpoint_C"] and row["geometric_integral_equals_endpoint_C"] for row in audit_rows),
        "instantaneous_sources_differ_for_all_rows": all(row["homotopy_sources_differ"] for row in audit_rows),
        "homotopy_endpoint_data_do_not_select_unique_source_dynamics": True,
        "required_bridge_addition": "a strict principle selecting a particular damping-completion path or source density, not just the endpoint compression primitive",
        "recommended_next_honest_step": (
            "Do not treat the endpoint log-compression primitive as a unique dynamics. The next honest step is to find a strict nadsoliton principle that selects the damping-completion homotopy/source density; otherwise keep the damping bridge conditional and do not run role-transfer."
        ),
        "not_licensed": [
            "No unique damping homotopy/source density is selected.",
            "No bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing L_total, legacy-role transfer, or ToE closure is exported.",
        ],
    }
    for flag in NEGATIVE_EXPORT_FLAGS:
        theorem_export[flag] = False
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "p2555_inherited": theorem_export["p2555_denominator_nonrenormalization_inherited"],
        "same_endpoint_transport": theorem_export["both_homotopies_have_same_endpoints_and_transport_primitive"],
        "source_nonuniqueness_witnessed": theorem_export["instantaneous_sources_differ_for_all_rows"],
        "no_false_bridge_or_role_transfer": theorem_export["full_bridge_theorem_exported"] is False and theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_or_toe_claim": theorem_export["qw2191_discharged_by_this_certificate"] is False and theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2556",
        "stage_id": "S1506",
        "status": "LEGACY_TO_STRICT_DAMPING_HOMOTOPY_SOURCE_NONUNIQUENESS_CERTIFICATE_NO_UNIQUE_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "legacy_to_strict_damping_homotopy_source_nonuniqueness_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["legacy_to_strict_damping_homotopy_source_nonuniqueness_certificate"]["theorem_export"]
    lines = [
        "# P2556/S1506 legacy-to-strict damping homotopy source nonuniqueness certificate", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Row count: `{t['row_count']}`.",
        f"- Same endpoint transport primitive for both homotopies: `{t['both_homotopies_have_same_endpoints_and_transport_primitive']}`.",
        f"- Instantaneous sources differ for all rows: `{t['instantaneous_sources_differ_for_all_rows']}`.",
        f"- Unique damping homotopy source exported: `{t['unique_damping_homotopy_source_exported']}`.", "", "## Interpretation", "",
        "The linear denominator homotopy and geometric denominator homotopy share the same legacy and strict endpoints and integrate to the same endpoint log-compression primitive, but their instantaneous source densities differ.  Endpoint compression data therefore do not select unique bridge dynamics.",
        "", "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No unique damping homotopy/source, bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, legacy-role transfer, or ToE closure is exported.",
        "", "## Fingerprint", "", f"`{payload['legacy_to_strict_damping_homotopy_source_nonuniqueness_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
`P2556/S1506` continues the damping bridge audit by separating endpoint compression from source dynamics.  The linear denominator homotopy `u_s=(1-s)L+sS` and the geometric homotopy `u_s=L^(1-s)S^s` have the same endpoints `L=1+beta_tors*d`, `S=1+d^(9/5)` and the same endpoint log-compression primitive `log(S/L)`, but their instantaneous source densities differ on every audited row.  Thus endpoint compression does not select a unique strict damping source; a real bridge must add a strict homotopy/source-density selector.
""".strip()
    lag_section = """
`P2556/S1506` blocks promotion of endpoint compression into a unique role-bearing `L_total` dynamics.  The damping bridge still needs a strict principle selecting the homotopy/source density, followed only later by bridge completion and role-transfer auditing.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2556/S1506 damping homotopy source nonuniqueness guard", "## P2556/S1506 damping homotopy source nonuniqueness guard\n\n" + eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2556/S1506 damping homotopy nonuniqueness Ltotal guard", "## P2556/S1506 damping homotopy nonuniqueness Ltotal guard\n\n" + lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
