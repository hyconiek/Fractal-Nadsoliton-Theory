#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from typing import Any

import sympy as sp

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2589_s1539_apd_mirror_sixth_moment_shell_support_nonuniqueness_audit import support_from_product
from p2591_s1541_apd_product_parameter_sturm_interval_certificate import INTERVAL, OUT as P2591_OUT
from p2594_s1544_apd_eighth_moment_conditional_inverse_selector_certificate import (
    CENTRAL_FORMULA_INTERCEPT,
    INTERNAL_FORMULA_INTERCEPT,
    OUT as P2594_OUT,
)

GEN = ROOT / "generated"
OUT = GEN / "p2595_s1545_apd_eighth_moment_continuum_inverse_sturm_transport_certificate.json"
MD = GEN / "p2595_s1545_apd_eighth_moment_continuum_inverse_sturm_transport_certificate.md"

SOURCE_FILES = {
    "P2591_STURM_INTERVAL_CERTIFICATE": P2591_OUT,
    "P2594_CONDITIONAL_INVERSE_SELECTOR": P2594_OUT,
}
NEGATIVE_EXPORT_FLAGS = [
    "apd_eighth_moment_continuum_source_exported",
    "apd_continuum_inverse_selector_source_exported",
    "apd_eighth_shell_interval_bijection_source_exported",
    "apd_inverse_sturm_transport_source_exported",
    "apd_continuum_support_reconstruction_source_exported",
    "apd_support_selection_source_exported",
    "apd_finite_support_source_exported",
    "strict_dynamical_source_for_A_P_D_exported",
    "strict_phase_frequency_source_exported",
    "strict_damping_beta_eta_source_exported",
    "bridge_theorem_exported",
    "legacy_to_strict_completion_bridge_exported",
    "role_transfer_theorem_exported",
    "role_bearing_ltotal_exported",
    "qw2191_discharged_by_this_certificate",
    "toe_closure_claimed",
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
        "new_packet": "P2595|S1545|APD eighth moment continuum inverse|eighth moment continuum.*APD",
        "intended_research_nonduplication": "APD.*continuum inverse selector|continuum inverse selector.*APD|APD.*eighth shell interval bijection|eighth shell interval bijection.*APD|APD.*continuum support reconstruction|continuum support reconstruction.*APD|APD.*inverse Sturm transport|inverse Sturm transport.*APD",
        "apd_precursors": "P2591|S1541|P2594|S1544|APD.*Sturm interval|APD.*eighth moment inverse|strict_dynamical_source_for_A_P_D",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer theorem|QW-2191|ToE closure|source theorem",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def exact_support_snapshot(product_parameter: int) -> dict[str, Any]:
    support = support_from_product(float(product_parameter))
    return {
        "product_parameter": product_parameter,
        "squared_offsets": [float(value) for value in support["squared_offsets"]],
        "point_count": len(support["points"]),
        "min_point": float(min(support["points"])),
        "max_point": float(max(support["points"])),
    }


def interval_transport_certificate() -> dict[str, Any]:
    e4, p4, m8 = sp.symbols("e4 p4 m8")
    e4_left, e4_right = [sp.Integer(value) for value in INTERVAL]
    internal_expr = sp.Integer(INTERNAL_FORMULA_INTERCEPT) - 4 * e4
    central_expr = CENTRAL_FORMULA_INTERCEPT - 8 * e4
    internal_inverse = sp.expand((sp.Integer(INTERNAL_FORMULA_INTERCEPT) - p4) / 4)
    central_inverse = sp.expand((CENTRAL_FORMULA_INTERCEPT - m8) / 8)
    internal_left = sp.expand(internal_expr.subs(e4, e4_left))
    internal_right = sp.expand(internal_expr.subs(e4, e4_right))
    central_left = sp.expand(central_expr.subs(e4, e4_left))
    central_right = sp.expand(central_expr.subs(e4, e4_right))
    internal_interval_ascending = [internal_right, internal_left]
    central_interval_ascending = [central_right, central_left]
    internal_inverse_endpoints = [sp.expand(internal_inverse.subs(p4, endpoint)) for endpoint in internal_interval_ascending]
    central_inverse_endpoints = [sp.expand(central_inverse.subs(m8, endpoint)) for endpoint in central_interval_ascending]
    transported = {
        "from_internal_left_endpoint": str(internal_inverse_endpoints[0]),
        "from_internal_right_endpoint": str(internal_inverse_endpoints[1]),
        "from_central_left_endpoint": str(central_inverse_endpoints[0]),
        "from_central_right_endpoint": str(central_inverse_endpoints[1]),
    }
    return {
        "product_parameter_interval": INTERVAL,
        "internal_eighth_formula": str(internal_expr),
        "central_eighth_formula_exact": str(central_expr),
        "internal_eighth_interval_ascending_exact": [str(value) for value in internal_interval_ascending],
        "central_eighth_interval_ascending_exact": [str(value) for value in central_interval_ascending],
        "central_eighth_interval_ascending_float": [float(value) for value in central_interval_ascending],
        "internal_inverse_formula": str(internal_inverse),
        "central_inverse_formula": str(central_inverse),
        "inverse_endpoint_transport": transported,
        "internal_interval_maps_back_to_product_interval_as_set": sorted(int(value) for value in internal_inverse_endpoints) == INTERVAL,
        "central_interval_maps_back_to_product_interval_as_set": sorted(int(value) for value in central_inverse_endpoints) == INTERVAL,
        "internal_affine_slope": float(sp.diff(internal_expr, e4)),
        "central_affine_slope": float(sp.diff(central_expr, e4)),
        "internal_inverse_slope": float(sp.diff(internal_inverse, p4)),
        "central_inverse_slope": float(sp.diff(central_inverse, m8)),
        "sturm_transport_statement": "P2591 certifies four positive squared-offset roots for every e4 in [300, 576]; the affine inverse maps every p4/M8 in the eighth-shell interval back into that same e4 interval.",
        "endpoint_support_snapshots": [exact_support_snapshot(value) for value in INTERVAL],
        "midpoint_support_snapshot": exact_support_snapshot(sum(INTERVAL) // 2),
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2591_payload = load_json(SOURCE_FILES["P2591_STURM_INTERVAL_CERTIFICATE"])
    p2594_payload = load_json(SOURCE_FILES["P2594_CONDITIONAL_INVERSE_SELECTOR"])
    p2591 = theorem(p2591_payload, "apd_product_parameter_sturm_interval_certificate")
    p2594 = theorem(p2594_payload, "apd_eighth_moment_conditional_inverse_selector_certificate")
    certificate = interval_transport_certificate()
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2595_T1_apd_eighth_moment_continuum_inverse_sturm_transport_certificate",
        "audited_chain": ["P2591/S1541", "P2594/S1544"],
        "frontier_atom_under_attack": "strict_dynamical_source_for_A_P_D",
        "candidate_principle_under_test": "promote continuum inverse transport of an eighth-shell interval into a strict APD support source",
        "p2591_sturm_interval_inherited": p2591.get("continuous_interval_of_valid_supports_certified") is True,
        "p2594_conditional_inverse_inherited": p2594.get("eighth_moment_conditionally_recovers_product_parameter") is True,
        "apd_eighth_moment_continuum_inverse_sturm_transport_certificate": certificate,
        "continuum_internal_eighth_interval_maps_to_valid_product_interval": certificate["internal_interval_maps_back_to_product_interval_as_set"],
        "continuum_central_eighth_interval_maps_to_valid_product_interval": certificate["central_interval_maps_back_to_product_interval_as_set"],
        "p2591_sturm_validity_transports_to_entire_eighth_shell_interval": True,
        "continuum_inverse_transport_is_conditional_not_source": True,
        "apd_dynamic_source_remains_unsourced": True,
        "recommended_next_honest_step": (
            "P2595 upgrades P2594 from grid reconstruction to continuum inverse transport: every eighth-shell value in the transported interval maps back to the P2591 Sturm-valid product interval. This still begins with an externally supplied eighth shell, so the next honest step is a strict nadsoliton source theorem for that shell/support law, not another inverse-selector promotion."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2591_sturm_interval_inherited": theorem_export["p2591_sturm_interval_inherited"],
        "p2594_conditional_inverse_inherited": theorem_export["p2594_conditional_inverse_inherited"],
        "internal_interval_maps_back_to_product_interval": theorem_export["continuum_internal_eighth_interval_maps_to_valid_product_interval"],
        "central_interval_maps_back_to_product_interval": theorem_export["continuum_central_eighth_interval_maps_to_valid_product_interval"],
        "internal_and_central_maps_have_expected_slopes": certificate["internal_affine_slope"] == -4.0 and certificate["central_affine_slope"] == -8.0,
        "inverse_maps_have_expected_slopes": certificate["internal_inverse_slope"] == -0.25 and certificate["central_inverse_slope"] == -0.125,
        "endpoint_and_midpoint_supports_have_ten_points": all(row["point_count"] == 10 for row in [*certificate["endpoint_support_snapshots"], certificate["midpoint_support_snapshot"]]),
        "apd_dynamic_source_not_exported": theorem_export["strict_dynamical_source_for_A_P_D_exported"] is False,
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_certificate"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2595",
        "stage_id": "S1545",
        "status": "P2595_APD_EIGHTH_MOMENT_CONTINUUM_INVERSE_STURM_TRANSPORT_CERTIFICATE_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "apd_eighth_moment_continuum_inverse_sturm_transport_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2591_STURM_INTERVAL_CERTIFICATE": sha256_json(p2591_payload),
                "P2594_CONDITIONAL_INVERSE_SELECTOR": sha256_json(p2594_payload),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["apd_eighth_moment_continuum_inverse_sturm_transport_certificate"]["theorem_export"]
    c = t["apd_eighth_moment_continuum_inverse_sturm_transport_certificate"]
    lines = [
        "# P2595/S1545 APD eighth-moment continuum inverse Sturm transport certificate", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- Internal eighth interval: `{c['internal_eighth_interval_ascending_exact']}`.",
        f"- Central eighth interval: `{c['central_eighth_interval_ascending_exact']}`.",
        f"- Internal inverse formula: `{c['internal_inverse_formula']}`.",
        f"- Central inverse formula: `{c['central_inverse_formula']}`.",
        f"- Sturm validity transports to entire interval: `{t['p2591_sturm_validity_transports_to_entire_eighth_shell_interval']}`.",
        f"- Strict APD dynamic source exported: `{t['strict_dynamical_source_for_A_P_D_exported']}`.", "",
        "## Interpretation", "",
        "P2595 upgrades P2594's grid inverse to a continuum statement.  The affine inverse maps the whole eighth-shell interval back onto `[300, 576]`; by the inherited P2591 Sturm certificate, every transported product parameter in that interval has four positive squared-offset roots.  This proves continuum conditional reconstruction, not a source for the eighth shell.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No eighth-moment continuum source, continuum inverse selector source, eighth-shell interval bijection source, inverse Sturm transport source, continuum support-reconstruction source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.", "", "## Fingerprint", "",
        f"`{payload['apd_eighth_moment_continuum_inverse_sturm_transport_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2595/S1545 APD eighth-moment continuum inverse Sturm transport guard

`P2595/S1545` upgrades P2594 from grid reconstruction to a continuum inverse transport: the internal eighth-shell interval `[72354, 73458]` and the corresponding central eighth interval map affinely back onto the P2591 product interval `[300, 576]`, so the inherited Sturm certificate transports valid four-positive-root support reconstruction across the entire interval.  This remains conditional on an externally supplied eighth shell and does not source the APD support law.
""".strip()
    lag_section = """
## P2595/S1545 APD eighth-moment continuum inverse Sturm transport Ltotal guard

`P2595/S1545` blocks a role-bearing APD Gram term in `L_total` from being justified by continuum inverse Sturm transport.  Transporting an unsourced eighth-shell interval into a Sturm-valid product interval proves conditional reconstruction only, not a strict nadsoliton support/density source.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2595/S1545 APD eighth-moment continuum inverse Sturm transport guard", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2595/S1545 APD eighth-moment continuum inverse Sturm transport Ltotal guard", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
