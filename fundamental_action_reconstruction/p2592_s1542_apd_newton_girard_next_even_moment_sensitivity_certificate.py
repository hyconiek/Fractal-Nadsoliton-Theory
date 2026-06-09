#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from typing import Any

import numpy as np
import sympy as sp

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2590_s1540_apd_finite_even_moment_shell_interval_nonuniqueness_audit import PRODUCT_PARAMETER_GRID
from p2589_s1539_apd_mirror_sixth_moment_shell_support_nonuniqueness_audit import (
    ELEMENTARY_SECOND_SUM,
    ELEMENTARY_THIRD_SUM,
    ENDPOINTS,
    EXPECTED_CARDINALITY,
    INTERNAL_FOURTH_SHELL,
    INTERNAL_SECOND_SHELL,
    INTERNAL_SIXTH_SHELL,
    MIRROR_CENTER,
    shell_support_witness,
    support_from_product,
)
from p2591_s1541_apd_product_parameter_sturm_interval_certificate import INTERVAL, OUT as P2591_OUT, theorem as p2590_theorem

GEN = ROOT / "generated"
OUT = GEN / "p2592_s1542_apd_newton_girard_next_even_moment_sensitivity_certificate.json"
MD = GEN / "p2592_s1542_apd_newton_girard_next_even_moment_sensitivity_certificate.md"

SOURCE_FILES = {
    "P2591_STURM_INTERVAL_CERTIFICATE": P2591_OUT,
}
ELEMENTARY_FIRST_SUM = INTERNAL_SECOND_SHELL
ENDPOINT_OFFSET = ENDPOINTS[1] - MIRROR_CENTER
NEGATIVE_EXPORT_FLAGS = [
    "apd_eighth_moment_shell_source_exported",
    "apd_next_even_moment_selector_source_exported",
    "apd_newton_girard_selector_source_exported",
    "apd_product_parameter_source_exported",
    "apd_support_selection_source_exported",
    "apd_finite_support_source_exported",
    "apd_positive_measure_source_exported",
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
        "new_packet": "P2592|S1542|APD Newton Girard next moment|Newton Girard.*APD",
        "intended_research_nonduplication": "APD.*next even moment shell|next even moment shell.*APD|APD.*eighth moment sensitivity|eighth moment sensitivity.*APD|APD.*central eighth support moment|central eighth support moment.*APD|APD.*moment prefix extension|moment prefix extension.*APD",
        "apd_precursors": "P2590|S1540|P2591|S1541|APD.*Sturm interval|APD.*finite even moment shell|strict_dynamical_source_for_A_P_D",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer theorem|QW-2191|ToE closure|source theorem",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def newton_girard_internal_p4(product_parameter: float) -> float:
    return float(
        ELEMENTARY_FIRST_SUM * INTERNAL_SIXTH_SHELL
        - ELEMENTARY_SECOND_SUM * INTERNAL_FOURTH_SHELL
        + ELEMENTARY_THIRD_SUM * INTERNAL_SECOND_SHELL
        - 4.0 * product_parameter
    )


def central_eighth_from_internal_p4(internal_p4: float) -> float:
    return float(2.0 * ENDPOINT_OFFSET ** 8 + 2.0 * internal_p4)


def numerical_snapshot(product_parameter: float) -> dict[str, Any]:
    witness = shell_support_witness(support_from_product(product_parameter))
    points = np.array(witness["points"], dtype=float)
    numerical_central_eighth = float(np.sum((points - MIRROR_CENTER) ** 8))
    formula_internal_p4 = newton_girard_internal_p4(product_parameter)
    formula_central_eighth = central_eighth_from_internal_p4(formula_internal_p4)
    return {
        "product_parameter": float(product_parameter),
        "squared_offsets": witness["squared_offsets"],
        "shared_internal_second_shell": witness["internal_second_shell"],
        "shared_internal_fourth_shell": witness["internal_fourth_shell"],
        "shared_internal_sixth_shell": witness["internal_sixth_shell"],
        "newton_girard_internal_eighth_shell_p4": formula_internal_p4,
        "central_eighth_moment_unweighted_formula": formula_central_eighth,
        "central_eighth_moment_unweighted_numeric": numerical_central_eighth,
        "formula_numeric_abs_error": abs(formula_central_eighth - numerical_central_eighth),
        "vandermonde_rank": witness["vandermonde_rank"],
        "has_expected_internal_second_fourth_sixth_shell": witness["has_expected_internal_second_fourth_sixth_shell"],
    }


def next_even_moment_certificate() -> dict[str, Any]:
    e4 = sp.symbols("e4")
    p4_expr = sp.expand(ELEMENTARY_FIRST_SUM * INTERNAL_SIXTH_SHELL - ELEMENTARY_SECOND_SUM * INTERNAL_FOURTH_SHELL + ELEMENTARY_THIRD_SUM * INTERNAL_SECOND_SHELL - 4 * e4)
    central_expr = sp.expand(2 * sp.Rational(21, 4) ** 8 + 2 * p4_expr)
    interval_endpoints = [numerical_snapshot(float(value)) for value in INTERVAL]
    grid_snapshots = [numerical_snapshot(float(value)) for value in PRODUCT_PARAMETER_GRID]
    central_values = [snapshot["central_eighth_moment_unweighted_formula"] for snapshot in grid_snapshots]
    internal_values = [snapshot["newton_girard_internal_eighth_shell_p4"] for snapshot in grid_snapshots]
    selector_recovery = [
        {
            "product_parameter": snapshot["product_parameter"],
            "recovered_product_from_internal_eighth_shell": float((float(p4_expr.subs(e4, 0)) - snapshot["newton_girard_internal_eighth_shell_p4"]) / 4.0),
            "recovered_product_from_central_eighth_moment": float((float(central_expr.subs(e4, 0)) - snapshot["central_eighth_moment_unweighted_formula"]) / 8.0),
        }
        for snapshot in grid_snapshots
    ]
    return {
        "numpy_version": np.__version__,
        "sympy_version": sp.__version__,
        "product_parameter_interval_inherited": INTERVAL,
        "elementary_sums_fixed": {
            "e1": ELEMENTARY_FIRST_SUM,
            "e2": ELEMENTARY_SECOND_SUM,
            "e3": ELEMENTARY_THIRD_SUM,
            "e4_free_product_parameter": "e4",
        },
        "fixed_lower_power_sums": {
            "p1_internal_second_shell": INTERNAL_SECOND_SHELL,
            "p2_internal_fourth_shell": INTERNAL_FOURTH_SHELL,
            "p3_internal_sixth_shell": INTERNAL_SIXTH_SHELL,
        },
        "newton_girard_formula_internal_eighth_shell_p4": str(p4_expr),
        "central_eighth_formula_with_endpoints": str(central_expr),
        "internal_eighth_shell_slope_by_product_parameter": float(sp.diff(p4_expr, e4)),
        "central_eighth_moment_slope_by_product_parameter": float(sp.diff(central_expr, e4)),
        "interval_endpoint_snapshots": interval_endpoints,
        "grid_snapshots": grid_snapshots,
        "grid_internal_eighth_shell_values": internal_values,
        "grid_central_eighth_moment_values": central_values,
        "grid_internal_eighth_shell_distinct_count": len({round(value, 8) for value in internal_values}),
        "grid_central_eighth_moment_distinct_count": len({round(value, 8) for value in central_values}),
        "selector_recovery_from_next_moment": selector_recovery,
        "max_formula_numeric_abs_error": max(snapshot["formula_numeric_abs_error"] for snapshot in grid_snapshots + interval_endpoints),
        "lower_shells_constant_but_next_shell_varies": (
            len({round(snapshot["shared_internal_second_shell"], 8) for snapshot in grid_snapshots}) == 1
            and len({round(snapshot["shared_internal_fourth_shell"], 8) for snapshot in grid_snapshots}) == 1
            and len({round(snapshot["shared_internal_sixth_shell"], 8) for snapshot in grid_snapshots}) == 1
            and len({round(value, 8) for value in internal_values}) == len(grid_snapshots)
        ),
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2591_payload = load_json(SOURCE_FILES["P2591_STURM_INTERVAL_CERTIFICATE"])
    p2591 = p2590_theorem(p2591_payload, "apd_product_parameter_sturm_interval_certificate")
    certificate = next_even_moment_certificate()
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2592_T1_apd_newton_girard_next_even_moment_sensitivity_certificate",
        "audited_chain": ["P2591/S1541"],
        "frontier_atom_under_attack": "strict_dynamical_source_for_A_P_D",
        "candidate_principle_under_test": "extend the finite even-moment shell prefix by the next even shell and treat that extension as an APD support selector/source",
        "p2591_sturm_interval_inherited": p2591.get("continuous_interval_of_valid_supports_certified") is True,
        "apd_newton_girard_next_even_moment_sensitivity_certificate": certificate,
        "lower_second_fourth_sixth_shells_do_not_select_next_shell": certificate["lower_shells_constant_but_next_shell_varies"],
        "next_even_moment_would_recover_product_parameter_if_supplied": True,
        "next_even_moment_is_selector_coordinate_not_source": True,
        "apd_dynamic_source_remains_unsourced": True,
        "recommended_next_honest_step": (
            "Do not treat the eighth shell as a strict APD source merely because Newton-Girard recovers the free product parameter from it. P2592 proves the next even moment is exactly the missing selector coordinate for the P2591 interval; the honest next step is to derive that next shell/support law from strict nadsoliton dynamics, or else keep the product-parameter freedom explicit."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2591_sturm_interval_inherited": theorem_export["p2591_sturm_interval_inherited"],
        "internal_eighth_shell_has_negative_four_slope": certificate["internal_eighth_shell_slope_by_product_parameter"] == -4.0,
        "central_eighth_moment_has_negative_eight_slope": certificate["central_eighth_moment_slope_by_product_parameter"] == -8.0,
        "grid_next_shell_values_distinct": certificate["grid_internal_eighth_shell_distinct_count"] == len(PRODUCT_PARAMETER_GRID),
        "formula_matches_numeric_snapshots": certificate["max_formula_numeric_abs_error"] <= 1.0e-6,
        "lower_shells_constant_but_next_shell_varies": certificate["lower_shells_constant_but_next_shell_varies"],
        "apd_dynamic_source_not_exported": theorem_export["strict_dynamical_source_for_A_P_D_exported"] is False,
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_certificate"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2592",
        "stage_id": "S1542",
        "status": "P2592_APD_NEWTON_GIRARD_NEXT_EVEN_MOMENT_SENSITIVITY_CERTIFICATE_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "apd_newton_girard_next_even_moment_sensitivity_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2591_STURM_INTERVAL_CERTIFICATE": sha256_json(p2591_payload),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["apd_newton_girard_next_even_moment_sensitivity_certificate"]["theorem_export"]
    certificate = t["apd_newton_girard_next_even_moment_sensitivity_certificate"]
    lines = [
        "# P2592/S1542 APD Newton-Girard next even-moment sensitivity certificate", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- Inherited product-parameter interval: `{certificate['product_parameter_interval_inherited']}`.",
        f"- Internal eighth-shell formula: `{certificate['newton_girard_formula_internal_eighth_shell_p4']}`.",
        f"- Central eighth-moment formula: `{certificate['central_eighth_formula_with_endpoints']}`.",
        f"- Lower shells constant but next shell varies: `{t['lower_second_fourth_sixth_shells_do_not_select_next_shell']}`.",
        f"- Strict APD dynamic source exported: `{t['strict_dynamical_source_for_A_P_D_exported']}`.", "",
        "## Interpretation", "",
        "P2592 identifies the exact next missing coordinate after the P2591 Sturm interval.  Newton-Girard gives `p4 = 74658 - 4*e4` for the internal squared-offset eighth shell, so the same fixed second/fourth/sixth shell prefix leaves the next even shell linearly free.  Supplying that eighth shell would recover the product parameter, but the certificate does not derive that shell from strict nadsoliton dynamics.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No eighth-moment shell source, next-even-moment selector source, Newton-Girard selector source, product-parameter source, support-selection source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.", "", "## Fingerprint", "",
        f"`{payload['apd_newton_girard_next_even_moment_sensitivity_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2592/S1542 APD Newton-Girard next even-moment sensitivity guard

`P2592/S1542` identifies the next missing coordinate in the P2591 APD product-parameter interval: Newton-Girard gives the internal eighth shell `p4 = 74658 - 4*e4`, hence the fixed second/fourth/sixth shell prefix leaves the next even shell linearly free.  Adding the eighth shell would select the product parameter only if that shell is strictly sourced; by itself it is not a strict APD support law.
""".strip()
    lag_section = """
## P2592/S1542 APD Newton-Girard next even-moment sensitivity Ltotal guard

`P2592/S1542` blocks a role-bearing APD Gram term in `L_total` from being justified by appending an unsourced eighth-moment shell to the finite moment prefix.  The next shell is a selector coordinate for the P2591 continuum, not a strict nadsoliton-derived density/support source.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2592/S1542 APD Newton-Girard next even-moment sensitivity guard", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2592/S1542 APD Newton-Girard next even-moment sensitivity Ltotal guard", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
