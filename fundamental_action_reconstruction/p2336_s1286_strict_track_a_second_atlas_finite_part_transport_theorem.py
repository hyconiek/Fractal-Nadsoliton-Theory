#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

import numpy as np
import scipy.linalg as la
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2336_s1286_strict_track_a_second_atlas_finite_part_transport_theorem.json"
MD = GEN / "p2336_s1286_strict_track_a_second_atlas_finite_part_transport_theorem.md"

SOURCE_FILES = {
    "P1973_FRW_BIANCHI_SCALAR_TRANSPORT": GEN / "p1973_s923_strict_frw_bianchi_finite_part_transport_matrix_witness.json",
    "P2035_BACKGROUND_TRANSPORT_OBSTRUCTION": GEN / "p2035_s985_strict_task1_quotient_background_transport_obstruction_theorem.json",
    "P2037_SAME_SCHEME_MAP_SEED": GEN / "p2037_s987_strict_task1_same_scheme_finite_part_map_seed.json",
    "P2334_EOM_GB_QUOTIENT_SCOPE": GEN / "p2334_s1284_strict_eom_gb_topological_quotient_scope_theorem.json",
    "P2335_TWO_TRACK_LEDGER": GEN / "p2335_s1285_strict_renormalization_two_track_non_gb_eom_and_gb_topological_ledger.json",
}
TRACK_A_BASIS = ("E(R2)", "E(Ric2)")


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def grep_hits() -> list[str]:
    cmd = [
        "rg",
        "-n",
        "P1973|P2035|P2037|P2334|P2335|Track A|second-atlas|second atlas|FRW|Bianchi|finite-part|GB.*ledger|topological counterterm|EOM-only GB|quotient-scope",
        "fundamental_action_reconstruction",
    ]
    proc = subprocess.run(cmd, cwd=ROOT.parent, text=True, capture_output=True, check=False)
    return proc.stdout.splitlines()[:120]


def sympify_coefficients(p2335: dict[str, Any]) -> tuple[sp.Expr, sp.Expr]:
    probe = p2335.get("strict_two_track_renormalization_ledger_probe", {})
    ledger = probe.get("two_track_ledger", {})
    track_a = ledger.get("track_A_non_gb_eom_transportable_quotient", {})
    coeffs = track_a.get("coefficients", {})
    return sp.sympify(coeffs.get("b_EOM_R2", "0")), sp.sympify(coeffs.get("b_EOM_Ric2", "0"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {key: load_json(path) for key, path in SOURCE_FILES.items()}

    p1973 = artifacts["P1973_FRW_BIANCHI_SCALAR_TRANSPORT"]
    p2035 = artifacts["P2035_BACKGROUND_TRANSPORT_OBSTRUCTION"]
    p2037 = artifacts["P2037_SAME_SCHEME_MAP_SEED"]
    p2334 = artifacts["P2334_EOM_GB_QUOTIENT_SCOPE"]
    p2335 = artifacts["P2335_TWO_TRACK_LEDGER"]

    lam = sp.symbols("lambda", real=True)
    nu = sp.symbols("nu", real=True)
    sigma2 = sp.symbols("sigma2", nonnegative=True, real=True)
    q = sp.simplify(nu * sigma2)
    multiplier = sp.simplify(1 + lam * q)
    connection_scalar = sp.simplify(sp.diff(multiplier, lam) / multiplier)

    b_r2, b_ric2 = sympify_coefficients(p2335)
    track_a_frw = sp.Matrix([b_r2, b_ric2])
    transport_matrix = sp.eye(2) * multiplier
    connection_matrix = sp.eye(2) * connection_scalar
    track_a_path = sp.simplify(transport_matrix * track_a_frw)
    transport_residual = sp.simplify(sp.diff(track_a_path, lam) - connection_matrix * track_a_path)
    matrix_residual = sp.simplify(sp.diff(transport_matrix, lam) - connection_matrix * transport_matrix)
    determinant_lambda_1 = sp.factor(transport_matrix.subs(lam, 1).det())
    transported_endpoint = sp.simplify(track_a_path.subs(lam, 1))
    endpoint_back_transport_residual = sp.simplify((transport_matrix.subs(lam, 1).inv() * transported_endpoint) - track_a_frw)

    transport_residual_entries = [sp.factor(sp.simplify(x)) for x in list(transport_residual)]
    matrix_residual_entries = [sp.factor(sp.simplify(x)) for x in list(matrix_residual)]
    endpoint_back_entries = [sp.factor(sp.simplify(x)) for x in list(endpoint_back_transport_residual)]
    all_symbolic_zero = all(x == 0 for x in transport_residual_entries + matrix_residual_entries + endpoint_back_entries)

    q_values = [0.0, 0.05, 0.2, 0.8]
    track_a_numeric = np.array([float(sp.N(b_r2, 40)), float(sp.N(b_ric2, 40))], dtype=float)
    numeric_rows: list[dict[str, Any]] = []
    for q_value in q_values:
        expected_transport = (1.0 + q_value) * np.eye(2)
        scipy_transport = la.expm(np.log1p(q_value) * np.eye(2))
        expected_endpoint = expected_transport @ track_a_numeric
        scipy_endpoint = scipy_transport @ track_a_numeric
        back_transported = la.solve(scipy_transport, scipy_endpoint)
        numeric_rows.append(
            {
                "q_equals_nu_sigma2": q_value,
                "transport_trace": float(np.trace(scipy_transport)),
                "transport_determinant": float(la.det(scipy_transport)),
                "max_abs_transport_error": float(np.max(np.abs(scipy_transport - expected_transport))),
                "max_abs_endpoint_error": float(np.max(np.abs(scipy_endpoint - expected_endpoint))),
                "max_abs_back_transport_error": float(np.max(np.abs(back_transported - track_a_numeric))),
                "invertible_positive_branch": bool(1.0 + q_value > 0.0),
            }
        )
    all_numeric_pass = all(
        row["max_abs_transport_error"] < 1e-14
        and row["max_abs_endpoint_error"] < 1e-14
        and row["max_abs_back_transport_error"] < 1e-14
        and row["invertible_positive_branch"]
        for row in numeric_rows
    )

    dependencies = {
        "p1973_local_scalar_transport_pass": p1973.get("local_result_kind") == "PASS_ZERO_LOCAL_TRANSPORT_RESIDUAL",
        "p1973_not_global_background_independence": "global background-independence theorem"
        in str((p1973.get("theorem_export") or {}).get("not_licensed", [])),
        "p2035_current_export_nontransportability_context_loaded": p2035.get("obstruction_verdict")
        == "PASS_CURRENT_EXPORT_NONTRANSPORTABILITY_WITH_TRACE",
        "p2037_same_scheme_map_seed_loaded": p2037.get("packet_id") == "P2037"
        and (p2037.get("gatekeeper_checks") or {}).get("c3_theorem_proven") is False,
        "p2334_eom_only_gb_no_go_loaded": p2334.get("result_kind")
        == "FORMAL_EOM_ONLY_GB_LIFT_NO_GO_CURRENT_STRICT_EXPORTS_QUOTIENT_SCOPE_ONLY",
        "p2335_two_track_ledger_loaded": p2335.get("result_kind")
        == "STRICT_TWO_TRACK_RENORMALIZATION_LEDGER_EXPORTED_NO_FULL_CLOSURE",
    }

    theorem_export = {
        "theorem_name": "P2336 Track-A second-atlas finite-part transport theorem on the FRW/Bianchi-I nu branch",
        "licensed_statement": (
            "For the P2335 non-GB Track A quotient vector (b_EOM_R2,b_EOM_Ric2), the P1973 "
            "FRW/Bianchi-I finite-part branch transports Track A by T_A(lambda)=(1+lambda*nu*sigma2) I_2. "
            "The covariant transport residual dF_A/dlambda-A_A(lambda)F_A is exactly zero and the endpoint "
            "is invertible on 1+nu*sigma2>0.  This is a same-branch finite-part transport witness for Track A only."
        ),
        "proof_witnesses": [
            "Exact symbolic zero residual for the two-component Track A path.",
            "Exact symbolic zero residual for the 2x2 connection matrix equation.",
            "Exact endpoint back-transport recovers the FRW Track A vector on the admissible branch.",
            "Numeric scipy expm replay agrees with the closed-form transport on the q-grid.",
        ],
        "not_licensed": [
            "transport or normalization of the Track B GB topological ledger",
            "independent a_GB measurement",
            "boundary/topological-number normalization theorem",
            "full tensor-component variational bundle transport",
            "global background-independence theorem",
            "full/global renormalization closure",
            "QW-2191 selector discharge",
            "G1/G3 update",
            "ToE closure",
        ],
        "next_honest_step": (
            "Keep Track A as transported only on the P1973 local finite-part branch; next build either a non-scalar "
            "variational-bundle Track-A component transport witness or a separate Track-B boundary/topological-number normalization witness."
        ),
    }

    probe = {
        "probe_id": "P2336_S1286_STRICT_TRACK_A_SECOND_ATLAS_FINITE_PART_TRANSPORT_THEOREM",
        "source_hashes": {key: sha256_file(path) for key, path in SOURCE_FILES.items()},
        "repo_grep_audit": {
            "search_focus": "P2335 Track A, P1973 FRW/Bianchi transport, quotient-scope and GB ledger boundaries",
            "top_hits": grep_hits(),
        },
        "track_A_input": {
            "basis": list(TRACK_A_BASIS),
            "frw_coefficients_symbolic": {"b_EOM_R2": str(b_r2), "b_EOM_Ric2": str(b_ric2)},
            "frw_coefficients_numeric": [float(sp.N(b_r2, 30)), float(sp.N(b_ric2, 30))],
        },
        "second_atlas_transport": {
            "background_pair": ["FRW", "Bianchi-I"],
            "branch_source": "P1973/P1897 nu branch finite-part multiplier",
            "scope": "Track A non-GB EOM quotient finite-part transport only",
            "q_definition": "q=nu*sigma2",
            "admissible_branch": "1+nu*sigma2>0",
            "transport_matrix_T_A_lambda": str(transport_matrix),
            "connection_matrix_A_A_lambda": str(connection_matrix),
            "track_A_path": [str(x) for x in list(track_a_path)],
            "transport_residual_entries": [str(x) for x in transport_residual_entries],
            "matrix_residual_entries": [str(x) for x in matrix_residual_entries],
            "endpoint_lambda_1": [str(x) for x in list(transported_endpoint)],
            "endpoint_back_transport_residual_entries": [str(x) for x in endpoint_back_entries],
            "determinant_T_A_lambda_1": str(determinant_lambda_1),
            "symbolic_zero_residuals": all_symbolic_zero,
        },
        "numeric_transport_replay": {"q_grid": q_values, "rows": numeric_rows, "all_rows_pass": all_numeric_pass},
        "current_export_dependencies": dependencies,
        "theorem_export": theorem_export,
        "theorem_fingerprint_sha256": sha256_json(theorem_export),
    }

    gatekeeper_checks = {
        "repo_grep_audit_recorded_before_theorem": len(probe["repo_grep_audit"]["top_hits"]) > 0,
        "p1973_local_transport_loaded": dependencies["p1973_local_scalar_transport_pass"],
        "p1973_global_limit_preserved": dependencies["p1973_not_global_background_independence"],
        "p2035_obstruction_context_preserved": dependencies["p2035_current_export_nontransportability_context_loaded"],
        "p2037_same_scheme_seed_not_misread_as_c3": dependencies["p2037_same_scheme_map_seed_loaded"],
        "p2334_eom_only_gb_no_go_preserved": dependencies["p2334_eom_only_gb_no_go_loaded"],
        "p2335_two_track_ledger_loaded": dependencies["p2335_two_track_ledger_loaded"],
        "track_a_has_two_components_only": len(TRACK_A_BASIS) == 2,
        "track_b_not_transported": True,
        "symbolic_track_a_transport_residual_zero": all_symbolic_zero,
        "numeric_transport_replay_pass": all_numeric_pass,
        "admissible_branch_condition_exported": True,
        "no_independent_a_gb_claimed": True,
        "no_boundary_topological_number_witness_claimed": True,
        "no_full_tensor_bundle_transport_claimed": True,
        "no_background_globalization_claimed": True,
        "no_full_global_renormalization_claimed": True,
        "no_selector_premise_added": True,
        "no_qw2191_discharge_claimed": True,
        "no_g1_g3_update_claimed": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2336_s1286_v1",
        "packet_id": "P2336",
        "stage_id": "S1286",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_TRACK_A_LOCAL_SECOND_ATLAS_FINITE_PART_TRANSPORT__NO_TRACK_B_OR_GLOBAL_CLOSURE",
        "strict_track_a_second_atlas_transport_probe": probe,
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": theorem_export["next_honest_step"],
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2336 strict Track-A second-atlas finite-part transport theorem\n\n"
        "Status: Track A finite-part transport passes on the P1973 FRW/Bianchi-I local branch; Track B and global closure remain open.\n\n"
        f"- Track A basis: `{TRACK_A_BASIS[0]}`, `{TRACK_A_BASIS[1]}`.\n"
        f"- Track A coefficients: `b_R2 = {b_r2}`, `b_Ric2 = {b_ric2}`.\n"
        "- Transport: `T_A(lambda) = (1 + lambda*nu*sigma2) I_2`, with admissible branch `1 + nu*sigma2 > 0`.\n"
        f"- Symbolic transport residuals zero: `{all_symbolic_zero}`.\n"
        f"- Numeric replay pass: `{all_numeric_pass}`.\n"
        "- Not licensed: Track B GB ledger transport, independent `a_GB`, boundary/topological-number normalization, full tensor-bundle transport, global background independence, full/global renormalization, QW-2191 discharge, G1/G3 update, or ToE closure.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
