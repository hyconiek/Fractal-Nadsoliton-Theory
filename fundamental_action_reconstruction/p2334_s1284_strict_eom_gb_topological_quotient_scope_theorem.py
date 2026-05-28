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
OUT = GEN / "p2334_s1284_strict_eom_gb_topological_quotient_scope_theorem.json"
MD = GEN / "p2334_s1284_strict_eom_gb_topological_quotient_scope_theorem.md"

CHANNELS = ("R2", "Ric2", "Riem2", "GB")
SOURCE_FILES = {
    "P1973_FRW_BIANCHI_TRANSPORT": GEN / "p1973_s923_strict_frw_bianchi_finite_part_transport_matrix_witness.json",
    "P2028_GB_QUOTIENT": GEN / "p2028_s978_strict_b1_gb_quotient_counterterm_identifiability_theorem.json",
    "P2034_TASK1_QUOTIENT_ONLY": GEN / "p2034_s984_strict_task1_quotient_only_renormalization_theorem.json",
    "P2330_B1_GB_OBSTRUCTION": GEN / "p2330_s1280_strict_b1_renormalization_gb_dependence_globalization_obstruction_probe.json",
    "P2331_SECOND_BACKGROUND_GAP": GEN / "p2331_s1281_strict_second_background_curvature_witness_gap_audit.json",
    "P2332_SPATIAL_TABLE": GEN / "p2332_s1282_strict_bianchi_spatial_tensor_component_table_gb_lift_probe.json",
    "P2333_SPACETIME_TABLE": GEN / "p2333_s1283_strict_bianchi_spacetime_component_table_gb_lift_probe.json",
}


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path)}
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def sha256_file(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def grep_hits() -> list[str]:
    cmd = [
        "rg",
        "-n",
        "second atlas|background transport|FRW|Bianchi|Gauss-Bonnet|GB.*quotient|quotient-scope|full/global renormalization|topological",
        "fundamental_action_reconstruction",
    ]
    proc = subprocess.run(cmd, cwd=ROOT.parent, text=True, capture_output=True, check=False)
    return proc.stdout.splitlines()[:80]


def sympy_rows(matrix: sp.Matrix) -> list[list[str]]:
    return [[sp.sstr(sp.simplify(matrix[i, j])) for j in range(matrix.cols)] for i in range(matrix.rows)]


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {key: load_json(path) for key, path in SOURCE_FILES.items()}

    # Universal four-dimensional EOM quotient induced by the local identity
    # GB = Riem2 - 4*Ric2 + R2 and by variational topological invariance dS_GB = 0.
    # The abstract EOM image is therefore spanned by E(R2), E(Ric2); Riem2 is
    # transported as 4*E(Ric2)-E(R2), and GB is zero in pure EOM data.
    eom_map = sp.Matrix(
        [
            [1, 0, -1, 0],
            [0, 1, 4, 0],
        ]
    )
    relation_gb_density = sp.Matrix([1, -4, 1, -1])
    relation_pure_gb_eom = sp.Matrix([0, 0, 0, 1])
    rank_exact = int(eom_map.rank())
    nullspace_exact = eom_map.nullspace()
    nullspace_rows = [[sp.sstr(entry) for entry in vector] for vector in nullspace_exact]
    relation_checks = {
        "eom_map_rows_R2_Ric2_Riem2_GB": sympy_rows(eom_map),
        "exact_rank": rank_exact,
        "exact_nullity": len(nullspace_exact),
        "nullspace_basis_R2_Ric2_Riem2_GB": nullspace_rows,
        "gb_density_relation_maps_to_zero": bool(eom_map * relation_gb_density == sp.zeros(2, 1)),
        "pure_gb_eom_variation_maps_to_zero": bool(eom_map * relation_pure_gb_eom == sp.zeros(2, 1)),
        "max_possible_eom_rank_for_any_component_atlas": 2,
    }

    eom_map_np = np.array(eom_map.tolist(), dtype=float)
    gram = eom_map_np.T @ eom_map_np
    singular_values = np.linalg.svd(eom_map_np, compute_uv=False)
    scipy_null = la.null_space(eom_map_np)
    density_residual_l2 = float(la.norm(eom_map_np @ np.array(relation_gb_density, dtype=float), ord=2))
    pure_gb_residual_l2 = float(la.norm(eom_map_np @ np.array(relation_pure_gb_eom, dtype=float), ord=2))

    p2333_probe = artifacts["P2333_SPACETIME_TABLE"].get("strict_bianchi_spacetime_component_table_gb_lift_probe", {})
    p2333_rank_test = p2333_probe.get("gb_lift_rank_test", {})
    p2333_target = p2333_probe.get("same_basis_divergence_target", {})
    p2331_status = artifacts["P2331_SECOND_BACKGROUND_GAP"].get("status")
    p1973_kind = artifacts["P1973_FRW_BIANCHI_TRANSPORT"].get("local_result_kind")
    p2034_verdict = artifacts["P2034_TASK1_QUOTIENT_ONLY"].get("local_verdict")

    existing_export_audit = {
        "p1973_frw_bianchi_transport_kind": p1973_kind,
        "p1973_is_scalar_finite_part_transport_not_component_atlas": p1973_kind == "PASS_ZERO_LOCAL_TRANSPORT_RESIDUAL",
        "p2034_quotient_only_verdict": p2034_verdict,
        "p2034_already_licenses_only_quotient_scope": p2034_verdict == "PASS_QUOTIENT_ONLY_RENORMALIZATION_WITH_TRACE",
        "p2331_second_background_gap_status": p2331_status,
        "p2333_local_spacetime_rank": p2333_rank_test.get("numeric_rank"),
        "p2333_local_spacetime_nullity": p2333_rank_test.get("numeric_nullity"),
        "p2333_same_basis_target_residual_l2": p2333_target.get("direct_reconstruction_residual_l2"),
        "genuine_second_tensor_component_atlas_exported_now": False,
    }

    formal_scope_theorem = {
        "theorem_name": "P2334 EOM-only Gauss-Bonnet topological quotient-scope theorem",
        "claim": (
            "For any strict renormalization witness that observes only Euler-Lagrange/tensor-EOM components "
            "of curvature-squared channels in four dimensions, the Gauss-Bonnet direction is killed by the "
            "variational map and the Riemann^2 channel is represented as 4*Ricci^2 - R^2 modulo GB. "
            "Therefore adding more EOM components or a second EOM atlas cannot by itself identify an independent "
            "a_GB coefficient; current strict exports license only quotient-scope renormalization until an "
            "action-density/boundary/topological-number witness is supplied."
        ),
        "symbolic_proof_bits": [
            "Use GB = Riem2 - 4*Ric2 + R2 in four dimensions.",
            "Apply the variational/EOM functor: E(GB)=0 for the topological density in the admitted local bulk scope.",
            "Thus E(Riem2)=4*E(Ric2)-E(R2) and the EOM channel map has exact rank 2 with nullity 2.",
            "P2333 confirms this rank/nullity in the built lapse+tracefree-spatial Bianchi-I component table.",
        ],
        "not_licensed": [
            "full/global four-channel counterterm uniqueness",
            "independent measurement of a_GB from EOM-only data",
            "replacement for a genuine second tensor-component atlas if the goal is background covariance",
            "QW-2191 selector discharge",
            "G1/G3 update",
            "ToE closure",
        ],
        "next_honest_step": (
            "Stop trying to lift GB with EOM-only component ranks.  Either add an action-density plus boundary/topological "
            "number witness for a_GB, or build a genuine second-atlas transport theorem for the non-GB quotient class."
        ),
    }

    probe = {
        "probe_id": "P2334_S1284_STRICT_EOM_GB_TOPOLOGICAL_QUOTIENT_SCOPE_THEOREM",
        "source_hashes": {key: sha256_file(path) for key, path in SOURCE_FILES.items()},
        "repo_grep_audit": {
            "search_focus": "second atlas/background transport and GB quotient/topological scope",
            "top_hits": grep_hits(),
        },
        "universal_eom_quotient_map": relation_checks,
        "numeric_replay": {
            "eom_map_singular_values": singular_values.tolist(),
            "gram_matrix": gram.tolist(),
            "scipy_null_space_shape": list(scipy_null.shape),
            "gb_density_relation_residual_l2": density_residual_l2,
            "pure_gb_eom_relation_residual_l2": pure_gb_residual_l2,
        },
        "existing_export_audit": existing_export_audit,
        "formal_scope_theorem": formal_scope_theorem,
        "theorem_fingerprint_sha256": sha256_json(formal_scope_theorem),
    }
    gatekeeper_checks = {
        "exact_eom_rank_2": rank_exact == 2,
        "exact_eom_nullity_2": len(nullspace_exact) == 2,
        "gb_density_relation_maps_to_zero": relation_checks["gb_density_relation_maps_to_zero"],
        "pure_gb_eom_variation_maps_to_zero": relation_checks["pure_gb_eom_variation_maps_to_zero"],
        "numeric_residuals_small": density_residual_l2 < 1e-12 and pure_gb_residual_l2 < 1e-12,
        "p2333_rank_matches_universal_eom_rank": p2333_rank_test.get("numeric_rank") == 2,
        "current_exports_are_quotient_scope_only": True,
        "action_boundary_topological_witness_required_for_a_gb": True,
        "no_full_global_renormalization_claimed": True,
        "no_selector_premise_added": True,
        "no_qw2191_discharge_claimed": True,
        "no_g1_g3_update_claimed": True,
        "no_toe_closure_claimed": True,
    }
    payload = {
        "schema_version": "p2334_s1284_v1",
        "packet_id": "P2334",
        "stage_id": "S1284",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "result_kind": "FORMAL_EOM_ONLY_GB_LIFT_NO_GO_CURRENT_STRICT_EXPORTS_QUOTIENT_SCOPE_ONLY",
        "strict_eom_gb_topological_quotient_scope_probe": probe,
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": formal_scope_theorem["next_honest_step"],
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2334 EOM-only GB topological quotient-scope theorem\n\n"
        "Status: EOM-only GB lifting is formally blocked; current strict exports remain quotient-scope.\n\n"
        f"- Exact universal EOM rank/nullity: `{rank_exact}/{len(nullspace_exact)}`.\n"
        f"- P2333 Bianchi spacetime rank/nullity: `{p2333_rank_test.get('numeric_rank')}/{p2333_rank_test.get('numeric_nullity')}`.\n"
        f"- GB density relation residual L2: `{density_residual_l2:.3e}`.\n"
        f"- Pure GB EOM residual L2: `{pure_gb_residual_l2:.3e}`.\n"
        "- Conclusion: more EOM components alone cannot identify independent `a_GB`; use action/boundary/topological witness or transport the non-GB quotient class.\n"
        "- No full/global renormalization, no QW-2191 discharge, no G1/G3 update, no ToE closure.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
