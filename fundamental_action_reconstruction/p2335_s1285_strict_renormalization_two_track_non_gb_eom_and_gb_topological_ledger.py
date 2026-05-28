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
OUT = GEN / "p2335_s1285_strict_renormalization_two_track_non_gb_eom_and_gb_topological_ledger.json"
MD = GEN / "p2335_s1285_strict_renormalization_two_track_non_gb_eom_and_gb_topological_ledger.md"

SOURCE_FILES = {
    "P1853_B1_COEFFICIENTS": GEN / "p1853_s803_strict_b1_symbolic_evaluation_and_scheme_trace_checkpoint.json",
    "P2028_GB_QUOTIENT": GEN / "p2028_s978_strict_b1_gb_quotient_counterterm_identifiability_theorem.json",
    "P2034_TASK1_QUOTIENT_ONLY": GEN / "p2034_s984_strict_task1_quotient_only_renormalization_theorem.json",
    "P2333_SPACETIME_TABLE": GEN / "p2333_s1283_strict_bianchi_spacetime_component_table_gb_lift_probe.json",
    "P2334_EOM_GB_QUOTIENT_SCOPE": GEN / "p2334_s1284_strict_eom_gb_topological_quotient_scope_theorem.json",
}
CHANNELS_4 = ("R2", "Ric2", "Riem2", "GB")
TWO_TRACK_CHANNELS = ("non_gb_E_R2", "non_gb_E_Ric2", "topological_GB_ledger")


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
        "two-track|topological counterterm|GB.*ledger|a_GB|action-density|boundary|topological-number|quotient-only|EOM-only GB|non-GB quotient|Riem2.*4\\*Ric2",
        "fundamental_action_reconstruction",
    ]
    proc = subprocess.run(cmd, cwd=ROOT.parent, text=True, capture_output=True, check=False)
    return proc.stdout.splitlines()[:80]


def coeff_exprs(p1853: dict[str, Any]) -> dict[str, sp.Expr]:
    coeffs = p1853.get("b1_symbolic_evaluation", {}).get("evaluated_coefficients", {})
    return {name: sp.sympify(coeffs.get(name, {}).get("symbolic", "0")) for name in ("a_R2", "a_Ric2", "a_Riem2", "a_GB")}


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {key: load_json(path) for key, path in SOURCE_FILES.items()}
    coeffs = coeff_exprs(artifacts["P1853_B1_COEFFICIENTS"])
    a_r2 = coeffs["a_R2"]
    a_ric2 = coeffs["a_Ric2"]
    a_riem2 = coeffs["a_Riem2"]
    a_gb = coeffs["a_GB"]

    # In 4D density algebra: GB = Riem2 - 4*Ric2 + R2, so
    # Riem2 = GB + 4*Ric2 - R2.  This gives an exact algebraic split of the
    # action-density counterterm into an EOM-transportable non-GB part and a
    # separate topological GB ledger coefficient.
    split = {
        "b_EOM_R2": sp.factor(a_r2 - a_riem2),
        "b_EOM_Ric2": sp.factor(a_ric2 + 4 * a_riem2),
        "b_GB_topological_ledger": sp.factor(a_gb + a_riem2),
    }
    # Reconstruction matrix from two-track coefficients (b_R2,b_Ric2,b_GB) to
    # original density channels (R2,Ric2,Riem2,GB_label) with GB_label kept as
    # a ledger label, not as an independent EOM force.
    # Original action is reconstructed after substituting GB_density as the
    # actual topological density, not as another independent EOM channel.
    reconstruction_from_split = {
        "R2_density_coefficient_after_GB_expansion": sp.factor(split["b_EOM_R2"] + split["b_GB_topological_ledger"]),
        "Ric2_density_coefficient_after_GB_expansion": sp.factor(split["b_EOM_Ric2"] - 4 * split["b_GB_topological_ledger"]),
        "Riem2_density_coefficient_after_GB_expansion": sp.factor(split["b_GB_topological_ledger"]),
        "GB_label_coefficient_after_expansion": sp.Integer(0),
    }
    reconstruction_residuals = {
        "R2": sp.simplify(reconstruction_from_split["R2_density_coefficient_after_GB_expansion"] - (a_r2 + a_gb)),
        "Ric2": sp.simplify(reconstruction_from_split["Ric2_density_coefficient_after_GB_expansion"] - (a_ric2 - 4 * a_gb)),
        "Riem2": sp.simplify(reconstruction_from_split["Riem2_density_coefficient_after_GB_expansion"] - (a_riem2 + a_gb)),
    }
    # The expanded representation intentionally folds the original explicit
    # a_GB label into the actual topological density coefficient a_Riem2+a_GB.
    action_density_identity_zero = all(value == 0 for value in reconstruction_residuals.values())

    eom_transport_vector = sp.Matrix([split["b_EOM_R2"], split["b_EOM_Ric2"]])
    topological_ledger_vector = sp.Matrix([split["b_GB_topological_ledger"]])
    two_track_vector = sp.Matrix([split["b_EOM_R2"], split["b_EOM_Ric2"], split["b_GB_topological_ledger"]])
    exact_split_rank = int(sp.Matrix([[1, 0, -1, 0], [0, 1, 4, 0], [0, 0, 1, 1]]).rank())

    numeric_two_track = np.array([float(sp.N(value, 50)) for value in two_track_vector], dtype=float)
    numeric_norms = {
        "non_gb_eom_vector_l2": float(la.norm(numeric_two_track[:2], ord=2)),
        "topological_ledger_abs": float(abs(numeric_two_track[2])),
        "two_track_vector_l2": float(la.norm(numeric_two_track, ord=2)),
    }

    p2034 = artifacts["P2034_TASK1_QUOTIENT_ONLY"]
    p2334 = artifacts["P2334_EOM_GB_QUOTIENT_SCOPE"]
    p2334_checks = p2334.get("gatekeeper_checks", {})
    p2334_probe = p2334.get("strict_eom_gb_topological_quotient_scope_probe", {})
    p2334_eom_map = p2334_probe.get("universal_eom_quotient_map", {})

    two_track_ledger = {
        "ledger_id": "P2335_two_track_strict_renormalization_v1",
        "track_A_non_gb_eom_transportable_quotient": {
            "basis": ["E(R2)", "E(Ric2)"],
            "coefficients": {"b_EOM_R2": str(split["b_EOM_R2"]), "b_EOM_Ric2": str(split["b_EOM_Ric2"])},
            "numeric_coefficients": [float(sp.N(split["b_EOM_R2"], 30)), float(sp.N(split["b_EOM_Ric2"], 30))],
            "transport_status": "eligible_for_second_atlas_EOM_transport_checks",
        },
        "track_B_gb_topological_counterterm_ledger": {
            "basis": ["GB_density = Riem2 - 4*Ric2 + R2"],
            "ledger_coefficient_b_GB_topological": str(split["b_GB_topological_ledger"]),
            "numeric_ledger_coefficient": float(sp.N(split["b_GB_topological_ledger"], 30)),
            "status": "separate_topological_counterterm_ledger_not_EOM_transportable",
            "requires_future_witness": "action-density plus boundary/topological-number normalization witness before independent a_GB claim",
        },
        "exact_action_density_split": {
            "rule": "After expanding the original explicit GB label as GB=Riem2-4*Ric2+R2: a_R2*R2 + a_Ric2*Ric2 + a_Riem2*Riem2 + a_GB*GB = b_R2*R2 + b_Ric2*Ric2 + b_GB*(Riem2 - 4*Ric2 + R2)",
            "b_R2_equals": "a_R2 - a_Riem2",
            "b_Ric2_equals": "a_Ric2 + 4*a_Riem2",
            "b_GB_equals": "a_GB + a_Riem2",
            "reconstruction_residuals": {key: str(value) for key, value in reconstruction_residuals.items()},
            "action_density_identity_zero": action_density_identity_zero,
            "exact_split_rank": exact_split_rank,
        },
        "numeric_norms": numeric_norms,
    }

    theorem_export = {
        "theorem_name": "P2335 strict renormalization two-track split: non-GB EOM quotient plus GB topological ledger",
        "claim": (
            "Current strict exports should be read as a two-track renormalization object: "
            "(A) an EOM-transportable non-GB quotient with coefficients b_R2=a_R2-a_Riem2 and "
            "b_Ric2=a_Ric2+4*a_Riem2, and (B) a separate topological GB ledger with coefficient "
            "b_GB=a_GB+a_Riem2.  This exactly reconstructs the action-density counterterm modulo "
            "the four-dimensional GB identity, while preserving P2334's result that EOM-only data cannot "
            "identify an independent a_GB."
        ),
        "proof_witnesses": [
            "Exact symbolic substitution Riem2 = GB + 4*Ric2 - R2.",
            "Zero symbolic reconstruction residuals for R2, Ric2, and the folded Riem2+GB topological density coefficient.",
            "P2334 gatekeepers show EOM-only GB lifting remains blocked, so the GB coefficient is separated into a ledger instead of promoted.",
        ],
        "not_licensed": [
            "independent a_GB measurement",
            "boundary/topological-number normalization theorem",
            "full/global renormalization closure",
            "QW-2191 selector discharge",
            "G1/G3 update",
            "ToE closure",
        ],
        "next_honest_step": (
            "Run second-atlas transport on Track A only, and separately build a boundary/topological-number "
            "normalization witness for Track B before any independent a_GB claim."
        ),
    }

    probe = {
        "probe_id": "P2335_S1285_STRICT_RENORMALIZATION_TWO_TRACK_LEDGER",
        "source_hashes": {key: sha256_file(path) for key, path in SOURCE_FILES.items()},
        "repo_grep_audit": {
            "search_focus": "action-density/boundary/topological witness and two-track GB counterterm ledger",
            "top_hits": grep_hits(),
        },
        "input_coefficients": {key: str(value) for key, value in coeffs.items()},
        "two_track_ledger": two_track_ledger,
        "current_export_dependencies": {
            "p2034_local_verdict": p2034.get("local_verdict"),
            "p2334_result_kind": p2334.get("result_kind"),
            "p2334_exact_eom_rank": p2334_eom_map.get("exact_rank"),
            "p2334_exact_eom_nullity": p2334_eom_map.get("exact_nullity"),
            "p2334_current_exports_quotient_scope_only": p2334_checks.get("current_exports_are_quotient_scope_only"),
        },
        "theorem_export": theorem_export,
        "theorem_fingerprint_sha256": sha256_json(theorem_export),
    }

    gatekeeper_checks = {
        "p1853_coefficients_loaded": artifacts["P1853_B1_COEFFICIENTS"].get("packet_id") == "P1853",
        "p2034_quotient_only_loaded": p2034.get("local_verdict") == "PASS_QUOTIENT_ONLY_RENORMALIZATION_WITH_TRACE",
        "p2334_eom_only_no_go_loaded": p2334.get("result_kind") == "FORMAL_EOM_ONLY_GB_LIFT_NO_GO_CURRENT_STRICT_EXPORTS_QUOTIENT_SCOPE_ONLY",
        "action_density_identity_zero": action_density_identity_zero,
        "non_gb_track_has_two_coefficients": len(two_track_ledger["track_A_non_gb_eom_transportable_quotient"]["numeric_coefficients"]) == 2,
        "gb_track_separated_as_topological_ledger": two_track_ledger["track_B_gb_topological_counterterm_ledger"]["status"] == "separate_topological_counterterm_ledger_not_EOM_transportable",
        "independent_a_gb_not_claimed": True,
        "boundary_topological_number_witness_still_required": True,
        "no_full_global_renormalization_claimed": True,
        "no_selector_premise_added": True,
        "no_qw2191_discharge_claimed": True,
        "no_g1_g3_update_claimed": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2335_s1285_v1",
        "packet_id": "P2335",
        "stage_id": "S1285",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "result_kind": "STRICT_TWO_TRACK_RENORMALIZATION_LEDGER_EXPORTED_NO_FULL_CLOSURE",
        "strict_two_track_renormalization_ledger_probe": probe,
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": theorem_export["next_honest_step"],
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2335 strict two-track renormalization ledger\n\n"
        "Status: non-GB EOM quotient separated from GB topological ledger; no full/global closure claimed.\n\n"
        f"- Track A coefficients: `b_R2 = {split['b_EOM_R2']}`, `b_Ric2 = {split['b_EOM_Ric2']}`.\n"
        f"- Track B topological ledger coefficient: `b_GB = {split['b_GB_topological_ledger']}`.\n"
        f"- Action-density reconstruction residuals zero: `{action_density_identity_zero}`.\n"
        "- Track A may be transported by a second-atlas EOM theorem; Track B needs boundary/topological-number normalization.\n"
        "- No independent `a_GB`, no full/global renormalization, no QW-2191 discharge, no G1/G3 update, no ToE closure.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
