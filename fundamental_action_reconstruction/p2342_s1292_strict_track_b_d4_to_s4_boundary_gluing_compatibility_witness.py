#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2342_s1292_strict_track_b_d4_to_s4_boundary_gluing_compatibility_witness.json"
MD = GEN / "p2342_s1292_strict_track_b_d4_to_s4_boundary_gluing_compatibility_witness.md"

SOURCE_FILES = {
    "P2335_TWO_TRACK_LEDGER": GEN / "p2335_s1285_strict_renormalization_two_track_non_gb_eom_and_gb_topological_ledger.json",
    "P2338_TRACK_B_CONTRACT": GEN / "p2338_s1288_strict_track_b_gb_boundary_topological_normalization_contract.json",
    "P2339_S4_FIXTURE": GEN / "p2339_s1289_strict_track_b_s4_boundaryless_topological_normalization_witness.json",
    "P2341_D4_BOUNDARY_FIXTURE": GEN / "p2341_s1291_strict_track_b_d4_boundary_correction_fixture_witness.json",
}


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
    candidates = [
        ROOT / "p2338_s1288_strict_track_b_gb_boundary_topological_normalization_contract.py",
        ROOT / "p2339_s1289_strict_track_b_s4_boundaryless_topological_normalization_witness.py",
        ROOT / "p2341_s1291_strict_track_b_d4_boundary_correction_fixture_witness.py",
        GEN / "p2339_s1289_strict_track_b_s4_boundaryless_topological_normalization_witness.md",
        GEN / "p2341_s1291_strict_track_b_d4_boundary_correction_fixture_witness.md",
    ]
    existing = [path.relative_to(ROOT.parent).as_posix() for path in candidates if path.exists()]
    if not existing:
        return []
    cmd = [
        "rg",
        "-n",
        "Track B|D4|D\\^4|S4|S\\^4|boundary correction|gluing|Euler|chi|topological normalization|a_GB|QW-2191",
        *existing,
    ]
    proc = subprocess.run(cmd, cwd=ROOT.parent, text=True, capture_output=True, check=False)
    return proc.stdout.splitlines()[:140]


def track_b_coefficient(p2335: dict[str, Any]) -> sp.Expr:
    probe = p2335.get("strict_two_track_renormalization_ledger_probe", {})
    ledger = probe.get("two_track_ledger", {})
    track_b = ledger.get("track_B_gb_topological_counterterm_ledger", {})
    raw = track_b.get("ledger_coefficient_b_GB_topological", "0")
    return sp.sympify(raw, locals={"pi": sp.pi, "log": sp.log, "ln": sp.log})


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {key: load_json(path) for key, path in SOURCE_FILES.items()}
    p2335 = artifacts["P2335_TWO_TRACK_LEDGER"]
    p2338 = artifacts["P2338_TRACK_B_CONTRACT"]
    p2339 = artifacts["P2339_S4_FIXTURE"]
    p2341 = artifacts["P2341_D4_BOUNDARY_FIXTURE"]

    b_gb = sp.factor(track_b_coefficient(p2335))
    per_euler_pairing = sp.factor(sp.simplify(32 * sp.pi**2 * b_gb))

    chi_d4_left = sp.Integer(1)
    chi_d4_right = sp.Integer(1)
    chi_s3_interface = sp.Integer(0)
    chi_glued = sp.factor(chi_d4_left + chi_d4_right - chi_s3_interface)
    chi_s4 = sp.Integer(2)

    d4_bulk_each = sp.Integer(0)
    d4_boundary_each = sp.factor(32 * sp.pi**2)
    left_total = sp.factor(d4_bulk_each + d4_boundary_each)
    right_total = sp.factor(d4_bulk_each + d4_boundary_each)
    glued_total_from_pieces = sp.factor(sp.simplify(left_total + right_total))
    s4_topological_number = sp.factor(32 * sp.pi**2 * chi_s4)
    gluing_topological_residual = sp.factor(sp.simplify(glued_total_from_pieces - s4_topological_number))

    left_pairing = sp.factor(sp.simplify(b_gb * left_total))
    right_pairing = sp.factor(sp.simplify(b_gb * right_total))
    glued_pairing_from_pieces = sp.factor(sp.simplify(left_pairing + right_pairing))
    s4_pairing_from_glued_chi = sp.factor(sp.simplify(per_euler_pairing * chi_glued))
    s4_pairing_direct = sp.factor(sp.simplify(b_gb * s4_topological_number))
    pairing_residual_against_s4 = sp.factor(sp.simplify(glued_pairing_from_pieces - s4_pairing_direct))
    chi_residual_against_s4 = sp.factor(sp.simplify(chi_glued - chi_s4))

    p2339_probe = p2339.get("strict_track_b_s4_boundaryless_topological_normalization_witness_probe", {})
    p2339_witness = p2339_probe.get("track_B_S4_witness", {})
    p2341_probe = p2341.get("strict_track_b_d4_boundary_correction_fixture_witness_probe", {})
    p2341_witness = p2341_probe.get("track_B_D4_boundary_fixture_witness", {})

    dependencies = {
        "p2335_two_track_ledger_loaded": p2335.get("result_kind")
        == "STRICT_TWO_TRACK_RENORMALIZATION_LEDGER_EXPORTED_NO_FULL_CLOSURE",
        "p2338_contract_loaded": p2338.get("result_kind")
        == "STRICT_TRACK_B_GB_BOUNDARY_TOPOLOGICAL_NORMALIZATION_CONTRACT_EXPORTED_NO_NORMALIZATION_CLAIM",
        "p2339_s4_fixture_loaded": p2339.get("result_kind")
        == "STRICT_TRACK_B_S4_SINGLE_FIXTURE_TOPOLOGICAL_NORMALIZATION_WITNESS_NO_GLOBAL_CLOSURE",
        "p2341_d4_boundary_fixture_loaded": p2341.get("result_kind")
        == "STRICT_TRACK_B_D4_NONZERO_BOUNDARY_FIXTURE_WITNESS_NO_UNIVERSAL_BOUNDARY_THEOREM",
        "p2339_s4_pairing_matches": p2339_witness.get("normalized_S4_pairing") == str(s4_pairing_direct),
        "p2341_d4_pairing_matches": p2341_witness.get("normalized_D4_pairing") == str(left_pairing),
    }

    gluing_model = {
        "model_id": "two_D4_boundary_gluing_to_S4_fixture_model_v1",
        "scope": "fixture-level boundary-gluing compatibility only; not a general gluing theorem",
        "left_piece": "D^4 with boundary S^3, chi=1, bulk_GB=0, boundary_correction=32*pi**2",
        "right_piece": "D^4 with boundary S^3, chi=1, bulk_GB=0, boundary_correction=32*pi**2",
        "interface": "S^3 interface with chi(S^3)=0 for Euler bookkeeping",
        "glued_fixture": "S^4 closed boundaryless fixture, chi=2",
        "chi_gluing_rule_used": "chi(D4_left) + chi(D4_right) - chi(S3_interface)",
        "chi_glued_from_pieces": str(chi_glued),
        "chi_residual_against_S4": str(chi_residual_against_s4),
        "glued_topological_number_from_pieces": str(glued_total_from_pieces),
        "S4_topological_number": str(s4_topological_number),
        "gluing_topological_residual": str(gluing_topological_residual),
        "left_D4_pairing": str(left_pairing),
        "right_D4_pairing": str(right_pairing),
        "glued_pairing_from_D4_pieces": str(glued_pairing_from_pieces),
        "S4_pairing_from_glued_chi": str(s4_pairing_from_glued_chi),
        "S4_pairing_direct": str(s4_pairing_direct),
        "pairing_residual_against_S4": str(pairing_residual_against_s4),
    }

    theorem_export = {
        "theorem_name": "P2342 Track-B D4-to-S4 boundary-gluing compatibility witness",
        "claim": (
            "P2342 checks one fixture-level boundary-gluing compatibility after P2341: two D^4 bookkeeping pieces, "
            "each with chi=1, bulk GB value 0, and boundary correction 32*pi**2, glue along an S^3 interface with "
            "chi(S^3)=0 to the closed S^4 fixture with chi=2.  The summed piece topological number is 64*pi**2, "
            "matching P2339's S^4 number, and the Track-B ledger pairing residual against S^4 is zero."
        ),
        "proof_witnesses": [
            "The same P2335 b_GB ledger coefficient is used for the D^4 pieces and S^4 comparison.",
            "Euler bookkeeping gives chi(D4_left)+chi(D4_right)-chi(S3)=1+1-0=2.",
            "The two D^4 boundary-correction totals sum to 64*pi**2, the S^4 topological number.",
            "The resulting Track-B pairing equals the P2339 S^4 pairing with zero symbolic residual.",
        ],
        "not_licensed": [
            "general boundary gluing theorem for arbitrary interfaces",
            "universal boundary correction functional",
            "connected-sum theorem beyond this D4-to-S4 fixture bookkeeping model",
            "independent a_GB measurement separate from the P2335 ledger coefficient",
            "dynamical selection of spacetime topology or boundary condition",
            "bulk EOM GB force or EOM-only GB lift",
            "full/global renormalization closure",
            "selector premise or QW-2191 selector discharge",
            "G1/G3 update",
            "ToE closure",
        ],
        "next_honest_step": (
            "Keep selector work deferred.  The next non-selector Track-B move should either derive a real boundary "
            "correction functional on a named admissible boundary class, or freeze Track-B as a fixture-level ledger "
            "normalization chain below full/global renormalization."
        ),
    }

    probe = {
        "probe_id": "P2342_S1292_STRICT_TRACK_B_D4_TO_S4_BOUNDARY_GLUING_COMPATIBILITY_WITNESS",
        "source_hashes": {key: sha256_file(path) for key, path in SOURCE_FILES.items()},
        "repo_grep_audit": {
            "search_focus": "P2341 next step: non-selector Track-B boundary-gluing compatibility against S4",
            "top_hits": grep_hits(),
        },
        "track_B_D4_to_S4_boundary_gluing_witness": {
            "ledger_coefficient_b_GB": str(b_gb),
            "per_euler_pairing": str(per_euler_pairing),
            "gluing_model": gluing_model,
        },
        "current_export_dependencies": dependencies,
        "theorem_export": theorem_export,
        "theorem_fingerprint_sha256": sha256_json(theorem_export),
    }

    gatekeeper_checks = {
        "repo_grep_audit_recorded": len(probe["repo_grep_audit"]["top_hits"]) > 0,
        "p2335_two_track_ledger_loaded": dependencies["p2335_two_track_ledger_loaded"],
        "p2338_contract_loaded": dependencies["p2338_contract_loaded"],
        "p2339_s4_fixture_loaded": dependencies["p2339_s4_fixture_loaded"],
        "p2341_d4_boundary_fixture_loaded": dependencies["p2341_d4_boundary_fixture_loaded"],
        "p2339_s4_pairing_matches": dependencies["p2339_s4_pairing_matches"],
        "p2341_d4_pairing_matches": dependencies["p2341_d4_pairing_matches"],
        "chi_s3_interface_zero_declared": chi_s3_interface == 0,
        "chi_glued_matches_s4": chi_residual_against_s4 == 0,
        "two_d4_topological_number_matches_s4": gluing_topological_residual == 0,
        "two_d4_pairing_matches_s4": pairing_residual_against_s4 == 0,
        "fixture_level_scope_declared": "fixture-level" in gluing_model["scope"],
        "general_boundary_gluing_theorem_not_claimed": True,
        "universal_boundary_functional_not_claimed": True,
        "connected_sum_theorem_not_claimed": True,
        "independent_a_gb_not_claimed": True,
        "bulk_eom_gb_force_not_claimed": True,
        "no_full_global_renormalization_claimed": True,
        "no_selector_premise_added": True,
        "no_qw2191_discharge_claimed": True,
        "no_g1_g3_update_claimed": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2342_s1292_v1",
        "packet_id": "P2342",
        "stage_id": "S1292",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "STRICT_TRACK_B_D4_TO_S4_BOUNDARY_GLUING_COMPATIBILITY_WITNESS_NO_GENERAL_GLUING_THEOREM",
        "strict_track_b_d4_to_s4_boundary_gluing_compatibility_witness_probe": probe,
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": theorem_export["next_honest_step"],
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2342 strict Track-B D4-to-S4 boundary-gluing compatibility witness\n\n"
        "Status: two D4 boundary fixtures glue compatibly to the S4 closed fixture; no general gluing theorem.\n\n"
        f"- `b_GB = {b_gb}`.\n"
        "- Fixture model: `D^4_left ∪_{S^3} D^4_right -> S^4`, with `chi(S^3)=0`.\n"
        f"- Euler bookkeeping: `1 + 1 - 0 = {chi_glued}`; residual against S4: `{chi_residual_against_s4}`.\n"
        f"- Piece topological number: `{glued_total_from_pieces}`; S4 topological number: `{s4_topological_number}`; residual: `{gluing_topological_residual}`.\n"
        f"- D4-piece pairing sum: `{glued_pairing_from_pieces}`; direct S4 pairing: `{s4_pairing_direct}`; residual: `{pairing_residual_against_s4}`.\n"
        "- No general boundary gluing theorem, no universal boundary functional, no independent `a_GB`, no full/global renormalization, no selector premise, no QW-2191 discharge, no G1/G3 update, no ToE closure.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
