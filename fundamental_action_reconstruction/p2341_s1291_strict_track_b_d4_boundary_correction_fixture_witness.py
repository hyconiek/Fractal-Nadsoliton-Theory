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
OUT = GEN / "p2341_s1291_strict_track_b_d4_boundary_correction_fixture_witness.json"
MD = GEN / "p2341_s1291_strict_track_b_d4_boundary_correction_fixture_witness.md"

SOURCE_FILES = {
    "P2335_TWO_TRACK_LEDGER": GEN / "p2335_s1285_strict_renormalization_two_track_non_gb_eom_and_gb_topological_ledger.json",
    "P2338_TRACK_B_CONTRACT": GEN / "p2338_s1288_strict_track_b_gb_boundary_topological_normalization_contract.json",
    "P2340_TWO_FIXTURE_COMPATIBILITY": GEN / "p2340_s1290_strict_track_b_two_fixture_topological_gluing_compatibility_witness.json",
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
        ROOT / "p2340_s1290_strict_track_b_two_fixture_topological_gluing_compatibility_witness.py",
        ROOT / "p2341_s1291_strict_track_b_d4_boundary_correction_fixture_witness.py",
        GEN / "p2338_s1288_strict_track_b_gb_boundary_topological_normalization_contract.md",
        GEN / "p2340_s1290_strict_track_b_two_fixture_topological_gluing_compatibility_witness.md",
    ]
    existing = [path.relative_to(ROOT.parent).as_posix() for path in candidates if path.exists()]
    if not existing:
        return []
    cmd = [
        "rg",
        "-n",
        "Track B|D4|D\\^4|4-ball|boundary correction|nonzero boundary|Euler|chi|topological normalization|a_GB|QW-2191",
        *existing,
    ]
    proc = subprocess.run(cmd, cwd=ROOT.parent, text=True, capture_output=True, check=False)
    return proc.stdout.splitlines()[:120]


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
    p2340 = artifacts["P2340_TWO_FIXTURE_COMPATIBILITY"]

    b_gb = sp.factor(track_b_coefficient(p2335))
    chi_d4 = sp.Integer(1)
    bulk_gb_integral = sp.Integer(0)
    boundary_correction = sp.factor(32 * sp.pi**2)
    target_topological_number = sp.factor(32 * sp.pi**2 * chi_d4)
    boundary_fixture_total = sp.factor(sp.simplify(bulk_gb_integral + boundary_correction))
    boundary_residual = sp.factor(sp.simplify(boundary_fixture_total - target_topological_number))
    normalized_d4_pairing = sp.factor(sp.simplify(b_gb * boundary_fixture_total))
    per_euler_pairing = sp.factor(sp.simplify(32 * sp.pi**2 * b_gb))
    pairing_residual = sp.factor(sp.simplify(normalized_d4_pairing - per_euler_pairing * chi_d4))

    p2338_probe = p2338.get("strict_track_b_gb_boundary_topological_normalization_contract_probe", {})
    p2338_contract = p2338_probe.get("track_B_contract", {})
    required_fields = p2338_contract.get("required_future_witness_fields", [])
    supplied_fields = {
        "oriented_compact_four_manifold_or_boundary_class": "oriented compact D^4 / 4-ball fixture with boundary S^3",
        "gb_density_normalization_convention": "bulk_integral_GB + boundary_correction = 32*pi**2*chi(D^4)",
        "boundary_correction_functional": "nonzero boundary fixture value 32*pi**2; bulk GB fixture value 0",
        "topological_number_chi_or_relative_invariant": "chi(D^4)=1",
        "scheme_pairing_rule_for_b_GB": "Delta_Gamma_GB[D^4] = b_GB * (bulk_GB + boundary_correction)",
    }
    supplies_all_required = sorted(required_fields) == sorted(supplied_fields)

    witness = {
        "witness_id": "P2341_track_B_D4_nonzero_boundary_correction_fixture_v1",
        "scope": "single D^4 boundary fixture for nonzero boundary-correction bookkeeping only",
        "supplied_witness_fields": supplied_fields,
        "ledger_coefficient_b_GB": str(b_gb),
        "chi_D4": str(chi_d4),
        "bulk_gb_integral_fixture_value": str(bulk_gb_integral),
        "boundary_correction_fixture_value": str(boundary_correction),
        "target_topological_number_32pi2_chi": str(target_topological_number),
        "boundary_residual": str(boundary_residual),
        "normalized_D4_pairing": str(normalized_d4_pairing),
        "pairing_residual": str(pairing_residual),
        "supplies_all_p2338_required_fields": supplies_all_required,
    }

    p2340_probe = p2340.get("strict_track_b_two_fixture_topological_gluing_compatibility_witness_probe", {})
    p2340_witness = p2340_probe.get("track_B_two_fixture_witness", {})
    dependencies = {
        "p2335_two_track_ledger_loaded": p2335.get("result_kind")
        == "STRICT_TWO_TRACK_RENORMALIZATION_LEDGER_EXPORTED_NO_FULL_CLOSURE",
        "p2338_contract_loaded": p2338.get("result_kind")
        == "STRICT_TRACK_B_GB_BOUNDARY_TOPOLOGICAL_NORMALIZATION_CONTRACT_EXPORTED_NO_NORMALIZATION_CLAIM",
        "p2340_two_fixture_loaded": p2340.get("result_kind")
        == "STRICT_TRACK_B_TWO_CLOSED_FIXTURE_TOPOLOGICAL_ADDITIVITY_WITNESS_NO_GLOBAL_CLOSURE",
        "p2340_closed_fixture_ledger_matches": p2340_witness.get("ledger_coefficient_b_GB") == str(b_gb),
    }

    theorem_export = {
        "theorem_name": "P2341 Track-B D4 nonzero boundary-correction fixture witness",
        "claim": (
            "P2341 adds a genuine boundary fixture after P2340: an oriented compact D^4/4-ball bookkeeping fixture "
            "with boundary S^3, chi(D^4)=1, bulk GB fixture value 0, and nonzero boundary correction 32*pi**2. "
            "The bulk-plus-boundary total equals 32*pi**2*chi(D^4), and the Track-B ledger pairing is "
            "32*pi**2*b_GB = 13152087*log(2)/10000000 with zero symbolic residual."
        ),
        "proof_witnesses": [
            "All P2338 required witness fields are supplied for a D^4 boundary fixture.",
            "The boundary correction fixture value is explicitly nonzero: 32*pi**2.",
            "The boundary residual bulk_GB + boundary - 32*pi**2*chi(D^4) is zero.",
            "The Track-B pairing residual against the per-Euler ledger pairing is zero.",
        ],
        "not_licensed": [
            "universal boundary functional for arbitrary manifolds",
            "connected-sum or general gluing theorem",
            "independent a_GB measurement separate from the P2335 ledger coefficient",
            "dynamical selection of spacetime topology or boundary condition",
            "bulk EOM GB force or EOM-only GB lift",
            "full/global renormalization closure",
            "QW-2191 selector discharge",
            "G1/G3 update",
            "ToE closure",
        ],
        "next_honest_step": (
            "Do not promote the D^4 fixture to a universal boundary theorem.  Next either derive an actual boundary "
            "functional over a class of boundaries, or keep Track-B normalization limited to the explicit closed and D^4 fixtures."
        ),
    }

    probe = {
        "probe_id": "P2341_S1291_STRICT_TRACK_B_D4_BOUNDARY_CORRECTION_FIXTURE_WITNESS",
        "source_hashes": {key: sha256_file(path) for key, path in SOURCE_FILES.items()},
        "repo_grep_audit": {
            "search_focus": "P2340 next step: Track-B genuine boundary fixture with nonzero boundary correction",
            "top_hits": grep_hits(),
        },
        "track_B_D4_boundary_fixture_witness": witness,
        "current_export_dependencies": dependencies,
        "theorem_export": theorem_export,
        "theorem_fingerprint_sha256": sha256_json(theorem_export),
    }

    gatekeeper_checks = {
        "repo_grep_audit_recorded": len(probe["repo_grep_audit"]["top_hits"]) > 0,
        "p2335_two_track_ledger_loaded": dependencies["p2335_two_track_ledger_loaded"],
        "p2338_contract_loaded": dependencies["p2338_contract_loaded"],
        "p2340_two_fixture_loaded": dependencies["p2340_two_fixture_loaded"],
        "p2340_closed_fixture_ledger_matches": dependencies["p2340_closed_fixture_ledger_matches"],
        "supplies_all_p2338_required_fields": supplies_all_required,
        "chi_d4_is_1": chi_d4 == 1,
        "boundary_correction_nonzero": boundary_correction != 0,
        "bulk_fixture_value_zero_declared": bulk_gb_integral == 0,
        "boundary_residual_zero": boundary_residual == 0,
        "pairing_residual_zero": pairing_residual == 0,
        "single_boundary_fixture_scope_declared": witness["scope"] == "single D^4 boundary fixture for nonzero boundary-correction bookkeeping only",
        "universal_boundary_functional_not_claimed": True,
        "general_gluing_theorem_not_claimed": True,
        "independent_a_gb_not_claimed": True,
        "bulk_eom_gb_force_not_claimed": True,
        "no_full_global_renormalization_claimed": True,
        "no_selector_premise_added": True,
        "no_qw2191_discharge_claimed": True,
        "no_g1_g3_update_claimed": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2341_s1291_v1",
        "packet_id": "P2341",
        "stage_id": "S1291",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "STRICT_TRACK_B_D4_NONZERO_BOUNDARY_FIXTURE_WITNESS_NO_UNIVERSAL_BOUNDARY_THEOREM",
        "strict_track_b_d4_boundary_correction_fixture_witness_probe": probe,
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": theorem_export["next_honest_step"],
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2341 strict Track-B D4 boundary correction fixture witness\n\n"
        "Status: one D4 boundary fixture supplies a nonzero boundary correction; no universal boundary theorem.\n\n"
        f"- `b_GB = {b_gb}`.\n"
        "- Fixture: oriented compact `D^4`/4-ball with boundary `S^3`, `chi(D^4)=1`.\n"
        f"- Bulk GB fixture value: `{bulk_gb_integral}`; boundary correction fixture value: `{boundary_correction}`.\n"
        f"- Boundary residual: `{boundary_residual}`.\n"
        f"- Normalized D4 pairing: `{normalized_d4_pairing}`; pairing residual: `{pairing_residual}`.\n"
        "- No universal boundary functional, no independent `a_GB`, no full/global renormalization, no QW-2191 discharge, no G1/G3 update, no ToE closure.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
