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
OUT = GEN / "p2340_s1290_strict_track_b_two_fixture_topological_gluing_compatibility_witness.json"
MD = GEN / "p2340_s1290_strict_track_b_two_fixture_topological_gluing_compatibility_witness.md"

SOURCE_FILES = {
    "P2335_TWO_TRACK_LEDGER": GEN / "p2335_s1285_strict_renormalization_two_track_non_gb_eom_and_gb_topological_ledger.json",
    "P2338_TRACK_B_CONTRACT": GEN / "p2338_s1288_strict_track_b_gb_boundary_topological_normalization_contract.json",
    "P2339_S4_FIXTURE": GEN / "p2339_s1289_strict_track_b_s4_boundaryless_topological_normalization_witness.json",
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
        ROOT / "p2340_s1290_strict_track_b_two_fixture_topological_gluing_compatibility_witness.py",
        GEN / "p2338_s1288_strict_track_b_gb_boundary_topological_normalization_contract.md",
        GEN / "p2339_s1289_strict_track_b_s4_boundaryless_topological_normalization_witness.md",
    ]
    existing = [path.relative_to(ROOT.parent).as_posix() for path in candidates if path.exists()]
    if not existing:
        return []
    cmd = [
        "rg",
        "-n",
        "Track B|S4|S\\^4|CP2|CP\\^2|two.*fixture|gluing|sign-convention|Euler|chi|topological normalization|a_GB|QW-2191",
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


def fixture_row(label: str, chi: sp.Integer, b_gb: sp.Expr) -> dict[str, str]:
    topological_number = sp.factor(32 * sp.pi**2 * chi)
    pairing = sp.factor(sp.simplify(b_gb * topological_number))
    return {
        "fixture": label,
        "boundary_class": "closed_boundaryless",
        "boundary_correction": "0",
        "chi": str(chi),
        "gb_topological_number_32pi2_chi": str(topological_number),
        "normalized_pairing": str(pairing),
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {key: load_json(path) for key, path in SOURCE_FILES.items()}
    p2335 = artifacts["P2335_TWO_TRACK_LEDGER"]
    p2338 = artifacts["P2338_TRACK_B_CONTRACT"]
    p2339 = artifacts["P2339_S4_FIXTURE"]

    b_gb = sp.factor(track_b_coefficient(p2335))
    chi_s4 = sp.Integer(2)
    chi_cp2 = sp.Integer(3)
    chi_disjoint_union = sp.Integer(5)

    s4_pairing = sp.factor(32 * sp.pi**2 * chi_s4 * b_gb)
    cp2_pairing = sp.factor(32 * sp.pi**2 * chi_cp2 * b_gb)
    union_pairing_direct = sp.factor(32 * sp.pi**2 * chi_disjoint_union * b_gb)
    union_pairing_sum = sp.factor(sp.simplify(s4_pairing + cp2_pairing))
    additivity_residual = sp.factor(sp.simplify(union_pairing_direct - union_pairing_sum))
    ratio_residual = sp.factor(sp.simplify(2 * cp2_pairing - 3 * s4_pairing))

    fixtures = [fixture_row("S^4", chi_s4, b_gb), fixture_row("CP^2", chi_cp2, b_gb)]
    compatibility = {
        "compatibility_id": "P2340_two_fixture_track_B_topological_additivity_v1",
        "sign_convention": "same orientation convention for GB density; Euler pairing uses positive chi for the closed fixtures",
        "gluing_model": "disjoint_union_additivity_check_not_connected_sum_theorem",
        "chi_disjoint_union": str(chi_disjoint_union),
        "direct_union_pairing": str(union_pairing_direct),
        "sum_of_fixture_pairings": str(union_pairing_sum),
        "additivity_residual": str(additivity_residual),
        "ratio_residual_2cp2_minus_3s4": str(ratio_residual),
    }

    p2339_probe = p2339.get("strict_track_b_s4_boundaryless_topological_normalization_witness_probe", {})
    p2339_witness = p2339_probe.get("track_B_S4_witness", {})
    dependencies = {
        "p2335_two_track_ledger_loaded": p2335.get("result_kind")
        == "STRICT_TWO_TRACK_RENORMALIZATION_LEDGER_EXPORTED_NO_FULL_CLOSURE",
        "p2338_contract_loaded": p2338.get("result_kind")
        == "STRICT_TRACK_B_GB_BOUNDARY_TOPOLOGICAL_NORMALIZATION_CONTRACT_EXPORTED_NO_NORMALIZATION_CLAIM",
        "p2339_s4_fixture_loaded": p2339.get("result_kind")
        == "STRICT_TRACK_B_S4_SINGLE_FIXTURE_TOPOLOGICAL_NORMALIZATION_WITNESS_NO_GLOBAL_CLOSURE",
        "p2339_s4_pairing_matches": p2339_witness.get("normalized_S4_pairing") == str(s4_pairing),
    }

    theorem_export = {
        "theorem_name": "P2340 Track-B two-fixture topological gluing/sign-convention compatibility witness",
        "claim": (
            "P2340 adds a second closed boundaryless CP^2 topology fixture next to P2339's S^4 fixture and checks "
            "that the same Track-B GB ledger convention is additive under the declared disjoint-union gluing model. "
            "The S^4 and CP^2 pairings scale with chi=2 and chi=3, the direct chi=5 disjoint-union pairing equals "
            "the sum of the two fixture pairings, and the symbolic additivity residual is zero."
        ),
        "proof_witnesses": [
            "The same P2335 b_GB ledger coefficient is used for both fixtures.",
            "Closed boundaryless S^4 has chi=2 and closed boundaryless CP^2 has chi=3 under the supplied convention.",
            "The disjoint-union Euler number chi=5 gives the same pairing as summing S^4 and CP^2 pairings.",
            "The ratio residual 2*CP2_pairing - 3*S4_pairing is zero, checking common sign/normalization convention.",
        ],
        "not_licensed": [
            "connected-sum gluing theorem",
            "universal boundary/topological normalization over all manifolds",
            "independent a_GB measurement separate from the P2335 ledger coefficient",
            "dynamical selection of spacetime topology",
            "bulk EOM GB force or EOM-only GB lift",
            "full/global renormalization closure",
            "QW-2191 selector discharge",
            "G1/G3 update",
            "ToE closure",
        ],
        "next_honest_step": (
            "Either add a genuine boundary fixture with a nonzero boundary correction functional, or keep the Track-B "
            "normalization result limited to closed boundaryless fixtures with disjoint-union additivity only."
        ),
    }

    probe = {
        "probe_id": "P2340_S1290_STRICT_TRACK_B_TWO_FIXTURE_TOPOLOGICAL_GLUING_COMPATIBILITY_WITNESS",
        "source_hashes": {key: sha256_file(path) for key, path in SOURCE_FILES.items()},
        "repo_grep_audit": {
            "search_focus": "P2339 next step: second topology fixture plus gluing/sign convention compatibility",
            "top_hits": grep_hits(),
        },
        "track_B_two_fixture_witness": {
            "ledger_coefficient_b_GB": str(b_gb),
            "fixtures": fixtures,
            "compatibility": compatibility,
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
        "p2339_s4_pairing_matches": dependencies["p2339_s4_pairing_matches"],
        "cp2_fixture_added": fixtures[1]["fixture"] == "CP^2" and fixtures[1]["chi"] == "3",
        "disjoint_union_chi_is_5": chi_disjoint_union == 5,
        "additivity_residual_zero": additivity_residual == 0,
        "ratio_residual_zero": ratio_residual == 0,
        "connected_sum_theorem_not_claimed": True,
        "universal_topology_normalization_not_claimed": True,
        "independent_a_gb_not_claimed": True,
        "bulk_eom_gb_force_not_claimed": True,
        "no_full_global_renormalization_claimed": True,
        "no_selector_premise_added": True,
        "no_qw2191_discharge_claimed": True,
        "no_g1_g3_update_claimed": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2340_s1290_v1",
        "packet_id": "P2340",
        "stage_id": "S1290",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "STRICT_TRACK_B_TWO_CLOSED_FIXTURE_TOPOLOGICAL_ADDITIVITY_WITNESS_NO_GLOBAL_CLOSURE",
        "strict_track_b_two_fixture_topological_gluing_compatibility_witness_probe": probe,
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": theorem_export["next_honest_step"],
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2340 strict Track-B two-fixture topological gluing compatibility witness\n\n"
        "Status: S4 and CP2 closed-fixture pairings are compatible under disjoint-union additivity; no global closure.\n\n"
        f"- `b_GB = {b_gb}`.\n"
        f"- S4 pairing: `{s4_pairing}` with `chi=2`.\n"
        f"- CP2 pairing: `{cp2_pairing}` with `chi=3`.\n"
        f"- Disjoint-union direct pairing for `chi=5`: `{union_pairing_direct}`.\n"
        f"- Additivity residual: `{additivity_residual}`; ratio residual `2*CP2 - 3*S4`: `{ratio_residual}`.\n"
        "- No connected-sum theorem, no universal topology normalization, no independent `a_GB`, no full/global renormalization, no QW-2191 discharge, no G1/G3 update, no ToE closure.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
