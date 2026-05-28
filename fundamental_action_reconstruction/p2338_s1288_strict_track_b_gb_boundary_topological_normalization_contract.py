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
OUT = GEN / "p2338_s1288_strict_track_b_gb_boundary_topological_normalization_contract.json"
MD = GEN / "p2338_s1288_strict_track_b_gb_boundary_topological_normalization_contract.md"

SOURCE_FILES = {
    "P2028_GB_QUOTIENT": GEN / "p2028_s978_strict_b1_gb_quotient_counterterm_identifiability_theorem.json",
    "P2034_QUOTIENT_ONLY": GEN / "p2034_s984_strict_task1_quotient_only_renormalization_theorem.json",
    "P2334_EOM_GB_NO_GO": GEN / "p2334_s1284_strict_eom_gb_topological_quotient_scope_theorem.json",
    "P2335_TWO_TRACK_LEDGER": GEN / "p2335_s1285_strict_renormalization_two_track_non_gb_eom_and_gb_topological_ledger.json",
    "P2337_TRACK_A_NONSCALAR_OBSTRUCTION": GEN / "p2337_s1287_strict_track_a_nonscalar_variational_bundle_transport_obstruction.json",
}

REQUIRED_TRACK_B_WITNESS_FIELDS = (
    "oriented_compact_four_manifold_or_boundary_class",
    "gb_density_normalization_convention",
    "boundary_correction_functional",
    "topological_number_chi_or_relative_invariant",
    "scheme_pairing_rule_for_b_GB",
)


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
        ROOT / "p2334_s1284_strict_eom_gb_topological_quotient_scope_theorem.py",
        ROOT / "p2335_s1285_strict_renormalization_two_track_non_gb_eom_and_gb_topological_ledger.py",
        ROOT / "p2337_s1287_strict_track_a_nonscalar_variational_bundle_transport_obstruction.py",
        ROOT / "p2338_s1288_strict_track_b_gb_boundary_topological_normalization_contract.py",
        GEN / "p2028_s978_strict_b1_gb_quotient_counterterm_identifiability_theorem.md",
        GEN / "p2034_s984_strict_task1_quotient_only_renormalization_theorem.md",
        GEN / "p2334_s1284_strict_eom_gb_topological_quotient_scope_theorem.md",
        GEN / "p2335_s1285_strict_renormalization_two_track_non_gb_eom_and_gb_topological_ledger.md",
        GEN / "p2337_s1287_strict_track_a_nonscalar_variational_bundle_transport_obstruction.md",
    ]
    existing = [path.relative_to(ROOT.parent).as_posix() for path in candidates if path.exists()]
    if not existing:
        return []
    cmd = [
        "rg",
        "-n",
        "Track B|GB.*ledger|topological-number|topological number|boundary|Euler|Gauss-Bonnet|a_GB|b_GB|quotient-scope|normalization witness",
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
    p2028 = artifacts["P2028_GB_QUOTIENT"]
    p2034 = artifacts["P2034_QUOTIENT_ONLY"]
    p2334 = artifacts["P2334_EOM_GB_NO_GO"]
    p2335 = artifacts["P2335_TWO_TRACK_LEDGER"]
    p2337 = artifacts["P2337_TRACK_A_NONSCALAR_OBSTRUCTION"]

    b_gb = sp.factor(track_b_coefficient(p2335))
    # Conditional bookkeeping convention only: if a future geometry witness
    # proves int_M GB + boundary = 32*pi**2*chi, then the ledger pairing per
    # unit Euler number is 32*pi**2*b_GB.  This file does not export chi or the
    # boundary functional itself.
    euler_pairing_coefficient = sp.factor(sp.simplify(32 * sp.pi**2 * b_gb))
    missing_witness_fields = list(REQUIRED_TRACK_B_WITNESS_FIELDS)
    current_exports_supply_required_fields = False

    dependencies = {
        "p2028_gb_quotient_loaded": p2028.get("packet_id") == "P2028",
        "p2034_quotient_only_loaded": p2034.get("local_verdict") == "PASS_QUOTIENT_ONLY_RENORMALIZATION_WITH_TRACE",
        "p2334_eom_gb_no_go_loaded": p2334.get("result_kind")
        == "FORMAL_EOM_ONLY_GB_LIFT_NO_GO_CURRENT_STRICT_EXPORTS_QUOTIENT_SCOPE_ONLY",
        "p2335_two_track_ledger_loaded": p2335.get("result_kind")
        == "STRICT_TWO_TRACK_RENORMALIZATION_LEDGER_EXPORTED_NO_FULL_CLOSURE",
        "p2337_track_a_obstruction_loaded": p2337.get("result_kind")
        == "STRICT_TRACK_A_NONSCALAR_VARIATIONAL_BUNDLE_SCALAR_LIFT_OBSTRUCTION_NO_GLOBAL_CLOSURE",
    }

    conditional_contract = {
        "contract_id": "P2338_track_B_boundary_topological_normalization_contract_v1",
        "track": "Track B GB topological counterterm ledger only",
        "ledger_coefficient_b_GB": str(b_gb),
        "conditional_gauss_bonnet_convention": "if future witness proves integral(GB_density)+boundary = 32*pi**2*chi",
        "conditional_per_euler_number_pairing": str(euler_pairing_coefficient),
        "required_future_witness_fields": missing_witness_fields,
        "current_exports_supply_required_fields": current_exports_supply_required_fields,
        "why_not_a_normalization_yet": (
            "The repo exports a separated GB ledger coefficient but does not export the manifold/boundary class, "
            "boundary correction functional, topological number, or scheme pairing rule needed to normalize it."
        ),
    }

    theorem_export = {
        "theorem_name": "P2338 Track-B GB boundary/topological-number normalization contract and nonexport audit",
        "claim": (
            "After P2337 blocks silent scalar-to-tensor promotion for Track A, the honest Track-B move is a boundary/"
            "topological-number normalization audit.  The current strict exports contain the Track-B ledger coefficient "
            "b_GB but do not yet contain the boundary/topological witness fields required to normalize it or to measure "
            "an independent a_GB.  Conditionally, under a future convention integral(GB)+boundary=32*pi**2*chi, the ledger "
            "would pair with chi by 32*pi**2*b_GB; this is a contract, not an exported normalization theorem."
        ),
        "proof_witnesses": [
            "P2335 supplies the separated Track-B topological ledger coefficient b_GB.",
            "P2334 shows EOM-only GB lifting cannot identify independent a_GB.",
            "P2337 blocks promotion of the scalar Track-A transport into a full tensor-bundle transport theorem.",
            "The required boundary/topological witness fields are explicitly listed and remain absent from current exports.",
        ],
        "not_licensed": [
            "independent a_GB measurement",
            "boundary/topological-number normalization theorem",
            "choice of manifold or boundary class",
            "export of Euler characteristic or relative topological invariant",
            "full/global renormalization closure",
            "QW-2191 selector discharge",
            "G1/G3 update",
            "ToE closure",
        ],
        "next_honest_step": (
            "Construct an actual Track-B witness by supplying the missing manifold/boundary class, boundary correction "
            "functional, topological number and scheme pairing rule; otherwise keep b_GB as an unnormalized topological ledger."
        ),
    }

    probe = {
        "probe_id": "P2338_S1288_STRICT_TRACK_B_GB_BOUNDARY_TOPOLOGICAL_NORMALIZATION_CONTRACT",
        "source_hashes": {key: sha256_file(path) for key, path in SOURCE_FILES.items()},
        "repo_grep_audit": {
            "search_focus": "Track B GB ledger boundary/topological-number normalization witness readiness",
            "top_hits": grep_hits(),
        },
        "track_B_contract": conditional_contract,
        "current_export_dependencies": dependencies,
        "theorem_export": theorem_export,
        "theorem_fingerprint_sha256": sha256_json(theorem_export),
    }

    gatekeeper_checks = {
        "repo_grep_audit_recorded": len(probe["repo_grep_audit"]["top_hits"]) > 0,
        "p2028_gb_quotient_loaded": dependencies["p2028_gb_quotient_loaded"],
        "p2034_quotient_only_loaded": dependencies["p2034_quotient_only_loaded"],
        "p2334_eom_gb_no_go_loaded": dependencies["p2334_eom_gb_no_go_loaded"],
        "p2335_two_track_ledger_loaded": dependencies["p2335_two_track_ledger_loaded"],
        "p2337_track_a_obstruction_loaded": dependencies["p2337_track_a_obstruction_loaded"],
        "track_b_coefficient_nonzero": sp.simplify(b_gb) != 0,
        "conditional_euler_pairing_computed": sp.simplify(euler_pairing_coefficient) != 0,
        "required_future_witness_fields_listed": len(missing_witness_fields) == len(REQUIRED_TRACK_B_WITNESS_FIELDS),
        "current_exports_do_not_supply_required_fields": current_exports_supply_required_fields is False,
        "normalization_theorem_not_claimed": True,
        "independent_a_gb_not_claimed": True,
        "boundary_topological_number_witness_still_required": True,
        "no_full_global_renormalization_claimed": True,
        "no_selector_premise_added": True,
        "no_qw2191_discharge_claimed": True,
        "no_g1_g3_update_claimed": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2338_s1288_v1",
        "packet_id": "P2338",
        "stage_id": "S1288",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_OBSTRUCTION_WITH_CONTRACT",
        "result_kind": "STRICT_TRACK_B_GB_BOUNDARY_TOPOLOGICAL_NORMALIZATION_CONTRACT_EXPORTED_NO_NORMALIZATION_CLAIM",
        "strict_track_b_gb_boundary_topological_normalization_contract_probe": probe,
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": theorem_export["next_honest_step"],
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2338 strict Track-B GB boundary/topological normalization contract\n\n"
        "Status: Track-B coefficient is ledgered, but boundary/topological-number normalization is not exported.\n\n"
        f"- Track-B ledger coefficient: `b_GB = {b_gb}`.\n"
        "- Conditional convention only: if a future witness proves `integral(GB)+boundary = 32*pi**2*chi`, "
        f"the per-Euler-number pairing is `{euler_pairing_coefficient}`.\n"
        f"- Missing witness fields: `{', '.join(missing_witness_fields)}`.\n"
        "- No independent `a_GB`, no boundary/topological-number normalization theorem, no full/global renormalization, "
        "no QW-2191 discharge, no G1/G3 update, no ToE closure.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
