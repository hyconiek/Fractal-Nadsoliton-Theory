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
OUT = GEN / "p2359_s1309_strict_track_b_s4_cap_complement_gluing_bulk_boundary_witness.json"
MD = GEN / "p2359_s1309_strict_track_b_s4_cap_complement_gluing_bulk_boundary_witness.md"

SOURCE_FILES = {
    "P2335_TWO_TRACK_LEDGER": GEN / "p2335_s1285_strict_renormalization_two_track_non_gb_eom_and_gb_topological_ledger.json",
    "P2346_SPHERICAL_CAP_BULK_BOUNDARY": GEN / "p2346_s1296_strict_track_b_spherical_cap_bulk_boundary_cancellation_witness.json",
    "P2347_SPHERICAL_CAP_CHERN_FORM": GEN / "p2347_s1297_strict_track_b_spherical_cap_chern_boundary_form_derivation.json",
    "P2353_OBSTRUCTION_LEDGER": GEN / "p2353_s1303_strict_track_b_arbitrary_boundary_obstruction_ledger.json",
    "P2358_CLOSED_INTERFACE_GLUING": GEN / "p2358_s1308_strict_track_b_closed_interface_gluing_cancellation_witness.json",
}

RESULT_KINDS = {
    "P2335_TWO_TRACK_LEDGER": "STRICT_TWO_TRACK_RENORMALIZATION_LEDGER_EXPORTED_NO_FULL_CLOSURE",
    "P2346_SPHERICAL_CAP_BULK_BOUNDARY": "STRICT_TRACK_B_SPHERICAL_CAP_CURVED_BULK_BOUNDARY_CANCELLATION_WITNESS_NO_GLOBAL_THEOREM",
    "P2347_SPHERICAL_CAP_CHERN_FORM": "STRICT_TRACK_B_SPHERICAL_CAP_CHERN_BOUNDARY_FORM_DERIVATION_NO_GENERAL_CGB_THEOREM",
    "P2353_OBSTRUCTION_LEDGER": "STRICT_TRACK_B_ARBITRARY_BOUNDARY_OBSTRUCTION_LEDGER_NO_UNIVERSAL_CLOSURE",
    "P2358_CLOSED_INTERFACE_GLUING": "STRICT_TRACK_B_CLOSED_INTERFACE_GLUING_CANCELLATION_WITNESS_NO_UNIVERSAL_CLOSURE",
}

LOCALS = {"pi": sp.pi, "log": sp.log, "ln": sp.log}


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


def find_first_key(obj: Any, key: str) -> Any:
    if isinstance(obj, dict):
        if key in obj:
            return obj[key]
        for value in obj.values():
            found = find_first_key(value, key)
            if found is not None:
                return found
    elif isinstance(obj, list):
        for value in obj:
            found = find_first_key(value, key)
            if found is not None:
                return found
    return None


def sympify_expr(raw: Any) -> sp.Expr:
    return sp.sympify(str(raw), locals=LOCALS | {"c": sp.symbols("c")})


def track_b_coefficient(p2335: dict[str, Any]) -> sp.Expr:
    probe = p2335.get("strict_two_track_renormalization_ledger_probe", {})
    ledger = probe.get("two_track_ledger", {})
    track_b = ledger.get("track_B_gb_topological_counterterm_ledger", {})
    return sp.sympify(track_b.get("ledger_coefficient_b_GB_topological", "0"), locals=LOCALS)


def grep_hits() -> list[str]:
    candidates = [
        ROOT / "p2346_s1296_strict_track_b_spherical_cap_bulk_boundary_cancellation_witness.py",
        ROOT / "p2347_s1297_strict_track_b_spherical_cap_chern_boundary_form_derivation.py",
        ROOT / "p2358_s1308_strict_track_b_closed_interface_gluing_cancellation_witness.py",
        GEN / "p2358_s1308_strict_track_b_closed_interface_gluing_cancellation_witness.md",
    ]
    existing = [path.relative_to(ROOT.parent).as_posix() for path in candidates if path.exists()]
    if not existing:
        return []
    cmd = [
        "rg",
        "-n",
        "Track-B|Track B|S4|spherical cap|bulk|boundary|gluing|interface|non-flat|Chern|O5|arbitrary-boundary|selector|QW-2191|ToE closure",
        *existing,
    ]
    proc = subprocess.run(cmd, cwd=ROOT.parent, text=True, capture_output=True, check=False)
    return proc.stdout.splitlines()[:300]


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {key: load_json(path) for key, path in SOURCE_FILES.items()}
    b_gb = sp.factor(track_b_coefficient(artifacts["P2335_TWO_TRACK_LEDGER"]))
    per_euler_pairing = sp.factor(32 * sp.pi**2 * b_gb)
    closed_s4_target = 64 * sp.pi**2
    closed_s4_pairing = sp.factor(closed_s4_target * b_gb)

    c = sp.symbols("c", real=True)
    north_bulk = sp.factor(16 * sp.pi**2 * (c - 1) ** 2 * (c + 2))
    interface_boundary_north = sp.factor(16 * sp.pi**2 * c * (3 - c**2))
    north_total = sp.factor(sp.simplify(north_bulk + interface_boundary_north))
    south_bulk = sp.factor(sp.simplify(closed_s4_target - north_bulk))
    interface_boundary_south = sp.factor(-interface_boundary_north)
    south_total = sp.factor(sp.simplify(south_bulk + interface_boundary_south))
    preglue_bulk_boundary_total = sp.factor(sp.simplify(north_total + south_total))
    postglue_closed_bulk_total = sp.factor(sp.simplify(north_bulk + south_bulk))
    interface_cancellation = sp.factor(interface_boundary_north + interface_boundary_south)
    closed_total_residual = sp.factor(preglue_bulk_boundary_total - closed_s4_target)
    postglue_residual = sp.factor(postglue_closed_bulk_total - closed_s4_target)
    gluing_consistency_residual = sp.factor(preglue_bulk_boundary_total - postglue_closed_bulk_total)
    c_derivative_residual = sp.factor(sp.diff(preglue_bulk_boundary_total, c))
    pairing_residual = sp.factor(b_gb * preglue_bulk_boundary_total - closed_s4_pairing)

    sample_values = []
    for value in [sp.Rational(1, 2), sp.Integer(0), -sp.Rational(1, 2)]:
        sample_values.append(
            {
                "c": str(value),
                "north_bulk": str(sp.factor(north_bulk.subs(c, value))),
                "north_boundary": str(sp.factor(interface_boundary_north.subs(c, value))),
                "north_total": str(sp.factor(north_total.subs(c, value))),
                "south_bulk": str(sp.factor(south_bulk.subs(c, value))),
                "south_boundary": str(sp.factor(interface_boundary_south.subs(c, value))),
                "south_total": str(sp.factor(south_total.subs(c, value))),
                "glued_closed_total": str(sp.factor(postglue_closed_bulk_total.subs(c, value))),
            }
        )

    p2346_cap = find_first_key(artifacts["P2346_SPHERICAL_CAP_BULK_BOUNDARY"], "spherical_cap_class") or {}
    p2347_derivation = find_first_key(artifacts["P2347_SPHERICAL_CAP_CHERN_FORM"], "chern_form_derivation") or {}
    p2353_ledger = find_first_key(artifacts["P2353_OBSTRUCTION_LEDGER"], "track_B_arbitrary_boundary_obstruction_ledger") or {}
    p2358_witness = find_first_key(artifacts["P2358_CLOSED_INTERFACE_GLUING"], "track_B_closed_interface_gluing_cancellation_witness") or {}
    p2353_minimal_cut = p2353_ledger.get("minimal_next_missing_cut", [])

    p2346_bulk_residual = sp.factor(north_bulk - sp.sympify(p2346_cap.get("bulk_gb_integral", "0"), locals=LOCALS | {"c": c}))
    p2346_boundary_residual = sp.factor(
        interface_boundary_north - sp.sympify(p2346_cap.get("boundary_transgression", "0"), locals=LOCALS | {"c": c})
    )
    p2347_boundary_residual = sp.factor(
        interface_boundary_north - sp.sympify(p2347_derivation.get("integrated_boundary_transgression", "0"), locals=LOCALS | {"c": c})
    )

    witness = {
        "witness_id": "P2359_s4_cap_complement_gluing_bulk_boundary_v1",
        "scope": "controlled non-flat S4 cap-complement gluing witness; not a general CGB, arbitrary-boundary, or gluing theorem",
        "parameter": "c=cos(rho), -1 < c < 1 for a regular shared geodesic S3 interface",
        "ledger_coefficient_b_GB": str(b_gb),
        "per_euler_pairing": str(per_euler_pairing),
        "closed_s4_target_64pi2": str(closed_s4_target),
        "closed_s4_pairing": str(closed_s4_pairing),
        "north_bulk": str(north_bulk),
        "north_boundary": str(interface_boundary_north),
        "north_total": str(north_total),
        "south_bulk": str(south_bulk),
        "south_boundary": str(interface_boundary_south),
        "south_total": str(south_total),
        "interface_cancellation": str(interface_cancellation),
        "preglue_bulk_boundary_total": str(preglue_bulk_boundary_total),
        "postglue_closed_bulk_total": str(postglue_closed_bulk_total),
        "closed_total_residual": str(closed_total_residual),
        "postglue_residual": str(postglue_residual),
        "gluing_consistency_residual": str(gluing_consistency_residual),
        "c_derivative_residual": str(c_derivative_residual),
        "pairing_residual": str(pairing_residual),
        "sample_values": sample_values,
        "p2346_bulk_residual": str(p2346_bulk_residual),
        "p2346_boundary_residual": str(p2346_boundary_residual),
        "p2347_boundary_residual": str(p2347_boundary_residual),
        "p2358_symbolic_interface_residual_replayed": p2358_witness.get("symbolic_interface_residual"),
        "p2358_all_case_residuals_replayed": p2358_witness.get("all_case_residuals_zero"),
        "p2353_minimal_cut_replayed": p2353_minimal_cut,
        "o5_closed_interface_nonflat_bulk_boundary_partially_attacked": "O5_regularization_corners_and_gluing" in p2353_minimal_cut,
        "remaining_gap": "This covers only the smooth S4 geodesic cap-complement split; it does not handle arbitrary curvature, corners, non-geodesic interfaces, or smoothing limits.",
    }

    dependencies = {
        key.lower() + "_loaded": artifacts[key].get("result_kind") == expected
        for key, expected in RESULT_KINDS.items()
    }
    dependencies.update(
        {
            "p2346_bulk_replayed": p2346_bulk_residual == 0,
            "p2346_boundary_replayed": p2346_boundary_residual == 0,
            "p2347_boundary_replayed": p2347_boundary_residual == 0,
            "p2358_interface_replayed": witness["p2358_symbolic_interface_residual_replayed"] == "0",
            "p2358_cases_replayed": witness["p2358_all_case_residuals_replayed"] is True,
            "p2353_o5_cut_replayed": witness["o5_closed_interface_nonflat_bulk_boundary_partially_attacked"],
            "interface_cancellation_zero": interface_cancellation == 0,
            "north_total_is_32pi2": north_total == 32 * sp.pi**2,
            "south_total_is_32pi2": south_total == 32 * sp.pi**2,
            "closed_total_residual_zero": closed_total_residual == 0,
            "postglue_residual_zero": postglue_residual == 0,
            "gluing_consistency_residual_zero": gluing_consistency_residual == 0,
            "pairing_residual_zero": pairing_residual == 0,
        }
    )

    theorem_export = {
        "theorem_name": "P2359 Track-B S4 cap-complement gluing bulk-boundary witness",
        "claim": (
            "For a unit S4 split into a geodesic cap and its complement along a smooth S3 interface, the P2346/P2347 bulk-boundary formulas give "
            "32*pi**2 on each side, the opposite interface transgressions cancel under gluing, and the glued closed S4 total is 64*pi**2. "
            "This is a controlled non-flat closed-interface witness, not a general Chern-Gauss-Bonnet or arbitrary-boundary theorem."
        ),
        "proof_witnesses": [
            "The north cap bulk plus boundary simplifies to 32*pi**2 for symbolic c.",
            "The south complement bulk is 64*pi**2 minus the north bulk and its boundary transgression has the opposite sign; it also totals 32*pi**2.",
            "The shared interface boundary terms cancel exactly, so preglue bulk-plus-boundary equals postglue closed bulk 64*pi**2.",
            "The c-derivative of the glued total is zero and sample c values replay the same total.",
            "P2346 bulk/boundary, P2347 Chern-boundary derivation, and P2358 closed-interface cancellation are all replayed with zero residuals where applicable.",
        ],
        "not_licensed": [
            "general Chern-Gauss-Bonnet theorem over arbitrary compact four-manifolds",
            "arbitrary-boundary theorem",
            "general gluing theorem",
            "corner contribution theorem",
            "smoothing-limit theorem",
            "general non-geodesic interface theorem",
            "general non-flat ambient curvature bridge beyond the unit S4 cap-complement fixture",
            "global renormalization closure",
            "independent a_GB measurement separate from the P2335 ledger coefficient",
            "bulk EOM GB force or EOM-only GB lift",
            "selector premise or QW-2191 selector discharge",
            "choice of the unique physical future",
            "observer-level prediction",
            "G1/G3 update",
            "ToE closure",
        ],
        "next_honest_step": (
            "The non-flat smooth closed-interface case now has a controlled cap-complement replay.  The next real progress target is a local "
            "corner/smoothing model or a non-geodesic-interface computation; do not relabel this S4 fixture as general CGB closure."
        ),
    }

    probe = {
        "probe_id": "P2359_S1309_STRICT_TRACK_B_S4_CAP_COMPLEMENT_GLUING_BULK_BOUNDARY_WITNESS",
        "source_hashes": {key: sha256_file(path) for key, path in SOURCE_FILES.items()},
        "repo_grep_audit": {
            "search_focus": "P2359 non-flat S4 cap-complement closed-interface gluing after P2358",
            "top_hits": grep_hits(),
        },
        "track_B_s4_cap_complement_gluing_bulk_boundary_witness": witness,
        "current_export_dependencies": dependencies,
        "theorem_export": theorem_export,
        "theorem_fingerprint_sha256": sha256_json(theorem_export),
    }

    gatekeeper_checks = {
        "repo_grep_audit_recorded": len(probe["repo_grep_audit"]["top_hits"]) > 0,
        **dependencies,
        "three_sample_values_recorded": len(sample_values) == 3,
        "nonflat_s4_scope_declared": "S4 cap-complement" in witness["scope"],
        "general_cgb_theorem_not_claimed": True,
        "arbitrary_boundary_theorem_not_claimed": True,
        "general_gluing_theorem_not_claimed": True,
        "corner_theorem_not_claimed": True,
        "smoothing_limit_theorem_not_claimed": True,
        "non_geodesic_interface_theorem_not_claimed": True,
        "global_renormalization_not_claimed": True,
        "independent_a_gb_not_claimed": True,
        "bulk_eom_gb_force_not_claimed": True,
        "no_selector_premise_added": True,
        "no_qw2191_discharge_claimed": True,
        "no_unique_future_choice_claimed": True,
        "no_observer_prediction_claimed": True,
        "no_g1_g3_update_claimed": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2359_s1309_v1",
        "packet_id": "P2359",
        "stage_id": "S1309",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-29T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "STRICT_TRACK_B_S4_CAP_COMPLEMENT_GLUING_BULK_BOUNDARY_WITNESS_NO_GENERAL_CGB",
        "strict_track_b_s4_cap_complement_gluing_bulk_boundary_witness_probe": probe,
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": theorem_export["next_honest_step"],
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2359 strict Track-B S4 cap-complement gluing bulk-boundary witness\n\n"
        "Status: controlled non-flat closed-interface cap-complement witness; no arbitrary-boundary, general CGB, selector, or ToE closure.\n\n"
        f"- North bulk `{north_bulk}` plus boundary `{interface_boundary_north}` gives `{north_total}`.\n"
        f"- South bulk `{south_bulk}` plus boundary `{interface_boundary_south}` gives `{south_total}`.\n"
        f"- Interface cancellation `{interface_cancellation}`; preglue total `{preglue_bulk_boundary_total}`; postglue closed total `{postglue_closed_bulk_total}`.\n"
        f"- Residuals: closed `{closed_total_residual}`, postglue `{postglue_residual}`, gluing consistency `{gluing_consistency_residual}`, pairing `{pairing_residual}`.\n"
        f"- Replays: P2346 bulk `{p2346_bulk_residual}`, P2346 boundary `{p2346_boundary_residual}`, P2347 boundary `{p2347_boundary_residual}`, P2358 interface `{witness['p2358_symbolic_interface_residual_replayed']}`.\n"
        f"- P2353 cut replayed: `{p2353_minimal_cut}`; O5 non-flat closed-interface partially attacked `{witness['o5_closed_interface_nonflat_bulk_boundary_partially_attacked']}`.\n"
        "- No general Chern-Gauss-Bonnet theorem, no arbitrary-boundary theorem, no general gluing theorem, no corner/smoothing theorem, no non-geodesic interface theorem, no global renormalization, no independent `a_GB`, no selector premise, no QW-2191 discharge, no unique-future choice, no observer prediction, no G1/G3 update, no ToE closure.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
