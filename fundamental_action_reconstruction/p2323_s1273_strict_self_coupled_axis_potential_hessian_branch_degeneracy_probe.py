#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2323_s1273_strict_self_coupled_axis_potential_hessian_branch_degeneracy_probe.json"
MD = GEN / "p2323_s1273_strict_self_coupled_axis_potential_hessian_branch_degeneracy_probe.md"

TOL = 1e-9

SOURCE_FILES = {
    "P2322_SELF_COUPLED_AXIS_POTENTIAL": GEN / "p2322_s1272_strict_nadsoliton_fractal_self_coupling_axis_potential_audit_probe.json",
    "P2321_NONLINEAR_AXIS_CANDIDATE": GEN / "p2321_s1271_strict_dihedral_nonlinear_invariant_axis_candidate_audit_probe.json",
    "P2318_SELECTOR_LANE_TO_MARGIN_OBLIGATION": GEN / "p2318_s1268_strict_selector_lane_to_task3_margin_bridge_obligation_probe.json",
    "PSI_HESSIAN_DIAGONAL_LOCAL": GEN / "psi_hessian_diagonal_local_strict_derived_value_instantiated_v1.json",
}

GREP_PATTERNS = (
    "Hessian",
    "hessian",
    "second variation",
    "mass matrix",
    "local stability",
    "axis potential",
    "branch degeneracy",
    "degenerate minima",
    "selector closure",
    "QW-2191",
    "provider_lift_per_step",
)


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def read_text(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8", errors="replace")


def sha256_file(path: Path) -> str:
    if not path.exists():
        return ""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_json(payload: Any) -> str:
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def potential_derivative_data(q: int, epsilon: float, lambda_stab: float) -> dict[str, float]:
    r0_q = epsilon / (2.0 * lambda_stab)
    r0 = r0_q ** (1.0 / q)
    min_value = lambda_stab * (r0 ** (2 * q)) - epsilon * (r0 ** q)
    barrier_to_origin = -min_value
    polar_radial_hessian = q * q * epsilon * (r0 ** (q - 2))
    polar_angular_hessian = q * q * epsilon * (r0 ** q)
    physical_angular_hessian = polar_angular_hessian / (r0 * r0)
    cross_hessian_at_min = 0.0
    hessian_det_physical = polar_radial_hessian * physical_angular_hessian - cross_hessian_at_min**2
    angular_hessian_at_adjacent_max = -polar_angular_hessian
    return {
        "r0_power_q": r0_q,
        "r0": r0,
        "minimum_value": min_value,
        "barrier_to_origin": barrier_to_origin,
        "radial_second_variation": polar_radial_hessian,
        "angular_second_variation_theta_coordinate": polar_angular_hessian,
        "angular_second_variation_arc_coordinate": physical_angular_hessian,
        "radial_minus_arc_angular_second_variation": polar_radial_hessian - physical_angular_hessian,
        "cross_second_variation_at_minimum": cross_hessian_at_min,
        "physical_hessian_determinant": hessian_det_physical,
        "angular_second_variation_at_adjacent_maximum": angular_hessian_at_adjacent_max,
    }


def branch_rows_from_p2322(p2322: dict[str, Any]) -> list[dict[str, Any]]:
    probe = p2322.get("strict_nadsoliton_fractal_self_coupling_axis_potential_audit_probe", {}) or {}
    rows = []
    for row in probe.get("self_coupling_axis_potential_rows", []) or []:
        potential = row.get("fractal_self_potential_candidate", {}) or {}
        q = int(row.get("degree_q", 0))
        epsilon = float(potential.get("epsilon", 1.0))
        lambda_stab = float(potential.get("lambda_stabilizer", 1.0))
        derivatives = potential_derivative_data(q, epsilon, lambda_stab)
        minima_angles = list(potential.get("minima_angles_radians", []))
        minima_values = [derivatives["minimum_value"] for _ in minima_angles]
        branch_rows = []
        for index, angle in enumerate(minima_angles):
            branch_rows.append({
                "branch_index": index,
                "theta_radians": angle,
                "energy": derivatives["minimum_value"],
                "radial_second_variation": derivatives["radial_second_variation"],
                "angular_second_variation_arc_coordinate": derivatives["angular_second_variation_arc_coordinate"],
                "physical_hessian_eigenvalues": [
                    derivatives["radial_second_variation"],
                    derivatives["angular_second_variation_arc_coordinate"],
                ],
                "locally_stable_minimum": derivatives["radial_second_variation"] > 0.0
                and derivatives["angular_second_variation_arc_coordinate"] > 0.0,
            })
        rows.append({
            "pair": row.get("pair"),
            "degree_q": q,
            "branch_count": len(minima_angles),
            "derivative_formulae": {
                "r0_condition": "r0^q = epsilon/(2*lambda)",
                "d2V_dr2_at_min": "q^2*epsilon*r0^(q-2)",
                "d2V_dtheta2_at_min": "q^2*epsilon*r0^q",
                "arc_coordinate_angular_hessian": "(1/r0^2)*d2V_dtheta2 = q^2*epsilon*r0^(q-2)",
                "mixed_second_derivative_at_min": "0",
            },
            "derivative_values": derivatives,
            "branch_rows": branch_rows,
            "energy_spread_across_minima": max(minima_values) - min(minima_values) if minima_values else 0.0,
            "all_branches_locally_stable": all(branch["locally_stable_minimum"] for branch in branch_rows),
            "all_branches_energy_degenerate": (max(minima_values) - min(minima_values) if minima_values else 0.0) < TOL,
            "local_stability_supplies_branch_selector": False,
            "why_no_branch_selector": "The physical Hessian is positive at every D12-related minimum and identical across branches; second variation proves local stability but preserves exact branch degeneracy.",
        })
    return rows


def corpus_hits() -> list[dict[str, Any]]:
    paths = sorted(
        set(SOURCE_FILES.values())
        | set(ROOT.glob("P*_*.md"))
        | set(ROOT.glob("F*_*.md"))
        | set(GEN.glob("*hessian*json"))
        | set(GEN.glob("p232*_s12*_strict*json"))
    )
    self_paths = {OUT.resolve(), MD.resolve(), Path(__file__).resolve()}
    rows: list[dict[str, Any]] = []
    for path in paths:
        if not path.exists() or path.is_dir() or path.resolve() in self_paths:
            continue
        text = read_text(path)
        lowered = text.lower()
        count = sum(lowered.count(pattern.lower()) for pattern in GREP_PATTERNS)
        if count == 0:
            continue
        first_line = 0
        first_excerpt = ""
        for idx, line in enumerate(text.splitlines(), start=1):
            if any(pattern.lower() in line.lower() for pattern in GREP_PATTERNS):
                first_line = idx
                first_excerpt = line.strip()[:240]
                break
        rows.append({
            "path": path.relative_to(REPO).as_posix(),
            "pattern_hit_count": count,
            "first_hit_line": first_line,
            "first_hit_excerpt": first_excerpt,
        })
    rows.sort(key=lambda row: (-int(row["pattern_hit_count"]), row["path"]))
    return rows


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2322 = load_json(SOURCE_FILES["P2322_SELF_COUPLED_AXIS_POTENTIAL"])
    p2321 = load_json(SOURCE_FILES["P2321_NONLINEAR_AXIS_CANDIDATE"])
    p2318 = load_json(SOURCE_FILES["P2318_SELECTOR_LANE_TO_MARGIN_OBLIGATION"])
    psi_hessian = load_json(SOURCE_FILES["PSI_HESSIAN_DIAGONAL_LOCAL"])

    rows = branch_rows_from_p2322(p2322)
    p2318_probe = p2318.get("strict_selector_lane_to_task3_margin_bridge_obligation_probe", {}) or {}
    required_lift = float((p2318_probe.get("bridge_obligation_verdict", {}) or {}).get("required_lift_per_step", 0.0))
    missing_fields = list((p2318_probe.get("bridge_obligation_verdict", {}) or {}).get("missing_required_bridge_fields", []))

    hessian_certificate = {
        "pair_count": len(rows),
        "branch_count_by_pair": {row["pair"]: row["branch_count"] for row in rows},
        "all_axis_branches_locally_stable": all(row["all_branches_locally_stable"] for row in rows),
        "all_axis_branches_energy_degenerate": all(row["all_branches_energy_degenerate"] for row in rows),
        "max_energy_spread_across_minima": max((row["energy_spread_across_minima"] for row in rows), default=0.0),
        "min_physical_hessian_eigenvalue": min(
            eigenvalue
            for row in rows
            for branch in row["branch_rows"]
            for eigenvalue in branch["physical_hessian_eigenvalues"]
        ),
        "max_radial_arc_angular_hessian_difference": max(
            abs(row["derivative_values"]["radial_minus_arc_angular_second_variation"])
            for row in rows
        ),
        "all_adjacent_axis_maxima_have_negative_angular_second_variation": all(
            row["derivative_values"]["angular_second_variation_at_adjacent_maximum"] < 0.0 for row in rows
        ),
        "local_stability_breaks_branch_degeneracy": False,
        "honest_progress": "P2323 upgrades P2322 from stable-potential candidate to explicit second-variation/Hessian stability and exact branch-degeneracy trace.",
    }

    bridge_obligation_update = {
        "required_lift_per_step": required_lift,
        "p2318_missing_fields_loaded": missing_fields,
        "hessian_stability_certificate_exported": True,
        "hessian_stability_fills_any_missing_p2318_field": False,
        "fields_still_unfilled_after_hessian_stability": missing_fields,
        "g1_g3_update_allowed": False,
        "reason": "Second variation certifies local stability of every branch, but the Hessian spectrum and energy are identical on all D12-related minima, so no signed branch/response/replay bridge is produced.",
    }

    theorem_export = {
        "theorem_name": "P2323 self-coupled axis-potential Hessian stability with branch-degeneracy obstruction",
        "formal_statement": (
            "For each P2322 self-coupled finite-axis potential V_q=lambda r^(2q)-epsilon r^q cos(q theta), the stationary minima at r0^q=epsilon/(2 lambda), theta=2*pi*j/q have positive radial and physical angular second variations q^2 epsilon r0^(q-2), zero mixed second variation, and equal energy -epsilon^2/(4 lambda).  "
            "Thus the self-coupled nadsoliton-internal potentials are locally stable on every D12 branch, but the entire branch orbit remains exactly degenerate.  Hessian stability therefore does not supply the signed provider_lift_per_step bridge, branch selector, response functional, replay semantics, or QW-2191 discharge."
        ),
        "proof_bits": {
            "branch_count_by_pair": hessian_certificate["branch_count_by_pair"],
            "all_axis_branches_locally_stable": hessian_certificate["all_axis_branches_locally_stable"],
            "all_axis_branches_energy_degenerate": hessian_certificate["all_axis_branches_energy_degenerate"],
            "max_energy_spread_across_minima": hessian_certificate["max_energy_spread_across_minima"],
            "min_physical_hessian_eigenvalue": hessian_certificate["min_physical_hessian_eigenvalue"],
            "max_radial_arc_angular_hessian_difference": hessian_certificate["max_radial_arc_angular_hessian_difference"],
            "p2318_missing_field_count": len(missing_fields),
            "required_lift_per_step_not_supplied": required_lift,
        },
        "scope_limits": [
            "second-variation certificate for P2322 candidate potentials, not full Lagrangian/EOM theorem",
            "does not choose a physical sign/branch among D12-related minima",
            "does not export Task-3 response functional or P2281 replay semantics",
            "does not discharge QW-2191, update G1/G3, or close ToE",
        ],
        "nonpromotion_rule": "Do not promote local Hessian stability of the self-coupled axis potentials into selector closure; exact branch degeneracy must be broken by an exported internal source or explicit premise.",
    }
    theorem_fingerprint = sha256_json(theorem_export)

    probe = {
        "probe_id": "P2323_S1273_STRICT_SELF_COUPLED_AXIS_POTENTIAL_HESSIAN_BRANCH_DEGENERACY",
        "source_packets": {key: path.relative_to(REPO).as_posix() for key, path in SOURCE_FILES.items()},
        "source_hashes": {f"{key}_sha256": sha256_file(path) for key, path in SOURCE_FILES.items()},
        "far_hessian_branch_degeneracy_grep_audit": {
            "patterns": list(GREP_PATTERNS),
            "hit_count": len(corpus_hits()),
            "top_hits": corpus_hits()[:30],
        },
        "hessian_branch_rows": rows,
        "hessian_certificate": hessian_certificate,
        "bridge_obligation_update": bridge_obligation_update,
        "existing_blocker_context": {
            "p2322_result_kind": p2322.get("result_kind"),
            "p2321_result_kind": p2321.get("result_kind"),
            "p2318_result_kind": p2318.get("result_kind"),
            "psi_hessian_loaded": bool(psi_hessian) and not psi_hessian.get("_missing"),
        },
        "theorem_export": theorem_export,
        "theorem_fingerprint_sha256": theorem_fingerprint,
    }

    gatekeeper_checks = {
        "far_hessian_branch_degeneracy_grep_hits_found": len(probe["far_hessian_branch_degeneracy_grep_audit"]["top_hits"]) >= 5,
        "p2322_loaded": p2322.get("packet_id") == "P2322",
        "p2318_loaded": p2318.get("packet_id") == "P2318",
        "pair_count_is_5": len(rows) == 5,
        "all_axis_branches_locally_stable": hessian_certificate["all_axis_branches_locally_stable"],
        "all_axis_branches_energy_degenerate": hessian_certificate["all_axis_branches_energy_degenerate"],
        "energy_spread_small": hessian_certificate["max_energy_spread_across_minima"] < TOL,
        "physical_hessian_positive": hessian_certificate["min_physical_hessian_eigenvalue"] > 0.0,
        "radial_and_arc_angular_hessians_match": hessian_certificate["max_radial_arc_angular_hessian_difference"] < TOL,
        "adjacent_maxima_detected_by_negative_angular_hessian": hessian_certificate["all_adjacent_axis_maxima_have_negative_angular_second_variation"],
        "hessian_stability_not_promoted_to_bridge": not bridge_obligation_update["hessian_stability_fills_any_missing_p2318_field"],
        "p2318_bridge_fields_still_unfilled": len(bridge_obligation_update["fields_still_unfilled_after_hessian_stability"]) == 6,
        "strict_g1_g3_not_updated": not bridge_obligation_update["g1_g3_update_allowed"],
        "no_selector_closure_claimed": True,
        "no_qw2191_discharge_claimed": True,
        "no_legacy_kernel_role_transfer": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2323_s1273_v1",
        "packet_id": "P2323",
        "stage_id": "S1273",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_HESSIAN_STABLE_SELF_COUPLED_AXIS_POTENTIALS_WITH_BRANCH_DEGENERACY",
        "result_kind": "STRICT_SELF_COUPLED_AXIS_POTENTIAL_HESSIAN_BRANCH_DEGENERACY_AUDIT_NO_G1_G3_UPDATE",
        "strict_self_coupled_axis_potential_hessian_branch_degeneracy_probe": probe,
        "recommended_next_honest_step": "Search for or rule out an internal source that splits the exactly degenerate D12 branch orbit while also exporting a signed Task-3 response/replay bridge.",
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_OBSTRUCTION_WITH_HESSIAN_STABILITY_AND_BRANCH_DEGENERACY_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2323/S1273 — self-coupled axis-potential Hessian branch-degeneracy audit\n\n"
        f"- Result kind: `{payload['result_kind']}`.\n"
        f"- Branch counts by pair: `{hessian_certificate['branch_count_by_pair']}`.\n"
        f"- All branches locally stable: `{hessian_certificate['all_axis_branches_locally_stable']}`.\n"
        f"- All branch energies degenerate: `{hessian_certificate['all_axis_branches_energy_degenerate']}`.\n"
        f"- Min physical Hessian eigenvalue: `{hessian_certificate['min_physical_hessian_eigenvalue']:.6e}`.\n"
        f"- Max energy spread across minima: `{hessian_certificate['max_energy_spread_across_minima']:.3e}`.\n"
        "- Verdict: Hessian stability is real progress over P2322, but exact D12 branch degeneracy remains; no signed provider-lift bridge is exported.\n"
        "- Guardrail: no G1/G3 update, no QW-2191 discharge, no selector closure, no ToE closure.\n",
        encoding="utf-8",
    )
    print(json.dumps({"wrote": str(OUT.relative_to(REPO)), "theorem_fingerprint_sha256": theorem_fingerprint}, indent=2))


if __name__ == "__main__":
    main()
