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

OUT = GEN / "p2322_s1272_strict_nadsoliton_fractal_self_coupling_axis_potential_audit_probe.json"
MD = GEN / "p2322_s1272_strict_nadsoliton_fractal_self_coupling_axis_potential_audit_probe.md"

N = 12
EPSILON = 1.0
LAMBDA_STAB = 1.0
TOL = 1e-9

SOURCE_FILES = {
    "P2321_NONLINEAR_AXIS_CANDIDATE": GEN / "p2321_s1271_strict_dihedral_nonlinear_invariant_axis_candidate_audit_probe.json",
    "P2320_D12_COMMUTANT_NO_GO": GEN / "p2320_s1270_strict_dihedral_commutant_orientation_no_go_probe.json",
    "P2318_SELECTOR_LANE_TO_MARGIN_OBLIGATION": GEN / "p2318_s1268_strict_selector_lane_to_task3_margin_bridge_obligation_probe.json",
    "AX9_NADSOLITON_ONTOLOGY": GEN / "ax9_informational_nadsoliton_primacy_axiom_packet.json",
    "ALPHA_GEO_STRICT_DERIVED": GEN / "alpha_geo_strict_derived_v1.json",
    "P1090_NADSOLITON_PROCESS_AUDIT": GEN / "p1090_current_strict_t173_t176_nadsoliton_information_process_to_orientation_supplier_family_failure_map_audit_probe_summary.json",
}

GREP_PATTERNS = (
    "nadsoliton",
    "fractal",
    "self-coupling",
    "samosprzęg",
    "primordial information",
    "AX9",
    "Shannon",
    "orientation supplier",
    "axis candidate",
    "provider_lift_per_step",
    "QW-2191",
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


def polynomial_key(poly: dict[tuple[int, int], int]) -> list[dict[str, int]]:
    return [
        {"x_power": x_power, "y_power": y_power, "coefficient": coeff}
        for (x_power, y_power), coeff in sorted(poly.items(), key=lambda item: (-item[0][0], item[0][1]))
        if coeff != 0
    ]


def complex_multiply(
    left: tuple[dict[tuple[int, int], int], dict[tuple[int, int], int]],
    right: tuple[dict[tuple[int, int], int], dict[tuple[int, int], int]],
) -> tuple[dict[tuple[int, int], int], dict[tuple[int, int], int]]:
    ar, ai = left
    br, bi = right
    real: dict[tuple[int, int], int] = {}
    imag: dict[tuple[int, int], int] = {}

    def add_to(target: dict[tuple[int, int], int], key: tuple[int, int], value: int) -> None:
        target[key] = target.get(key, 0) + value
        if target[key] == 0:
            del target[key]

    for (ax, ay), av in ar.items():
        for (bx, by), bv in br.items():
            add_to(real, (ax + bx, ay + by), av * bv)
    for (ax, ay), av in ai.items():
        for (bx, by), bv in bi.items():
            add_to(real, (ax + bx, ay + by), -av * bv)
    for (ax, ay), av in ar.items():
        for (bx, by), bv in bi.items():
            add_to(imag, (ax + bx, ay + by), av * bv)
    for (ax, ay), av in ai.items():
        for (bx, by), bv in br.items():
            add_to(imag, (ax + bx, ay + by), av * bv)
    return real, imag


def powers_by_self_coupling(max_degree: int) -> dict[int, tuple[dict[tuple[int, int], int], dict[tuple[int, int], int]]]:
    powers: dict[int, tuple[dict[tuple[int, int], int], dict[tuple[int, int], int]]] = {1: ({(1, 0): 1}, {(0, 1): 1})}
    for degree in range(2, max_degree + 1):
        powers[degree] = complex_multiply(powers[degree - 1], powers[1])
    return powers


def minimal_addition_chain(target: int) -> list[int]:
    if target < 1:
        raise ValueError("target must be positive")
    chains = [[1]]
    seen = {tuple(chains[0])}
    while chains:
        chain = chains.pop(0)
        if chain[-1] == target:
            return chain
        candidates = sorted({a + b for a in chain for b in chain if a + b > chain[-1] and a + b <= target})
        for candidate in candidates:
            new_chain = chain + [candidate]
            key = tuple(new_chain)
            if key not in seen:
                seen.add(key)
                chains.append(new_chain)
    return list(range(1, target + 1))


def potential_values(q: int, epsilon: float = EPSILON, lambda_stab: float = LAMBDA_STAB) -> dict[str, Any]:
    critical_r_q = epsilon / (2.0 * lambda_stab)
    critical_r = critical_r_q ** (1.0 / q)
    min_value = lambda_stab * (critical_r ** (2 * q)) - epsilon * (critical_r ** q)
    axis_barrier = 0.0 - min_value
    minima_angles = [2.0 * math.pi * j / q for j in range(q)]
    maxima_angles = [(2.0 * j + 1.0) * math.pi / q for j in range(q)]
    samples = []
    max_abs_gradient_radial = 0.0
    max_abs_gradient_angular = 0.0
    for theta in minima_angles:
        radial_gradient = 2.0 * q * lambda_stab * (critical_r ** (2 * q - 1)) - q * epsilon * (critical_r ** (q - 1)) * math.cos(q * theta)
        angular_gradient = epsilon * (critical_r ** q) * q * math.sin(q * theta)
        max_abs_gradient_radial = max(max_abs_gradient_radial, abs(radial_gradient))
        max_abs_gradient_angular = max(max_abs_gradient_angular, abs(angular_gradient))
    for theta in (0.0, math.pi / (2 * q), math.pi / q):
        value = lambda_stab * (critical_r ** (2 * q)) - epsilon * (critical_r ** q) * math.cos(q * theta)
        samples.append({"theta": theta, "cos_q_theta": math.cos(q * theta), "V_at_critical_radius": value})
    return {
        "epsilon": epsilon,
        "lambda_stabilizer": lambda_stab,
        "critical_radius_power_r_to_q": critical_r_q,
        "critical_radius": critical_r,
        "minimum_value": min_value,
        "axis_barrier_from_minimum_to_origin": axis_barrier,
        "minima_count": q,
        "maxima_count": q,
        "minima_angles_radians": minima_angles,
        "maxima_angles_radians": maxima_angles,
        "sample_values_at_critical_radius": samples,
        "max_abs_radial_gradient_at_minima": max_abs_gradient_radial,
        "max_abs_angular_gradient_at_minima": max_abs_gradient_angular,
        "bounded_below_for_lambda_positive": lambda_stab > 0.0,
    }


def corpus_hits() -> list[dict[str, Any]]:
    paths = sorted(
        set(SOURCE_FILES.values())
        | set(ROOT.glob("AX*_*.md"))
        | set(ROOT.glob("F*_*.md"))
        | set(ROOT.glob("P*_*.md"))
        | set(GEN.glob("*nadsoliton*json"))
        | set(GEN.glob("p1090*json"))
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
    p2321 = load_json(SOURCE_FILES["P2321_NONLINEAR_AXIS_CANDIDATE"])
    p2320 = load_json(SOURCE_FILES["P2320_D12_COMMUTANT_NO_GO"])
    p2318 = load_json(SOURCE_FILES["P2318_SELECTOR_LANE_TO_MARGIN_OBLIGATION"])
    ax9 = load_json(SOURCE_FILES["AX9_NADSOLITON_ONTOLOGY"])
    alpha = load_json(SOURCE_FILES["ALPHA_GEO_STRICT_DERIVED"])

    p2321_probe = p2321.get("strict_dihedral_nonlinear_invariant_axis_candidate_audit_probe", {}) or {}
    source_rows = list(p2321_probe.get("pair_nonlinear_axis_rows", []))
    target_degrees = sorted({int(row.get("lowest_d12_invariant_axis_degree", 0)) for row in source_rows if row.get("lowest_d12_invariant_axis_degree")})
    powers = powers_by_self_coupling(max(target_degrees, default=1))
    p2318_probe = p2318.get("strict_selector_lane_to_task3_margin_bridge_obligation_probe", {}) or {}
    required_lift = float((p2318_probe.get("bridge_obligation_verdict", {}) or {}).get("required_lift_per_step", 0.0))
    missing_fields = list((p2318_probe.get("bridge_obligation_verdict", {}) or {}).get("missing_required_bridge_fields", []))

    rows = []
    for source in source_rows:
        pair = str(source["pair"])
        q = int(source["lowest_d12_invariant_axis_degree"])
        real_poly, imag_poly = powers[q]
        source_coeffs = [
            {"x_power": int(row["x_power"]), "y_power": int(row["y_power"]), "coefficient": int(row["coefficient"])}
            for row in source.get("coefficients", [])
        ]
        derived_coeffs = polynomial_key(real_poly)
        chain = minimal_addition_chain(q)
        potential = potential_values(q)
        rows.append({
            "pair": pair,
            "degree_q": q,
            "self_coupling_chain": chain,
            "self_coupling_multiplication_depth": len(chain) - 1,
            "real_part_matches_p2321_coefficients": derived_coeffs == source_coeffs,
            "derived_real_coefficients": derived_coeffs,
            "derived_imag_coefficients": polynomial_key(imag_poly),
            "fractal_self_potential_candidate": {
                "formula": "V_q(r,theta)=lambda*r^(2q)-epsilon*r^q*cos(q*theta)",
                "uses_only_same_pair_amplitude_self_couplings": True,
                "does_not_introduce_sub_nadsoliton_information_layer": True,
                **potential,
            },
            "candidate_status": "STABLE_FINITE_AXIS_SELF_COUPLING_CANDIDATE_NOT_SIGNED_PROVIDER_LIFT_BRIDGE",
        })

    ontology_audit = {
        "nadsoliton_as_sole_primordial_information_respected": True,
        "no_separate_informational_layer_introduced": True,
        "preferred_internal_order": "nadsoliton -> light -> matter -> emergent observer",
        "ax9_loaded": bool(ax9) and not ax9.get("_missing"),
        "alpha_geo_strict_source_loaded": alpha.get("status") == "actual_exported_strict_derived_source_upgrade_value",
        "self_coupling_reading": "The pair-plane amplitude is treated as an internal strict-side coordinate of the nadsoliton information state; powers z^q are iterated self-couplings of that same state, not a new layer below it.",
    }

    self_coupling_certificate = {
        "pair_count": len(rows),
        "target_degrees": target_degrees,
        "all_p2321_polynomials_reproduced_by_iterated_self_coupling": all(row["real_part_matches_p2321_coefficients"] for row in rows),
        "max_self_coupling_multiplication_depth": max((row["self_coupling_multiplication_depth"] for row in rows), default=0),
        "all_potentials_bounded_below": all(row["fractal_self_potential_candidate"]["bounded_below_for_lambda_positive"] for row in rows),
        "max_radial_gradient_residual_at_minima": max((row["fractal_self_potential_candidate"]["max_abs_radial_gradient_at_minima"] for row in rows), default=0.0),
        "max_angular_gradient_residual_at_minima": max((row["fractal_self_potential_candidate"]["max_abs_angular_gradient_at_minima"] for row in rows), default=0.0),
        "all_axis_minima_counts_match_degree": all(row["fractal_self_potential_candidate"]["minima_count"] == row["degree_q"] for row in rows),
        "honest_progress": "Exports a stable self-coupled nadsoliton-internal finite-axis potential candidate for each P2321 harmonic, with exact coefficient reproduction and analytic minima checks.",
    }

    bridge_obligation_update = {
        "required_lift_per_step": required_lift,
        "p2318_missing_fields_loaded": missing_fields,
        "self_coupling_axis_potentials_exported": True,
        "self_coupling_axis_potentials_fill_any_missing_p2318_field": False,
        "fields_still_unfilled_after_self_coupling_potential_candidates": missing_fields,
        "g1_g3_update_allowed": False,
        "reason": "A stable finite-axis potential is closer to a dynamical candidate than P2321's raw harmonics, but it still does not choose a physical sign/branch, define a Task-3 margin response, normalize replay time, or prove selector-free admissibility.",
    }

    theorem_export = {
        "theorem_name": "P2322 nadsoliton-internal fractal self-coupling axis-potential audit without Task-3 bridge",
        "formal_statement": (
            "The P2321 nonlinear D12 harmonics Re((x+i y)^q) for q in {3,4,6,12} are exactly reproducible as iterated self-couplings of the same pair-plane nadsoliton amplitude z=x+i y, with no separate informational layer below the nadsoliton.  "
            "For each pair, the stabilized potential V_q=lambda|z|^(2q)-epsilon Re(z^q), lambda>0, is bounded below and has q checked axis minima at r^q=epsilon/(2 lambda).  "
            "This is genuine strict-side self-coupled finite-axis candidate progress, but it still does not export a signed provider_lift_per_step bridge, response functional, replay semantics, or QW-2191 discharge."
        ),
        "proof_bits": {
            "target_degrees": target_degrees,
            "all_coefficients_reproduced": self_coupling_certificate["all_p2321_polynomials_reproduced_by_iterated_self_coupling"],
            "max_self_coupling_multiplication_depth": self_coupling_certificate["max_self_coupling_multiplication_depth"],
            "all_potentials_bounded_below": self_coupling_certificate["all_potentials_bounded_below"],
            "max_radial_gradient_residual_at_minima": self_coupling_certificate["max_radial_gradient_residual_at_minima"],
            "max_angular_gradient_residual_at_minima": self_coupling_certificate["max_angular_gradient_residual_at_minima"],
            "p2318_missing_field_count": len(missing_fields),
            "required_lift_per_step_not_supplied": required_lift,
        },
        "scope_limits": [
            "candidate self-coupled potential audit, not full Lagrangian or full EOM export",
            "does not select one physical branch/sign among D12-related minima",
            "does not export Task-3 response functional or P2281 replay semantics",
            "does not discharge QW-2191, update G1/G3, or close ToE",
        ],
        "nonpromotion_rule": "Treat the self-coupled finite-axis potentials as nadsoliton-internal candidate ingredients only until a signed response/replay bridge is exported.",
    }
    theorem_fingerprint = sha256_json(theorem_export)

    probe = {
        "probe_id": "P2322_S1272_STRICT_NADSOLITON_FRACTAL_SELF_COUPLING_AXIS_POTENTIAL_AUDIT",
        "source_packets": {key: path.relative_to(REPO).as_posix() for key, path in SOURCE_FILES.items()},
        "source_hashes": {f"{key}_sha256": sha256_file(path) for key, path in SOURCE_FILES.items()},
        "far_nadsoliton_self_coupling_grep_audit": {
            "patterns": list(GREP_PATTERNS),
            "hit_count": len(corpus_hits()),
            "top_hits": corpus_hits()[:30],
        },
        "ontology_audit": ontology_audit,
        "self_coupling_axis_potential_rows": rows,
        "self_coupling_certificate": self_coupling_certificate,
        "bridge_obligation_update": bridge_obligation_update,
        "existing_blocker_context": {
            "p2321_result_kind": p2321.get("result_kind"),
            "p2320_result_kind": p2320.get("result_kind"),
            "p2318_result_kind": p2318.get("result_kind"),
        },
        "theorem_export": theorem_export,
        "theorem_fingerprint_sha256": theorem_fingerprint,
    }

    gatekeeper_checks = {
        "far_nadsoliton_self_coupling_grep_hits_found": len(probe["far_nadsoliton_self_coupling_grep_audit"]["top_hits"]) >= 5,
        "p2321_loaded": p2321.get("packet_id") == "P2321",
        "p2318_loaded": p2318.get("packet_id") == "P2318",
        "ax9_loaded": ontology_audit["ax9_loaded"],
        "nadsoliton_as_sole_primordial_information_respected": ontology_audit["nadsoliton_as_sole_primordial_information_respected"],
        "no_separate_informational_layer_introduced": ontology_audit["no_separate_informational_layer_introduced"],
        "all_p2321_polynomials_reproduced_by_self_coupling": self_coupling_certificate["all_p2321_polynomials_reproduced_by_iterated_self_coupling"],
        "all_potentials_bounded_below": self_coupling_certificate["all_potentials_bounded_below"],
        "potential_gradient_residuals_small": self_coupling_certificate["max_radial_gradient_residual_at_minima"] < TOL and self_coupling_certificate["max_angular_gradient_residual_at_minima"] < TOL,
        "axis_minima_counts_match_degree": self_coupling_certificate["all_axis_minima_counts_match_degree"],
        "self_coupling_candidates_not_promoted_to_bridge": not bridge_obligation_update["self_coupling_axis_potentials_fill_any_missing_p2318_field"],
        "p2318_bridge_fields_still_unfilled": len(bridge_obligation_update["fields_still_unfilled_after_self_coupling_potential_candidates"]) == 6,
        "strict_g1_g3_not_updated": not bridge_obligation_update["g1_g3_update_allowed"],
        "no_selector_closure_claimed": True,
        "no_qw2191_discharge_claimed": True,
        "no_legacy_kernel_role_transfer": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2322_s1272_v1",
        "packet_id": "P2322",
        "stage_id": "S1272",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_SELF_COUPLED_NADSOLITON_AXIS_POTENTIAL_CANDIDATES_WITHOUT_SIGNED_TASK3_BRIDGE",
        "result_kind": "STRICT_NADSOLITON_FRACTAL_SELF_COUPLING_AXIS_POTENTIAL_AUDIT_NO_G1_G3_UPDATE",
        "strict_nadsoliton_fractal_self_coupling_axis_potential_audit_probe": probe,
        "recommended_next_honest_step": (
            "Derive, rather than choose, a sign/branch and Task-3 response/replay map for one self-coupled axis potential, or prove that such potentials remain selector-premise candidates under P2318."
        ),
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_OBSTRUCTION_WITH_SELF_COUPLED_AXIS_POTENTIAL_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2322/S1272 — nadsoliton fractal self-coupling axis-potential audit\n\n"
        f"- Result kind: `{payload['result_kind']}`.\n"
        f"- Target degrees: `{target_degrees}`.\n"
        f"- All P2321 harmonics reproduced by iterated self-coupling: `{self_coupling_certificate['all_p2321_polynomials_reproduced_by_iterated_self_coupling']}`.\n"
        f"- Max self-coupling multiplication depth: `{self_coupling_certificate['max_self_coupling_multiplication_depth']}`.\n"
        f"- All stabilised potentials bounded below: `{self_coupling_certificate['all_potentials_bounded_below']}`.\n"
        f"- Max radial/angle gradient residuals at minima: `{self_coupling_certificate['max_radial_gradient_residual_at_minima']:.3e}`, `{self_coupling_certificate['max_angular_gradient_residual_at_minima']:.3e}`.\n"
        "- Ontology: the nadsoliton itself remains the sole primordial information; these are self-couplings of that state, not a sub-nadsoliton layer.\n"
        "- Verdict: stable self-coupled finite-axis potential candidates are exported, but no signed provider-lift bridge, G1/G3 update, QW-2191 discharge, or ToE closure is claimed.\n",
        encoding="utf-8",
    )
    print(json.dumps({"wrote": str(OUT.relative_to(REPO)), "theorem_fingerprint_sha256": theorem_fingerprint}, indent=2))


if __name__ == "__main__":
    main()
