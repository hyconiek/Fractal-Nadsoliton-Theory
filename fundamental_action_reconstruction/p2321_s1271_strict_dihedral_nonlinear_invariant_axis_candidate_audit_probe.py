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

OUT = GEN / "p2321_s1271_strict_dihedral_nonlinear_invariant_axis_candidate_audit_probe.json"
MD = GEN / "p2321_s1271_strict_dihedral_nonlinear_invariant_axis_candidate_audit_probe.md"

N = 12
TOL = 1e-9

SOURCE_FILES = {
    "P2320_D12_COMMUTANT_NO_GO": GEN / "p2320_s1270_strict_dihedral_commutant_orientation_no_go_probe.json",
    "P2319_D12_RESPONSE_NO_GO": GEN / "p2319_s1269_strict_dihedral_orientation_response_no_go_probe.json",
    "P2318_SELECTOR_LANE_TO_MARGIN_OBLIGATION": GEN / "p2318_s1268_strict_selector_lane_to_task3_margin_bridge_obligation_probe.json",
    "P2317_FOURIER_DEGENERACY_LANE_AUDIT": GEN / "p2317_s1267_strict_fourier_degeneracy_existing_lane_audit_probe.json",
    "H29_WAVE_RETARDATION_PROXY_SELECTOR_REDUCTION_AUDIT": ROOT / "H29_WAVE_RETARDATION_PROXY_SELECTOR_REDUCTION_AUDIT.md",
    "P719_LOW_COMPLEXITY_ODD_POLYNOMIAL_AUDIT": GEN / "p719_current_strict_t176_low_complexity_odd_polynomial_two_readout_provider_class_audit_probe.json",
}

GREP_PATTERNS = (
    "nonlinear invariant",
    "invariant ring",
    "higher-degree",
    "polynomial",
    "angular harmonic",
    "orientation axis",
    "D12",
    "QW-2191",
    "provider_lift_per_step",
    "selector closure",
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


def group_maps() -> list[dict[str, Any]]:
    maps: list[dict[str, Any]] = []
    for rotation in range(N):
        for reflected in (False, True):
            sign = -1 if reflected else 1
            image = [((sign * index) + rotation) % N for index in range(N)]
            maps.append({"rotation": rotation, "reflected": reflected, "image": image})
    return maps


D12 = group_maps()


def transform_vector(vector: list[float], image: list[int]) -> list[float]:
    out = [0.0 for _ in range(N)]
    for source, target in enumerate(image):
        out[target] = vector[source]
    return out


def dot(a: list[float], b: list[float]) -> float:
    return sum(x * y for x, y in zip(a, b))


def real_fourier_pair(m: int) -> tuple[list[float], list[float]]:
    c = [math.sqrt(2.0 / N) * math.cos(2.0 * math.pi * m * x / N) for x in range(N)]
    s = [-math.sqrt(2.0 / N) * math.sin(2.0 * math.pi * m * x / N) for x in range(N)]
    return c, s


def pair_action_matrix(m: int, image: list[int]) -> list[list[float]]:
    c, s = real_fourier_pair(m)
    tc = transform_vector(c, image)
    ts = transform_vector(s, image)
    return [[dot(c, tc), dot(c, ts)], [dot(s, tc), dot(s, ts)]]


def apply_pair_action(matrix: list[list[float]], x: float, y: float) -> tuple[float, float]:
    return (matrix[0][0] * x + matrix[0][1] * y, matrix[1][0] * x + matrix[1][1] * y)


def real_z_power_coefficients(q: int) -> list[dict[str, Any]]:
    coeffs: list[dict[str, Any]] = []
    for k in range(0, q + 1, 2):
        coeff = math.comb(q, k) * ((-1) ** (k // 2))
        coeffs.append({"x_power": q - k, "y_power": k, "coefficient": coeff})
    return coeffs


def eval_real_z_power(coeffs: list[dict[str, Any]], x: float, y: float) -> float:
    return sum(float(row["coefficient"]) * (x ** int(row["x_power"])) * (y ** int(row["y_power"])) for row in coeffs)


def polynomial_string(coeffs: list[dict[str, Any]]) -> str:
    terms = []
    for row in coeffs:
        coeff = int(row["coefficient"])
        x_power = int(row["x_power"])
        y_power = int(row["y_power"])
        monomial = []
        if x_power:
            monomial.append("x" if x_power == 1 else f"x^{x_power}")
        if y_power:
            monomial.append("y" if y_power == 1 else f"y^{y_power}")
        body = "*".join(monomial) if monomial else "1"
        terms.append(f"{coeff:+d}*{body}")
    return " ".join(terms).lstrip("+")


def invariant_residual_for_pair(m: int, q: int, coeffs: list[dict[str, Any]]) -> float:
    sample_points = []
    for radius in (0.35, 0.8, 1.0, 1.4):
        for idx in range(24):
            theta = 2.0 * math.pi * idx / 24.0 + 0.017 * (m + 1)
            sample_points.append((radius * math.cos(theta), radius * math.sin(theta)))
    max_residual = 0.0
    for row in D12:
        action = pair_action_matrix(m, row["image"])
        for x, y in sample_points:
            tx, ty = apply_pair_action(action, x, y)
            residual = abs(eval_real_z_power(coeffs, tx, ty) - eval_real_z_power(coeffs, x, y))
            max_residual = max(max_residual, residual)
    return max_residual


def angular_variation(q: int) -> dict[str, float]:
    values = [math.cos(q * 2.0 * math.pi * idx / 720.0) for idx in range(720)]
    return {"min": min(values), "max": max(values), "spread": max(values) - min(values)}


def pair_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for m in range(1, N // 2):
        q = N // math.gcd(N, m)
        coeffs = real_z_power_coefficients(q)
        variation = angular_variation(q)
        rows.append({
            "pair": f"pair{m}",
            "m": m,
            "partner_modes": [m, N - m],
            "rotation_order_q": q,
            "lowest_d12_invariant_axis_degree": q,
            "polynomial_id": f"Re((x+i*y)^{q})",
            "polynomial_expansion": polynomial_string(coeffs),
            "coefficients": coeffs,
            "d12_sample_invariance_residual": invariant_residual_for_pair(m, q, coeffs),
            "unit_circle_angular_form": f"cos({q}*theta)",
            "unit_circle_angular_min": variation["min"],
            "unit_circle_angular_max": variation["max"],
            "unit_circle_angular_spread": variation["spread"],
            "directed_extremal_ray_count": 2 * q,
            "maxima_count": q,
            "minima_count": q,
            "breaks_continuous_o2_to_finite_axis_set": True,
            "is_linear_or_quadratic": q <= 2,
            "supplies_signed_policy_orientation": False,
            "why_not_signed_policy_orientation": "D12-invariant harmonic fixes only a finite axis/extrema set after a nonlinear degree-q premise; coefficient sign, preferred extremum branch, response functional, and replay semantics remain unexported.",
        })
    return rows


def corpus_hits() -> list[dict[str, Any]]:
    paths = sorted(
        set(SOURCE_FILES.values())
        | set(ROOT.glob("H*_*.md"))
        | set(GEN.glob("p7*_current_strict*polynomial*json"))
        | set(GEN.glob("p231*_s12*_strict*dihedral*json"))
        | set(GEN.glob("p232*_s12*_strict*dihedral*json"))
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
    p2320 = load_json(SOURCE_FILES["P2320_D12_COMMUTANT_NO_GO"])
    p2319 = load_json(SOURCE_FILES["P2319_D12_RESPONSE_NO_GO"])
    p2318 = load_json(SOURCE_FILES["P2318_SELECTOR_LANE_TO_MARGIN_OBLIGATION"])
    p2317 = load_json(SOURCE_FILES["P2317_FOURIER_DEGENERACY_LANE_AUDIT"])

    rows = pair_rows()
    p2318_probe = p2318.get("strict_selector_lane_to_task3_margin_bridge_obligation_probe", {}) or {}
    required_lift = float((p2318_probe.get("bridge_obligation_verdict", {}) or {}).get("required_lift_per_step", 0.0))
    missing_fields = list((p2318_probe.get("bridge_obligation_verdict", {}) or {}).get("missing_required_bridge_fields", []))

    nonlinear_axis_certificate = {
        "pair_count": len(rows),
        "all_pairs_have_nonlinear_d12_axis_harmonic": all(row["lowest_d12_invariant_axis_degree"] >= 3 for row in rows),
        "degree_sequence_by_pair": {row["pair"]: row["lowest_d12_invariant_axis_degree"] for row in rows},
        "max_d12_sample_invariance_residual": max(row["d12_sample_invariance_residual"] for row in rows),
        "min_angular_spread": min(row["unit_circle_angular_spread"] for row in rows),
        "all_axis_harmonics_nonradial": all(row["unit_circle_angular_spread"] > 1.9 for row in rows),
        "all_axis_harmonics_outside_p2319_linear_quadratic_scope": all(not row["is_linear_or_quadratic"] for row in rows),
        "all_axis_harmonics_fail_signed_policy_orientation": all(not row["supplies_signed_policy_orientation"] for row in rows),
        "honest_upgrade_over_p2319_p2320": "P2319/P2320 no-go remains true for linear/quadratic and full linear-commutant classes, but D12-invariant nonlinear harmonics Re(z^q) do exist and can cut continuous O(2) to finite axis sets.  This is axis-candidate evidence, not a signed provider-lift bridge.",
    }

    bridge_obligation_update = {
        "required_lift_per_step": required_lift,
        "p2318_missing_fields_loaded": missing_fields,
        "nonlinear_axis_candidates_exported": True,
        "nonlinear_axis_candidates_fill_any_missing_p2318_field": False,
        "fields_still_unfilled_after_nonlinear_axis_candidates": missing_fields,
        "g1_g3_update_allowed": False,
        "reason": "The new harmonics supply a D12-invariant finite-axis object but no signed scalar lift, no response map, no preferred sign/branch, no P2281 replay clock, and no admissibility theorem without selector premise.",
    }

    theorem_export = {
        "theorem_name": "P2321 D12 nonlinear invariant axis-candidate audit without Task-3 bridge",
        "formal_statement": (
            "For each real Fourier pair m on strict Z12, the lowest-degree nonradial D12-invariant polynomial on the pair plane is Re((x+i y)^q) with q=12/gcd(12,m). "
            "These degrees are pair1=12, pair2=6, pair3=4, pair4=3, pair5=12.  The harmonics are D12-invariant and nonradial, so they are genuine nonlinear finite-axis candidates beyond the P2319/P2320 linear/quadratic no-go scope.  "
            "However, they do not define a signed policy orientation or provider_lift_per_step bridge: coefficient sign, preferred extremum branch, response functional, replay semantics, and selector-free admissibility remain unexported."
        ),
        "proof_bits": {
            "degree_sequence_by_pair": nonlinear_axis_certificate["degree_sequence_by_pair"],
            "max_d12_sample_invariance_residual": nonlinear_axis_certificate["max_d12_sample_invariance_residual"],
            "min_angular_spread": nonlinear_axis_certificate["min_angular_spread"],
            "all_axis_harmonics_outside_p2319_linear_quadratic_scope": nonlinear_axis_certificate["all_axis_harmonics_outside_p2319_linear_quadratic_scope"],
            "p2318_missing_field_count": len(missing_fields),
            "required_lift_per_step_not_supplied": required_lift,
        },
        "scope_limits": [
            "constructive nonlinear axis-candidate audit, not Task-3 closure",
            "does not choose a coefficient sign or preferred extremum branch",
            "does not export margin response functional or replay semantics",
            "does not discharge QW-2191",
            "does not update G1/G3 or close ToE",
        ],
        "nonpromotion_rule": "Nonlinear D12 finite-axis harmonics may be tracked as candidate selector ingredients, but they must not be promoted to strict signed provider_lift_per_step without an exported sign/response/replay bridge.",
    }
    theorem_fingerprint = sha256_json(theorem_export)

    probe = {
        "probe_id": "P2321_S1271_STRICT_DIHEDRAL_NONLINEAR_INVARIANT_AXIS_CANDIDATE_AUDIT",
        "source_packets": {key: path.relative_to(REPO).as_posix() for key, path in SOURCE_FILES.items()},
        "source_hashes": {f"{key}_sha256": sha256_file(path) for key, path in SOURCE_FILES.items()},
        "far_nonlinear_invariant_grep_audit": {
            "patterns": list(GREP_PATTERNS),
            "hit_count": len(corpus_hits()),
            "top_hits": corpus_hits()[:30],
        },
        "nonlinear_invariant_class_contract": {
            "class_id": "PAIR_PLANE_LOWEST_DEGREE_D12_NONRADIAL_POLYNOMIAL_INVARIANTS",
            "allowed_inputs": ["strict Z12 pair-plane representation", "D12 group action", "polynomial invariant Re((x+i*y)^q)"],
            "excluded_inputs": ["coefficient sign as physical selector", "preferred extremum branch", "Task-3 response functional", "P2281 replay semantics", "legacy kernel role transfer"],
        },
        "pair_nonlinear_axis_rows": rows,
        "nonlinear_axis_certificate": nonlinear_axis_certificate,
        "bridge_obligation_update": bridge_obligation_update,
        "existing_blocker_context": {
            "p2320_result_kind": p2320.get("result_kind"),
            "p2319_result_kind": p2319.get("result_kind"),
            "p2318_result_kind": p2318.get("result_kind"),
            "p2317_result_kind": p2317.get("result_kind"),
        },
        "theorem_export": theorem_export,
        "theorem_fingerprint_sha256": theorem_fingerprint,
    }

    gatekeeper_checks = {
        "far_nonlinear_invariant_grep_hits_found": len(probe["far_nonlinear_invariant_grep_audit"]["top_hits"]) >= 5,
        "p2320_loaded": p2320.get("packet_id") == "P2320",
        "p2318_loaded": p2318.get("packet_id") == "P2318",
        "pair_count_is_5": len(rows) == 5,
        "degree_sequence_matches_d12_pair_orders": nonlinear_axis_certificate["degree_sequence_by_pair"] == {"pair1": 12, "pair2": 6, "pair3": 4, "pair4": 3, "pair5": 12},
        "d12_invariance_residual_small": nonlinear_axis_certificate["max_d12_sample_invariance_residual"] < TOL,
        "all_axis_harmonics_nonradial": nonlinear_axis_certificate["all_axis_harmonics_nonradial"],
        "outside_p2319_linear_quadratic_scope": nonlinear_axis_certificate["all_axis_harmonics_outside_p2319_linear_quadratic_scope"],
        "nonlinear_axis_candidates_not_promoted_to_bridge": not bridge_obligation_update["nonlinear_axis_candidates_fill_any_missing_p2318_field"],
        "p2318_bridge_fields_still_unfilled": len(bridge_obligation_update["fields_still_unfilled_after_nonlinear_axis_candidates"]) == 6,
        "strict_g1_g3_not_updated": not bridge_obligation_update["g1_g3_update_allowed"],
        "no_selector_closure_claimed": True,
        "no_qw2191_discharge_claimed": True,
        "no_legacy_kernel_role_transfer": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2321_s1271_v1",
        "packet_id": "P2321",
        "stage_id": "S1271",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_NONLINEAR_D12_AXIS_CANDIDATES_WITHOUT_SIGNED_TASK3_BRIDGE",
        "result_kind": "STRICT_DIHEDRAL_NONLINEAR_INVARIANT_AXIS_CANDIDATE_AUDIT_NO_G1_G3_UPDATE",
        "strict_dihedral_nonlinear_invariant_axis_candidate_audit_probe": probe,
        "recommended_next_honest_step": (
            "Either derive a strict sign/branch and response/replay bridge for one nonlinear D12 axis harmonic, or prove that such harmonics remain selector-premise-only under the current P2318 interface."
        ),
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_OBSTRUCTION_WITH_NEW_NONLINEAR_AXIS_CANDIDATE_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2321/S1271 — nonlinear D12 invariant axis-candidate audit\n\n"
        f"- Result kind: `{payload['result_kind']}`.\n"
        f"- Degree sequence: `{nonlinear_axis_certificate['degree_sequence_by_pair']}`.\n"
        f"- Max sampled D12 invariance residual: `{nonlinear_axis_certificate['max_d12_sample_invariance_residual']:.3e}`.\n"
        f"- Min angular spread on unit circle: `{nonlinear_axis_certificate['min_angular_spread']:.3e}`.\n"
        "- Verdict: nonlinear D12 harmonics `Re((x+i*y)^q)` are genuine finite-axis candidates beyond linear/quadratic no-go scope, but they do not supply signed Task-3 provider-lift bridge fields.\n"
        "- Guardrail: no G1/G3 update, no QW-2191 discharge, no selector closure, no ToE closure.\n",
        encoding="utf-8",
    )
    print(json.dumps({"wrote": str(OUT.relative_to(REPO)), "theorem_fingerprint_sha256": theorem_fingerprint}, indent=2))


if __name__ == "__main__":
    main()
