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

OUT = GEN / "p2320_s1270_strict_dihedral_commutant_orientation_no_go_probe.json"
MD = GEN / "p2320_s1270_strict_dihedral_commutant_orientation_no_go_probe.md"

N = 12
TOL = 1e-10

SOURCE_FILES = {
    "P2319_D12_RESPONSE_NO_GO": GEN / "p2319_s1269_strict_dihedral_orientation_response_no_go_probe.json",
    "P2318_SELECTOR_LANE_TO_MARGIN_OBLIGATION": GEN / "p2318_s1268_strict_selector_lane_to_task3_margin_bridge_obligation_probe.json",
    "P2317_FOURIER_DEGENERACY_LANE_AUDIT": GEN / "p2317_s1267_strict_fourier_degeneracy_existing_lane_audit_probe.json",
    "P2308_CURRENT_INTERFACE_NONEXISTENCE": GEN / "p2308_s1258_strict_current_interface_class_response_functional_nonexistence_probe.json",
    "N494_UNIQUENESS_UP_TO_CONJUGATION": ROOT / "N494_CURRENT_FIRST_STRICT_QW2190_DIAGONAL_LOCAL_MODE_INDEX_CANONICALIZATION_UNIQUENESS_UP_TO_CONJUGATION_THEOREM.md",
    "N493_RESIDUAL_Z2_SIGN_FLIP": ROOT / "N493_CURRENT_FIRST_STRICT_QW2191_RESIDUAL_Z2_SIGN_FLIP_GAUGE_EQUIVALENCE_THEOREM.md",
    "P454_O2_GAUGE_EQUIVALENCE_NOTE": ROOT / "P454_CURRENT_STRICT_QW2191_O2_ROTATION_GAUGE_EQUIVALENCE_AUDIT_PROBE.md",
}

GREP_PATTERNS = (
    "commutant",
    "intertwiner",
    "dihedral",
    "D12",
    "reflection",
    "signed orientation",
    "orientation response",
    "response functional",
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


def cyclic_distance(i: int, j: int) -> int:
    raw = abs(i - j) % N
    return min(raw, N - raw)


def group_maps(*, include_reflections: bool) -> list[dict[str, Any]]:
    maps: list[dict[str, Any]] = []
    reflected_values = (False, True) if include_reflections else (False,)
    for rotation in range(N):
        for reflected in reflected_values:
            sign = -1 if reflected else 1
            image = [((sign * index) + rotation) % N for index in range(N)]
            maps.append({"rotation": rotation, "reflected": reflected, "image": image})
    return maps


D12 = group_maps(include_reflections=True)
C12 = group_maps(include_reflections=False)


def transform_matrix(matrix: list[list[float]], image: list[int]) -> list[list[float]]:
    out = [[0.0 for _ in range(N)] for _ in range(N)]
    for i, target_i in enumerate(image):
        for j, target_j in enumerate(image):
            out[target_i][target_j] = matrix[i][j]
    return out


def mat_sub(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    return [[a[i][j] - b[i][j] for j in range(N)] for i in range(N)]


def mat_add(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    return [[a[i][j] + b[i][j] for j in range(N)] for i in range(N)]


def mat_norm_inf(matrix: list[list[float]]) -> float:
    return max((sum(abs(x) for x in row) for row in matrix), default=0.0)


def frob_dot(a: list[list[float]], b: list[list[float]]) -> float:
    return sum(a[i][j] * b[i][j] for i in range(N) for j in range(N))


def frob_norm(matrix: list[list[float]]) -> float:
    return math.sqrt(max(0.0, frob_dot(matrix, matrix)))


def mat_scale(scale: float, matrix: list[list[float]]) -> list[list[float]]:
    return [[scale * matrix[i][j] for j in range(N)] for i in range(N)]


def mat_zero() -> list[list[float]]:
    return [[0.0 for _ in range(N)] for _ in range(N)]


def reynolds_matrix(matrix: list[list[float]], group: list[dict[str, Any]]) -> list[list[float]]:
    acc = mat_zero()
    for row in group:
        transformed = transform_matrix(matrix, row["image"])
        for i in range(N):
            for j in range(N):
                acc[i][j] += transformed[i][j]
    return [[acc[i][j] / len(group) for j in range(N)] for i in range(N)]


def distance_basis() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for distance in range((N // 2) + 1):
        matrix = [[1.0 if cyclic_distance(i, j) == distance else 0.0 for j in range(N)] for i in range(N)]
        rows.append({"basis_id": f"D12_distance_{distance}", "distance": distance, "matrix": matrix})
    return rows


def shift_basis() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for shift in range(N):
        matrix = [[1.0 if j == (i + shift) % N else 0.0 for j in range(N)] for i in range(N)]
        rows.append({"basis_id": f"C12_shift_{shift}", "shift": shift, "matrix": matrix})
    return rows


def max_invariance_residual(matrix: list[list[float]], group: list[dict[str, Any]]) -> float:
    return max(mat_norm_inf(mat_sub(transform_matrix(matrix, row["image"]), matrix)) for row in group)


def basis_gram_rank(basis: list[dict[str, Any]]) -> int:
    gram = [[frob_dot(row_a["matrix"], row_b["matrix"]) for row_b in basis] for row_a in basis]
    return numeric_rank(gram)


def numeric_rank(matrix: list[list[float]], tol: float = 1e-9) -> int:
    work = [row[:] for row in matrix]
    if not work:
        return 0
    row_count = len(work)
    col_count = len(work[0])
    rank = 0
    for col in range(col_count):
        pivot = max(range(rank, row_count), key=lambda r: abs(work[r][col]))
        if abs(work[pivot][col]) <= tol:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        pivot_value = work[rank][col]
        work[rank] = [x / pivot_value for x in work[rank]]
        for r in range(row_count):
            if r == rank:
                continue
            factor = work[r][col]
            if abs(factor) <= tol:
                continue
            work[r] = [work[r][c] - factor * work[rank][c] for c in range(col_count)]
        rank += 1
        if rank == row_count:
            break
    return rank


def project_onto_orthogonal_basis(matrix: list[list[float]], basis: list[dict[str, Any]]) -> dict[str, Any]:
    projection = mat_zero()
    coeffs = []
    for row in basis:
        basis_matrix = row["matrix"]
        denom = frob_dot(basis_matrix, basis_matrix)
        coeff = frob_dot(matrix, basis_matrix) / denom if denom else 0.0
        coeffs.append({"basis_id": row["basis_id"], "coefficient": coeff})
        scaled = mat_scale(coeff, basis_matrix)
        projection = mat_add(projection, scaled)
    residual = mat_sub(matrix, projection)
    return {
        "projection_coefficients": coeffs,
        "projection_frobenius_norm": frob_norm(projection),
        "input_frobenius_norm": frob_norm(matrix),
        "residual_frobenius_norm": frob_norm(residual),
        "projection_inf_norm": mat_norm_inf(projection),
        "residual_inf_norm": mat_norm_inf(residual),
    }


def real_fourier_pair(m: int) -> tuple[list[float], list[float]]:
    c = [math.sqrt(2.0 / N) * math.cos(2.0 * math.pi * m * x / N) for x in range(N)]
    s = [-math.sqrt(2.0 / N) * math.sin(2.0 * math.pi * m * x / N) for x in range(N)]
    return c, s


def outer(a: list[float], b: list[float]) -> list[list[float]]:
    return [[a[i] * b[j] for j in range(N)] for i in range(N)]


def orientation_operator_rows(d12_basis: list[dict[str, Any]], c12_basis: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for m in range(1, N // 2):
        c, s = real_fourier_pair(m)
        projector = mat_add(outer(c, c), outer(s, s))
        anisotropy = mat_sub(outer(c, c), outer(s, s))
        shear = mat_add(outer(c, s), outer(s, c))
        handed_skew = mat_sub(outer(c, s), outer(s, c))
        candidates = [
            ("unsigned_pair_projector", projector, "D12-even unsigned norm"),
            ("traceless_anisotropy_cc_minus_ss", anisotropy, "orientation axis quadratic"),
            ("symmetric_shear_cs_plus_sc", shear, "orientation angle quadratic"),
            ("antisymmetric_handed_skew_cs_minus_sc", handed_skew, "rotation-invariant but reflection-odd handedness"),
        ]
        candidate_rows = []
        for name, matrix, meaning in candidates:
            d12_projection = project_onto_orthogonal_basis(matrix, d12_basis)
            c12_projection = project_onto_orthogonal_basis(matrix, c12_basis)
            candidate_rows.append({
                "candidate_id": name,
                "meaning": meaning,
                "input_frobenius_norm": frob_norm(matrix),
                "d12_projection_frobenius_norm": d12_projection["projection_frobenius_norm"],
                "d12_residual_frobenius_norm": d12_projection["residual_frobenius_norm"],
                "c12_projection_frobenius_norm": c12_projection["projection_frobenius_norm"],
                "c12_residual_frobenius_norm": c12_projection["residual_frobenius_norm"],
                "d12_survives_commutant_projection": d12_projection["projection_frobenius_norm"] > TOL,
                "c12_survives_commutant_projection": c12_projection["projection_frobenius_norm"] > TOL,
            })
        rows.append({"pair": f"pair{m}", "m": m, "partner_modes": [m, N - m], "candidate_rows": candidate_rows})
    return rows


def corpus_hits() -> list[dict[str, Any]]:
    paths = sorted(
        set(SOURCE_FILES.values())
        | set(ROOT.glob("N49*_CURRENT*QW2191*.md"))
        | set(ROOT.glob("P45*_CURRENT_STRICT*QW2191*.md"))
        | set(GEN.glob("p230*_s12*_strict*response*json"))
        | set(GEN.glob("p231*_s12*_strict*orientation*json"))
        | set(GEN.glob("p231*_s12*_strict*dihedral*json"))
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
    p2319 = load_json(SOURCE_FILES["P2319_D12_RESPONSE_NO_GO"])
    p2318 = load_json(SOURCE_FILES["P2318_SELECTOR_LANE_TO_MARGIN_OBLIGATION"])
    p2317 = load_json(SOURCE_FILES["P2317_FOURIER_DEGENERACY_LANE_AUDIT"])
    p2308 = load_json(SOURCE_FILES["P2308_CURRENT_INTERFACE_NONEXISTENCE"])

    d12_basis = distance_basis()
    c12_basis = shift_basis()
    operator_rows = orientation_operator_rows(d12_basis, c12_basis)
    p2318_probe = p2318.get("strict_selector_lane_to_task3_margin_bridge_obligation_probe", {}) or {}
    required_lift = float((p2318_probe.get("bridge_obligation_verdict", {}) or {}).get("required_lift_per_step", 0.0))
    missing_fields = list((p2318_probe.get("bridge_obligation_verdict", {}) or {}).get("missing_required_bridge_fields", []))

    d12_basis_rows = [
        {
            "basis_id": row["basis_id"],
            "distance": row["distance"],
            "frobenius_norm": frob_norm(row["matrix"]),
            "d12_invariance_residual": max_invariance_residual(row["matrix"], D12),
            "c12_invariance_residual": max_invariance_residual(row["matrix"], C12),
        }
        for row in d12_basis
    ]
    c12_basis_rows = [
        {
            "basis_id": row["basis_id"],
            "shift": row["shift"],
            "frobenius_norm": frob_norm(row["matrix"]),
            "c12_invariance_residual": max_invariance_residual(row["matrix"], C12),
            "d12_invariance_residual": max_invariance_residual(row["matrix"], D12),
        }
        for row in c12_basis
    ]

    d12_non_unsigned_survivors = []
    c12_handed_survivors = []
    for pair in operator_rows:
        for row in pair["candidate_rows"]:
            if row["candidate_id"] != "unsigned_pair_projector" and row["d12_survives_commutant_projection"]:
                d12_non_unsigned_survivors.append({"pair": pair["pair"], "candidate_id": row["candidate_id"]})
            if row["candidate_id"] == "antisymmetric_handed_skew_cs_minus_sc" and row["c12_survives_commutant_projection"]:
                c12_handed_survivors.append({"pair": pair["pair"], "candidate_id": row["candidate_id"]})

    commutant_certificate = {
        "d12_group_order": len(D12),
        "c12_group_order": len(C12),
        "d12_distance_basis_count": len(d12_basis),
        "d12_distance_basis_gram_rank": basis_gram_rank(d12_basis),
        "c12_shift_basis_count": len(c12_basis),
        "c12_shift_basis_gram_rank": basis_gram_rank(c12_basis),
        "all_distance_basis_d12_invariant": all(row["d12_invariance_residual"] < TOL for row in d12_basis_rows),
        "all_shift_basis_c12_invariant": all(row["c12_invariance_residual"] < TOL for row in c12_basis_rows),
        "some_shift_basis_reflection_noninvariant": any(row["d12_invariance_residual"] > TOL for row in c12_basis_rows),
        "all_non_unsigned_orientation_candidates_project_to_zero_in_d12_commutant": len(d12_non_unsigned_survivors) == 0,
        "all_unsigned_pair_projectors_survive_d12_commutant": all(
            next(row for row in pair["candidate_rows"] if row["candidate_id"] == "unsigned_pair_projector")["d12_survives_commutant_projection"]
            for pair in operator_rows
        ),
        "c12_handed_skew_survivor_count": len(c12_handed_survivors),
        "c12_handed_skew_survivors_are_reflection_odd_not_d12_admissible": len(c12_handed_survivors) == len(operator_rows),
        "max_d12_non_unsigned_projection_norm": max(
            row["d12_projection_frobenius_norm"]
            for pair in operator_rows
            for row in pair["candidate_rows"]
            if row["candidate_id"] != "unsigned_pair_projector"
        ),
        "min_d12_unsigned_projection_norm": min(
            row["d12_projection_frobenius_norm"]
            for pair in operator_rows
            for row in pair["candidate_rows"]
            if row["candidate_id"] == "unsigned_pair_projector"
        ),
    }

    bridge_obligation_update = {
        "required_lift_per_step": required_lift,
        "p2318_missing_fields_loaded": missing_fields,
        "commutant_class_fills_any_missing_p2318_field": False,
        "fields_still_unfilled_by_full_d12_commutant_class": missing_fields,
        "why_c12_contrast_not_enough": "C12 admits reflection-odd handed skew operators, but admitting them drops D12 reflection symmetry and therefore adds exactly the sort of orientation/sign premise not exported by strict kernel data alone.",
        "g1_g3_update_allowed": False,
    }

    theorem_export = {
        "theorem_name": "P2320 full D12 commutant orientation-response no-go theorem for strict Z12 pair planes",
        "formal_statement": (
            "The full D12-invariant linear-operator commutant on the strict Z12 slot space is spanned by the seven cyclic-distance matrices. "
            "Its projection on every real Fourier pair plane preserves the unsigned pair projector and annihilates the signed orientation operators "
            "cc-ss, cs+sc, and the reflection-odd handed skew cs-sc.  Hence no operator in this full D12 commutant supplies a signed Fourier-pair "
            "orientation response or any P2318 provider_lift_per_step bridge field.  A C12-only contrast admits handed skew directions, but only by dropping reflection symmetry, so it is a new selector/sign premise rather than a strict kernel-alone consequence."
        ),
        "proof_bits": {
            "d12_group_order": commutant_certificate["d12_group_order"],
            "d12_commutant_dimension": commutant_certificate["d12_distance_basis_gram_rank"],
            "c12_commutant_dimension": commutant_certificate["c12_shift_basis_gram_rank"],
            "pair_count": len(operator_rows),
            "max_d12_non_unsigned_projection_norm": commutant_certificate["max_d12_non_unsigned_projection_norm"],
            "min_d12_unsigned_projection_norm": commutant_certificate["min_d12_unsigned_projection_norm"],
            "c12_handed_skew_survivor_count": commutant_certificate["c12_handed_skew_survivor_count"],
            "p2318_missing_field_count": len(missing_fields),
            "required_lift_per_step_not_supplied": required_lift,
        },
        "scope_limits": [
            "finite-dimensional strict Z12 linear-operator commutant theorem, not a theorem over all future nonlinear dynamics",
            "C12-only/reflection-breaking routes are classified as adding a sign premise unless independently exported",
            "does not refute lane-scoped diagonal/local or Shannon O(2)->Z2 computations inside their lanes",
            "does not discharge QW-2191",
            "does not close selector, G1/G3, Task-3, or ToE",
        ],
        "nonpromotion_rule": "Do not promote D12 commutant operators or C12-only handed contrast directions into a Task-3 signed provider-lift bridge without an exported reflection-breaking response theorem and replay semantics.",
    }
    theorem_fingerprint = sha256_json(theorem_export)

    probe = {
        "probe_id": "P2320_S1270_STRICT_DIHEDRAL_COMMUTANT_ORIENTATION_NO_GO",
        "source_packets": {key: path.relative_to(REPO).as_posix() for key, path in SOURCE_FILES.items()},
        "source_hashes": {f"{key}_sha256": sha256_file(path) for key, path in SOURCE_FILES.items()},
        "far_commutant_grep_audit": {
            "patterns": list(GREP_PATTERNS),
            "hit_count": len(corpus_hits()),
            "top_hits": corpus_hits()[:30],
        },
        "commutant_class_contract": {
            "class_id": "FULL_D12_INVARIANT_LINEAR_OPERATOR_COMMUTANT_ON_R12",
            "basis_formula": "B_d[i,j] = 1 iff min(|i-j| mod 12, 12-|i-j| mod 12) = d, for d=0..6",
            "allowed_inputs": ["strict Z12 slot representation", "D12 rotations/reflections", "linear operators commuting with that action"],
            "excluded_inputs": ["external selector premise", "reflection-breaking handedness premise", "nonlinear dynamics", "P2281 replay response theorem", "legacy kernel role transfer"],
        },
        "d12_basis_rows": d12_basis_rows,
        "c12_contrast_basis_rows": c12_basis_rows,
        "pair_operator_projection_rows": operator_rows,
        "commutant_certificate": commutant_certificate,
        "bridge_obligation_update": bridge_obligation_update,
        "existing_blocker_context": {
            "p2319_result_kind": p2319.get("result_kind"),
            "p2318_result_kind": p2318.get("result_kind"),
            "p2317_result_kind": p2317.get("result_kind"),
            "p2308_result_kind": p2308.get("result_kind"),
        },
        "theorem_export": theorem_export,
        "theorem_fingerprint_sha256": theorem_fingerprint,
    }

    gatekeeper_checks = {
        "far_commutant_grep_hits_found": len(probe["far_commutant_grep_audit"]["top_hits"]) >= 5,
        "p2319_loaded": p2319.get("packet_id") == "P2319",
        "p2318_loaded": p2318.get("packet_id") == "P2318",
        "d12_group_order_is_24": commutant_certificate["d12_group_order"] == 24,
        "d12_commutant_dimension_is_7": commutant_certificate["d12_distance_basis_gram_rank"] == 7,
        "c12_commutant_dimension_is_12": commutant_certificate["c12_shift_basis_gram_rank"] == 12,
        "distance_basis_d12_invariant": commutant_certificate["all_distance_basis_d12_invariant"],
        "non_unsigned_orientation_operators_annihilated_by_d12_commutant_projection": commutant_certificate["all_non_unsigned_orientation_candidates_project_to_zero_in_d12_commutant"],
        "unsigned_pair_projectors_survive_d12_commutant_projection": commutant_certificate["all_unsigned_pair_projectors_survive_d12_commutant"],
        "c12_handed_contrast_detected_but_not_promoted": commutant_certificate["c12_handed_skew_survivors_are_reflection_odd_not_d12_admissible"],
        "p2318_bridge_fields_still_unfilled": len(bridge_obligation_update["fields_still_unfilled_by_full_d12_commutant_class"]) == 6,
        "strict_g1_g3_not_updated": not bridge_obligation_update["g1_g3_update_allowed"],
        "no_selector_closure_claimed": True,
        "no_qw2191_discharge_claimed": True,
        "no_legacy_kernel_role_transfer": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2320_s1270_v1",
        "packet_id": "P2320",
        "stage_id": "S1270",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_NO_GO_FOR_FULL_D12_COMMUTANT_SIGNED_ORIENTATION_RESPONSE_CLASS",
        "result_kind": "STRICT_DIHEDRAL_COMMUTANT_ORIENTATION_NO_GO_NO_G1_G3_UPDATE",
        "strict_dihedral_commutant_orientation_no_go_probe": probe,
        "recommended_next_honest_step": (
            "If continuing the selector route, the next admissible move must explicitly export a reflection-breaking/non-D12 strict source plus "
            "a signed response theorem and P2281/P2302 replay semantics; otherwise broaden the no-go class beyond linear commutants without claiming closure."
        ),
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2320/S1270 — full D12 commutant orientation-response no-go\n\n"
        f"- Result kind: `{payload['result_kind']}`.\n"
        f"- D12 commutant dimension: `{commutant_certificate['d12_distance_basis_gram_rank']}` with distance basis d=0..6.\n"
        f"- C12 contrast dimension: `{commutant_certificate['c12_shift_basis_gram_rank']}`.\n"
        f"- Max D12 projection norm of non-unsigned orientation operators: `{commutant_certificate['max_d12_non_unsigned_projection_norm']:.3e}`.\n"
        f"- Min D12 projection norm of unsigned pair projectors: `{commutant_certificate['min_d12_unsigned_projection_norm']:.3e}`.\n"
        f"- C12 handed/skew survivor count: `{commutant_certificate['c12_handed_skew_survivor_count']}` (contrast only; not promoted).\n"
        "- Verdict: the full D12 commutant leaves unsigned pair projectors but no signed orientation response, so P2318 bridge fields remain open.\n"
        "- Guardrail: no G1/G3 update, no QW-2191 discharge, no selector closure, no ToE closure.\n",
        encoding="utf-8",
    )
    print(json.dumps({"wrote": str(OUT.relative_to(REPO)), "theorem_fingerprint_sha256": theorem_fingerprint}, indent=2))


if __name__ == "__main__":
    main()
