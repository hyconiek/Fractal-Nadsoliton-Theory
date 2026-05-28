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

OUT = GEN / "p2319_s1269_strict_dihedral_orientation_response_no_go_probe.json"
MD = GEN / "p2319_s1269_strict_dihedral_orientation_response_no_go_probe.md"

N = 12
OMEGA = 0.18575
PHI = 0.16250
BETA = 1.0
ETA = 1.8
ALPHA_STRICT = 4.0 * math.log(2.0)
TOL = 1e-10

SOURCE_FILES = {
    "P2317_FOURIER_DEGENERACY_LANE_AUDIT": GEN / "p2317_s1267_strict_fourier_degeneracy_existing_lane_audit_probe.json",
    "P2318_SELECTOR_LANE_TO_MARGIN_OBLIGATION": GEN / "p2318_s1268_strict_selector_lane_to_task3_margin_bridge_obligation_probe.json",
    "P2308_CURRENT_INTERFACE_NONEXISTENCE": GEN / "p2308_s1258_strict_current_interface_class_response_functional_nonexistence_probe.json",
    "P2306_RESPONSE_ORIENTATION_INTERFACE": GEN / "p2306_s1256_strict_response_orientation_functional_interface_audit_probe.json",
    "P454_O2_GAUGE_EQUIVALENCE_NOTE": ROOT / "P454_CURRENT_STRICT_QW2191_O2_ROTATION_GAUGE_EQUIVALENCE_AUDIT_PROBE.md",
    "N494_UNIQUENESS_UP_TO_CONJUGATION": ROOT / "N494_CURRENT_FIRST_STRICT_QW2190_DIAGONAL_LOCAL_MODE_INDEX_CANONICALIZATION_UNIQUENESS_UP_TO_CONJUGATION_THEOREM.md",
    "H37_SIGN_DISTINCTION_STATE_AUDIT": ROOT / "H37_SIGN_DISTINCTION_STATE_AUDIT.md",
}

GREP_PATTERNS = (
    "dihedral",
    "D12",
    "reflection",
    "orientation",
    "O(2)",
    "Z2",
    "response functional",
    "signed orientation",
    "QW-2191",
    "selector closure",
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


def strict_kernel(distance: int) -> float:
    if distance == 0:
        return 0.0
    return math.cos(OMEGA * distance + PHI) / (1.0 + BETA * (distance ** ETA))


def cyclic_distance(i: int, j: int) -> int:
    raw = abs(i - j) % N
    return min(raw, N - raw)


def kernel_matrix() -> list[list[float]]:
    return [[strict_kernel(cyclic_distance(i, j)) if i != j else 0.0 for j in range(N)] for i in range(N)]


def dot(a: list[float], b: list[float]) -> float:
    return sum(x * y for x, y in zip(a, b))


def mat_vec(matrix: list[list[float]], vector: list[float]) -> list[float]:
    return [dot(row, vector) for row in matrix]


def mat_norm_inf(matrix: list[list[float]]) -> float:
    return max((sum(abs(x) for x in row) for row in matrix), default=0.0)


def vec_norm2(vector: list[float]) -> float:
    return math.sqrt(sum(x * x for x in vector))


def mat_sub(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    return [[a[i][j] - b[i][j] for j in range(N)] for i in range(N)]


def outer(a: list[float], b: list[float]) -> list[list[float]]:
    return [[a[i] * b[j] for j in range(N)] for i in range(N)]


def mat_add(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    return [[a[i][j] + b[i][j] for j in range(N)] for i in range(N)]


def mat_scale(scale: float, a: list[list[float]]) -> list[list[float]]:
    return [[scale * a[i][j] for j in range(N)] for i in range(N)]


def real_fourier_pair(m: int) -> tuple[list[float], list[float]]:
    c = [math.sqrt(2.0 / N) * math.cos(2.0 * math.pi * m * x / N) for x in range(N)]
    s = [-math.sqrt(2.0 / N) * math.sin(2.0 * math.pi * m * x / N) for x in range(N)]
    return c, s


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


def transform_matrix(matrix: list[list[float]], image: list[int]) -> list[list[float]]:
    out = [[0.0 for _ in range(N)] for _ in range(N)]
    for i, target_i in enumerate(image):
        for j, target_j in enumerate(image):
            out[target_i][target_j] = matrix[i][j]
    return out


def reynolds_vector(vector: list[float]) -> list[float]:
    acc = [0.0 for _ in range(N)]
    for row in D12:
        transformed = transform_vector(vector, row["image"])
        for i in range(N):
            acc[i] += transformed[i]
    return [x / len(D12) for x in acc]


def reynolds_matrix(matrix: list[list[float]]) -> list[list[float]]:
    acc = [[0.0 for _ in range(N)] for _ in range(N)]
    for row in D12:
        transformed = transform_matrix(matrix, row["image"])
        for i in range(N):
            for j in range(N):
                acc[i][j] += transformed[i][j]
    return [[acc[i][j] / len(D12) for j in range(N)] for i in range(N)]


def max_group_invariance_residual_matrix(matrix: list[list[float]]) -> float:
    return max(mat_norm_inf(mat_sub(transform_matrix(matrix, row["image"]), matrix)) for row in D12)


def corpus_hits() -> list[dict[str, Any]]:
    paths = sorted(
        set(SOURCE_FILES.values())
        | set(ROOT.glob("P45*_CURRENT_STRICT*QW2191*.md"))
        | set(ROOT.glob("N49*_CURRENT*QW2191*.md"))
        | set(GEN.glob("p230*_s12*_strict*response*json"))
        | set(GEN.glob("p231*_s12*_strict*orientation*json"))
        | set(GEN.glob("p231*_s12*_strict*margin*json"))
        | set(GEN.glob("p231*_s12*_strict*fourier*json"))
    )
    rows: list[dict[str, Any]] = []
    self_paths = {OUT.resolve(), MD.resolve(), Path(__file__).resolve()}
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


def pair_no_go_rows(kernel: list[list[float]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for m in range(1, N // 2):
        c, s = real_fourier_pair(m)
        cc = outer(c, c)
        ss = outer(s, s)
        cs = outer(c, s)
        sc = outer(s, c)
        projector = mat_add(cc, ss)
        anisotropy = mat_sub(cc, ss)
        shear = mat_add(cs, sc)
        handed_skew = mat_sub(cs, sc)
        averaged_c = reynolds_vector(c)
        averaged_s = reynolds_vector(s)
        averaged_anisotropy = reynolds_matrix(anisotropy)
        averaged_shear = reynolds_matrix(shear)
        averaged_handed_skew = reynolds_matrix(handed_skew)
        averaged_projector = reynolds_matrix(projector)
        kc = mat_vec(kernel, c)
        ks = mat_vec(kernel, s)
        block = [[dot(c, kc), dot(c, ks)], [dot(s, kc), dot(s, ks)]]
        lambda_scalar = 0.5 * (block[0][0] + block[1][1])
        scalar_residual = max(
            abs(block[0][0] - lambda_scalar),
            abs(block[1][1] - lambda_scalar),
            abs(block[0][1]),
            abs(block[1][0]),
        )
        rows.append({
            "pair": f"pair{m}",
            "m": m,
            "partner_modes": [m, N - m],
            "kernel_block": block,
            "kernel_scalar_identity_residual": scalar_residual,
            "reynolds_linear_cos_norm": vec_norm2(averaged_c),
            "reynolds_linear_sin_norm": vec_norm2(averaged_s),
            "reynolds_traceless_anisotropy_inf_norm": mat_norm_inf(averaged_anisotropy),
            "reynolds_shear_inf_norm": mat_norm_inf(averaged_shear),
            "reynolds_handed_skew_inf_norm": mat_norm_inf(averaged_handed_skew),
            "projector_reynolds_preservation_inf_residual": mat_norm_inf(mat_sub(averaged_projector, projector)),
            "orientation_odd_linear_survives_d12_average": vec_norm2(averaged_c) > TOL or vec_norm2(averaged_s) > TOL,
            "orientation_anisotropy_survives_d12_average": mat_norm_inf(averaged_anisotropy) > TOL or mat_norm_inf(averaged_shear) > TOL,
            "handed_pseudoscalar_survives_d12_average": mat_norm_inf(averaged_handed_skew) > TOL,
            "only_unsigned_pair_norm_survives": mat_norm_inf(mat_sub(averaged_projector, projector)) < TOL,
        })
    return rows


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2317 = load_json(SOURCE_FILES["P2317_FOURIER_DEGENERACY_LANE_AUDIT"])
    p2318 = load_json(SOURCE_FILES["P2318_SELECTOR_LANE_TO_MARGIN_OBLIGATION"])
    p2308 = load_json(SOURCE_FILES["P2308_CURRENT_INTERFACE_NONEXISTENCE"])
    p2306 = load_json(SOURCE_FILES["P2306_RESPONSE_ORIENTATION_INTERFACE"])
    kernel = kernel_matrix()
    alpha_kernel = mat_scale(ALPHA_STRICT, kernel)
    rows = pair_no_go_rows(kernel)

    p2318_probe = p2318.get("strict_selector_lane_to_task3_margin_bridge_obligation_probe", {}) or {}
    required_lift = float((p2318_probe.get("bridge_obligation_verdict", {}) or {}).get("required_lift_per_step", 0.0))
    missing_fields = list((p2318_probe.get("bridge_obligation_verdict", {}) or {}).get("missing_required_bridge_fields", []))
    p2308_probe = p2308.get("strict_current_interface_class_response_functional_nonexistence_probe", {}) or {}
    p2306_probe = p2306.get("strict_response_orientation_functional_interface_audit_probe", {}) or {}

    symmetry_class_contract = {
        "class_id": "D12_INVARIANT_STRICT_KERNEL_SCALAR_LINEAR_QUADRATIC_RESPONSE_CLASS",
        "allowed_inputs": [
            "strict gate kernel K_strict_gate on Z12 with frozen QW-2049 parameters",
            "strict scalar alpha_geo_strict_derived_v1 = 4 ln 2 used only as a scalar multiplier",
            "identity / pair norm projectors obtained by D12-equivariant linear algebra",
            "linear or quadratic response candidates after D12 Reynolds projection",
        ],
        "excluded_inputs": [
            "external selector premise",
            "lane-specific diagonal/local or Shannon element-order sign convention promoted without bridge",
            "non-D12 source term, boundary condition, replay response map, or time-clock normalization",
            "legacy kernel role transfer",
        ],
        "why_this_is_the_next_honest_class": "P2317 established kernel-alone pair degeneracy; P2318 showed lane defects can be scaled but lack signed response/margin fields.  This packet tests whether the strict D12-symmetric kernel/scalar class itself can supply the missing signed orientation response.",
    }

    no_go_certificate = {
        "kernel_d12_invariance_residual": max_group_invariance_residual_matrix(kernel),
        "alpha_scaled_kernel_d12_invariance_residual": max_group_invariance_residual_matrix(alpha_kernel),
        "pair_count": len(rows),
        "all_kernel_pair_blocks_scalar_identity": all(row["kernel_scalar_identity_residual"] < TOL for row in rows),
        "all_orientation_odd_linear_reynolds_projections_zero": all(
            not row["orientation_odd_linear_survives_d12_average"] for row in rows
        ),
        "all_traceless_quadratic_orientation_projections_zero": all(
            not row["orientation_anisotropy_survives_d12_average"] for row in rows
        ),
        "all_handed_skew_projections_zero_under_reflection": all(
            not row["handed_pseudoscalar_survives_d12_average"] for row in rows
        ),
        "all_unsigned_pair_norm_projectors_survive": all(row["only_unsigned_pair_norm_survives"] for row in rows),
        "max_linear_projection_norm": max(max(row["reynolds_linear_cos_norm"], row["reynolds_linear_sin_norm"]) for row in rows),
        "max_traceless_quadratic_projection_norm": max(
            max(row["reynolds_traceless_anisotropy_inf_norm"], row["reynolds_shear_inf_norm"]) for row in rows
        ),
        "max_handed_skew_projection_norm": max(row["reynolds_handed_skew_inf_norm"] for row in rows),
        "max_unsigned_projector_preservation_residual": max(row["projector_reynolds_preservation_inf_residual"] for row in rows),
        "consequence": "The D12-invariant strict kernel/scalar class leaves only unsigned pair norms.  It cannot export a signed orientation-to-policy response functional or provider_lift_per_step bridge.",
    }

    bridge_obligation_update = {
        "required_lift_per_step": required_lift,
        "p2318_missing_fields_loaded": missing_fields,
        "fields_still_unfilled_by_d12_kernel_scalar_class": missing_fields,
        "d12_class_supplies_signed_scalar_lift_per_step": False,
        "d12_class_supplies_time_step_normalization": False,
        "d12_class_supplies_margin_response_functional": False,
        "d12_class_supplies_orientation_to_policy_sign_rule": False,
        "d12_class_supplies_p2281_replay_semantics": False,
        "d12_class_supplies_admissibility_theorem_no_selector_premise": False,
        "g1_g3_update_allowed": False,
    }

    theorem_export = {
        "theorem_name": "P2319 D12-invariant strict-kernel orientation-response no-go theorem for the current linear/quadratic class",
        "formal_statement": (
            "Within the current strict class generated only by the Z12 strict gate kernel, the scalar 4 ln 2, identity/pair projectors, "
            "and D12-invariant linear or quadratic Reynolds-projected responses, no signed Fourier-pair orientation response survives.  "
            "Linear pair functionals average to zero, traceless orientation quadratics average to zero, reflection-odd handed operators average to zero, "
            "and only unsigned pair norms/projectors remain.  Therefore this class cannot fill the P2318 signed response/margin bridge fields or update G1/G3."
        ),
        "proof_bits": {
            "d12_group_order": len(D12),
            "pair_count": len(rows),
            "kernel_d12_invariance_residual": no_go_certificate["kernel_d12_invariance_residual"],
            "alpha_scaled_kernel_d12_invariance_residual": no_go_certificate["alpha_scaled_kernel_d12_invariance_residual"],
            "max_linear_projection_norm": no_go_certificate["max_linear_projection_norm"],
            "max_traceless_quadratic_projection_norm": no_go_certificate["max_traceless_quadratic_projection_norm"],
            "max_handed_skew_projection_norm": no_go_certificate["max_handed_skew_projection_norm"],
            "max_unsigned_projector_preservation_residual": no_go_certificate["max_unsigned_projector_preservation_residual"],
            "required_lift_per_step_not_supplied": required_lift,
            "p2318_missing_field_count": len(missing_fields),
        },
        "scope_limits": [
            "not a universal theorem over future non-D12 dynamics or future response maps",
            "not a rejection of lane-scoped diagonal/local or Shannon O(2)->Z2 cuts inside their own lanes",
            "not a QW-2191 discharge",
            "not selector closure",
            "not ToE closure",
        ],
        "nonpromotion_rule": "Do not use D12-invariant strict kernel/scalar data alone as a signed selector or Task-3 provider_lift_per_step bridge.",
    }

    theorem_fingerprint = sha256_json(theorem_export)

    probe = {
        "probe_id": "P2319_S1269_STRICT_DIHEDRAL_ORIENTATION_RESPONSE_NO_GO",
        "source_packets": {key: path.relative_to(REPO).as_posix() for key, path in SOURCE_FILES.items()},
        "source_hashes": {f"{key}_sha256": sha256_file(path) for key, path in SOURCE_FILES.items()},
        "far_symmetry_response_grep_audit": {
            "patterns": list(GREP_PATTERNS),
            "hit_count": len(corpus_hits()),
            "top_hits": corpus_hits()[:30],
        },
        "symmetry_class_contract": symmetry_class_contract,
        "strict_parameters": {"n": N, "omega": OMEGA, "phi": PHI, "beta": BETA, "eta": ETA, "alpha_strict": "4 ln 2"},
        "d12_group_order": len(D12),
        "pair_orientation_reynolds_rows": rows,
        "no_go_certificate": no_go_certificate,
        "bridge_obligation_update": bridge_obligation_update,
        "existing_blocker_context": {
            "p2317_result_kind": p2317.get("result_kind"),
            "p2318_result_kind": p2318.get("result_kind"),
            "p2308_nonexistence_current_class": (p2308_probe.get("strict_nonexistence_verdict", {}) or {}).get("nonexistence_proven_for_current_interface_class"),
            "p2306_candidate_count": len(p2306_probe.get("candidate_response_functionals", []) or []),
        },
        "theorem_export": theorem_export,
        "theorem_fingerprint_sha256": theorem_fingerprint,
    }

    gatekeeper_checks = {
        "far_symmetry_response_grep_hits_found": len(probe["far_symmetry_response_grep_audit"]["top_hits"]) >= 5,
        "p2317_loaded": p2317.get("packet_id") == "P2317",
        "p2318_loaded": p2318.get("packet_id") == "P2318",
        "p2308_current_class_nonexistence_loaded": probe["existing_blocker_context"]["p2308_nonexistence_current_class"] is True,
        "d12_group_order_is_24": len(D12) == 24,
        "kernel_d12_invariance_verified": no_go_certificate["kernel_d12_invariance_residual"] < TOL,
        "alpha_scaled_kernel_d12_invariance_verified": no_go_certificate["alpha_scaled_kernel_d12_invariance_residual"] < TOL,
        "all_kernel_pair_blocks_scalar_identity": no_go_certificate["all_kernel_pair_blocks_scalar_identity"],
        "orientation_odd_linear_responses_project_to_zero": no_go_certificate["all_orientation_odd_linear_reynolds_projections_zero"],
        "traceless_quadratic_orientation_responses_project_to_zero": no_go_certificate["all_traceless_quadratic_orientation_projections_zero"],
        "handed_reflection_odd_responses_project_to_zero": no_go_certificate["all_handed_skew_projections_zero_under_reflection"],
        "only_unsigned_pair_norms_survive": no_go_certificate["all_unsigned_pair_norm_projectors_survive"],
        "p2318_bridge_fields_still_unfilled": len(bridge_obligation_update["fields_still_unfilled_by_d12_kernel_scalar_class"]) == 6,
        "strict_g1_g3_not_updated": not bridge_obligation_update["g1_g3_update_allowed"],
        "no_selector_closure_claimed": True,
        "no_qw2191_discharge_claimed": True,
        "no_legacy_kernel_role_transfer": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2319_s1269_v1",
        "packet_id": "P2319",
        "stage_id": "S1269",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_NO_GO_FOR_D12_INVARIANT_STRICT_KERNEL_SCALAR_ORIENTATION_RESPONSE_CLASS",
        "result_kind": "STRICT_DIHEDRAL_ORIENTATION_RESPONSE_NO_GO_NO_G1_G3_UPDATE",
        "strict_dihedral_orientation_response_no_go_probe": probe,
        "recommended_next_honest_step": (
            "Either introduce and justify a non-D12 strict source/response map with replay semantics, or prove a broader nonexistence theorem "
            "for larger response classes; do not reuse kernel/scalar symmetry data alone as a selector."
        ),
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "# P2319/S1269 — D12-invariant strict-kernel orientation-response no-go\n\n"
        f"- Result kind: `{payload['result_kind']}`.\n"
        f"- D12 group order: `{len(D12)}`.\n"
        f"- Kernel D12 invariance residual: `{no_go_certificate['kernel_d12_invariance_residual']:.3e}`.\n"
        f"- Max linear Reynolds projection norm on pair planes: `{no_go_certificate['max_linear_projection_norm']:.3e}`.\n"
        f"- Max traceless quadratic Reynolds projection norm: `{no_go_certificate['max_traceless_quadratic_projection_norm']:.3e}`.\n"
        f"- Max handed/skew Reynolds projection norm: `{no_go_certificate['max_handed_skew_projection_norm']:.3e}`.\n"
        f"- Unsigned pair projectors preserved with max residual: `{no_go_certificate['max_unsigned_projector_preservation_residual']:.3e}`.\n"
        "- Verdict: the strict D12-invariant kernel/scalar class leaves only unsigned pair norms and cannot export a signed Task-3 provider-lift response.\n"
        "- Guardrail: no G1/G3 update, no QW-2191 discharge, no selector closure, no ToE closure.\n",
        encoding="utf-8",
    )
    print(json.dumps({"wrote": str(OUT.relative_to(REPO)), "theorem_fingerprint_sha256": theorem_fingerprint}, indent=2))


if __name__ == "__main__":
    main()
