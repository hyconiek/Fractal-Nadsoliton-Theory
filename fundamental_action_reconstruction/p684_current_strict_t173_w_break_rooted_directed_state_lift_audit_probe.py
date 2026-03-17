#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

AS_OF = "2026-03-17"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F647 = GENERATED / "f647_strict_witness_provider_export_packet_for_seed_v1_realization_attempt.json"

IN_A = {
    1: GENERATED / "a_1_pair1_orientation_projector_operator_strict_core_v1.json",
    2: GENERATED / "a_2_pair2_orientation_projector_operator_strict_core_v1.json",
    3: GENERATED / "a_3_pair3_orientation_projector_operator_strict_core_v1.json",
    4: GENERATED / "a_4_pair4_orientation_projector_operator_strict_core_v1.json",
    5: GENERATED / "a_5_pair5_orientation_projector_operator_strict_core_v1.json",
}

IN_O_1M = {
    2: GENERATED / "o12_pair1_pair2_selector_chart_transport_operator_strict_derived_from_sigma_int_alpha12_v1.json",
    3: GENERATED / "o13_pair1_pair3_selector_chart_transport_operator_axis_only_alpha13_mod_pi_strict_core_v1.json",
    4: GENERATED / "o14_pair1_pair4_selector_chart_transport_operator_axis_only_alpha14_mod_pi_strict_core_v1.json",
    5: GENERATED / "o15_pair1_pair5_selector_chart_transport_operator_axis_only_alpha15_mod_pi_strict_core_v1.json",
}

OUT_JSON = GENERATED / "p684_current_strict_t173_w_break_rooted_directed_state_lift_audit_probe.json"
OUT_SUMMARY = GENERATED / "p684_current_strict_t173_w_break_rooted_directed_state_lift_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def is_numeric_list_len(obj: Any, n: int) -> bool:
    return (
        isinstance(obj, list)
        and len(obj) == n
        and all(isinstance(v, (int, float)) and math.isfinite(float(v)) for v in obj)
    )


def sign(x: float) -> int:
    return 1 if x > 0 else -1


def dot(a: list[float], b: list[float]) -> float:
    return sum(float(x) * float(y) for x, y in zip(a, b))


def matvec12(m: list[list[float]], v: list[float]) -> list[float]:
    return [sum(float(m[i][j]) * float(v[j]) for j in range(12)) for i in range(12)]


def outer2(v: list[float]) -> list[list[float]]:
    return [[float(v[i]) * float(v[j]) for j in range(2)] for i in range(2)]


def max_abs_diff_mat2(a: list[list[float]], b: list[list[float]]) -> float:
    return max(abs(float(a[i][j] - b[i][j])) for i in range(2) for j in range(2))


def max_abs_diff_vec(a: list[float], b: list[float]) -> float:
    return max(abs(float(x - y)) for x, y in zip(a, b))


def load_transport_matrix(path: Path) -> list[list[float]]:
    obj = load_json(path)
    outs = obj.get("outputs") or {}
    for key in ("O12", "O13", "O14", "O15", "O"):
        M = outs.get(key)
        if isinstance(M, list) and len(M) == 12 and all(is_numeric_list_len(row, 12) for row in M):
            return [[float(x) for x in row] for row in M]
    raise ValueError(f"{path} outputs does not contain a 12x12 transport matrix under keys O12/O13/O14/O15/O")


def load_pair_data(m: int) -> tuple[list[float], list[float], list[list[float]]]:
    a = load_json(IN_A[m])
    data = a.get("data") or {}
    u_key = "u_1" if m == 1 else f"u_{m}"
    coords_key = "u_1_coords_in_c1_s1" if m == 1 else f"u_{m}_coords_in_c{m}_s{m}"
    A_key = "A_1_pair1_matrix_in_c1_s1" if m == 1 else f"A_{m}_pair{m}_matrix_in_c{m}_s{m}"
    u_raw = data.get(u_key)
    coords_raw = data.get(coords_key)
    A_raw = data.get(A_key)
    if not is_numeric_list_len(u_raw, 12):
        raise ValueError(f"Invalid {IN_A[m].relative_to(REPO)}: data.{u_key} must be length-12 numeric list")
    if not is_numeric_list_len(coords_raw, 2):
        raise ValueError(f"Invalid {IN_A[m].relative_to(REPO)}: data.{coords_key} must be length-2 numeric list")
    if not (isinstance(A_raw, list) and len(A_raw) == 2 and all(is_numeric_list_len(row, 2) for row in A_raw)):
        raise ValueError(f"Invalid {IN_A[m].relative_to(REPO)}: data.{A_key} must be a 2x2 numeric list")
    u = [float(v) for v in u_raw]
    coords = [float(v) for v in coords_raw]
    A = [[float(A_raw[0][0]), float(A_raw[0][1])], [float(A_raw[1][0]), float(A_raw[1][1])]]
    return u, coords, A


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_F647] + [IN_A[m] for m in sorted(IN_A)] + [IN_O_1M[m] for m in sorted(IN_O_1M)]
    missing = [str(p.relative_to(REPO)) for p in prereq if not p.exists()]
    if missing:
        artifact = {
            "stage": "P684",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_JSON)
        return

    f647 = load_json(IN_F647)
    payload = (f647.get("constructed_source_object") or {}).get("exported_payload") or {}
    w_break_raw = payload.get("w_break_by_x")
    if not is_numeric_list_len(w_break_raw, 12):
        raise SystemExit("Invalid F647: constructed_source_object.exported_payload.w_break_by_x must be length-12 numeric list")
    w_break = [float(v) for v in w_break_raw]

    tol_dot_nonzero = 1e-12
    tol_projector = 1e-12
    tol_transport = 1e-9

    u: dict[int, list[float]] = {}
    coords: dict[int, list[float]] = {}
    A: dict[int, list[list[float]]] = {}
    for m in sorted(IN_A):
        uu, cc, AA = load_pair_data(m)
        u[m] = uu
        coords[m] = cc
        A[m] = AA

    dot_root = dot(w_break, u[1])
    if abs(dot_root) <= tol_dot_nonzero:
        artifact = {
            "stage": "P684",
            "status": "FAIL_NO_ROOT_SIGN_LIFT_FROM_W_BREAK_ON_PAIR1",
            "as_of": AS_OF,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "root": {"pair": "pair1", "dot_w_break_u_1": float(dot_root), "dot_nonzero_abs_tol": tol_dot_nonzero},
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_JSON)
        return

    s: dict[int, int] = {1: sign(dot_root)}
    u1_star = [float(s[1]) * x for x in u[1]]

    transport_alignment: dict[str, Any] = {}
    for m in (2, 3, 4, 5):
        O = load_transport_matrix(IN_O_1M[m])
        v = matvec12(O, u1_star)
        d = dot(v, u[m])
        if abs(d) <= tol_dot_nonzero:
            raise SystemExit(f"FAIL: dot(O_1{m} u1*, u_{m}) is numerically zero; cannot propagate sign")
        s[m] = sign(d)
        u_m_star = [float(s[m]) * x for x in u[m]]
        max_abs = max_abs_diff_vec(v, u_m_star)
        transport_alignment[f"pair{m}"] = {"dot": float(d), "sign_s_m": int(s[m]), "max_abs_diff": float(max_abs)}

    projector_match: dict[str, Any] = {}
    for m in range(1, 6):
        coords_star = [float(s[m]) * coords[m][0], float(s[m]) * coords[m][1]]
        A_star = outer2(coords_star)
        max_abs = max_abs_diff_mat2(A[m], A_star)
        projector_match[f"pair{m}"] = {"max_abs_diff": float(max_abs), "matches": bool(max_abs <= tol_projector)}

    all_projectors_match = all(bool(projector_match[f"pair{m}"]["matches"]) for m in range(1, 6))
    all_transports_align = all(
        bool(transport_alignment[f"pair{m}"]["max_abs_diff"] <= tol_transport) for m in (2, 3, 4, 5)
    )

    status = "PASS_W_BREAK_ROOTED_DIRECTED_STATE_LIFT_DESCENDS_TO_PROJECTIVE_AND_TRANSPORT_ALIGNS"
    if not all_projectors_match:
        status = "FAIL_DIRECTED_SECTION_DOES_NOT_DESCEND_TO_PROJECTIVE_PROJECTORS"
    elif not all_transports_align:
        status = "FAIL_ROOTED_TRANSPORT_ALIGNMENT_EXCEEDS_TOLERANCE"

    artifact = {
        "stage": "P684",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "goal": "audit_w_break_rooted_residual_sign_fix_on_pair1_and_rooted_transport_directed_vector_section_on_C_v1",
        "inputs": {
            "seed_w_break_ref": str(IN_F647.relative_to(REPO)),
            "a_m_refs": {str(m): str(IN_A[m].relative_to(REPO)) for m in sorted(IN_A)},
            "rooted_transport_refs": {str(m): str(IN_O_1M[m].relative_to(REPO)) for m in sorted(IN_O_1M)},
        },
        "root": {"pair": "pair1", "dot_w_break_u_1": float(dot_root), "root_sign_s_1": int(s[1])},
        "sign_lift": {
            "rule": "rooted_transport_from_w_break_on_pair1",
            "signs_by_pair": {f"pair{m}": int(s[m]) for m in range(1, 6)},
            "counts_as_strict_physical_orientation_datum": False,
        },
        "compatibility": {
            "descends_to_projective_projectors": {
                "max_abs_diff_by_pair": {f"pair{m}": float(projector_match[f'pair{m}']['max_abs_diff']) for m in range(1, 6)},
                "tol": tol_projector,
                "all_pairs_match": all_projectors_match,
            },
            "rooted_transport_alignment": {
                "by_pair": transport_alignment,
                "tol": tol_transport,
                "all_pairs_align": all_transports_align,
                "note": "Alignment is checked only along rooted edges O_1m; some edges are axis-only representatives (see N512 boundary).",
            },
        },
        "hard_limits": [
            "Convention/section choice only; does not claim a directed/sign-sensitive physical orientation datum in strict core.",
            "Does not claim Aut(Z_12)-invariant sign canonicity (rooted sign uses w_break and a chosen chart embedding).",
            "Does not claim strict-core selector closure / admissible S_sel_int.",
            "Does not claim kernel-alone/global QW-2191 discharge.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P684",
        "status": status,
        "root_sign_s_1": int(s[1]),
        "signs_by_pair": {f"pair{m}": int(s[m]) for m in range(1, 6)},
        "descends_to_projective": all_projectors_match,
        "rooted_transport_alignment_ok": all_transports_align,
        "counts_as_strict_physical_orientation_datum": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()

