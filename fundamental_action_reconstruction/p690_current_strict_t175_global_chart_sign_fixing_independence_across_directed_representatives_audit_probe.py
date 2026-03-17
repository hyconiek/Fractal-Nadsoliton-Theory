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

IN_F690_OBJECT = (
    GENERATED
    / "selector_state_global_c_v1_directed_sign_fixed_from_strict_core_payload_weights_strict_convention_v1.json"
)
IN_F684_STATE = (
    GENERATED
    / "selector_state_global_c_v1_directed_rooted_transport_from_s_sel_int_w_break_strict_convention_v1.json"
)
IN_F647 = GENERATED / "f647_strict_witness_provider_export_packet_for_seed_v1_realization_attempt.json"

OUT = (
    GENERATED
    / "p690_current_strict_t175_global_chart_sign_fixing_independence_across_directed_representatives_audit_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p690_current_strict_t175_global_chart_sign_fixing_independence_across_directed_representatives_audit_probe_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def is_numeric_list_len(obj: Any, n: int) -> bool:
    return (
        isinstance(obj, list)
        and len(obj) == n
        and all(isinstance(v, (int, float)) and math.isfinite(float(v)) for v in obj)
    )


def dot(a: list[float], b: list[float]) -> float:
    if len(a) != len(b):
        raise ValueError("dot: mismatched lengths")
    return float(sum(float(x) * float(y) for x, y in zip(a, b)))


def max_abs_diff(a: list[float], b: list[float]) -> float:
    return max(abs(float(x) - float(y)) for x, y in zip(a, b))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_F690_OBJECT, IN_F684_STATE, IN_F647]
    missing = [str(p.relative_to(REPO)) for p in prereq if not p.exists()]
    if missing:
        artifact = {
            "stage": "P690",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    tol_match = 1e-9
    tol_positive = 0.0

    f690 = load_json(IN_F690_OBJECT)
    state_b = load_json(IN_F684_STATE)
    f647 = load_json(IN_F647)

    weight_choice = ((f690.get("construction") or {}).get("weight_choice_by_pair") or {})
    if not isinstance(weight_choice, dict):
        artifact = {
            "stage": "P690",
            "status": "INVALID_F690_OBJECT_SHAPE",
            "as_of": AS_OF,
            "error": "F690 object must contain construction.weight_choice_by_pair dict",
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    u_fixed_a = (((f690.get("outputs") or {}).get("u_vectors_directed_sign_fixed")) or {})
    if not isinstance(u_fixed_a, dict):
        artifact = {
            "stage": "P690",
            "status": "INVALID_F690_OBJECT_SHAPE",
            "as_of": AS_OF,
            "error": "F690 object must contain outputs.u_vectors_directed_sign_fixed dict",
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    payload = ((f647.get("constructed_source_object") or {}).get("exported_payload") or {})
    if not isinstance(payload, dict):
        raise SystemExit("Invalid F647 shape: expected constructed_source_object.exported_payload dict")

    w_break = payload.get("w_break_by_x")
    w_ref = payload.get("w_ref_unnormalized_by_x")
    if not is_numeric_list_len(w_break, 12) or not is_numeric_list_len(w_ref, 12):
        raise SystemExit("Invalid F647 payload weights: expected length-12 arrays w_break_by_x and w_ref_unnormalized_by_x")

    w_by_key = {
        "w_break_by_x": [float(v) for v in w_break],
        "w_ref_unnormalized_by_x": [float(v) for v in w_ref],
    }

    u_outs_b = ((state_b.get("outputs") or {}).get("u_vectors_directed") or {})
    if not isinstance(u_outs_b, dict):
        raise SystemExit("Invalid F684 directed state: expected outputs.u_vectors_directed dict")

    # Apply the same sign-fixing rule to the F684 directed representative.
    per_pair: dict[str, Any] = {}
    errors: list[dict[str, Any]] = []
    max_diff_overall = 0.0
    all_ok = True

    for m in range(1, 6):
        pair = f"pair{m}"
        u_key = f"u_{m}"
        u_b = u_outs_b.get(u_key)
        u_a = u_fixed_a.get(u_key)
        if not is_numeric_list_len(u_b, 12) or not is_numeric_list_len(u_a, 12):
            errors.append({"pair": pair, "error": "missing or invalid u vector (expected length-12 numeric list)"})
            all_ok = False
            continue

        weight_key = weight_choice.get(pair)
        if not isinstance(weight_key, str) or weight_key not in w_by_key:
            errors.append({"pair": pair, "error": f"invalid weight key in F690: {weight_key!r}"})
            all_ok = False
            continue

        w = w_by_key[weight_key]
        d_raw = dot(w, [float(v) for v in u_b])
        if d_raw == 0.0:
            errors.append({"pair": pair, "error": "dot(w,u)=0; cannot fix sign deterministically", "weight_key": weight_key})
            all_ok = False
            continue

        t = 1 if d_raw > 0.0 else -1
        u_b_fixed = [float(t) * float(v) for v in u_b]
        d_after = dot(w, u_b_fixed)

        diff = float(max_abs_diff([float(v) for v in u_a], u_b_fixed))
        max_diff_overall = max(max_diff_overall, diff)

        ok_match = diff <= tol_match
        ok_positive = d_after > tol_positive
        ok = bool(ok_match and ok_positive)
        all_ok = all_ok and ok

        per_pair[pair] = {
            "pair": pair,
            "weight_key": weight_key,
            "t_on_input_state": int(t),
            "dot_raw": float(d_raw),
            "dot_after_fix": float(d_after),
            "max_abs_diff_vs_F690": diff,
            "ok_match": bool(ok_match),
            "ok_dot_positive": bool(ok_positive),
            "ok": ok,
        }

    status = (
        "PASS_SIGN_FIXED_DIRECTED_REPRESENTATIVE_INDEPENDENT_OF_STARTING_DIRECTED_STATE"
        if all_ok and not errors
        else "FAIL_SIGN_FIXED_DIRECTED_REPRESENTATIVE_NOT_INDEPENDENT_OR_INVALID"
    )

    artifact = {
        "stage": "P690",
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "goal": "audit_that_F690_chart_sign_fixing_yields_the_same_sign_fixed_directed_representative_when_applied_to_multiple_exported_directed_states__no_false_pass",
        "inputs": {
            "f690_sign_fixed_state_ref": str(IN_F690_OBJECT.relative_to(REPO)),
            "f684_directed_state_ref": str(IN_F684_STATE.relative_to(REPO)),
            "f647_weights_payload_ref": str(IN_F647.relative_to(REPO)),
        },
        "status": status,
        "tolerances": {"tol_match_max_abs_diff": tol_match, "dot_after_fix_must_be_gt": tol_positive},
        "pair_audit": {
            "per_pair": per_pair,
            "errors": errors,
            "max_abs_diff_overall": max_diff_overall,
        },
        "result": {"independence_ok": all_ok and not errors},
        "counts_as_strict_physical_orientation_datum": False,
        "hard_limits": [
            "Convention-layer gauge fixing only; does not claim a strict physical sign/orientation datum.",
            "Does not claim Aut(Z_12)-invariant sign canonicity (N462 discipline).",
            "Does not claim kernel-alone/global QW-2191 discharge.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P690",
        "status": status,
        "independence_ok": bool(all_ok and not errors),
        "max_abs_diff_overall": max_diff_overall,
        "counts_as_strict_physical_orientation_datum": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

