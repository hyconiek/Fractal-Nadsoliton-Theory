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
IN_STATE_PREMISE = GENERATED / "selector_state_global_c_v1_directed_strict_v1.json"

OUT_OBJECT = (
    GENERATED
    / "selector_state_global_c_v1_directed_sign_fixed_from_strict_core_payload_weights_strict_convention_v1.json"
)
OUT_SUMMARY = GENERATED / "f690_current_strict_t175_global_chart_sign_fixing_from_strict_core_payload_weights_export_packet_summary.json"


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


def sign_nonzero(x: float) -> int:
    if x > 0.0:
        return 1
    if x < 0.0:
        return -1
    raise ValueError("sign_nonzero: x=0 (cannot fix sign)")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_F647, IN_STATE_PREMISE]
    missing = [str(p.relative_to(REPO)) for p in prereq if not p.exists()]
    if missing:
        summary = {
            "stage": "F690",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT_OBJECT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_OBJECT)
        return

    f647 = load_json(IN_F647)
    state_premise = load_json(IN_STATE_PREMISE)

    payload = ((f647.get("constructed_source_object") or {}).get("exported_payload") or {})
    if not isinstance(payload, dict):
        raise SystemExit("Invalid F647 shape: expected constructed_source_object.exported_payload dict")

    w_break = payload.get("w_break_by_x")
    w_ref = payload.get("w_ref_unnormalized_by_x")
    if not is_numeric_list_len(w_break, 12) or not is_numeric_list_len(w_ref, 12):
        raise SystemExit("Invalid F647 payload weights: expected length-12 arrays w_break_by_x and w_ref_unnormalized_by_x")

    u_outs = ((state_premise.get("outputs") or {}).get("u_vectors_directed") or {})
    if not isinstance(u_outs, dict):
        raise SystemExit("Invalid premise-based state: expected outputs.u_vectors_directed dict")

    # Deterministic chartwise sign-fixing rule:
    # - pair1, pair5: use w_break_by_x (reflection-breaking seed weight)
    # - pair2..pair4: use w_ref_unnormalized_by_x (strict-core reference weight)
    weight_key_by_pair = {
        "pair1": "w_break_by_x",
        "pair2": "w_ref_unnormalized_by_x",
        "pair3": "w_ref_unnormalized_by_x",
        "pair4": "w_ref_unnormalized_by_x",
        "pair5": "w_break_by_x",
    }

    w_by_key = {
        "w_break_by_x": [float(v) for v in w_break],
        "w_ref_unnormalized_by_x": [float(v) for v in w_ref],
    }

    t_by_pair: dict[str, int] = {}
    dot_raw_by_pair: dict[str, float] = {}
    dot_by_pair: dict[str, float] = {}
    u_fixed: dict[str, list[float]] = {}
    errors: list[dict[str, Any]] = []

    for m in range(1, 6):
        pair = f"pair{m}"
        u_key = f"u_{m}"
        vec = u_outs.get(u_key)
        if not is_numeric_list_len(vec, 12):
            errors.append({"pair": pair, "error": f"missing or invalid {u_key} (expected length-12 numeric list)"})
            continue

        weight_key = weight_key_by_pair[pair]
        w = w_by_key[weight_key]
        d_raw = dot(w, [float(v) for v in vec])
        try:
            t = sign_nonzero(d_raw)
        except Exception as exc:
            errors.append(
                {
                    "pair": pair,
                    "error": f"cannot_fix_sign: dot={d_raw} ({exc})",
                    "weight_key": weight_key,
                }
            )
            continue

        # Flip u if needed to force dot(w,u_fixed) > 0.
        u_vec = [float(v) for v in vec]
        if t == -1:
            u_vec = [-x for x in u_vec]
        d_abs = abs(float(d_raw))

        t_by_pair[pair] = int(t)
        dot_raw_by_pair[pair] = float(d_raw)
        dot_by_pair[pair] = float(d_abs)
        u_fixed[u_key] = u_vec

    if errors or len(u_fixed) != 5:
        status = "NOT_COMPUTABLE_SIGN_FIXING_FAILED_ON_SOME_CHARTS"
    else:
        status = "PASS_EXPORTED_GLOBAL_SIGN_FIXED_DIRECTED_STATE_OBJECT"

    obj = {
        "object": "SelectorState_global_C_v1_directed_sign_fixed_from_strict_core_payload_weights_strict_convention_v1",
        "stage": "F690",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "intent": (
            "Export an explicit chart-level sign-fixing (Z2 0-cochain) datum and the resulting sign-fixed directed representative "
            "on C_v1, constructed deterministically from already exported strict-core payload weights. This is a strict_convention "
            "gauge-fixing layer: it does not claim a strict physical sign/orientation datum, does not claim Aut(Z_12)-invariant sign "
            "canonicity, and does not imply kernel-alone/global QW-2191 discharge."
        ),
        "domain": {"configuration_space": "C_v1", "charts": ["pair1", "pair2", "pair3", "pair4", "pair5"]},
        "depends_on": {
            "weights_payload_ref": str(IN_F647.relative_to(REPO)),
            "premise_based_directed_state_ref": str(IN_STATE_PREMISE.relative_to(REPO)),
        },
        "construction": {
            "rule": "u_i_fixed := sign(<w_i,u_i>)*u_i with <w_i,u_i_fixed> > 0",
            "weight_choice_by_pair": weight_key_by_pair,
            "weight_keys_available": sorted(w_by_key.keys()),
            "notes": [
                "This is an explicit convention-layer gauge fixing; it is not claimed as a strict physical sign datum.",
                "Different exported directed representatives may differ by chart sign gauge; this rule aims to pick one canonical representative.",
            ],
        },
        "outputs": {
            "t_by_pair": {k: int(v) for k, v in sorted(t_by_pair.items())},
            "dot_raw_by_pair": {k: float(v) for k, v in sorted(dot_raw_by_pair.items())},
            "dot_by_pair": {k: float(v) for k, v in sorted(dot_by_pair.items())},
            "u_vectors_directed_sign_fixed": u_fixed,
            "errors": errors,
        },
        "counts_as_strict_physical_orientation_datum": False,
        "hard_limits": [
            "Convention layer only; does not claim a strict-core physical sign/orientation datum.",
            "Does not claim Aut(Z_12)-invariant sign canonicity (N462 discipline).",
            "Does not claim kernel-alone/global QW-2191 discharge.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F690",
        "status": status,
        "exported_object_ref": str(OUT_OBJECT.relative_to(REPO)),
        "t_by_pair": {k: int(v) for k, v in sorted(t_by_pair.items())},
        "counts_as_strict_physical_orientation_datum": False,
        "no_false_pass": True,
    }

    OUT_OBJECT.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_OBJECT)


if __name__ == "__main__":
    main()
