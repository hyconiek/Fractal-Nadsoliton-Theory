#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

AS_OF = "2026-03-17"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F677_DIRECTED_CLOSURE = GENERATED / "selector_closure_global_c_v1_directed_strict_v1.json"
IN_F692_SIGN_FIXED_DIRECTED_CLOSURE = (
    GENERATED
    / "selector_closure_global_c_v1_directed_sign_fixed_from_strict_core_payload_weights_strict_convention_v1.json"
)
IN_F690_SIGN_FIXED_STATE = (
    GENERATED
    / "selector_state_global_c_v1_directed_sign_fixed_from_strict_core_payload_weights_strict_convention_v1.json"
)

OUT = GENERATED / "p693_current_strict_t173_output_sign_lift_gauge_covariance_audit_probe.json"
OUT_SUMMARY = GENERATED / "p693_current_strict_t173_output_sign_lift_gauge_covariance_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def require_sign_dict(obj: Any, *, where: str) -> dict[str, int]:
    if not isinstance(obj, dict):
        raise ValueError(f"{where} must be an object mapping pair->±1")
    out: dict[str, int] = {}
    for k, v in obj.items():
        if v not in (-1, 1):
            raise ValueError(f"{where}.{k} must be ±1 (got {v!r})")
        out[str(k)] = int(v)
    return out


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_F677_DIRECTED_CLOSURE, IN_F692_SIGN_FIXED_DIRECTED_CLOSURE, IN_F690_SIGN_FIXED_STATE]
    missing = [str(p.relative_to(REPO)) for p in prereq if not p.exists()]
    if missing:
        artifact = {
            "stage": "P693",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    f677 = load_json(IN_F677_DIRECTED_CLOSURE)
    f692 = load_json(IN_F692_SIGN_FIXED_DIRECTED_CLOSURE)
    f690 = load_json(IN_F690_SIGN_FIXED_STATE)

    try:
        s_prem = require_sign_dict(
            ((f677.get("sign_lift") or {}).get("signs_by_pair")),
            where="F677.sign_lift.signs_by_pair",
        )
        s_fix = require_sign_dict(
            ((f692.get("output_sign_lift") or {}).get("signs_by_pair")),
            where="F692.output_sign_lift.signs_by_pair",
        )
        t = require_sign_dict(
            ((f690.get("outputs") or {}).get("t_by_pair")),
            where="F690.outputs.t_by_pair",
        )
    except Exception as exc:
        artifact = {
            "stage": "P693",
            "status": "INVALID_INPUT_SHAPE",
            "as_of": AS_OF,
            "error": str(exc),
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    expected_pairs = [f"pair{i}" for i in range(1, 6)]
    missing_pairs = [p for p in expected_pairs if p not in s_prem or p not in s_fix or p not in t]
    if missing_pairs:
        status = "NOT_COMPUTABLE_MISSING_PAIR_KEYS"
        artifact = {
            "stage": "P693",
            "status": status,
            "as_of": AS_OF,
            "missing_pairs": missing_pairs,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    checks: list[dict[str, Any]] = []
    mismatches: list[str] = []
    implied: dict[str, int] = {}
    for p in expected_pairs:
        implied[p] = int(t[p] * s_prem[p])
        ok = implied[p] == int(s_fix[p])
        checks.append(
            {
                "pair": p,
                "t_by_pair": int(t[p]),
                "s_out_premise": int(s_prem[p]),
                "implied_s_out_sign_fixed": int(implied[p]),
                "s_out_sign_fixed": int(s_fix[p]),
                "pass": bool(ok),
            }
        )
        if not ok:
            mismatches.append(p)

    discharged = len(mismatches) == 0
    status = (
        "PASS_OUTPUT_SIGN_LIFT_GAUGE_COVARIANT_UNDER_CHART_SIGN_RELIFT"
        if discharged
        else "FAIL_OUTPUT_SIGN_LIFT_NOT_GAUGE_COVARIANT_UNDER_CHART_SIGN_RELIFT"
    )

    artifact = {
        "stage": "P693",
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "goal": "audit_output_sign_lift_gauge_covariance_between_premise_based_and_sign_fixed_directed_closure_exports__no_false_pass",
        "inputs": {
            "premise_based_directed_closure_ref": str(IN_F677_DIRECTED_CLOSURE.relative_to(REPO)),
            "sign_fixed_directed_closure_ref": str(IN_F692_SIGN_FIXED_DIRECTED_CLOSURE.relative_to(REPO)),
            "sign_fixed_state_ref": str(IN_F690_SIGN_FIXED_STATE.relative_to(REPO)),
        },
        "identity_tested": "s_out_sign_fixed(pair) == t_by_pair(pair) * s_out_premise(pair)",
        "signs": {
            "s_out_premise_by_pair": {p: int(s_prem[p]) for p in expected_pairs},
            "t_by_pair": {p: int(t[p]) for p in expected_pairs},
            "s_out_sign_fixed_by_pair": {p: int(s_fix[p]) for p in expected_pairs},
            "implied_s_out_sign_fixed_by_pair": {p: int(implied[p]) for p in expected_pairs},
        },
        "checks_by_pair": checks,
        "result": {
            "gauge_covariant": bool(discharged),
            "mismatched_pairs": mismatches,
        },
        "hard_limits": [
            "Does not claim a directed/sign-sensitive physical orientation datum in strict core.",
            "Does not claim kernel-alone/global QW-2191 discharge.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
        "status": status,
    }

    summary = {
        "stage": "P693",
        "status": status,
        "gauge_covariant": bool(discharged),
        "mismatched_pair_count": len(mismatches),
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

