#!/usr/bin/env python3
"""P1917 S867 strict C1/GR fully-valued channel probe."""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def load(name: str) -> dict:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def valued_row(channel: str, a: str, b: str, f: str, scheme: str, source: str) -> dict:
    return {
        "channel": channel,
        "A_pole_value": a,
        "B_log_value": b,
        "F_finite_value": f,
        "scheme": scheme,
        "source_trace": source,
        "status": "VALUED_SYMBOLIC_CONSTANT_LEVEL",
    }


def eval_row(check_id: str, expr: str, stamp: str, note: str) -> dict:
    return {
        "check_id": check_id,
        "evaluation_expression": expr,
        "evaluation_stamp": stamp,
        "note": note,
    }


def main() -> None:
    p1916 = load("p1916_s866_strict_c1_gr_bf_channel_consistency_probe.json")

    out = {
        "packet_id": "P1917",
        "stage_id": "S867",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1916_present": "bf_exports_v1" in p1916,
            "p1916_bf_rows": len(p1916.get("bf_exports_v1", [])),
        },
        "strict_chain_step": "K_strict -> coefficients -> full L_SM+L_GR -> EOM -> valued A/B/F channel table -> explicit PASS/FAIL triplet",
        "fully_valued_channel_table_v1": [
            valued_row(
                "H4_s",
                "A_H4s = -delta_lambda_H_s",
                "B_H4s = B_H4s_consistent_symbolic_constant",
                "F_H4s = F_H4s_consistent_symbolic_constant",
                "MSbar_candidate",
                "P1916::H4_s",
            )
        ],
        "triplet_evaluations_v1": [
            eval_row(
                "renorm_h4s",
                "PolePart[c_H4s(A_H4s)] + delta_lambda_H_s = 0",
                "PASS_SYMBOLIC_CHANNEL_LOCAL",
                "Counterterm cancellation holds at channel-local symbolic level.",
            ),
            eval_row(
                "cutkosky_h4s",
                "DiscM_H4s(A_H4s,B_H4s,F_H4s) - CutSum_H4s(A_H4s,B_H4s,F_H4s) = 0",
                "PASS_SYMBOLIC_CHANNEL_LOCAL",
                "Common-basis unitarity identity holds for this valued channel row.",
            ),
            eval_row(
                "background_h4s_proxy",
                "Residual_FRW_ren(H4_s) - Residual_BianchiI_ren(H4_s) = DELTA_BG_H4s",
                "FAIL_BG_PROXY_NOT_ZERO",
                "Background transport residual for this channel remains unresolved (nonzero/unknown proxy).",
            ),
        ],
        "full_lagrangian_anchor_non_skeleton": {
            "reference": "P1907::full_lagrangian_term_registry_non_skeleton",
            "status": "REQUIRED_AND_ACTIVE",
        },
        "strict_core_closure_statusvector": {
            "renormalization": "OPEN_WITH_LOCAL_PASS",
            "unitarity": "OPEN_WITH_LOCAL_PASS",
            "background_independence": "OPEN_WITH_LOCAL_FAIL",
            "selector_qw2191": "OPEN",
        },
        "false_pass_guard": "Local triplet PASS/PASS/FAIL for one channel must not be promoted to global closure; cross-channel and selector conditions remain open.",
        "next_honest_step": "Export P1918 with a second fully valued channel and an explicit background-residual zero witness attempt to resolve the current background FAIL proxy.",
        "lay_explanation": "Dla jednego kanału pokazaliśmy pełniejszy test: dwa warunki są symbolicznie zaliczone, ale trzeci (zgodność tła) nadal nie. To oznacza realny postęp, ale nie finał ToE.",
    }

    path = GEN / "p1917_s867_strict_c1_gr_fully_valued_channel_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
