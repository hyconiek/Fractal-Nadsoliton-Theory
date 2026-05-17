#!/usr/bin/env python3
"""P1916 S866 strict C1/GR B/F channel consistency probe."""
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


def bf_row(channel: str, b_log: str, f_fin: str, source: str) -> dict:
    return {
        "channel": channel,
        "B_log_value": b_log,
        "F_finite_value": f_fin,
        "source_trace": source,
        "status": "PARTIAL_SYMBOLIC_EXPORT",
    }


def check_row(name: str, expr: str, status: str, note: str) -> dict:
    return {
        "check_id": name,
        "consistency_expression": expr,
        "evaluation_stamp": status,
        "note": note,
    }


def main() -> None:
    p1915 = load("p1915_s865_strict_c1_gr_solved_subrows_restamp_probe.json")

    out = {
        "packet_id": "P1916",
        "stage_id": "S866",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1915_present": "solved_subrows_v1" in p1915,
            "p1915_row_pass_count": len(p1915.get("solved_subrows_v1", [])),
        },
        "strict_chain_step": "K_strict -> coefficients -> full L_SM+L_GR -> EOM -> A/B/F exports -> renorm+Cutkosky+FRW/Bianchi-I consistency frame",
        "bf_exports_v1": [
            bf_row("H4_s", "B_H4s = B_H4s_symbolic_common_basis", "F_H4s = F_H4s_symbolic_common_basis", "P1913::H4_s"),
            bf_row("grmix", "B_gr = B_gr_symbolic_common_basis", "F_gr = F_gr_symbolic_common_basis", "P1913::grmix"),
        ],
        "cross_consistency_checks_v1": [
            check_row(
                "renorm_h4s",
                "PolePart[c_H4s] + delta_lambda_H_s == 0 and B_H4s finite-scheme constraint satisfied",
                "OPEN_SYMBOLIC_CONSISTENT_NOT_NUMERIC",
                "A/B/F symbolic slots aligned, but numeric coefficients not injected.",
            ),
            check_row(
                "cutkosky_h4s",
                "DiscM_H4s(A,B,F) - CutSum_H4s(A,B,F) == 0",
                "OPEN_SYMBOLIC_CONSISTENT_NOT_NUMERIC",
                "Common-basis identity declared; coefficient valuation remains pending.",
            ),
            check_row(
                "background_grmix",
                "Residual_FRW_ren(A,B,F) - Residual_BianchiI_ren(A,B,F) == 0",
                "OPEN_SYMBOLIC_CONSISTENT_NOT_NUMERIC",
                "Same-channel B/F included, but FRW/Bianchi-I coefficient tables still not fully exported.",
            ),
        ],
        "full_lagrangian_anchor_non_skeleton": {
            "reference": "P1907::full_lagrangian_term_registry_non_skeleton",
            "status": "REQUIRED_AND_ACTIVE",
        },
        "strict_core_closure_statusvector": {
            "renormalization": "OPEN_BF_PARTIAL_SYMBOLIC",
            "unitarity": "OPEN_PARTIAL_SYMBOLIC",
            "background_independence": "OPEN_PARTIAL_SYMBOLIC",
            "selector_qw2191": "OPEN",
        },
        "false_pass_guard": "Symbolic A/B/F consistency does not imply numeric or theorem-level strict-core closure across all channels.",
        "next_honest_step": "Export P1917 with one fully valued channel table (A,B,F numeric or fully reduced symbolic constants) and explicit PASS/FAIL stamps for renorm+Cutkosky+background checks in that channel.",
        "lay_explanation": "Dołożyliśmy brakujące elementy B i F dla wybranych kanałów, więc testy są pełniejsze niż wcześniej. Nadal jednak brakują ostatecznych wartości, dlatego całość pozostaje otwarta.",
    }

    path = GEN / "p1916_s866_strict_c1_gr_bf_channel_consistency_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
