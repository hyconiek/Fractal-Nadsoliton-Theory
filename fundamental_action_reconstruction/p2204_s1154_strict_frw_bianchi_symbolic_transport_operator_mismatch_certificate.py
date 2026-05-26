#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2203 = GEN / "p2203_s1153_strict_frw_bianchi_transport_residual_map_under_shared_majorant.json"
OUT = GEN / "p2204_s1154_strict_frw_bianchi_symbolic_transport_operator_mismatch_certificate.json"
MD = GEN / "p2204_s1154_strict_frw_bianchi_symbolic_transport_operator_mismatch_certificate.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2203 = load(IN_2203)
    packet = p2203.get("strict_frw_bianchi_transport_residual_map_under_shared_majorant", {}) or {}

    d, omega0, phi, beta, eta, m = sp.symbols("d omega0 phi beta eta m", positive=True, finite=True)
    x = sp.symbols("x", real=True)

    k_frw = sp.cos(omega0 * d + phi) / (1 + beta * d**eta)
    k_bi = sp.cos((omega0 * m) * d + phi) / (1 + beta * d**eta)
    delta_sq = sp.expand_trig(sp.simplify(k_frw**2 - k_bi**2))

    # exact sufficient condition candidate for zero mismatch
    sufficient_condition = sp.Eq(m, 1)
    delta_under_condition = sp.simplify(delta_sq.subs({m: 1}))

    # first-order around m=1 (for local theorem-gap diagnostics)
    first_dm_at_1 = sp.simplify(sp.diff(delta_sq, m).subs({m: 1}))

    rows = packet.get("residual_map_rows", []) or []
    max_residual = float(packet.get("summary", {}).get("max_residual", 0.0) or 0.0)
    near_one_rows = [r for r in rows if abs(float(r.get("omega_mult_probe", 0.0)) - 1.0) < 1e-12]
    at_one_residual = float(near_one_rows[0]["transport_residual_l1"]) if near_one_rows else None

    theorem_gap_statement = {
        "task3_requirement": "global FRW<->Bianchi transport closure for nu branch",
        "symbolic_zero_mismatch_condition": str(sufficient_condition),
        "symbolic_zero_mismatch_verified": bool(delta_under_condition == 0),
        "first_order_mismatch_term_at_m_eq_1": str(first_dm_at_1),
        "numeric_at_m_eq_1_residual": at_one_residual,
        "global_closure_proven": False,
        "reason": "symbolic sufficient condition and local expansion exported; no global transport theorem proved",
    }

    payload = {
        "schema_version": "p2204_s1154_v1",
        "packet_id": "P2204",
        "stage_id": "S1154",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_FRW_BIANCHI_SYMBOLIC_TRANSPORT_OPERATOR_MISMATCH_CERTIFICATE",
        "strict_frw_bianchi_symbolic_transport_operator_mismatch_certificate": {
            "certificate_id": "STRICT_FRW_BIANCHI_SYMBOLIC_TRANSPORT_OPERATOR_MISMATCH_CERTIFICATE_V1",
            "source_packet": str(IN_2203.relative_to(ROOT)),
            "symbolic_objects": {
                "delta_sq_kernel": str(delta_sq),
                "delta_sq_under_m_eq_1": str(delta_under_condition),
                "first_dm_at_m_eq_1": str(first_dm_at_1),
            },
            "residual_context": {
                "max_residual_from_p2203": max_residual,
                "residual_map_size": len(rows),
            },
            "theorem_gap_statement": theorem_gap_statement,
        },
        "recommended_next_honest_step": {
            "id": "P2205_candidate",
            "goal": "connect symbolic mismatch term to explicit nu-branch transport operator and derive computable sufficient vanishing constraints"
        },
        "gatekeeper_checks": {
            "symbolic_mismatch_certificate_exported": True,
            "symbolic_zero_mismatch_at_m_eq_1": bool(delta_under_condition == 0),
            "numeric_residual_context_present": at_one_residual is not None,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_cutkosky_closure_proven": False,
            "full_d3_covariance_transport_proven": False,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2204 S1154: symbolic transport-operator mismatch certificate",
            "",
            f"- symbolic zero condition: `{sufficient_condition}`",
            f"- delta under m=1: `{delta_under_condition}`",
            f"- first derivative term at m=1 exported: `{first_dm_at_1}`",
            f"- max residual context (P2203): `{max_residual:.12e}`",
            "",
            "Symbolic sufficient condition exported; no global Task-3 closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
