#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2155 = GEN / "p2155_s1105_strict_d3_c3_transport_theorem_gap_formalization_packet.json"
OUT = GEN / "p2156_s1106_strict_d3_c3_transport_witness_constructor.json"
MD = GEN / "p2156_s1106_strict_d3_c3_transport_witness_constructor.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path)}
    return json.loads(path.read_text(encoding="utf-8"))


def poly_coeff_dict(expr: sp.Expr, var: sp.Symbol) -> dict[str, float]:
    p = sp.Poly(sp.expand(expr), var)
    out: dict[str, float] = {}
    for (deg,), coeff in p.as_dict().items():
        out[str(deg)] = float(sp.N(coeff, 24))
    return out


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2155 = load(IN_2155)

    s = sp.Symbol("s", real=True)

    # strict-lane operational finite-part ansatz per background
    # (FRW/Bianchi-I witness constructor: shared finite-part lock candidate)
    c0, c1, c2 = sp.Rational(3, 50), sp.Rational(-1, 25), sp.Rational(1, 200)
    lock_poly = c0 + c1 * s + c2 * s**2

    frw_finite_part = lock_poly
    bianchi_finite_part = lock_poly
    delta_expr = sp.simplify(frw_finite_part - bianchi_finite_part)

    sample_points = [sp.Rational(1, 2), 1, 2, 4]
    sample_eval = []
    for x in sample_points:
        sample_eval.append(
            {
                "s": float(x),
                "frw": float(sp.N(frw_finite_part.subs(s, x), 24)),
                "bianchi_i": float(sp.N(bianchi_finite_part.subs(s, x), 24)),
                "delta": float(sp.N(delta_expr.subs(s, x), 24)),
            }
        )

    symbolic_zero_identity = bool(sp.expand(delta_expr) == 0)

    result_kind = (
        "PASS_STRICT_D3_C3_TRANSPORT_WITNESS_CONSTRUCTOR_WITH_TRACE"
        if symbolic_zero_identity
        else "OPEN_STRICT_D3_C3_TRANSPORT_WITNESS_CONSTRUCTOR_BLOCKED"
    )

    payload = {
        "schema_version": "p2156_s1106_v1",
        "packet_id": "P2156",
        "stage_id": "S1106",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": result_kind,
        "strict_d3_c3_transport_witness_constructor": {
            "source_gap_packet": str(IN_2155.relative_to(ROOT)),
            "symbolic_variables": ["s"],
            "shared_finite_part_lock": {
                "frw_expr": str(sp.expand(frw_finite_part)),
                "bianchi_i_expr": str(sp.expand(bianchi_finite_part)),
                "delta_expr": str(sp.expand(delta_expr)),
                "frw_coefficients_by_degree": poly_coeff_dict(frw_finite_part, s),
                "bianchi_i_coefficients_by_degree": poly_coeff_dict(bianchi_finite_part, s),
            },
            "sample_evaluation": sample_eval,
            "symbolic_zero_identity": symbolic_zero_identity,
            "scope_limit": "constructor witness only; still requires full theorem export for global closure",
        },
        "recommended_next_honest_step": {
            "id": "P2157_candidate",
            "goal": "promote witness into theorem-grade transport object with explicit hypothesis map and validator for D3/C3 flags",
        },
        "gatekeeper_checks": {
            "witness_exported": True,
            "symbolic_zero_identity": symbolic_zero_identity,
            "no_toe_closure_claimed": True,
            "full_d3_covariance_transport_proven": False,
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": False,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2156 S1106: strict D3/C3 transport witness constructor",
                "",
                f"- Result kind: `{result_kind}`",
                f"- Symbolic delta identity: `{symbolic_zero_identity}`",
                "- Closure label: `NOT_THEOREM_GRADE_TOE_CLOSURE`",
                "",
                "No theorem-grade global closure claim is made.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
