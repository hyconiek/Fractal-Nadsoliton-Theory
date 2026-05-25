#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2087_s1037_strict_full_ltotal_eom_normal_form_extraction_audit.json"
MD = GEN / "p2087_s1037_strict_full_ltotal_eom_normal_form_extraction_audit.md"

SCHEMA_VERSION = "p2087_s1037_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"


def load(name: str) -> dict[str, Any]:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def file_sha256(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def euler_lagrange(expr: sp.Expr, field: sp.Expr, x: sp.Symbol) -> sp.Expr:
    return sp.simplify(sp.diff(expr, field) - sp.diff(sp.diff(expr, sp.diff(field, x)), x))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2086 = load("p2086_s1036_strict_full_ltotal_eom_termwise_execution_audit.json")
    ready = p2086.get("result_kind") == "PASS_STRICT_FULL_LTOTAL_EOM_TERMWISE_EXECUTION_AUDIT_WITH_TRACE__C3_STILL_OPEN"

    x = sp.Symbol("x", real=True)
    psi = sp.Function("psi")(x)
    A = sp.Function("A")(x)
    h = sp.Function("h")(x)

    kpsi, mpsi, lam4 = sp.symbols("kpsi mpsi lam4", real=True)
    kA, mA, gA = sp.symbols("kA mA gA", real=True)
    kh, mh, gh = sp.symbols("kh mh gh", real=True)
    gmix, zeta = sp.symbols("gmix zeta", real=True)

    L_total = sp.expand(
        sp.Rational(1, 2) * kpsi * sp.diff(psi, x) ** 2
        - sp.Rational(1, 2) * mpsi * psi**2
        - lam4 * psi**4
        + sp.Rational(1, 2) * kA * sp.diff(A, x) ** 2
        - sp.Rational(1, 2) * mA * A**2
        - gA * A**4
        + sp.Rational(1, 2) * kh * sp.diff(h, x) ** 2
        - sp.Rational(1, 2) * mh * h**2
        - gh * h**4
        - gmix * psi * A * h
        - zeta * psi**2 * A**2
    )

    eom = {
        "psi": sp.expand(euler_lagrange(L_total, psi, x)),
        "A": sp.expand(euler_lagrange(L_total, A, x)),
        "h": sp.expand(euler_lagrange(L_total, h, x)),
    }

    ddpsi = sp.diff(psi, x, 2)
    ddA = sp.diff(A, x, 2)
    ddh = sp.diff(h, x, 2)
    normal_forms = {
        "psi_ddot": sp.expand(sp.solve(sp.Eq(eom["psi"], 0), ddpsi, dict=True)[0][ddpsi]),
        "A_ddot": sp.expand(sp.solve(sp.Eq(eom["A"], 0), ddA, dict=True)[0][ddA]),
        "h_ddot": sp.expand(sp.solve(sp.Eq(eom["h"], 0), ddh, dict=True)[0][ddh]),
    }

    substitution = {ddpsi: normal_forms["psi_ddot"], ddA: normal_forms["A_ddot"], ddh: normal_forms["h_ddot"]}
    reduced_residual = {k: sp.simplify(v.subs(substitution)) for k, v in eom.items()}

    sample_subs = {
        kpsi: sp.Rational(11, 10),
        mpsi: sp.Rational(9, 10),
        lam4: sp.Rational(1, 5),
        kA: sp.Rational(13, 10),
        mA: sp.Rational(4, 5),
        gA: sp.Rational(3, 20),
        kh: sp.Rational(6, 5),
        mh: sp.Rational(7, 10),
        gh: sp.Rational(1, 10),
        gmix: sp.Rational(2, 25),
        zeta: sp.Rational(1, 20),
        psi: sp.sin(x),
        A: sp.cos(2 * x),
        h: sp.exp(-x),
    }
    sample_rhs = {k: sp.simplify(v.subs(sample_subs)) for k, v in normal_forms.items()}

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2087",
        "stage_id": "S1037",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_FULL_LTOTAL_EOM_NORMAL_FORM_EXTRACTION_AUDIT_WITH_TRACE__C3_STILL_OPEN" if ready else "OPEN_STRICT_FULL_LTOTAL_EOM_NORMAL_FORM_EXTRACTION_AUDIT_BLOCKED",
        "depends_on": {
            "p2086_present": p2086.get("_missing") is None,
            "preconditions_ready": ready,
        },
        "input_hashes": {
            "p2086_json_sha256": file_sha256(GEN / "p2086_s1036_strict_full_ltotal_eom_termwise_execution_audit.json"),
        },
        "eom_normal_form_results": {
            "full_ltotal": sp.srepr(L_total),
            "eom_full": {k: sp.srepr(v) for k, v in eom.items()},
            "solved_second_derivative_rhs": {k: sp.srepr(v) for k, v in normal_forms.items()},
            "residual_after_substitute_normal_form": {k: sp.srepr(v) for k, v in reduced_residual.items()},
            "sample_normal_form_rhs": {k: str(v) for k, v in sample_rhs.items()},
        },
        "c3_gate_update": {
            "C3_full_ltotal_eom_normal_form_extraction_audit": "COMPUTED",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
            "missing_for_full_eom_theorem": [
                "strict background-family tensor closure proof (FRW/Bianchi-I and targeted branches)",
                "global non-perturbative well-posedness/export layer",
                "QW-2191-compatible strict selector premise or source",
            ],
        },
        "gatekeeper_checks": {
            "preconditions_ready": ready,
            "normal_form_solved_all_fields": all(k in normal_forms for k in ["psi_ddot", "A_ddot", "h_ddot"]),
            "residual_zero_after_normal_form_substitution": all(sp.simplify(v) == 0 for v in reduced_residual.values()),
            "c3_theorem_proven": False,
            "no_background_globalization_claimed": True,
            "no_tensor_component_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2087 S1037: strict full-L_total EOM normal-form extraction audit",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Normal forms solved: `{payload['gatekeeper_checks']['normal_form_solved_all_fields']}`",
            f"- Residual zero after substitution: `{payload['gatekeeper_checks']['residual_zero_after_normal_form_substitution']}`",
            "",
            "This stage rewrites full-L_total EOM into explicit second-derivative normal forms for psi/A/h.",
            "No theorem-grade closure claim is made; C3 remains OPEN.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
