#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2086_s1036_strict_full_ltotal_eom_termwise_execution_audit.json"
MD = GEN / "p2086_s1036_strict_full_ltotal_eom_termwise_execution_audit.md"

SCHEMA_VERSION = "p2086_s1036_v1"
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
    p2073 = load("p2073_s1023_strict_full_ltotal_eom_derivation_scaffold_audit.json")
    ready = p2073.get("result_kind") == "PASS_STRICT_FULL_LTOTAL_EOM_DERIVATION_SCAFFOLD_AUDIT_WITH_TRACE__C3_STILL_OPEN"

    x = sp.Symbol("x", real=True)
    psi = sp.Function("psi")(x)
    A = sp.Function("A")(x)
    h = sp.Function("h")(x)

    kpsi, mpsi, lam4 = sp.symbols("kpsi mpsi lam4", real=True)
    kA, mA, gA = sp.symbols("kA mA gA", real=True)
    kh, mh, gh = sp.symbols("kh mh gh", real=True)
    gmix, zeta = sp.symbols("gmix zeta", real=True)

    L_terms = {
        "L_psi_kin": sp.Rational(1, 2) * kpsi * sp.diff(psi, x) ** 2,
        "L_psi_mass": -sp.Rational(1, 2) * mpsi * psi**2,
        "L_psi_self": -lam4 * psi**4,
        "L_A_kin": sp.Rational(1, 2) * kA * sp.diff(A, x) ** 2,
        "L_A_mass": -sp.Rational(1, 2) * mA * A**2,
        "L_A_self": -gA * A**4,
        "L_h_kin": sp.Rational(1, 2) * kh * sp.diff(h, x) ** 2,
        "L_h_mass": -sp.Rational(1, 2) * mh * h**2,
        "L_h_self": -gh * h**4,
        "L_mix_trilinear": -gmix * psi * A * h,
        "L_mix_biquadratic": -zeta * psi**2 * A**2,
    }
    L_total = sp.expand(sum(L_terms.values()))

    fields = {"psi": psi, "A": A, "h": h}
    termwise = {fname: {} for fname in fields}
    for tname, texpr in L_terms.items():
        for fname, fexpr in fields.items():
            termwise[fname][tname] = sp.expand(euler_lagrange(texpr, fexpr, x))

    eom_full = {fname: sp.expand(euler_lagrange(L_total, fexpr, x)) for fname, fexpr in fields.items()}
    eom_sum_of_terms = {
        fname: sp.expand(sum(termwise[fname].values())) for fname in fields
    }
    residuals = {fname: sp.expand(eom_full[fname] - eom_sum_of_terms[fname]) for fname in fields}

    subs_sample = {
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
    }
    sample_forms = {
        psi: sp.sin(x),
        A: sp.cos(2 * x),
        h: sp.exp(-x),
    }
    residual_samples = {
        fname: sp.simplify(residuals[fname].subs(subs_sample).subs(sample_forms)) for fname in fields
    }

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2086",
        "stage_id": "S1036",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_FULL_LTOTAL_EOM_TERMWISE_EXECUTION_AUDIT_WITH_TRACE__C3_STILL_OPEN" if ready else "OPEN_STRICT_FULL_LTOTAL_EOM_TERMWISE_EXECUTION_AUDIT_BLOCKED",
        "depends_on": {
            "p2073_present": p2073.get("_missing") is None,
            "preconditions_ready": ready,
        },
        "input_hashes": {
            "p2073_json_sha256": file_sha256(GEN / "p2073_s1023_strict_full_ltotal_eom_derivation_scaffold_audit.json"),
        },
        "eom_execution_protocol": {
            "coordinate": "x",
            "fields": ["psi(x)", "A(x)", "h(x)"],
            "term_count": len(L_terms),
            "strategy": "explicit termwise Euler-Lagrange decomposition + symbolic recomposition residual check",
        },
        "eom_execution_results": {
            "lagrangian_terms": {k: sp.srepr(v) for k, v in L_terms.items()},
            "full_ltotal": sp.srepr(L_total),
            "termwise_variation_map": {
                fname: {tname: sp.srepr(texpr) for tname, texpr in fmap.items()} for fname, fmap in termwise.items()
            },
            "eom_full": {fname: sp.srepr(expr) for fname, expr in eom_full.items()},
            "eom_sum_of_terms": {fname: sp.srepr(expr) for fname, expr in eom_sum_of_terms.items()},
            "symbolic_recomposition_residual": {fname: sp.srepr(expr) for fname, expr in residuals.items()},
            "numeric_probe_residual": {fname: str(expr) for fname, expr in residual_samples.items()},
            "termwise_execution_status": "TERMWISE_EXECUTED_NON_THEOREM_GRADE",
        },
        "kernel_split_context": p2073.get("kernel_split_context", {}),
        "c3_gate_update": {
            "C3_full_ltotal_eom_termwise_execution_audit": "COMPUTED",
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
            "all_symbolic_residual_zero": all(sp.simplify(v) == 0 for v in residuals.values()),
            "all_numeric_probe_residual_zero": all(sp.simplify(v) == 0 for v in residual_samples.values()),
            "term_count_expected": len(L_terms) == 11,
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
            "# P2086 S1036: strict full-L_total EOM termwise execution audit",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Term count: `{len(L_terms)}`",
            f"- Symbolic recomposition zero residual: `{payload['gatekeeper_checks']['all_symbolic_residual_zero']}`",
            f"- Numeric probe zero residual: `{payload['gatekeeper_checks']['all_numeric_probe_residual_zero']}`",
            "",
            "This stage executes full-L_total termwise symbolic variation and verifies exact recomposition.",
            "No theorem-grade closure claim is made; C3 remains OPEN.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
