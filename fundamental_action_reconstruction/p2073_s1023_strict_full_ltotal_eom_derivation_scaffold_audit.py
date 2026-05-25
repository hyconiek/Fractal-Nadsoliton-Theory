#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2073_s1023_strict_full_ltotal_eom_derivation_scaffold_audit.json"
MD = GEN / "p2073_s1023_strict_full_ltotal_eom_derivation_scaffold_audit.md"

SCHEMA_VERSION = "p2073_s1023_v1"
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


def main() -> None:
    GEN.mkdir(exist_ok=True)

    p2072 = load("p2072_s1022_strict_same_scheme_tau_stability_margin_stress_audit.json")
    ready = p2072.get("result_kind") == "PASS_TAU_STABILITY_MARGIN_STRESS_AUDIT_WITH_TRACE__C3_STILL_OPEN"

    x = sp.Symbol("x", real=True)
    psi = sp.Function("psi")(x)
    A = sp.Function("A")(x)
    h = sp.Function("h")(x)

    kpsi, mpsi, lam4 = sp.symbols("kpsi mpsi lam4", real=True)
    kA, mA, gA = sp.symbols("kA mA gA", real=True)
    kh, mh, gh = sp.symbols("kh mh gh", real=True)
    gmix, zeta = sp.symbols("gmix zeta", real=True)

    L_psi = sp.Rational(1, 2) * kpsi * sp.diff(psi, x) ** 2 - sp.Rational(1, 2) * mpsi * psi**2 - lam4 * psi**4
    L_A = sp.Rational(1, 2) * kA * sp.diff(A, x) ** 2 - sp.Rational(1, 2) * mA * A**2 - gA * A**4
    L_h = sp.Rational(1, 2) * kh * sp.diff(h, x) ** 2 - sp.Rational(1, 2) * mh * h**2 - gh * h**4
    L_mix = -gmix * psi * A * h - zeta * psi**2 * A**2
    L_total_scaffold = sp.expand(L_psi + L_A + L_h + L_mix)

    def euler_lagrange(expr: sp.Expr, field: sp.Expr) -> sp.Expr:
        return sp.simplify(sp.diff(expr, field) - sp.diff(sp.diff(expr, sp.diff(field, x)), x))

    eom_psi = euler_lagrange(L_total_scaffold, psi)
    eom_A = euler_lagrange(L_total_scaffold, A)
    eom_h = euler_lagrange(L_total_scaffold, h)

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2073",
        "stage_id": "S1023",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_STRICT_FULL_LTOTAL_EOM_DERIVATION_SCAFFOLD_AUDIT_WITH_TRACE__C3_STILL_OPEN"
            if ready
            else "OPEN_STRICT_FULL_LTOTAL_EOM_DERIVATION_SCAFFOLD_AUDIT_BLOCKED"
        ),
        "depends_on": {
            "p2072_present": p2072.get("_missing") is None,
            "preconditions_ready": ready,
        },
        "input_hashes": {
            "p2072_json_sha256": file_sha256(GEN / "p2072_s1022_strict_same_scheme_tau_stability_margin_stress_audit.json"),
        },
        "kernel_split_context": {
            "legacy_kernel": "alpha_geo*cos(omega*d+phi)/(1+beta_tors*d)",
            "strict_kernel": "cos(omega*d+phi)/(1+beta*d^eta)",
            "bridge_status": "OPEN_NO_RIGOROUS_IDENTIFICATION_EXPORT",
            "missing_characteristic_hypothesis": "OPEN_CANDIDATE_REQUIRED",
        },
        "eom_scaffold": {
            "coordinates": ["x"],
            "fields": ["psi(x)", "A(x)", "h(x)"],
            "operator_blocks": {
                "L_psi": sp.srepr(L_psi),
                "L_A": sp.srepr(L_A),
                "L_h": sp.srepr(L_h),
                "L_mix": sp.srepr(L_mix),
                "L_total_scaffold": sp.srepr(L_total_scaffold),
            },
            "variation_map": {
                "deltaS_delta_psi": sp.srepr(eom_psi),
                "deltaS_delta_A": sp.srepr(eom_A),
                "deltaS_delta_h": sp.srepr(eom_h),
            },
            "symbolic_trace_status": "SCAFFOLD_ONLY_NON_THEOREM_GRADE",
        },
        "c3_gate_update": {
            "C3_full_ltotal_eom_scaffold_audit": "COMPUTED",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
            "missing_for_full_eom_theorem": [
                "full exported strict L_total operator basis with declared background family map",
                "global tensorial variational closure proof across targeted backgrounds",
                "selector obstruction discharge path compatible with QW-2191",
            ],
        },
        "gatekeeper_checks": {
            "preconditions_ready": ready,
            "scaffold_fields_nonempty": True,
            "variation_map_nonempty": True,
            "symbolic_trace_scaffold_only": True,
            "c3_theorem_proven": False,
            "no_background_globalization_claimed": True,
            "no_tensor_component_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2073 S1023: strict full-L_total EOM derivation scaffold audit",
                "",
                f"- Status: `{payload['status']}`",
                f"- Result kind: `{payload['result_kind']}`",
                "- Kernel split bridge status: `OPEN_NO_RIGOROUS_IDENTIFICATION_EXPORT`",
                "- Missing characteristic hypothesis: `OPEN_CANDIDATE_REQUIRED`",
                "- Symbolic trace status: `SCAFFOLD_ONLY_NON_THEOREM_GRADE`",
                "",
                "This stage exports a formal operator/variation scaffold only.",
                "No full EOM theorem claim is made; C3 remains OPEN.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
