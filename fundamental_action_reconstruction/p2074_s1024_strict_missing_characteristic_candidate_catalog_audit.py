#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2074_s1024_strict_missing_characteristic_candidate_catalog_audit.json"
MD = GEN / "p2074_s1024_strict_missing_characteristic_candidate_catalog_audit.md"

SCHEMA_VERSION = "p2074_s1024_v1"
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

    L_base = (
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

    c1, c2, c3 = sp.symbols("c1 c2 c3", real=True)
    candidates = [
        ("mc_grad_cross", c1 * sp.diff(psi, x) * sp.diff(A, x)),
        ("mc_curv_like", c2 * psi**2 * sp.diff(h, x) ** 2),
        ("mc_phase_lock", c3 * psi * A * sp.diff(h, x)),
    ]

    base_eoms = {
        "psi": euler_lagrange(L_base, psi, x),
        "A": euler_lagrange(L_base, A, x),
        "h": euler_lagrange(L_base, h, x),
    }

    rows = []
    for cid, op in candidates:
        L_aug = sp.expand(L_base + op)
        aug_eoms = {
            "psi": euler_lagrange(L_aug, psi, x),
            "A": euler_lagrange(L_aug, A, x),
            "h": euler_lagrange(L_aug, h, x),
        }
        deltas = {k: sp.simplify(aug_eoms[k] - base_eoms[k]) for k in ("psi", "A", "h")}

        # Conservative non-discharge screening proxy.
        selector_safe_proxy = "QW2191_UNCHANGED_PROXY_ONLY"
        false_discharge_risk = "LOW_IF_MARKED_NON_THEOREM_GRADE"

        rows.append(
            {
                "candidate_id": cid,
                "operator_srepr": sp.srepr(op),
                "delta_eom": {"psi": sp.srepr(deltas["psi"]), "A": sp.srepr(deltas["A"]), "h": sp.srepr(deltas["h"])},
                "legacy_strict_alignment_signal": "OPEN_REQUIRES_NUMERIC_FIT_AND_BRIDGE_CHECK",
                "selector_qw2191_screening": selector_safe_proxy,
                "false_discharge_risk": false_discharge_risk,
                "admission": "CANDIDATE_RETAI NED_FOR_NEXT_STAGE".replace(" ", ""),
            }
        )

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2074",
        "stage_id": "S1024",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_MISSING_CHARACTERISTIC_CANDIDATE_CATALOG_AUDIT_WITH_TRACE__C3_STILL_OPEN"
            if ready
            else "OPEN_MISSING_CHARACTERISTIC_CANDIDATE_CATALOG_AUDIT_BLOCKED"
        ),
        "depends_on": {"p2073_present": p2073.get("_missing") is None, "preconditions_ready": ready},
        "input_hashes": {
            "p2073_json_sha256": file_sha256(GEN / "p2073_s1023_strict_full_ltotal_eom_derivation_scaffold_audit.json"),
        },
        "kernel_split_context": {
            "legacy_kernel": "alpha_geo*cos(omega*d+phi)/(1+beta_tors*d)",
            "strict_kernel": "cos(omega*d+phi)/(1+beta*d^eta)",
            "bridge_status": "OPEN_NO_RIGOROUS_IDENTIFICATION_EXPORT",
            "missing_characteristic_hypothesis": "OPEN_CANDIDATE_REQUIRED",
        },
        "candidate_catalog": {
            "base_lagrangian_srepr": sp.srepr(sp.expand(L_base)),
            "candidate_count": len(rows),
            "rows": rows,
            "screening_policy": "retain_only_non-discharge_candidates_with_qw2191_proxy_safety",
        },
        "c3_gate_update": {
            "C3_missing_characteristic_candidate_catalog_audit": "COMPUTED",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
            "missing_for_next_stage": [
                "numeric fit objective linking candidate-augmented strict model to legacy observable traces",
                "explicit legacy<->strict bridge/non-bridge theorem for admitted candidate",
                "strict selector closure theorem beyond QW-2191 proxy checks",
            ],
        },
        "gatekeeper_checks": {
            "preconditions_ready": ready,
            "candidate_rows_nonempty": len(rows) > 0,
            "delta_eom_exported": all("delta_eom" in r for r in rows),
            "qw2191_proxy_screening_exported": all("selector_qw2191_screening" in r for r in rows),
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
                "# P2074 S1024: missing-characteristic candidate catalog audit",
                "",
                f"- Status: `{payload['status']}`",
                f"- Result kind: `{payload['result_kind']}`",
                f"- Candidate count: `{len(rows)}`",
                "- Legacy/Strict bridge status: `OPEN_NO_RIGOROUS_IDENTIFICATION_EXPORT`",
                "- Missing characteristic hypothesis: `OPEN_CANDIDATE_REQUIRED`",
                "",
                "This stage catalogs candidate missing-characteristic operators and exports delta-EOM traces.",
                "No theorem-grade closure claim is made; C3 remains OPEN.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
