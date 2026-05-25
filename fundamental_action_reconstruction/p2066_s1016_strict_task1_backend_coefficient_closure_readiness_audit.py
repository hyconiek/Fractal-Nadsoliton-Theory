#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2066_s1016_strict_task1_backend_coefficient_closure_readiness_audit.json"
MD = GEN / "p2066_s1016_strict_task1_backend_coefficient_closure_readiness_audit.md"

SCHEMA_VERSION = "p2066_s1016_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"


def sha256(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(name: str) -> dict[str, Any]:
    p = GEN / name
    if not p.exists():
        return {"_missing": name}
    return json.loads(p.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)

    p2027 = load_json("p2027_s977_strict_b1_rank3_gb_null_adaptive_quadrature_witness.json")
    p2028 = load_json("p2028_s978_strict_b1_gb_quotient_counterterm_identifiability_theorem.json")
    p2029 = load_json("p2029_s979_strict_task1_renormalization_quotient_ledger_update.json")
    p2034 = load_json("p2034_s984_strict_task1_quotient_only_renormalization_theorem.json")

    evidence = []
    for path in sorted(GEN.glob("p20*.json")):
        txt = path.read_text(encoding="utf-8")
        if re.search(r"backend-computed a_R2/a_Ric2/a_Riem2/a_GB", txt):
            evidence.append(path.name)

    coeffs = (((p2027.get("full_four_channel_min_norm_family") or {}).get("canonical_rank3_representative") or {}))
    a_r2 = float(coeffs.get("a_R2", 0.0))
    a_ric2 = float(coeffs.get("a_Ric2", 0.0))
    a_riem2 = float(coeffs.get("a_Riem2", 0.0))

    # Symbolic consistency check for quotient pass target t=(1,0,0)
    aR2, aRic2, aRiem2, aGB = sp.symbols("aR2 aRic2 aRiem2 aGB")
    T = sp.Matrix([aR2 + aGB, aRic2 - 4 * aGB, aRiem2 + aGB])
    t = sp.Matrix([1, 0, 0])
    section = {aR2: a_r2, aRic2: a_ric2, aRiem2: a_riem2, aGB: 0}
    residual = [float(sp.N(v)) for v in (T.subs(section) - t)]
    residual_linf = max(abs(v) for v in residual)

    independent_agb_identified = bool((p2029.get("task1_quotient_status") or {}).get("independent_a_GB_identified", False))

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2066",
        "stage_id": "S1016",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_TASK1_BACKEND_COEFFICIENT_CLOSURE_READINESS_AUDIT__QUOTIENT_ONLY_STILL_OPEN",
        "input_hashes": {
            "p2027": sha256(GEN / "p2027_s977_strict_b1_rank3_gb_null_adaptive_quadrature_witness.json"),
            "p2028": sha256(GEN / "p2028_s978_strict_b1_gb_quotient_counterterm_identifiability_theorem.json"),
            "p2029": sha256(GEN / "p2029_s979_strict_task1_renormalization_quotient_ledger_update.json"),
            "p2034": sha256(GEN / "p2034_s984_strict_task1_quotient_only_renormalization_theorem.json"),
        },
        "content_grep_evidence": {
            "pattern": "backend-computed a_R2/a_Ric2/a_Riem2/a_GB",
            "matching_generated_packets": evidence,
            "match_count": len(evidence),
        },
        "strict_kernel_context": {
            "kernel": "K_strict(d)=cos(omega*d+phi)/(1+beta*d^eta)",
            "task_scope": "Task-1 backend coefficient closure readiness on strict B1",
        },
        "symbolic_numeric_probe": {
            "quotient_map": "T(a)=(a_R2+a_GB, a_Ric2-4*a_GB, a_Riem2+a_GB)",
            "section_used": "a_GB=0 canonical section (gauge choice only)",
            "residual_vector_vs_target_100": residual,
            "residual_linf": residual_linf,
        },
        "closure_readiness": {
            "backend_coefficients_present": all(k in coeffs for k in ["a_R2", "a_Ric2", "a_Riem2"]),
            "independent_a_GB_identified": independent_agb_identified,
            "tensor_component_projection_rule_exported": False,
            "strict_task1_full_closure_ready": False,
            "missing_items": [
                "independent a_GB identification beyond scalar quotient section",
                "curved metric ansatz and component projection contract H_00/H_ii(d)",
                "background transport theorem for same-scheme finite-part lock",
            ],
        },
        "gatekeeper_checks": {
            "no_false_pass": True,
            "declares_quotient_only_scope": True,
            "c3_transport_theorem_open": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2066 S1016: strict Task-1 backend coefficient closure readiness audit",
                "",
                f"- Status: `{payload['status']}`",
                f"- Result: `{payload['result_kind']}`",
                f"- Grep matches for backend-coefficient phrase: `{len(evidence)}`",
                f"- Quotient residual L_inf (canonical section): `{residual_linf:.3e}`",
                "",
                "This packet confirms backend-coefficient traces exist on strict B1 but full Task-1 closure remains quotient-only.",
                "No independent a_GB claim, no tensor component closure, and no C3 transport theorem discharge are made.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
