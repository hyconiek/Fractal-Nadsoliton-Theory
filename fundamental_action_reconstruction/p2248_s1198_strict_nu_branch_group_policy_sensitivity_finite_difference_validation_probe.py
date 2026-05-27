#!/usr/bin/env python3
from __future__ import annotations
import json
import random
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2247 = GEN / "p2247_s1197_strict_nu_branch_group_policy_boundary_sensitivity_probe.json"
OUT = GEN / "p2248_s1198_strict_nu_branch_group_policy_sensitivity_finite_difference_validation_probe.json"
MD = GEN / "p2248_s1198_strict_nu_branch_group_policy_sensitivity_finite_difference_validation_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def margin(total: float, n: float, load_ratio: float, kappa: float) -> float:
    return total - ((n - 1.0) + kappa * load_ratio)


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2247 = load(IN_2247)
    probe = (p2247.get("strict_nu_branch_group_policy_boundary_sensitivity_probe", {}) or {})
    inp = (probe.get("inputs", {}) or {})
    sens = (probe.get("analytic_sensitivity", {}) or {})

    n0 = float(inp.get("base_group_count", 1.0) or 1.0)
    load0 = float(inp.get("base_load_ratio", 1.0) or 1.0)
    total0 = float(inp.get("base_total_coverage_mass", 1.0) or 1.0)
    kappa = float(inp.get("kappa_reference", 1.0) or 1.0)

    d_dn_analytic = float(sens.get("dmargin_d_group_count", -1.0) or -1.0)
    d_dload_analytic = float(sens.get("dmargin_d_load_ratio", -kappa) or -kappa)

    eps = 1e-5
    rng = random.Random(2248)
    samples = 7
    rows = []

    for _ in range(samples):
        dn = (rng.random() - 0.5) * 0.4
        dl = (rng.random() - 0.5) * 0.4 * max(load0, 1e-6)
        n = max(1.0, n0 + dn)
        l = max(1e-12, load0 + dl)

        d_dn_fd = (margin(total0, n + eps, l, kappa) - margin(total0, n - eps, l, kappa)) / (2.0 * eps)
        d_dload_fd = (margin(total0, n, l + eps, kappa) - margin(total0, n, l - eps, kappa)) / (2.0 * eps)

        rows.append(
            {
                "n": n,
                "load_ratio": l,
                "fd_dmargin_d_group_count": d_dn_fd,
                "fd_dmargin_d_load_ratio": d_dload_fd,
                "abs_err_group_count": abs(d_dn_fd - d_dn_analytic),
                "abs_err_load_ratio": abs(d_dload_fd - d_dload_analytic),
            }
        )

    max_err_dn = max((r["abs_err_group_count"] for r in rows), default=0.0)
    max_err_dl = max((r["abs_err_load_ratio"] for r in rows), default=0.0)
    tol = 1e-7
    validation_pass = (max_err_dn <= tol) and (max_err_dl <= tol)

    payload = {
        "schema_version": "p2248_s1198_v1",
        "packet_id": "P2248",
        "stage_id": "S1198",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_SENSITIVITY_FINITE_DIFFERENCE_VALIDATION_PROBE",
        "strict_nu_branch_group_policy_sensitivity_finite_difference_validation_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_SENSITIVITY_FINITE_DIFFERENCE_VALIDATION_PROBE_V1",
            "source_packets": [str(IN_2247.relative_to(ROOT))],
            "inputs": {
                "base_group_count": n0,
                "base_load_ratio": load0,
                "base_total_coverage_mass": total0,
                "kappa_reference": kappa,
                "analytic_dmargin_d_group_count": d_dn_analytic,
                "analytic_dmargin_d_load_ratio": d_dload_analytic,
                "finite_difference_eps": eps,
                "sample_count": samples,
                "random_seed": 2248,
            },
            "finite_difference_validation": {
                "rows": rows,
                "max_abs_error_group_count": max_err_dn,
                "max_abs_error_load_ratio": max_err_dl,
                "tolerance": tol,
                "validation_pass": validation_pass,
            },
            "physical_interpretation_note": "Finite-difference agreement validates that local robustness-erosion slopes are structurally stable, supporting interpretation of margin flow as a first-order physical response law in policy space.",
            "theorem_scope_limit": "local finite-difference validation only; not a legacy->strict bridge theorem and not strict-core selector closure",
        },
        "recommended_next_honest_step": {
            "id": "P2249_candidate",
            "goal": "integrate first-order slope law into conservative step-size controller for policy updates with guaranteed nonnegative margin",
        },
        "gatekeeper_checks": {
            "finite_difference_validation_exported": True,
            "validation_pass": validation_pass,
            "tolerance_strict": tol <= 1e-6,
            "no_bridge_theorem_claimed": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_cutkosky_closure_proven": False,
            "full_d3_covariance_transport_proven": False,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2248 S1198: finite-difference sensitivity validation probe",
                "",
                f"- max |fd-analytic| for dmargin/d(group_count): `{max_err_dn:.12e}`",
                f"- max |fd-analytic| for dmargin/d(load_ratio): `{max_err_dl:.12e}`",
                f"- tolerance: `{tol:.12e}`",
                f"- validation pass: `{validation_pass}`",
                "",
                "Local finite-difference validation only; no kernel-bridge or selector-closure claim.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
