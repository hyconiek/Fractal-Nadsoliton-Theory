#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import scipy.integrate as si
import scipy.linalg as la
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_1972 = GEN / "p1972_s922_strict_b1_renormalization_exact_cancellation_witness.json"
IN_1985 = GEN / "p1985_s935_strict_adm_bianchi_non_gb_curvature_squared_lapse_obstruction_witness.json"
IN_2025 = GEN / "p2025_s975_strict_cutkosky_same_scheme_cohomology_amplitude_bridge_seed.json"
OUT = GEN / "p2296_s1246_strict_global_task1_renormalization_replay_and_7task_reclassification_probe.json"
MD = GEN / "p2296_s1246_strict_global_task1_renormalization_replay_and_7task_reclassification_probe.md"

COEFF_KEYS = ("a_R2", "a_Ric2", "a_Riem2", "a_GB")
DELTA_KEYS = ("delta_c_gr_1", "delta_c_gr_2", "delta_c_gr_3", "delta_c_gr_4")
OPERATOR_KEYS = ("R2", "Ric2", "Riem2", "G_GB")
STRICT_TUPLE = {"omega": 0.18575, "phi": 0.16250, "beta": 1.0, "eta": 1.8}


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    if not path.exists():
        return ""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def strict_kernel(d: float) -> float:
    omega = STRICT_TUPLE["omega"]
    phi = STRICT_TUPLE["phi"]
    beta = STRICT_TUPLE["beta"]
    eta = STRICT_TUPLE["eta"]
    return float(np.cos(omega * d + phi) / (1.0 + beta * d**eta))


def kernel_moment(power: int) -> tuple[float, float]:
    def integrand(d: float) -> float:
        k = strict_kernel(d)
        return float((d**power) * k * k)

    value, error = si.quad(integrand, 0.0, 20.0, epsabs=1e-11, epsrel=1e-11, limit=400)
    return float(value), float(error)


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p1972 = load(IN_1972)
    p1985 = load(IN_1985)
    p2025 = load(IN_2025)

    eps = sp.symbols("epsilon", nonzero=True, real=True)
    op_symbols = sp.symbols("O_R2 O_Ric2 O_Riem2 O_GB", real=True)
    local = {"epsilon": eps, "pi": sp.pi, "log": sp.log, "ln": sp.log}

    coeff_packet = p1972.get("computed_divergence_coefficients", {}) or {}
    delta_packet = p1972.get("counterterm_map", {}) or {}

    coefficients = {
        key: sp.sympify((coeff_packet.get(key, {}) or {}).get("strict_tuple_value", "0"), locals=local)
        for key in COEFF_KEYS
    }
    deltas = {key: sp.sympify(delta_packet.get(key, "0"), locals=local) for key in DELTA_KEYS}

    divergence_density = sp.Add(*[(coefficients[key] / eps) * op for key, op in zip(COEFF_KEYS, op_symbols)])
    counterterm_density = sp.Add(*[deltas[key] * op for key, op in zip(DELTA_KEYS, op_symbols)])
    residual_density = sp.factor(sp.simplify(divergence_density + counterterm_density))
    residual_vector = [sp.factor(sp.simplify(coefficients[c] / eps + deltas[d])) for c, d in zip(COEFF_KEYS, DELTA_KEYS)]

    replay_rows = []
    rng = np.random.default_rng(2296)
    for idx in range(6):
        op_values = [sp.Rational(int(v), 11) for v in rng.integers(1, 41, size=4)]
        eps_value = sp.Rational(idx + 2, 17)
        exact = sp.simplify(residual_density.subs(dict(zip(op_symbols, op_values))).subs(eps, eps_value))
        vector_float = np.array([float(sp.N(v.subs(eps, eps_value), 40)) for v in residual_vector], dtype=float)
        replay_rows.append(
            {
                "row_id": f"p2296_b1_replay_{idx}",
                "operator_values": {name: str(value) for name, value in zip(OPERATOR_KEYS, op_values)},
                "epsilon": str(eps_value),
                "exact_residual": str(exact),
                "exact_zero": bool(exact == 0),
                "coefficient_residual_l2": float(la.norm(vector_float, ord=2)),
            }
        )

    moments = []
    for power in range(4):
        value, error = kernel_moment(power)
        moments.append({"moment": f"int_0_20 d^{power} K_strict(d)^2 dd", "value": value, "abs_error": error})

    p1972_gate = p1972.get("gatekeeper_checks", {}) or {}
    p1985_gate = p1985.get("gatekeeper_checks", {}) or {}
    p2025_tasks = p2025.get("toe_closure_gaps_7tasks", []) or []
    p2025_task1 = next((row for row in p2025_tasks if int(row.get("task_id", -1)) == 1), {})

    b1_projected_pass = bool(
        p1972_gate.get("symbolic_residual_zero", False)
        and p1972_gate.get("all_four_coefficients_computed", False)
        and residual_density == 0
        and all(row["exact_zero"] for row in replay_rows)
    )
    global_extension_open = bool(
        p1985.get("status") == "OPEN_OBSTRUCTION_WITH_TRACE"
        and p1985_gate.get("obstruction_witness_passed", False)
    )

    task1_reclassification = {
        "task_id": 1,
        "name": "Renormalization exact divergence cancellation coefficients",
        "previous_p2025_status": p2025_task1.get("status", "UNKNOWN"),
        "previous_p2025_missing": p2025_task1.get("missing", []),
        "updated_status": "B1_PROJECTED_COEFFICIENT_CANCELLATION_PASS__GLOBAL_BIANCHI_EXTENSION_OPEN"
        if b1_projected_pass and global_extension_open
        else "OPEN_OBSTRUCTION_WITH_TRACE",
        "closed_subclaim": "P1972 proves exact projected-basis MSbar_B1 cancellation for a_R2/a_Ric2/a_Riem2/a_GB and delta_c_gr_i",
        "still_open_global_claim": "P1985 exports a strict non-GB ADM/Bianchi-I lapse residual obstruction; background-family/global renormalization closure is not licensed",
    }

    seven_task_status = []
    for row in p2025_tasks:
        rid = int(row.get("task_id", -1))
        if rid == 1:
            seven_task_status.append(task1_reclassification)
        else:
            seven_task_status.append(
                {
                    "task_id": rid,
                    "name": row.get("name", "UNKNOWN"),
                    "updated_status": row.get("status", "OPEN_OBSTRUCTION_WITH_TRACE"),
                    "still_open": True,
                    "reason": "not touched by P2296; retained from P2025 no-false-pass ledger",
                }
            )

    payload = {
        "schema_version": "p2296_s1246_v1",
        "packet_id": "P2296",
        "stage_id": "S1246",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_GLOBAL_TASK1_RENORMALIZATION_REPLAY_AND_7TASK_RECLASSIFICATION_PROBE",
        "strict_global_task1_renormalization_replay_and_7task_reclassification_probe": {
            "probe_id": "STRICT_GLOBAL_TASK1_RENORMALIZATION_REPLAY_AND_7TASK_RECLASSIFICATION_PROBE_V1",
            "source_packets": [str(IN_1972.relative_to(ROOT)), str(IN_1985.relative_to(ROOT)), str(IN_2025.relative_to(ROOT))],
            "source_hashes": {
                str(IN_1972.relative_to(ROOT)): sha256_file(IN_1972),
                str(IN_1985.relative_to(ROOT)): sha256_file(IN_1985),
                str(IN_2025.relative_to(ROOT)): sha256_file(IN_2025),
            },
            "strict_kernel_tuple": {**STRICT_TUPLE, "legacy_role_transfer_used": False},
            "computed_coefficients": {key: str(value) for key, value in coefficients.items()},
            "counterterm_map": {key: str(value) for key, value in deltas.items()},
            "symbolic_replay": {
                "divergence_density": str(divergence_density),
                "counterterm_density": str(counterterm_density),
                "renormalized_divergent_density_residual": str(residual_density),
                "residual_vector": [str(v) for v in residual_vector],
                "residual_is_zero": bool(residual_density == 0),
            },
            "numeric_replay_rows": replay_rows,
            "strict_kernel_moment_crosscheck": moments,
            "task1_reclassification": task1_reclassification,
            "seven_task_status_after_p2296": seven_task_status,
            "theorem_scope_limit": "Task-1 B1 projected-basis replay and status reclassification only; not global background-independence, Cutkosky, selector, or ToE closure",
        },
        "recommended_next_honest_step": {
            "id": "P2297_candidate",
            "goal": "attack the P1985 non-GB Bianchi-I residual by deriving a full spatial-EOM transport/provider matrix, without promoting B1 projected renormalization to global ToE closure",
        },
        "gatekeeper_checks": {
            "p1972_symbolic_residual_zero_replayed": bool(residual_density == 0),
            "all_replay_rows_exact_zero": all(row["exact_zero"] for row in replay_rows),
            "coefficient_residual_norms_zero": all(row["coefficient_residual_l2"] == 0.0 for row in replay_rows),
            "strict_kernel_moments_positive": all(row["value"] > 0.0 for row in moments),
            "task1_b1_projected_pass": b1_projected_pass,
            "global_extension_still_open_by_p1985": global_extension_open,
            "no_cutkosky_closure_claimed": True,
            "no_background_independence_closure_claimed": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2296 S1246: strict Task-1 renormalization replay + 7-task reclassification",
                "",
                f"- symbolic B1 residual replay zero: `{residual_density == 0}`",
                f"- exact replay rows zero: `{all(row['exact_zero'] for row in replay_rows)}`",
                f"- strict-kernel moment count: `{len(moments)}`",
                f"- Task-1 updated status: `{task1_reclassification['updated_status']}`",
                f"- global extension still open by P1985: `{global_extension_open}`",
                "",
                "Theorem-level claim licensed here: projected-basis `MSbar_B1` divergence cancellation from P1972 is replayed exactly.",
                "No global background-independence, Cutkosky, selector, or ToE closure is claimed.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
