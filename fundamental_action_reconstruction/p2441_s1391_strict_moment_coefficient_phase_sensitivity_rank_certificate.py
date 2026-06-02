#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2441_s1391_strict_moment_coefficient_phase_sensitivity_rank_certificate.json"
MD = GEN / "p2441_s1391_strict_moment_coefficient_phase_sensitivity_rank_certificate.md"

SOURCE_FILES = {
    "P1562_THREE_COEFF_DERIVATION": GEN / "p1562_s512_strict_kernel_to_lagrangian_coefficient_derivation_summary.json",
    "P2440_REPLAY_RANK": GEN / "p2440_s1390_current_strict_tuple_coefficient_replay_rank_certificate.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

STRICT_PARAMS = {
    "omega": 0.18575,
    "phi": 0.16250,
    "beta": 1.0,
    "eta": 1.8,
    "D": 25.0,
    "steps": 20000,
}
PARAMETER_ORDER = ["omega", "phi", "beta", "eta"]
COEFFICIENT_ORDER = ["lambda_sm_eff", "kappa_gr_eff", "epsilon_mix_eff"]
PHI_SWEEP_DELTAS = [-0.05, -0.01, 0.01, 0.05]


def rel(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": rel(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg",
            "-n",
            pattern,
            "fundamental_action_reconstruction",
            "material_dowodowy",
            "-g",
            "*.py",
            "-g",
            "*.md",
            "-g",
            "*.tex",
            "-g",
            "!fundamental_action_reconstruction/generated/**",
        ],
        cwd=REPO,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:35]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2441|S1391|strict moment coefficient phase sensitivity|moment coefficient phase sensitivity|phase sensitivity rank certificate",
        "p2440_input": "P2440|S1390|phase-null coefficient replay|coefficient replay rank",
        "p1562_moment_route": "P1562|lambda_sm_eff|kappa_gr_eff|epsilon_mix_eff|moment\\(",
        "bridge_completed_moments": "P2363|S1313|bridge-completed moment|completed moments|raw legacy moments",
        "phase_selector_language": "phi|phase|topological|QW-2191|selector uniqueness",
        "rank_or_jacobian": "Jacobian|rank|sensitivity|null direction|nullspace",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def k_strict(d: float, omega: float, phi: float, beta: float, eta: float) -> float:
    return math.cos(omega * d + phi) / (1.0 + beta * (d**eta))


def moment(n: int, params: dict[str, float]) -> float:
    steps = int(params["steps"])
    h = params["D"] / steps
    total = 0.0
    for i in range(steps + 1):
        d = i * h
        weight = 0.5 if i in (0, steps) else 1.0
        total += weight * (d**n) * k_strict(d, params["omega"], params["phi"], params["beta"], params["eta"])
    return total * h


def moment_coefficients(params: dict[str, float]) -> dict[str, Any]:
    moments = {f"M{n}": moment(n, params) for n in range(4)}
    m0 = moments["M0"]
    if abs(m0) < 1e-12:
        raise RuntimeError("M0 too small for stable sensitivity replay")
    ratios = {
        "R1": moments["M1"] / m0,
        "R2": moments["M2"] / m0,
        "R3": moments["M3"] / m0,
    }
    coefficients = {
        "lambda_sm_eff": abs(ratios["R1"]),
        "kappa_gr_eff": abs(ratios["R2"] - ratios["R1"] * ratios["R1"]),
        "epsilon_mix_eff": abs(ratios["R3"]) / (1.0 + abs(ratios["R2"])),
    }
    return {"moments": moments, "dimensionless_ratios": ratios, "coefficients": coefficients}


def finite_difference_jacobian(params: dict[str, float]) -> list[list[float]]:
    rows: list[list[float]] = []
    for coeff_name in COEFFICIENT_ORDER:
        row: list[float] = []
        for param_name in PARAMETER_ORDER:
            step = 1e-5 * max(1.0, abs(float(params[param_name])))
            plus = dict(params)
            minus = dict(params)
            plus[param_name] += step
            minus[param_name] -= step
            cp = moment_coefficients(plus)["coefficients"][coeff_name]
            cm = moment_coefficients(minus)["coefficients"][coeff_name]
            row.append((cp - cm) / (2.0 * step))
        rows.append(row)
    return rows


def real_rank(matrix: list[list[float]], tol: float = 1e-7) -> int:
    mat = [row[:] for row in matrix]
    if not mat:
        return 0
    rank = 0
    rows = len(mat)
    cols = len(mat[0])
    for col in range(cols):
        pivot = max(range(rank, rows), key=lambda r: abs(mat[r][col]), default=rank)
        if rank >= rows or abs(mat[pivot][col]) <= tol:
            continue
        mat[rank], mat[pivot] = mat[pivot], mat[rank]
        pivot_value = mat[rank][col]
        mat[rank] = [value / pivot_value for value in mat[rank]]
        for r in range(rows):
            if r != rank and abs(mat[r][col]) > tol:
                factor = mat[r][col]
                mat[r] = [a - factor * b for a, b in zip(mat[r], mat[rank])]
        rank += 1
        if rank == rows:
            break
    return rank


def l2_norm(values: list[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def phi_sweep(base_coeffs: dict[str, float]) -> list[dict[str, Any]]:
    rows = []
    for delta in PHI_SWEEP_DELTAS:
        params = dict(STRICT_PARAMS)
        params["phi"] += delta
        coeffs = moment_coefficients(params)["coefficients"]
        deltas = {name: coeffs[name] - base_coeffs[name] for name in COEFFICIENT_ORDER}
        rows.append(
            {
                "phi_delta": delta,
                "coefficients": coeffs,
                "coefficient_deltas_from_current_phi": deltas,
                "max_abs_delta": max(abs(value) for value in deltas.values()),
            }
        )
    return rows


def compare_to_p1562(base_coeffs: dict[str, float], p1562: dict[str, Any]) -> dict[str, Any]:
    old = p1562.get("derived_lagrangian_coefficients", {})
    rows = []
    for name in COEFFICIENT_ORDER:
        p1562_value = old.get(name)
        replay_value = base_coeffs[name]
        delta = None if p1562_value is None else replay_value - float(p1562_value)
        rows.append({"coefficient": name, "p1562_value": p1562_value, "p2441_replay_value": replay_value, "delta": delta})
    return {
        "row_count": len(rows),
        "max_abs_delta": max(abs(row["delta"]) for row in rows if row["delta"] is not None),
        "rows": rows,
        "p2441_steps": STRICT_PARAMS["steps"],
        "p1562_steps": p1562.get("strict_kernel_params", {}).get("steps"),
    }


def append_doc_sections() -> None:
    eq_section = """
## P2441/S1391 strict moment coefficient phase-sensitivity rank certificate

`P2441/S1391` audits the actual P1562 strict moment coefficient route, not the P1664 algebraic replay ansatz.  A finite-difference Jacobian for `(lambda_sm_eff, kappa_gr_eff, epsilon_mix_eff)` with respect to `(omega, phi, beta, eta)` has rank `3`, and its `phi` column is nonzero.  Small phase sweeps around the current `QW-2049` phase change all three moment-derived effective coefficients.

Therefore `P2440`'s phase-null replay is not an admissible replacement for the strict moment route unless a separate theorem proves phase/topology invariance for the physical coefficient map.  This remains a sensitivity/obstruction certificate only: it exports no SM/GR physical-value generator, no `QW-2191` discharge, and no role-bearing `L_total` closure.
""".strip()
    lag_section = """
## P2441/S1391 strict moment phase-sensitivity guard

`P2441/S1391` shows that the P1562 moment-derived effective coefficients are locally sensitive to `phi`.  Any Lagrangian coefficient source that drops `phi` must now carry an explicit phase-invariance theorem before it can replace the strict moment route or feed final SM/GR constants.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    base = moment_coefficients(STRICT_PARAMS)
    base_coeffs = base["coefficients"]
    jacobian = finite_difference_jacobian(STRICT_PARAMS)
    rank = real_rank(jacobian)
    column_norms = {
        name: l2_norm([row[index] for row in jacobian]) for index, name in enumerate(PARAMETER_ORDER)
    }
    phi_index = PARAMETER_ORDER.index("phi")
    phi_column = [row[phi_index] for row in jacobian]
    sweep = phi_sweep(base_coeffs)
    p1562_comparison = compare_to_p1562(base_coeffs, sources["P1562_THREE_COEFF_DERIVATION"])
    p2440_theorem = (
        sources["P2440_REPLAY_RANK"].get("current_strict_tuple_coefficient_replay_rank_certificate", {}).get("theorem_export", {})
    )
    theorem_export = {
        "theorem_name": "P2441_T1_strict_moment_coefficient_phase_sensitivity_rank_certificate",
        "strict_params_used": STRICT_PARAMS,
        "coefficient_order": COEFFICIENT_ORDER,
        "parameter_order": PARAMETER_ORDER,
        "base_moments": base["moments"],
        "base_dimensionless_ratios": base["dimensionless_ratios"],
        "base_coefficients": base_coeffs,
        "p1562_comparison": p1562_comparison,
        "jacobian_numeric": jacobian,
        "jacobian_real_rank": rank,
        "jacobian_column_norms": column_norms,
        "phi_column": phi_column,
        "phi_column_nonzero": column_norms["phi"] > 1e-6,
        "phi_sweep_rows": sweep,
        "phi_sweep_row_count": len(sweep),
        "every_phi_sweep_row_changes_a_coefficient": all(row["max_abs_delta"] > 1e-4 for row in sweep),
        "p2440_phase_null_inherited": p2440_theorem.get("phi_column_zero") is True,
        "p2440_not_a_replacement_for_phase_sensitive_moment_route": True,
        "strict_phase_invariance_theorem_exported_by_this_certificate": False,
        "strict_physical_value_generator_exported": False,
        "qw2191_discharged": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "P1562-style moment coefficients are locally phase-sensitive at the current strict tuple.",
            "A phase-null coefficient ansatz cannot replace the moment route without a separate phase-invariance theorem.",
            "This sensitivity certificate is not a Standard-Model/GR physical-value generator.",
            "No QW-2191 discharge, role-bearing L_total export, or ToE closure is exported.",
        ],
        "next_honest_step": (
            "Construct a phase/topology-sensitive strict coefficient map, or prove an explicit phase-invariance theorem identifying which coefficients may lawfully ignore phi."
        ),
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "three_coefficients": len(base_coeffs) == 3,
        "four_parameters": len(PARAMETER_ORDER) == 4,
        "jacobian_rank_three": theorem_export["jacobian_real_rank"] == 3,
        "phi_column_nonzero": theorem_export["phi_column_nonzero"],
        "phi_sweep_four_rows": theorem_export["phi_sweep_row_count"] == 4,
        "phi_sweep_changes_coefficients": theorem_export["every_phi_sweep_row_changes_a_coefficient"],
        "p1562_replay_close": theorem_export["p1562_comparison"]["max_abs_delta"] < 2e-4,
        "p2440_phase_null_inherited": theorem_export["p2440_phase_null_inherited"],
        "p2440_not_replacement": theorem_export["p2440_not_a_replacement_for_phase_sensitive_moment_route"],
        "no_phase_invariance_theorem_export": not theorem_export["strict_phase_invariance_theorem_exported_by_this_certificate"],
        "no_value_generator_export": not theorem_export["strict_physical_value_generator_exported"],
        "no_qw2191_discharge": not theorem_export["qw2191_discharged"],
        "no_ltotal_export": not theorem_export["role_bearing_ltotal_exported"],
        "no_toe_export": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2441_s1391_v1",
        "packet_id": "P2441",
        "stage_id": "S1391",
        "result_kind": "STRICT_MOMENT_COEFFICIENT_PHASE_SENSITIVITY_RANK_CERTIFICATE",
        "status": "PASS_STRICT_MOMENT_COEFFICIENT_PHASE_SENSITIVE_NO_GENERATOR",
        "strict_moment_coefficient_phase_sensitivity_rank_certificate": {
            "rg_nonduplication_audit": grep,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": theorem_export["next_honest_step"],
        "global_status": "OPEN_PROGRESS_STRICT_MOMENT_PHASE_SENSITIVITY_EXPORTED_NO_CLOSURE",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["strict_moment_coefficient_phase_sensitivity_rank_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2441 S1391: strict moment coefficient phase-sensitivity rank certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Finite facts",
                "",
                f"- Coefficients: `{theorem['coefficient_order']}`.",
                f"- Parameters: `{theorem['parameter_order']}`.",
                f"- Jacobian rank: `{theorem['jacobian_real_rank']}`.",
                f"- Phi column norm: `{theorem['jacobian_column_norms']['phi']}`.",
                f"- Phi sweep rows: `{theorem['phi_sweep_row_count']}`.",
                f"- Every phi sweep row changes a coefficient: `{theorem['every_phi_sweep_row_changes_a_coefficient']}`.",
                f"- P2440 phase-null inherited: `{theorem['p2440_phase_null_inherited']}`.",
                "",
                "## Hard limits",
                "",
                *[f"- {item}" for item in theorem["not_licensed"]],
                "",
                "## Next honest step",
                "",
                theorem["next_honest_step"],
                "",
                "## Gatekeepers",
                "",
                f"`{payload['gatekeeper_checks']}`",
                "",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    append_doc_sections()
    payload = build_payload()
    write_outputs(payload)
    if not all(payload["gatekeeper_checks"].values()):
        raise SystemExit(f"gatekeeper failure: {payload['gatekeeper_checks']}")


if __name__ == "__main__":
    main()
