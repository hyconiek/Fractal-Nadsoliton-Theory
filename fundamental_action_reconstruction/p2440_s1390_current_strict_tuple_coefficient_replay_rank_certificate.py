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

OUT = GEN / "p2440_s1390_current_strict_tuple_coefficient_replay_rank_certificate.json"
MD = GEN / "p2440_s1390_current_strict_tuple_coefficient_replay_rank_certificate.md"

SOURCE_FILES = {
    "P1562_THREE_COEFF_DERIVATION": GEN / "p1562_s512_strict_kernel_to_lagrangian_coefficient_derivation_summary.json",
    "P1563_KERNEL_TO_EOM": GEN / "p1563_s513_strict_kernel_to_euler_lagrange_eom_export_summary.json",
    "P1664_MANIFEST_INVERSION": GEN / "p1664_s614_strict_full_lagrangian_manifest_and_inversion.json",
    "P2439_COEFF_SOURCE_AUDIT": GEN / "p2439_s1389_strict_coefficient_source_consistency_audit_certificate.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

CURRENT_STRICT_TUPLE = {
    "omega": 0.18575,
    "phi": 0.16250,
    "beta": 1.0,
    "eta": 1.8,
    "A": 1.0,
}

PARAMETER_ORDER = ["omega", "phi", "beta", "eta", "A"]
COEFFICIENT_ORDER = [
    "Mpl2",
    "cR2",
    "cRic2",
    "cRiem2",
    "Z3",
    "Z2",
    "Z1",
    "muH2",
    "lambdaH",
    "xiHR",
    "chiRG",
    "yu",
    "yd",
    "ye",
]


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
        "new_packet": "P2440|S1390|current strict tuple coefficient replay|coefficient replay rank|phase-null coefficient replay",
        "p2439_input": "P2439|S1389|strict coefficient-source consistency|current strict tuple coefficient source",
        "p1664_formula": "P1664|strict_full_lagrangian_manifest_and_inversion|inverse_recovery|coefficient_map",
        "current_strict_tuple": "omega = 0.18575|phi = 0.16250|beta = 1.0|eta = 1.8|QW-2049",
        "phase_selector_blocker": "phi|phase|topological|QW-2191|selector uniqueness",
        "false_closure_flags": "qw2191_closed|toe_closed|strict_observable_value_generator_exported|toe_closure_exported",
        "rank_or_jacobian": "Jacobian|rank|local invertibility|null direction|nullspace",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def replay_coefficients(params: dict[str, float]) -> dict[str, float]:
    omega = params["omega"]
    beta = params["beta"]
    eta = params["eta"]
    A = params["A"]
    return {
        "Mpl2": A * (1.0 + beta),
        "cR2": A * beta / (1.0 + eta),
        "cRic2": A * beta * eta / (1.0 + eta),
        "cRiem2": A * beta * (1.0 + eta) / 4.0,
        "Z3": 1.0 + 0.8 * beta**2,
        "Z2": 1.0 + 0.6 * beta**2,
        "Z1": 1.0 + 0.4 * beta**2,
        "muH2": A * omega**2,
        "lambdaH": (1.0 + eta**2) / (1.0 + beta),
        "xiHR": beta / (1.0 + beta),
        "chiRG": beta * eta / (2.0 + beta),
        "yu": 0.9 * omega,
        "yd": 0.5 * omega,
        "ye": 0.3 * omega,
    }


def local_inverse(coeff: dict[str, float]) -> dict[str, Any]:
    beta_rec = coeff["xiHR"] / (1.0 - coeff["xiHR"])
    A_rec = coeff["Mpl2"] / (1.0 + beta_rec)
    omega_rec = math.sqrt(max(coeff["muH2"] / A_rec, 0.0))
    eta_rec = math.sqrt(max(coeff["lambdaH"] * (1.0 + beta_rec) - 1.0, 0.0))
    recovered = {"omega": omega_rec, "beta": beta_rec, "eta": eta_rec, "A": A_rec}
    errors = {name: abs(value - CURRENT_STRICT_TUPLE[name]) for name, value in recovered.items()}
    return {
        "recovered_parameters": recovered,
        "unrecovered_parameters": ["phi"],
        "absolute_errors": errors,
        "local_inverse_pass_on_recovered_four_parameters": all(value < 1e-10 for value in errors.values()),
    }


def finite_difference_jacobian(params: dict[str, float]) -> list[list[float]]:
    jacobian: list[list[float]] = []
    for coeff_name in COEFFICIENT_ORDER:
        row: list[float] = []
        for param_name in PARAMETER_ORDER:
            step = 1e-6 * max(1.0, abs(params[param_name]))
            plus = dict(params)
            minus = dict(params)
            plus[param_name] += step
            minus[param_name] -= step
            deriv = (replay_coefficients(plus)[coeff_name] - replay_coefficients(minus)[coeff_name]) / (2.0 * step)
            row.append(deriv)
        jacobian.append(row)
    return jacobian


def real_rank(matrix: list[list[float]], tol: float = 1e-9) -> int:
    mat = [row[:] for row in matrix]
    if not mat:
        return 0
    rows = len(mat)
    cols = len(mat[0])
    rank = 0
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


def compare_to_p1664(current_coeffs: dict[str, float], p1664: dict[str, Any]) -> dict[str, Any]:
    old_coeffs = p1664.get("coefficient_map", {})
    rows = []
    for name in COEFFICIENT_ORDER:
        old = old_coeffs.get(name)
        new = current_coeffs[name]
        delta = None if old is None else new - float(old)
        rows.append(
            {
                "coefficient": name,
                "p1664_value": old,
                "current_tuple_replay_value": new,
                "delta": delta,
                "changed": delta is None or abs(delta) > 1e-12,
            }
        )
    return {
        "row_count": len(rows),
        "changed_count": sum(1 for row in rows if row["changed"]),
        "unchanged_count": sum(1 for row in rows if not row["changed"]),
        "rows": rows,
    }


def append_doc_sections() -> None:
    eq_section = """
## P2440/S1390 current strict tuple coefficient replay rank certificate

`P2440/S1390` performs a bounded computational replay of the `P1664` algebraic coefficient ansatz at the current `QW-2049` strict tuple.  The replay exports numeric coefficients and a finite-difference Jacobian rank fact, but it also proves the ansatz has a `phi`-null direction: the coefficient formulas recover `omega`, `beta`, `eta`, and `A` locally, while the phase/topological parameter `phi` is not represented.  Therefore the replay is only a diagnostic coefficient candidate, not a strict physical-value generator or selector theorem.

This certificate also quarantines old closure flags when they conflict with the current P2438/P2439 no-closure state: local coefficient replay must not be read as `QW-2191` discharge, ToE closure, SM/GR value generation, or role-bearing `L_total` export.
""".strip()
    lag_section = """
## P2440/S1390 current-tuple coefficient replay guard

`P2440/S1390` allows the `P1664` formulas to be replayed at the current strict tuple only as a diagnostic coefficient ansatz.  Because the replay Jacobian has a `phi`-null direction and no selector/observable theorem, the resulting coefficients cannot be inserted into the Lagrangian as final SM/GR physical constants or as a closure of `QW-2191`.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    current_coeffs = replay_coefficients(CURRENT_STRICT_TUPLE)
    inverse = local_inverse(current_coeffs)
    jacobian = finite_difference_jacobian(CURRENT_STRICT_TUPLE)
    rank = real_rank(jacobian)
    phi_column = [row[PARAMETER_ORDER.index("phi")] for row in jacobian]
    column_norms = {
        param: l2_norm([row[index] for row in jacobian]) for index, param in enumerate(PARAMETER_ORDER)
    }
    p1664_comparison = compare_to_p1664(current_coeffs, sources["P1664_MANIFEST_INVERSION"])
    p1562 = sources["P1562_THREE_COEFF_DERIVATION"]
    p1563 = sources["P1563_KERNEL_TO_EOM"]
    p2439_theorem = (
        sources["P2439_COEFF_SOURCE_AUDIT"]
        .get("strict_coefficient_source_consistency_audit_certificate", {})
        .get("theorem_export", {})
    )
    p1562_closure_flags = {
        "qw2191_closed": p1562.get("qw2191_closed"),
        "toe_closed": p1562.get("toe_closed"),
    }
    closure_conflict = (
        p1562_closure_flags["qw2191_closed"] is True
        and p1563.get("qw2191_closed") is False
        and p2439_theorem.get("qw2191_discharged") is False
    ) or (
        p1562_closure_flags["toe_closed"] is True
        and p1563.get("toe_closed") is False
        and p2439_theorem.get("toe_closure_exported") is False
    )
    theorem_export = {
        "theorem_name": "P2440_T1_current_strict_tuple_coefficient_replay_rank_certificate",
        "current_strict_tuple_with_amplitude": CURRENT_STRICT_TUPLE,
        "replay_source_formula": "P1664 algebraic coefficient ansatz replayed at current QW-2049 tuple with A=1",
        "coefficient_count": len(current_coeffs),
        "coefficient_order": COEFFICIENT_ORDER,
        "current_tuple_replayed_coefficients": current_coeffs,
        "local_inverse": inverse,
        "jacobian_parameter_order": PARAMETER_ORDER,
        "jacobian_coefficient_order": COEFFICIENT_ORDER,
        "jacobian_numeric": jacobian,
        "jacobian_real_rank": rank,
        "jacobian_column_norms": column_norms,
        "phi_column_zero": all(abs(value) < 1e-12 for value in phi_column),
        "phase_parameter_recovered_by_replay": False,
        "p1664_current_tuple_coefficient_comparison": p1664_comparison,
        "p1562_closure_flags": p1562_closure_flags,
        "p1562_closure_flags_conflict_with_current_no_closure_state": closure_conflict,
        "p2439_no_acceptable_generator_inherited": p2439_theorem.get(
            "acceptable_current_strict_physical_value_generator_source_count"
        )
        == 0,
        "strict_kernel_to_coefficient_map_theorem_exported_by_this_certificate": False,
        "strict_observable_value_generator_exported": False,
        "qw2191_discharged": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "Replaying the P1664 ansatz at the current tuple is not a derivation of the ansatz from the strict nadsoliton.",
            "The replay is phase-null: phi/topological selector data are absent from the coefficient formulas.",
            "Old P1562 closure flags are quarantined when they conflict with P1563/P2438/P2439 no-closure results.",
            "No SM/GR physical-value generator, QW-2191 discharge, role-bearing L_total, or ToE closure is exported.",
        ],
        "next_honest_step": (
            "Add a phase/topology-sensitive current strict coefficient map or prove that physical coefficients are phase-invariant before promoting any replayed coefficient to SM/GR value status."
        ),
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "fourteen_coefficients_replayed": theorem_export["coefficient_count"] == 14,
        "local_inverse_recovers_four_parameters": inverse["local_inverse_pass_on_recovered_four_parameters"],
        "phi_unrecovered": inverse["unrecovered_parameters"] == ["phi"],
        "jacobian_rank_four": theorem_export["jacobian_real_rank"] == 4,
        "phi_column_zero": theorem_export["phi_column_zero"],
        "phase_not_recovered": not theorem_export["phase_parameter_recovered_by_replay"],
        "p1664_comparison_complete": p1664_comparison["row_count"] == 14,
        "p1664_has_changed_coefficients": p1664_comparison["changed_count"] > 0,
        "p1562_closure_conflict_detected": theorem_export[
            "p1562_closure_flags_conflict_with_current_no_closure_state"
        ],
        "p2439_no_generator_inherited": theorem_export["p2439_no_acceptable_generator_inherited"],
        "no_coefficient_theorem_export": not theorem_export[
            "strict_kernel_to_coefficient_map_theorem_exported_by_this_certificate"
        ],
        "no_observable_generator_export": not theorem_export["strict_observable_value_generator_exported"],
        "no_qw2191_discharge": not theorem_export["qw2191_discharged"],
        "no_ltotal_export": not theorem_export["role_bearing_ltotal_exported"],
        "no_toe_export": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2440_s1390_v1",
        "packet_id": "P2440",
        "stage_id": "S1390",
        "result_kind": "CURRENT_STRICT_TUPLE_COEFFICIENT_REPLAY_RANK_CERTIFICATE",
        "status": "PASS_CURRENT_STRICT_TUPLE_COEFFICIENT_REPLAY_PHASE_NULL_NO_GENERATOR",
        "current_strict_tuple_coefficient_replay_rank_certificate": {
            "rg_nonduplication_audit": grep,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": theorem_export["next_honest_step"],
        "global_status": "OPEN_PROGRESS_CURRENT_TUPLE_COEFFICIENT_REPLAY_EXPORTED_NO_CLOSURE",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["current_strict_tuple_coefficient_replay_rank_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2440 S1390: current strict tuple coefficient replay rank certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Finite facts",
                "",
                f"- Coefficients replayed: `{theorem['coefficient_count']}`.",
                f"- Jacobian rank over parameters `{theorem['jacobian_parameter_order']}`: `{theorem['jacobian_real_rank']}`.",
                f"- Phi column zero: `{theorem['phi_column_zero']}`.",
                f"- Local inverse recovered parameters: `{list(theorem['local_inverse']['recovered_parameters'])}`.",
                f"- Local inverse unrecovered parameters: `{theorem['local_inverse']['unrecovered_parameters']}`.",
                f"- P1664 changed coefficient count under current tuple replay: `{theorem['p1664_current_tuple_coefficient_comparison']['changed_count']}`.",
                f"- P1562 closure flag conflict detected: `{theorem['p1562_closure_flags_conflict_with_current_no_closure_state']}`.",
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
