#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_R11 = GENERATED / "r11_symmetry_certified_declared_control_transport_packet_for_psi_block_route.json"
IN_R12 = GENERATED / "r12_explicit_canonical_psi_block_export_packet_with_kernel_mixing_for_host_route.json"

OUT_OBJECT = GENERATED / "m_control_canonical_psi_declared_coefficient_filled_v1.json"
OUT_JSON = (
    GENERATED
    / "p476_current_strict_r11_r12_declared_control_pullback_m_control_coefficient_filled_export_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p476_current_strict_r11_r12_declared_control_pullback_m_control_coefficient_filled_export_probe_summary.json"
)


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _p(path: Path) -> str:
    try:
        return str(path.relative_to(REPO))
    except ValueError:
        return str(path)


def _clean_scalar(value: float) -> float:
    if abs(value) < 1e-15:
        return 0.0
    return round(float(value), 15)


def _build_linear_expr(terms: list[tuple[float, str]]) -> tuple[str, list[dict[str, Any]]]:
    pieces: list[str] = []
    term_list: list[dict[str, Any]] = []
    for coeff, term in terms:
        c = _clean_scalar(coeff)
        if c == 0.0:
            continue
        term_list.append({"coeff": c, "term": term})
        if c == 1.0:
            pieces.append(f"({term})")
        elif c == -1.0:
            pieces.append(f"-({term})")
        else:
            pieces.append(f"({c})*({term})")
    if not pieces:
        return "0.0", []
    expr = pieces[0]
    for piece in pieces[1:]:
        if piece.startswith("-"):
            expr += " - " + piece[1:]
        else:
            expr += " + " + piece
    return expr, term_list


def _require_matrix_rows(name: str, value: Any, n_rows: int, n_cols: int) -> list[list[Any]]:
    if not (isinstance(value, list) and len(value) == n_rows):
        raise ValueError(f"expected {name} to be a {n_rows}x{n_cols} list")
    for i, row in enumerate(value):
        if not (isinstance(row, list) and len(row) == n_cols):
            raise ValueError(f"expected {name}[{i}] to have length {n_cols}")
    return value


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    required = {"R11": IN_R11, "R12": IN_R12}
    missing = [k for k, p in required.items() if not p.is_file()]
    if missing:
        payload = {
            "stage": "P476",
            "date_utc": datetime.now(timezone.utc).date().isoformat(),
            "goal": "export_coefficient_filled_declared_control_pullback_m_control_from_canonical_psi_block__r11_r12",
            "status": "FAIL_MISSING_REQUIRED_INPUTS",
            "missing_required_inputs": missing,
            "required_paths": {k: _p(v) for k, v in required.items()},
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(
            json.dumps(
                {
                    "stage": "P476",
                    "status": payload["status"],
                    "missing_required_inputs": payload["missing_required_inputs"],
                    "no_false_pass": True,
                },
                indent=2,
                ensure_ascii=True,
            )
            + "\n",
            encoding="ascii",
        )
        print(OUT_SUMMARY)
        return

    r11 = _read_json(IN_R11)
    r12 = _read_json(IN_R12)

    T = _require_matrix_rows("R11.transport_packet.matrix_rows", r11["transport_packet"]["matrix_rows"], 12, 4)
    H = _require_matrix_rows("R12.canonical_psi_block_export.matrix_rows", r12["canonical_psi_block_export"]["matrix_rows"], 12, 12)

    control_basis = r11["transport_packet"]["domain_basis"]
    psi_basis = r11["transport_packet"]["target_basis"]
    if control_basis != ["c1", "s1", "c2", "s2"]:
        raise ValueError("expected R11 transport domain_basis = ['c1','s1','c2','s2']")
    if psi_basis != [f"psi{i}" for i in range(12)]:
        raise ValueError("expected R11 transport target_basis = ['psi0'..'psi11']")

    diag_terms = [str(H[i][i]) for i in range(12)]
    ksym_terms: dict[tuple[int, int], str] = {}
    for i in range(12):
        for j in range(i + 1, 12):
            ksym_terms[(i, j)] = str(H[i][j])

    matrix_rows: list[list[str]] = [["0.0"] * 4 for _ in range(4)]
    entry_terms: list[list[list[dict[str, Any]]]] = [[[] for _ in range(4)] for _ in range(4)]

    # Compute only the upper triangle and mirror it to enforce symmetry and avoid float-drift asymmetry.
    for a in range(4):
        for b in range(a, 4):
            terms: list[tuple[float, str]] = []
            for i in range(12):
                terms.append((float(T[i][a]) * float(T[i][b]), diag_terms[i]))
            for i in range(12):
                for j in range(i + 1, 12):
                    coeff = float(T[i][a]) * float(T[j][b]) + float(T[j][a]) * float(T[i][b])
                    terms.append((coeff, ksym_terms[(i, j)]))
            expr, kept = _build_linear_expr(terms)
            matrix_rows[a][b] = expr
            matrix_rows[b][a] = expr
            entry_terms[a][b] = kept
            entry_terms[b][a] = kept

    pair1_block = [row[:2] for row in matrix_rows[:2]]
    pair2_block = [row[2:] for row in matrix_rows[2:]]

    obj = {
        "object": "M_control_canonical_psi_declared_coefficient_filled_v1",
        "status": "EXPORTED_DECLARED_COEFFICIENT_FILLED_CONTROL_PULLBACK_FROM_CANONICAL_PSI_BLOCK",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "intent": (
            "Export a coefficient-filled declared control pullback M_control := T_control^T H_PsiPsi T_control using "
            "the explicit declared control transport matrix T_control from R11 and the explicit coefficient-filled canonical "
            "Psi x Psi block H_PsiPsi from R12 (full declared transport support). This is a declared-control export only: "
            "no physical canonicalization, host matching, factorization, selector closure, or QW-2191 discharge is claimed."
        ),
        "inputs": {
            "r11_transport_packet": _p(IN_R11),
            "r12_canonical_psi_block": _p(IN_R12),
        },
        "domain_basis": control_basis,
        "target_basis": psi_basis,
        "assembly_formula": "M_control = T_control^T H_PsiPsi_declared_full_support T_control",
        "computed": {
            "m_control_matrix_rows": matrix_rows,
            "m_control_entry_terms": entry_terms,
            "pair1_block_rows": pair1_block,
            "pair2_block_rows": pair2_block,
        },
        "hard_limits": [
            "Declared-control export only (uses R11 transport without selector-relevant physical canonicalization).",
            "Does not identify the QW-2186 existing-feedback host operator with the canonical Psi block (C10_B1 remains open).",
            "Does not claim an intertwiner/equality witness to the computed current-pair H3 block (P15/P16 remain open).",
            "Does not claim strict-core selector closure / admissible S_sel_int.",
            "Does not claim global discharge of QW-2191.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    payload = {
        "stage": "P476",
        "date": datetime.now(timezone.utc).date().isoformat(),
        "goal": "export_coefficient_filled_declared_control_pullback_m_control_from_canonical_psi_block__r11_r12",
        "status": "PASS_EXPORTED_COEFFICIENT_FILLED_DECLARED_CONTROL_PULLBACK_M_CONTROL",
        "exported_object": _p(OUT_OBJECT),
        "notes": [
            "This is a declared-control export: it uses the explicit declared R11 transport and the explicit R12 canonical Psi block.",
            "It is not a host-matching or factorization witness for existing kernel feedback.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": payload["stage"],
        "status": payload["status"],
        "exported_object": payload["exported_object"],
        "no_false_pass": True,
    }

    OUT_OBJECT.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_JSON.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

