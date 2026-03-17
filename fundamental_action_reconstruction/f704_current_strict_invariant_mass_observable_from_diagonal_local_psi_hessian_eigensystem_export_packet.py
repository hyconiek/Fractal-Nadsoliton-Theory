#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

AS_OF = "2026-03-17"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F459_OBJECT = GENERATED / "psi_hessian_diagonal_local_strict_derived_value_instantiated_v1.json"

OUT_OBJECT = GENERATED / "mass_observable_diagonal_local_strict_derived_v1.json"
OUT_SUMMARY = (
    GENERATED
    / "f704_current_strict_invariant_mass_observable_from_diagonal_local_psi_hessian_eigensystem_export_packet_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def is_finite_number(x: Any) -> bool:
    return isinstance(x, (int, float)) and math.isfinite(float(x))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    if not IN_F459_OBJECT.exists():
        artifact = {
            "stage": "F704",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": [str(IN_F459_OBJECT.relative_to(REPO))],
            "no_false_pass": True,
        }
        OUT_OBJECT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_OBJECT)
        return

    f459 = load_json(IN_F459_OBJECT)
    out = (f459.get("outputs") or {}) if isinstance(f459, dict) else {}
    evals = out.get("eigenvalues_ascending")

    if not (
        isinstance(evals, list)
        and len(evals) == 12
        and all(is_finite_number(v) and float(v) > 0.0 for v in evals)
    ):
        artifact = {
            "stage": "F704",
            "status": "INVALID_F459_EIGENVALUES_SHAPE",
            "as_of": AS_OF,
            "error": "Expected F459.outputs.eigenvalues_ascending as length-12 list[positive finite number].",
            "input_ref": str(IN_F459_OBJECT.relative_to(REPO)),
            "no_false_pass": True,
        }
        OUT_OBJECT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_OBJECT)
        return

    evals_m2 = [float(v) for v in evals]
    masses_m = [float(math.sqrt(v)) for v in evals_m2]

    spectrum = [
        {"mode_index": i, "m2_proxy": float(evals_m2[i]), "m_proxy": float(masses_m[i])}
        for i in range(12)
    ]

    obj = {
        "object": "MassObservable_diagonal_local_strict_derived_v1",
        "stage": "F704",
        "status": "PASS_EXPORTED_STRICT_INVARIANT_MASS_OBSERVABLE_OBJECT",
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "goal": (
            "Export a named strict object defining the observer-limit mass-like readout as the basis-invariant eigenvalue spectrum "
            "of the diagonal/local Psi-sector Hessian H_psi at the exported value instantiation (dimensionless proxy only)."
        ),
        "inputs": {
            "f459_h_psi_eigensystem_object": str(IN_F459_OBJECT.relative_to(REPO)),
        },
        "definition": {
            "observable": "eigenvalues of H_psi (m^2 proxy), and their square roots (m proxy)",
            "basis_invariant": True,
            "normalization": "strict_repo_internal_dimensionless",
            "notes": [
                "This object is a proxy readout only; it does not assign physical units.",
                "See N703 for strict scope/meaning discipline: these are quadratic coefficients in the strict normalization.",
            ],
        },
        "outputs": {
            "n": 12,
            "eigenvalues_m2_proxy_ascending": evals_m2,
            "m_proxy_ascending": masses_m,
            "spectrum": spectrum,
            "min_m2_proxy": float(min(evals_m2)),
            "max_m2_proxy": float(max(evals_m2)),
        },
        "hard_limits": [
            "Does not assign any physical mass units or calibration (no GeV/SI claim).",
            "Does not identify Standard Model particles or claim host matching.",
            "Does not imply kernel-alone/global QW-2191 discharge.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F704",
        "status": obj["status"],
        "as_of": AS_OF,
        "min_m2_proxy": obj["outputs"]["min_m2_proxy"],
        "max_m2_proxy": obj["outputs"]["max_m2_proxy"],
        "exported_object_ref": str(OUT_OBJECT.relative_to(REPO)),
        "no_false_pass": True,
    }

    OUT_OBJECT.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_OBJECT)


if __name__ == "__main__":
    main()

