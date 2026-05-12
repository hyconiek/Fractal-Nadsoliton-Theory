#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

FAR = Path(__file__).resolve().parent
GEN = FAR / "generated"


def _read_rows(path: Path) -> list[dict[str, float | str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        rdr = csv.DictReader(f)
        cols = ["observable", "predicted", "observed", "sigma"]
        if rdr.fieldnames != cols:
            raise SystemExit(f"CSV must have columns {cols}: {path}")
        rows = []
        for r in rdr:
            rows.append(
                {
                    "observable": r["observable"],
                    "predicted": float(r["predicted"]),
                    "observed": float(r["observed"]),
                    "sigma": float(r["sigma"]),
                }
            )
        return rows


def main() -> None:
    ap = argparse.ArgumentParser(description="P1357: diagnose whether current trial values are theory-predictive or template-calibrated.")
    ap.add_argument("--inputs", type=Path, default=GEN / "p1353_r81_nb_qw2191_trial_inputs_template.csv")
    ap.add_argument("--residual", type=Path, default=GEN / "p1353_r81_nb_qw2191_residual_table_summary.json")
    ap.add_argument("--out", type=Path, default=GEN / "p1357_r81_physical_value_diagnostic_summary.json")
    a = ap.parse_args()

    rows = _read_rows(a.inputs)
    residual = json.loads(a.residual.read_text(encoding="utf-8"))

    exact_matches = 0
    nearly_matches = 0
    n = len(rows)
    tol = 1e-12
    for r in rows:
        d = abs(float(r["predicted"]) - float(r["observed"]))
        if d <= tol:
            exact_matches += 1
        if d <= float(r["sigma"]):
            nearly_matches += 1

    pass_status = residual.get("pass_fail_status_v1")
    diagnostic = {
        "packet": "P1357",
        "as_of": "2026-05-12",
        "input_file": str(a.inputs),
        "residual_file": str(a.residual),
        "row_count": n,
        "exact_prediction_equals_observed_count": exact_matches,
        "within_one_sigma_count": nearly_matches,
        "residual_pass_status": pass_status,
        "physical_value_verdict": (
            "NOT_YET_THEORY_PREDICTIVE_VALUES"
            if exact_matches == n
            else "PARTIAL_SIGNAL_REQUIRES_PROVENANCE_AND_NONTRIVIAL_PREDICTION_CHECK"
        ),
        "reason": (
            "Current template has predicted=observed for all rows, so PASS is pipeline sanity check, not independent physical prediction from kernel."
            if exact_matches == n
            else "Some nontrivial deltas exist; still requires provenance proving values were generated from strict kernel pipeline rather than fitted placeholders."
        ),
        "next_priority": "BUILD_PROVENANCE_LOCKED_KERNEL_TO_VALUE_GENERATOR_P1358",
    }

    a.out.write_text(json.dumps(diagnostic, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1357] wrote {a.out}; verdict={diagnostic['physical_value_verdict']}")


if __name__ == "__main__":
    main()
