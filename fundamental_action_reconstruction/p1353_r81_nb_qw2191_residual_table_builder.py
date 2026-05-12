#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

FAR = Path(__file__).resolve().parent
GEN = FAR / "generated"


def _load_contract(path: Path) -> dict:
    if not path.exists():
        raise SystemExit(f"Missing contract file: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def _read_inputs(path: Path) -> list[dict]:
    if not path.exists():
        raise SystemExit(f"Missing input CSV: {path}")
    rows: list[dict] = []
    with path.open("r", encoding="utf-8", newline="") as f:
        rdr = csv.DictReader(f)
        required = {"observable", "predicted", "observed", "sigma"}
        if set(rdr.fieldnames or []) != required:
            raise SystemExit(
                "Input CSV columns must be exactly: observable,predicted,observed,sigma"
            )
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
    ap = argparse.ArgumentParser(
        description="P1353: build residual table and PASS/FAIL decision from preregistered NB R8.1 trial inputs."
    )
    ap.add_argument(
        "--contract",
        type=Path,
        default=GEN / "p1352_qw2191_nb_r81_single_scale_extraction_trial_contract_summary.json",
    )
    ap.add_argument(
        "--inputs",
        type=Path,
        default=GEN / "p1353_r81_nb_qw2191_trial_inputs_template.csv",
    )
    ap.add_argument("--z-threshold", type=float, default=5.0)
    ap.add_argument(
        "--out",
        type=Path,
        default=GEN / "p1353_r81_nb_qw2191_residual_table_summary.json",
    )
    a = ap.parse_args()

    contract = _load_contract(a.contract)
    if contract.get("trial") != "NB-QW2191-R8.1-SINGLE-SCALE-EXTRACTION-TRIAL-V1":
        raise SystemExit("Unsupported contract trial id.")

    rows = _read_inputs(a.inputs)
    residual_rows = []
    all_pass = True
    for row in rows:
        sigma = row["sigma"]
        if sigma <= 0:
            raise SystemExit(f"Sigma must be positive for {row['observable']}")
        residual = row["predicted"] - row["observed"]
        z_score = abs(residual) / sigma
        pass_flag = z_score <= a.z_threshold
        all_pass = all_pass and pass_flag
        residual_rows.append(
            {
                "observable": row["observable"],
                "predicted": row["predicted"],
                "observed": row["observed"],
                "sigma": sigma,
                "residual": residual,
                "abs_z": z_score,
                "pass": pass_flag,
            }
        )

    out = {
        "packet": "P1353",
        "as_of": "2026-05-12",
        "trial": contract["trial"],
        "scale": contract["contract"]["scale"],
        "scheme": contract["contract"]["scheme"],
        "threshold": {"abs_z_max": a.z_threshold},
        "rows": residual_rows,
        "pass_fail_status_v1": "PASS" if all_pass else "FAIL",
        "incident_log_v1": [] if all_pass else ["At least one observable exceeded abs_z threshold."],
        "next_priority": (
            "R8_1_EXTERNAL_BLIND_AUDIT_PREPARE_PUBLICATION"
            if all_pass
            else "R8_1_ROLLBACK_AND_CORRECTIVE_PACKET"
        ),
    }

    a.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1353] wrote {a.out}; status={out['pass_fail_status_v1']}")


if __name__ == "__main__":
    main()
