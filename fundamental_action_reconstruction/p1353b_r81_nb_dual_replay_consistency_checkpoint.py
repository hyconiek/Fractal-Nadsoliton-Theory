#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

FAR = Path(__file__).resolve().parent
GEN = FAR / "generated"


REQ_COLS = ["observable", "predicted", "observed", "sigma"]


def _read_csv(path: Path) -> dict[str, dict[str, float]]:
    if not path.exists():
        raise SystemExit(f"Missing CSV: {path}")
    out: dict[str, dict[str, float]] = {}
    with path.open("r", encoding="utf-8", newline="") as f:
        rdr = csv.DictReader(f)
        if rdr.fieldnames != REQ_COLS:
            raise SystemExit(f"CSV {path} must have columns: {','.join(REQ_COLS)}")
        for r in rdr:
            obs = r["observable"]
            out[obs] = {
                "predicted": float(r["predicted"]),
                "observed": float(r["observed"]),
                "sigma": float(r["sigma"]),
            }
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="P1353b dual-replay consistency checkpoint for NB-QW2191-R8.1")
    ap.add_argument("--contract", type=Path, default=GEN / "p1352_qw2191_nb_r81_single_scale_extraction_trial_contract_summary.json")
    ap.add_argument("--team-a", type=Path, default=GEN / "p1353b_team_a_inputs.csv")
    ap.add_argument("--team-b", type=Path, default=GEN / "p1353b_team_b_inputs.csv")
    ap.add_argument("--delta-threshold", type=float, default=2.0, help="max allowed cross-team delta in sigma units")
    ap.add_argument("--out", type=Path, default=GEN / "p1353b_r81_nb_dual_replay_consistency_summary.json")
    a = ap.parse_args()

    contract = json.loads(a.contract.read_text(encoding="utf-8"))
    if contract.get("trial") != "NB-QW2191-R8.1-SINGLE-SCALE-EXTRACTION-TRIAL-V1":
        raise SystemExit("Unsupported contract trial id for P1353b")

    A = _read_csv(a.team_a)
    B = _read_csv(a.team_b)

    if set(A) != set(B):
        raise SystemExit("Team A and Team B observable sets differ")

    rows = []
    consistent = True
    for obs in sorted(A):
        arow, brow = A[obs], B[obs]
        sigma_ref = max(arow["sigma"], brow["sigma"])
        if sigma_ref <= 0:
            raise SystemExit(f"Non-positive sigma in {obs}")

        pred_delta = arow["predicted"] - brow["predicted"]
        obs_delta = arow["observed"] - brow["observed"]
        pred_delta_z = abs(pred_delta) / sigma_ref
        obs_delta_z = abs(obs_delta) / sigma_ref
        row_ok = pred_delta_z <= a.delta_threshold and obs_delta_z <= a.delta_threshold
        consistent = consistent and row_ok

        rows.append(
            {
                "observable": obs,
                "team_a": arow,
                "team_b": brow,
                "predicted_delta": pred_delta,
                "observed_delta": obs_delta,
                "predicted_delta_abs_z": pred_delta_z,
                "observed_delta_abs_z": obs_delta_z,
                "within_threshold": row_ok,
            }
        )

    out = {
        "packet": "P1353b",
        "as_of": "2026-05-12",
        "trial": contract["trial"],
        "delta_threshold_abs_z": a.delta_threshold,
        "cross_implementation_consistency": "PASS" if consistent else "FAIL",
        "rows": rows,
        "next_priority": "R8_1_EXTERNAL_BLIND_AUDIT_PREPARE_PUBLICATION" if consistent else "R8_1_INTERNAL_MISMATCH_CORRECTIVE_PACKET",
    }
    a.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1353b] wrote {a.out}; consistency={out['cross_implementation_consistency']}")


if __name__ == "__main__":
    main()
