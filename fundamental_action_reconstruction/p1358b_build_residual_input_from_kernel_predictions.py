#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path

FAR = Path(__file__).resolve().parent
GEN = FAR / "generated"


def _read_pred(path: Path) -> dict[str, float]:
    out = {}
    with path.open("r", encoding="utf-8", newline="") as f:
        rdr = csv.DictReader(f)
        if rdr.fieldnames != ["observable", "predicted"]:
            raise SystemExit("predicted CSV must have columns observable,predicted")
        for r in rdr:
            out[r["observable"]] = float(r["predicted"])
    return out


def _read_ref(path: Path) -> dict[str, tuple[float, float]]:
    out = {}
    with path.open("r", encoding="utf-8", newline="") as f:
        rdr = csv.DictReader(f)
        if rdr.fieldnames != ["observable", "observed", "sigma"]:
            raise SystemExit("reference CSV must have columns observable,observed,sigma")
        for r in rdr:
            out[r["observable"]] = (float(r["observed"]), float(r["sigma"]))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="P1358b: join kernel predictions with reference data for residual evaluation")
    ap.add_argument("--pred", type=Path, default=GEN / "p1358_kernel_predicted_values.csv")
    ap.add_argument("--ref", type=Path, default=GEN / "p1358_reference_observed_values.csv")
    ap.add_argument("--out", type=Path, default=GEN / "p1358_residual_input.csv")
    a = ap.parse_args()

    pred = _read_pred(a.pred)
    ref = _read_ref(a.ref)
    if set(pred) != set(ref):
        raise SystemExit("observable sets differ between prediction and reference")

    with a.out.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["observable", "predicted", "observed", "sigma"])
        for obs in sorted(pred):
            observed, sigma = ref[obs]
            w.writerow([obs, f"{pred[obs]:.12f}", f"{observed:.12f}", f"{sigma:.12f}"])

    print(f"[P1358b] wrote {a.out}")


if __name__ == "__main__":
    main()
