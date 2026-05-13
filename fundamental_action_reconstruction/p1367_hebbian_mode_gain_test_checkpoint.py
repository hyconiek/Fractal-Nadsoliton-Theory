#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

FAR = Path(__file__).resolve().parent
GEN = FAR / "generated"


def k_strict(d: float, omega: float, phi: float, beta: float, eta: float) -> float:
    return math.cos(omega * d + phi) / (1.0 + beta * (d ** eta))


def moment(order: int, n: int, omega: float, phi: float, beta: float, eta: float) -> float:
    acc = 0.0
    for i in range(1, n + 1):
        d = i / n
        acc += (d ** order) * k_strict(d, omega, phi, beta, eta)
    return acc / n


def map_values(m0: float, m1: float, m2: float, m3: float) -> dict[str, float]:
    return {
        "g1": abs(m0),
        "g2": abs(m1) * 3.0,
        "g3": abs(m2) * 12.0,
        "gravity_effective_observable_1": abs(m3) * 40.0,
    }


def z_report(pred: dict[str, float], ref: dict[str, tuple[float, float]]) -> dict:
    rows = []
    mx = 0.0
    for k, p in pred.items():
        obs, sig = ref[k]
        z = abs(p - obs) / sig
        mx = max(mx, z)
        rows.append({"observable": k, "predicted": p, "observed": obs, "sigma": sig, "abs_z": z})
    return {"rows": rows, "max_abs_z": mx}


def load_ref(path: Path) -> dict[str, tuple[float, float]]:
    out = {}
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f)
        if r.fieldnames != ["observable", "observed", "sigma"]:
            raise SystemExit("reference CSV schema mismatch")
        for row in r:
            out[row["observable"]] = (float(row["observed"]), float(row["sigma"]))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="P1367 Hebbian-like mode gain A/B test")
    ap.add_argument("--n-grid", type=int, default=4096)
    ap.add_argument("--omega", type=float, default=0.18575)
    ap.add_argument("--phi", type=float, default=0.16250)
    ap.add_argument("--beta", type=float, default=1.0)
    ap.add_argument("--eta", type=float, default=1.8)
    ap.add_argument("--gain", type=float, default=1.15)
    ap.add_argument("--suppress", type=float, default=0.85)
    ap.add_argument("--ref", type=Path, default=GEN / "p1358_reference_observed_values.csv")
    ap.add_argument("--out", type=Path, default=GEN / "p1367_hebbian_mode_gain_test_summary.json")
    a = ap.parse_args()

    m0, m1, m2, m3 = [moment(i, a.n_grid, a.omega, a.phi, a.beta, a.eta) for i in range(4)]

    base = map_values(m0, m1, m2, m3)
    coherent = map_values(m0 * a.gain, m1 * a.gain, m2 * a.suppress, m3 * a.suppress)
    incoherent = map_values(m0 * a.suppress, m1 * a.suppress, m2 * a.gain, m3 * a.gain)

    ref = load_ref(a.ref)
    r_base = z_report(base, ref)
    r_coh = z_report(coherent, ref)
    r_inc = z_report(incoherent, ref)

    out = {
        "packet": "P1367",
        "as_of": "2026-05-12",
        "params": {"gain": a.gain, "suppress": a.suppress},
        "base": r_base,
        "coherent_ab": r_coh,
        "incoherent_ab": r_inc,
        "coherent_beats_incoherent": r_coh["max_abs_z"] < r_inc["max_abs_z"],
        "next_priority": "P1367B_LOOP_P1364_P1365_TO_P1362_P1363",
    }
    a.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1367] wrote {a.out}; coherent_beats_incoherent={out['coherent_beats_incoherent']}")


if __name__ == "__main__":
    main()
