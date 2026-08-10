#!/usr/bin/env python3
"""Standard-library replay for the P528, P533 and P534 exported ledgers."""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def hex_fraction(text: str) -> Fraction:
    return Fraction.from_float(float.fromhex(text))


def rational(text: str | int | float) -> Fraction:
    return Fraction(str(text))


def check_p528() -> dict:
    data = json.loads((ROOT / "FIN_P528_Stability_Replay_Certificate.json").read_text(encoding="utf-8"))
    charts, bridges = data["charts"], data["bridges"]
    assert len(charts) == 208 and len(bridges) == 207
    for left, right in zip(charts[:-1], charts[1:]):
        assert hex_fraction(left["omega_interval"][1]) == hex_fraction(right["omega_interval"][0])
    for row in charts:
        assert hex_fraction(row["state_inclusion_margin_lower"]) > 0
        assert hex_fraction(row["state_defect_upper"]) < 1
        assert hex_fraction(row["lminus_second_lower"]) > 0
        assert hex_fraction(row["lplus_negative_upper"]) < 0
        assert hex_fraction(row["lplus_first_positive_lower"]) > 0
        assert hex_fraction(row["dP_lower"]) <= hex_fraction(row["dP_upper"])
        assert hex_fraction(row["d2P_lower"]) <= hex_fraction(row["d2P_upper"])
    for row in bridges:
        assert hex_fraction(row["inclusion_margin_lower"]) > 0
        assert hex_fraction(row["defect_upper"]) < 1
        assert hex_fraction(row["nesting_margin_lower"]) >= 0
    return {"charts": len(charts), "bridges": len(bridges), "accepted": True}


def check_p533() -> dict:
    data = json.loads((ROOT / "FIN_P533_Interval_FFT_Certificate.json").read_text(encoding="utf-8"))
    assert "2e-9" not in data["arithmetic_model"]
    for block in (data["coarse"], data["anchor"]):
        assert rational(block["maximum_fft_disk_radius"]) > 0
        assert rational(block["minimum_response_denominator"]) > 0
        for lo, hi in block["sigma_intervals"] + block["fingerprint_intervals"]:
            assert rational(lo) <= rational(hi)
            assert rational(lo) >= 0
        reference = block["fingerprint_intervals"][3]
        assert rational(reference[0]) <= 1 <= rational(reference[1])
    for text in data["combined_C384_to_limit_upper"]:
        assert 0 <= hex_fraction(text) < Fraction(7, 100)
    return {"u_count": len(data["u_values"]), "accepted": True}


def check_p534() -> dict:
    data = json.loads((ROOT / "FIN_P534_Rational_PSD_Certificate.json").read_text(encoding="utf-8"))
    rows = data["rows"]
    assert len(rows) == 573
    previous = Fraction(0)
    for index, row in enumerate(rows):
        lo, hi = (rational(value) for value in row["interval"])
        assert lo == previous and lo <= hi
        previous = hi
        bounds = [(rational(a), rational(b)) for a, b in row["bounds"]]
        assert all(a <= b for a, b in bounds)
        if row["rational_class"] == "PSD":
            assert all(a > 0 for a, _ in bounds)
        elif row["rational_class"] == "NONPSD":
            assert any(b < 0 for _, b in bounds)
        else:
            assert row["rational_class"] == "UNRESOLVED"
        assert row["rational_class"] == row["source_class"]
    assert previous == 1
    return {"boxes": len(rows), "accepted": True}


def main() -> None:
    output = {"P528": check_p528(), "P533": check_p533(), "P534": check_p534()}
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
