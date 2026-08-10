#!/usr/bin/env python3
"""Standard-library exact replay for FIN P518 and P524 exported ledgers.

Hexadecimal binary64 endpoints are parsed exactly and converted to Fraction.
The checker does not regenerate transcendental interval bounds; that explicit
trust boundary is part of both certificate formats.
"""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def q(hexadecimal: str) -> Fraction:
    return Fraction.from_float(float.fromhex(hexadecimal))


def check_p518() -> dict:
    payload = json.loads((ROOT / "FIN_P518_Krawczyk_Replay_Certificate.json").read_text(encoding="utf-8"))
    charts = payload["parameter_charts"]
    bridges = payload["bridges"]
    assert len(charts) == 401
    assert len(bridges) == 400
    for row in charts + bridges:
        assert q(row["inclusion_margin_lower"]) > 0
        assert q(row["defect_upper"]) < 1
    for row in bridges:
        assert q(row["nesting_margin_lower"]) >= 0
    assert q(charts[0]["klo"]) == 0
    assert q(charts[-1]["khi"]) == 1
    for left, right in zip(charts[:-1], charts[1:]):
        assert q(left["khi"]) == q(right["klo"])
    return {"charts": len(charts), "bridges": len(bridges), "inclusions": len(charts) + len(bridges), "accepted": True}


def check_p524() -> dict:
    payload = json.loads((ROOT / "FIN_P524_PSD_Replay_Certificate.json").read_text(encoding="utf-8"))
    rows = payload["rows"]
    assert rows
    assert q(rows[0]["interval"][0]) == 0
    assert q(rows[-1]["interval"][1]) == 1
    unresolved = []
    for index, row in enumerate(rows):
        lo, hi = map(q, row["interval"])
        assert lo < hi
        if index:
            assert q(rows[index - 1]["interval"][1]) == lo
        bounds = [(q(pair[0]), q(pair[1])) for pair in row["eigenvalue_bounds"]]
        assert len(bounds) == 11
        assert all(lower <= upper for lower, upper in bounds)
        if row["class"] == "PSD":
            assert all(lower > 0 for lower, _upper in bounds)
        elif row["class"] == "NONPSD":
            assert any(upper < 0 for _lower, upper in bounds)
        elif row["class"] == "UNRESOLVED":
            unresolved.append([lo, hi])
        else:
            raise AssertionError(f"unknown class {row['class']}")
    components = []
    for lo, hi in unresolved:
        if not components or components[-1][1] != lo:
            components.append([lo, hi])
        else:
            components[-1][1] = hi
    assert len(components) == 1
    return {"boxes": len(rows), "unresolved_boxes": len(unresolved), "unresolved_component": [str(x) for x in components[0]], "accepted": True}


def main() -> None:
    print(json.dumps({"P518": check_p518(), "P524": check_p524()}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
