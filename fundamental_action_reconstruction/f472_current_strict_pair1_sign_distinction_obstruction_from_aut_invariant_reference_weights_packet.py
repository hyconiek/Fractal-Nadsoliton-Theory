#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

AS_OF = "2026-03-16"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_A1 = GENERATED / "a_1_pair1_orientation_projector_operator_strict_core_v1.json"
IN_R_ORD = GENERATED / "r_ord_z12_v1_reference_distribution.json"

OUT_OBJ = GENERATED / "sign_distinction_obstruction_pair1_aut_invariant_reference_weights_strict_v1.json"
OUT_SUMMARY = (
    GENERATED
    / "f472_current_strict_pair1_sign_distinction_obstruction_from_aut_invariant_reference_weights_packet_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def units_mod_n(n: int) -> list[int]:
    out: list[int] = []
    for a in range(n):
        if math.gcd(a, n) == 1:
            out.append(a)
    return out


def aut_orbits_mul_units(n: int) -> list[list[int]]:
    units = units_mod_n(n)
    seen: set[int] = set()
    orbits: list[list[int]] = []
    for x in range(n):
        if x in seen:
            continue
        orb = sorted({(a * x) % n for a in units})
        for y in orb:
            seen.add(y)
        orbits.append(orb)
    return orbits


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    missing = [str(p.relative_to(REPO)) for p in (IN_A1,) if not p.exists()]
    if missing:
        raise SystemExit(
            json.dumps(
                {
                    "stage": "F472",
                    "status": "NOT_COMPUTABLE_MISSING_INPUTS",
                    "missing": missing,
                    "expected": ["F456 exports A_1(pair1) with u_1"],
                    "no_false_pass": True,
                },
                ensure_ascii=True,
            )
        )

    a1 = load_json(IN_A1)
    try:
        u1 = (a1.get("data") or {}).get("u_1")
        if not (isinstance(u1, list) and len(u1) == 12 and all(isinstance(v, (int, float)) for v in u1)):
            raise ValueError("A1.data.u_1 must be a length-12 numeric list")
        u1 = [float(v) for v in u1]
        coords = (a1.get("data") or {}).get("u_1_coords_in_c1_s1")
        coords_ok = isinstance(coords, list) and len(coords) == 2 and all(isinstance(v, (int, float)) for v in coords)
        u1_coords = [float(v) for v in coords] if coords_ok else None
    except Exception as exc:
        raise SystemExit(json.dumps({"stage": "F472", "status": "INVALID_A1_SHAPE", "error": str(exc)}, ensure_ascii=True))

    n = 12

    def reflect(x: int) -> int:
        return (-x) % n

    odd_u1_max_abs = max(abs(u1[x] + u1[reflect(x)]) for x in range(n))
    odd_u1 = bool(odd_u1_max_abs <= 1e-15)

    units = units_mod_n(n)
    orbits = aut_orbits_mul_units(n)

    # Verify orbits partition Z_n
    orbit_union = sorted({x for orb in orbits for x in orb})
    orbits_disjoint = len(orbit_union) == n and len(orbit_union) == len({x for x in orbit_union})

    orbit_sums = []
    all_orbit_sums_zero = True
    for orb in orbits:
        s = sum(u1[x] for x in orb)
        orbit_sums.append({"orbit": orb, "sum_u1_over_orbit": float(s)})
        if abs(float(s)) > 1e-15:
            all_orbit_sums_zero = False

    # Optionally cross-check that the current ord-reference weights are Aut-invariant (N479) and therefore even.
    ord_by_x: list[int] | None = None
    ord_aut_invariant: bool | None = None
    r_ord_aut_invariant: bool | None = None
    dot_u1_with_ord: float | None = None
    dot_u1_with_r_ord: float | None = None
    if IN_R_ORD.exists():
        r_ord = load_json(IN_R_ORD)
        try:
            ord_by_x_raw = r_ord.get("ord_z12_by_x")
            if not (
                isinstance(ord_by_x_raw, list)
                and len(ord_by_x_raw) == 12
                and all(isinstance(v, int) for v in ord_by_x_raw)
            ):
                raise ValueError("r_ord.ord_z12_by_x must be a length-12 int list")
            ord_by_x = [int(v) for v in ord_by_x_raw]
        except Exception as exc:
            raise SystemExit(
                json.dumps({"stage": "F472", "status": "INVALID_R_ORD_SHAPE", "error": str(exc)}, ensure_ascii=True)
            )

        ord_aut_invariant = True
        for a in units:
            for x in range(n):
                if ord_by_x[(a * x) % n] != ord_by_x[x]:
                    ord_aut_invariant = False
                    break
            if not ord_aut_invariant:
                break

        alpha_geo = 4.0 * math.log(2.0)
        r_weights = [math.exp(-alpha_geo * float(o)) for o in ord_by_x]
        r_ord_aut_invariant = True
        for a in units:
            for x in range(n):
                if abs(r_weights[(a * x) % n] - r_weights[x]) > 0.0:
                    r_ord_aut_invariant = False
                    break
            if not r_ord_aut_invariant:
                break

        dot_u1_with_ord = float(sum(float(o) * v for o, v in zip(ord_by_x, u1)))
        dot_u1_with_r_ord = float(sum(w * v for w, v in zip(r_weights, u1)))

    generated_utc = datetime.now(timezone.utc).isoformat()

    obj = {
        "object": "SignDistinctionObstruction_pair1_aut_invariant_reference_weights_strict_v1",
        "stage": "F472",
        "status": "actual_exported_pair1_sign_distinction_obstruction_for_aut_invariant_reference_weights_on_current_exported_instance__no_false_pass",
        "as_of": AS_OF,
        "generated_utc": generated_utc,
        "intent": (
            "Record a strict obstruction that any Aut(Z_12)-invariant (direction-free) reference weight family w(x) cannot distinguish "
            "u_1 from -u_1 on the current exported pair1 sine axis via a sign-sensitive scalar of the form Σ_x w(x) u_1(x). "
            "Reason: Aut(Z_12) contains inversion x↦-x, so Aut-invariant weights are even; while the current exported u_1 is odd under reflection, "
            "so the weighted sum cancels."
        ),
        "inputs": {
            "A1_pair1_ref": str(IN_A1.relative_to(REPO)),
            "r_ord_z12_ref_if_present": (str(IN_R_ORD.relative_to(REPO)) if IN_R_ORD.exists() else None),
        },
        "data": {
            "n": n,
            "aut_units_mod_n": units,
            "aut_orbits_under_units_multiplication": orbits,
        },
        "checks": {
            "u1_odd_under_reflection_max_abs": float(odd_u1_max_abs),
            "u1_odd_under_reflection": bool(odd_u1),
            "u1_coords_in_c1_s1": u1_coords,
            "aut_orbits_partition_Zn": bool(orbits_disjoint),
            "sum_u1_over_each_aut_orbit": orbit_sums,
            "all_aut_orbit_sums_zero": bool(all_orbit_sums_zero),
            "ord_z12_by_x_aut_invariant_if_present": ord_aut_invariant,
            "r_ord_weights_aut_invariant_if_present": r_ord_aut_invariant,
        },
        "computed_if_present": {
            "dot_u1_with_ord_z12_by_x": dot_u1_with_ord,
            "dot_u1_with_r_ord_weights": dot_u1_with_r_ord,
        },
        "meaning": [
            "Any Aut(Z_12)-invariant weight w(x) is even under reflection x↦-x (because -1∈Aut(Z_12)).",
            "For any odd u(x) (u(-x)=-u(x)), the scalar Σ_x w(x)u(x) cancels, so it cannot distinguish u from -u.",
            "Therefore direction-free (Aut-invariant) reference weights cannot supply a strict sign-distinction observable on pair1 of the form Σ_x w(x)u_1(x) on this exported instance.",
        ],
        "hard_limits": [
            "Obstruction record only; does not discharge H37 and does not export a sign-sensitive physical orientation datum.",
            "Does not claim strict-core selector closure / admissible S_sel_int.",
            "Does not claim global discharge of QW-2191.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    OUT_OBJ.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")

    summary = {
        "stage": "F472",
        "status": "PASS_EXPORTED_PAIR1_SIGN_DISTINCTION_OBSTRUCTION_OBJECT",
        "exported": {"obstruction_object": str(OUT_OBJ.relative_to(REPO))},
        "u1_odd_under_reflection_max_abs": float(odd_u1_max_abs),
        "all_aut_orbit_sums_zero": bool(all_orbit_sums_zero),
        "no_false_pass": True,
    }
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")

    print(OUT_OBJ)


if __name__ == "__main__":
    main()

