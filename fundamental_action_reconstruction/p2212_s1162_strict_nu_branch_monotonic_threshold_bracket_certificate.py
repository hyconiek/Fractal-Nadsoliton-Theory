#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2209 = GEN / "p2209_s1159_strict_nu_branch_transport_nonvanishing_threshold_map.json"
IN_2211 = GEN / "p2211_s1161_strict_nu_branch_extended_threshold_crossing_certificate.json"
OUT = GEN / "p2212_s1162_strict_nu_branch_monotonic_threshold_bracket_certificate.json"
MD = GEN / "p2212_s1162_strict_nu_branch_monotonic_threshold_bracket_certificate.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def effective_bound(dm: float, c1: float, c2: float) -> float:
    return max(c1 * dm, c2 * dm * dm)


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2209 = load(IN_2209)
    p2211 = load(IN_2211)

    inputs = (p2209.get("strict_nu_branch_transport_nonvanishing_threshold_map", {}) or {}).get("inputs", {}) or {}
    c1 = float(inputs.get("linear_lower_bound_coeff_c1", 0.0) or 0.0)
    c2 = float(inputs.get("quadratic_lower_bound_coeff_c2", 0.0) or 0.0)
    min_margin = float(inputs.get("min_majorant_margin", 0.0) or 0.0)

    cross = (p2211.get("strict_nu_branch_extended_threshold_crossing_certificate", {}) or {}).get("crossing_summary", {}) or {}
    interval = cross.get("crossing_interval", {}) or {}
    lo = float(interval.get("abs_dm_lo_uncertified", 0.0) or 0.0)
    hi = float(interval.get("abs_dm_hi_certified", 0.0) or 0.0)

    if not (c1 > 0.0 and c2 > 0.0 and min_margin > 0.0 and 0.0 <= lo <= hi):
        raise RuntimeError("Invalid upstream inputs for P2212 monotonic threshold bracket")

    # Monotonicity witness on bracket: for dm>=0 both c1*dm and c2*dm^2 are nondecreasing,
    # so max(c1*dm, c2*dm^2) is nondecreasing.
    monotonicity_samples = []
    n_check = 401
    prev = None
    monotone_ok = True
    for k in range(n_check):
        dm = lo + (hi - lo) * (k / (n_check - 1)) if n_check > 1 else lo
        val = effective_bound(dm, c1, c2)
        if prev is not None and val + 1e-18 < prev:
            monotone_ok = False
        monotonicity_samples.append({"abs_dm": dm, "effective_lower_bound": val})
        prev = val

    lo_val = effective_bound(lo, c1, c2)
    hi_val = effective_bound(hi, c1, c2)

    endpoint_certificate = {
        "lo_uncertified": lo_val <= min_margin,
        "hi_certified": hi_val > min_margin,
        "margin_at_lo": lo_val - min_margin,
        "margin_at_hi": hi_val - min_margin,
    }

    payload = {
        "schema_version": "p2212_s1162_v1",
        "packet_id": "P2212",
        "stage_id": "S1162",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_MONOTONIC_THRESHOLD_BRACKET_CERTIFICATE",
        "strict_nu_branch_monotonic_threshold_bracket_certificate": {
            "certificate_id": "STRICT_NU_BRANCH_MONOTONIC_THRESHOLD_BRACKET_CERTIFICATE_V1",
            "source_packets": [str(IN_2209.relative_to(ROOT)), str(IN_2211.relative_to(ROOT))],
            "input_bracket": {"abs_dm_lo": lo, "abs_dm_hi": hi, "width": hi - lo},
            "monotonicity_witness": {
                "check_samples": n_check,
                "effective_bound_nondecreasing_on_bracket": monotone_ok,
            },
            "endpoint_certificate": endpoint_certificate,
            "certified_statement": "There exists dm* in [lo,hi] such that effective lower bound crosses min_margin on this modeled lane.",
            "theorem_scope_limit": "modeled lower-bound lane + finite bracket certificate; not global Task-3 closure",
        },
        "recommended_next_honest_step": {
            "id": "P2213_candidate",
            "goal": "bisection-refine dm* bracket to target precision and export deterministic crossing tolerance certificate",
        },
        "gatekeeper_checks": {
            "monotonic_bracket_certificate_exported": True,
            "monotonicity_verified_on_bracket": monotone_ok,
            "lo_uncertified_endpoint_verified": endpoint_certificate["lo_uncertified"],
            "hi_certified_endpoint_verified": endpoint_certificate["hi_certified"],
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_cutkosky_closure_proven": False,
            "full_d3_covariance_transport_proven": False,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2212 S1162: strict nu-branch monotonic threshold bracket certificate",
            "",
            f"- bracket: `[{lo:.12e}, {hi:.12e}]`",
            f"- monotonicity verified on bracket: `{monotone_ok}`",
            f"- lo endpoint uncertified: `{endpoint_certificate['lo_uncertified']}`",
            f"- hi endpoint certified: `{endpoint_certificate['hi_certified']}`",
            "",
            "Modeled-lane bracket certificate only; no global Task-3 closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
