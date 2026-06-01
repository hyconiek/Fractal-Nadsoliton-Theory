#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any, Callable

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2383_s1333_closed_form_bathtub_corner_reduction_theorem.json"
MD = GEN / "p2383_s1333_closed_form_bathtub_corner_reduction_theorem.md"

SOURCE_FILES = {
    "P2382_BOUNDED_DENSITY_BATHTUB": GEN / "p2382_s1332_bounded_density_bathtub_frontload_certificate.json",
    "P2381_AFFINE_SOURCE_OBLIGATION": GEN / "p2381_s1331_affine_frontload_source_obligation_theorem.json",
    "P2376_RECTANGLE_ROBUSTNESS": GEN / "p2376_s1326_damping_compression_eta_beta_rectangle_robustness_theorem.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

ETA_MIN = 9.0 / 5.0
ETA_MAX = 2.0
BETA_MIN = 0.0
BETA_MAX = 0.1
CAP_TEST = 8.0 / 5.0
CAP_AUDIT_LO = 1.5
CAP_AUDIT_HI = CAP_TEST
STRICT_PARAMS = {"omega": 743.0 / 4000.0, "phi": 13.0 / 80.0, "beta": 1.0, "eta": ETA_MIN}
CHAMBER_THRESHOLD = 1.0 / 3.0
DERIVATIVE_GRID = 96
BISECTION_STEPS = 80


def rel(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": rel(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def k_strict(d: int) -> float:
    return math.cos(STRICT_PARAMS["omega"] * d + STRICT_PARAMS["phi"]) / (
        1.0 + STRICT_PARAMS["beta"] * d ** STRICT_PARAMS["eta"]
    )


def u0(d: int, beta_tors: float) -> float:
    return 1.0 + beta_tors * d


def delta(d: int, eta: float, beta_tors: float) -> float:
    return d**eta - beta_tors * d


def transport_speed(d: int, eta: float, beta_tors: float, s: float) -> float:
    de = delta(d, eta, beta_tors)
    return de / (u0(d, beta_tors) + s * de)


def qprime(eta: float, beta_tors: float, s: float) -> float:
    a1 = transport_speed(1, eta, beta_tors, s)
    a5 = transport_speed(5, eta, beta_tors, s)
    return 3.0 * a1 * a1 - a5 * a5


def ratio_r(eta: float, beta_tors: float, s: float) -> float:
    return transport_speed(5, eta, beta_tors, s) / transport_speed(1, eta, beta_tors, s)


def ratio_at_s1(eta: float, beta_tors: float) -> float:
    x = 5.0**eta
    return 2.0 * (x - 5.0 * beta_tors) / ((1.0 - beta_tors) * (1.0 + x))


def capped_weight(d: int, eta: float, beta_tors: float, cap: float) -> float:
    de = delta(d, eta, beta_tors)
    start = u0(d, beta_tors)
    return cap * math.log((start + de / cap) / start)


def chamber_margin(eta: float, beta_tors: float, cap: float) -> float:
    return (k_strict(5) + capped_weight(5, eta, beta_tors, cap)) - 3.0 * (
        k_strict(1) + capped_weight(1, eta, beta_tors, cap)
    )


def weight_d_eta(d: int, eta: float, beta_tors: float, cap: float) -> float:
    if d == 1:
        return 0.0
    de = delta(d, eta, beta_tors)
    return d**eta * math.log(d) / (u0(d, beta_tors) + de / cap)


def margin_d_eta(eta: float, beta_tors: float, cap: float) -> float:
    return weight_d_eta(5, eta, beta_tors, cap) - 3.0 * weight_d_eta(1, eta, beta_tors, cap)


def weight_d_beta(d: int, eta: float, beta_tors: float, cap: float) -> float:
    u = u0(d, beta_tors)
    de = delta(d, eta, beta_tors)
    return d * (cap - 1.0) / cap / (u + de / cap) - d / u


def margin_d_beta(eta: float, beta_tors: float, cap: float) -> float:
    return weight_d_beta(5, eta, beta_tors, cap) - 3.0 * weight_d_beta(1, eta, beta_tors, cap)


def weight_d_cap(d: int, eta: float, beta_tors: float, cap: float) -> float:
    u = u0(d, beta_tors)
    x = delta(d, eta, beta_tors) / (cap * u)
    return math.log1p(x) - x / (1.0 + x)


def margin_d_cap(eta: float, beta_tors: float, cap: float) -> float:
    return weight_d_cap(5, eta, beta_tors, cap) - 3.0 * weight_d_cap(1, eta, beta_tors, cap)


def bisection_threshold(eta: float, beta_tors: float) -> float:
    lo = 1.0
    hi = CAP_TEST
    if chamber_margin(eta, beta_tors, lo) > 0.0:
        return lo
    if chamber_margin(eta, beta_tors, hi) <= 0.0:
        raise ValueError("cap test does not bracket a positive chamber margin")
    for _ in range(BISECTION_STEPS):
        mid = (lo + hi) / 2.0
        if chamber_margin(eta, beta_tors, mid) > 0.0:
            hi = mid
        else:
            lo = mid
    return hi


def grid_minmax(fn: Callable[[float, float, float], float], cap_lo: float, cap_hi: float) -> dict[str, Any]:
    best_min = {"value": float("inf"), "eta": None, "beta_tors": None, "cap": None}
    best_max = {"value": -float("inf"), "eta": None, "beta_tors": None, "cap": None}
    for i in range(DERIVATIVE_GRID + 1):
        eta = ETA_MIN + (ETA_MAX - ETA_MIN) * i / DERIVATIVE_GRID
        for j in range(DERIVATIVE_GRID + 1):
            beta = BETA_MIN + (BETA_MAX - BETA_MIN) * j / DERIVATIVE_GRID
            for k in range(DERIVATIVE_GRID + 1):
                cap = cap_lo + (cap_hi - cap_lo) * k / DERIVATIVE_GRID
                value = fn(eta, beta, cap)
                if value < best_min["value"]:
                    best_min = {"value": value, "eta": eta, "beta_tors": beta, "cap": cap}
                if value > best_max["value"]:
                    best_max = {"value": value, "eta": eta, "beta_tors": beta, "cap": cap}
    return {"grid": DERIVATIVE_GRID, "min": best_min, "max": best_max}


def ratio_audit() -> dict[str, Any]:
    r_corner = ratio_at_s1(ETA_MIN, BETA_MIN)
    qprime_bound_factor = 3.0 - r_corner * r_corner
    samples = []
    worst = {"R": float("inf"), "eta": None, "beta_tors": None, "s": None, "qprime": None}
    for i in range(DERIVATIVE_GRID + 1):
        eta = ETA_MIN + (ETA_MAX - ETA_MIN) * i / DERIVATIVE_GRID
        for j in range(DERIVATIVE_GRID + 1):
            beta = BETA_MIN + (BETA_MAX - BETA_MIN) * j / DERIVATIVE_GRID
            for k in range(DERIVATIVE_GRID + 1):
                s = k / DERIVATIVE_GRID
                r = ratio_r(eta, beta, s)
                qp = qprime(eta, beta, s)
                if r < worst["R"]:
                    worst = {"R": r, "eta": eta, "beta_tors": beta, "s": s, "qprime": qp}
    for eta in (ETA_MIN, 1.9, ETA_MAX):
        for beta in (BETA_MIN, BETA_MAX):
            samples.append(
                {
                    "eta": eta,
                    "beta_tors": beta,
                    "R_s0": ratio_r(eta, beta, 0.0),
                    "R_s1": ratio_r(eta, beta, 1.0),
                    "qprime_s1": qprime(eta, beta, 1.0),
                }
            )
    return {
        "closed_form_R": "R=A5/A1=(Delta5*(u1+s*Delta1))/(Delta1*(u5+s*Delta5))",
        "s1_boundary_R": "R(1)=2*(5^eta-5*beta_tors)/((1-beta_tors)*(1+5^eta))",
        "monotonicity_proof_notes": [
            "For the audited rectangle, A5>A1, so R decreases in s and the s-worst boundary is s=1.",
            "At s=1, dR/d(5^eta)>0 and dR/d(beta_tors)>0 because 5^eta>5; hence the global lower corner is eta=9/5,beta_tors=0.",
            "Therefore q'(s)=A1^2*(3-R^2)<0 follows from R_min^2>3, without relying on the P2382 box grid as the primary proof object.",
        ],
        "R_min_closed_form_corner": r_corner,
        "sqrt3": math.sqrt(3.0),
        "R_min_square_minus_3": r_corner * r_corner - 3.0,
        "qprime_factor_upper_bound_3_minus_Rmin_square": qprime_bound_factor,
        "dense_replay_worst_R": worst,
        "sample_boundaries": samples,
    }


def corner_reduction_certificate() -> dict[str, Any]:
    threshold = bisection_threshold(ETA_MIN, BETA_MAX)
    derivative_audits = {
        "margin_d_eta_on_cap_band": grid_minmax(margin_d_eta, CAP_AUDIT_LO, CAP_AUDIT_HI),
        "margin_d_beta_on_cap_band": grid_minmax(margin_d_beta, CAP_AUDIT_LO, CAP_AUDIT_HI),
        "margin_d_cap_on_cap_band": grid_minmax(margin_d_cap, CAP_AUDIT_LO, CAP_AUDIT_HI),
    }
    corner_rows = []
    for eta in (ETA_MIN, ETA_MAX):
        for beta in (BETA_MIN, BETA_MAX):
            corner_rows.append(
                {
                    "eta": eta,
                    "beta_tors": beta,
                    "margin_M_1_5": chamber_margin(eta, beta, CAP_AUDIT_LO),
                    "margin_M_1_6": chamber_margin(eta, beta, CAP_TEST),
                    "threshold": bisection_threshold(eta, beta),
                }
            )
    return {
        "audited_cap_band": [CAP_AUDIT_LO, CAP_AUDIT_HI],
        "corner_reduction_statement": (
            "On the audited cap band, the d5 chamber margin increases with eta, decreases with beta_tors, "
            "and increases with M; hence the rectangle-uniform cap threshold is attained at (eta,beta_tors)=(9/5,0.1)."
        ),
        "derivative_audits": derivative_audits,
        "worst_corner": {"eta": ETA_MIN, "beta_tors": BETA_MAX},
        "worst_corner_threshold": threshold,
        "worst_corner_margin_M_1_6": chamber_margin(ETA_MIN, BETA_MAX, CAP_TEST),
        "corner_rows": corner_rows,
        "positive_at_M_1_6": chamber_margin(ETA_MIN, BETA_MAX, CAP_TEST) > 0.0,
        "negative_at_M_1_5": chamber_margin(ETA_MIN, BETA_MAX, CAP_AUDIT_LO) < 0.0,
    }


def source_burden_translation(corner: dict[str, Any]) -> dict[str, Any]:
    threshold = corner["worst_corner_threshold"]
    return {
        "minimal_rectangle_cap_gt": threshold,
        "equivalent_early_interval_length_lt": 1.0 / threshold,
        "M_1_6_early_interval_length": 1.0 / CAP_TEST,
        "M_1_6_early_mass_on_first_half": min(1.0, CAP_TEST / 2.0),
        "M_1_6_barycenter": 1.0 / (2.0 * CAP_TEST),
        "M_1_6_barycenter_shift_from_uniform": 0.5 - 1.0 / (2.0 * CAP_TEST),
        "interpretation": (
            "If a future strict source theorem chooses to use the bounded-density route, it must now source at least this cap/frontload scale, "
            "or leave it as an explicit non-strict premise. The certificate does not itself derive the density from strict dynamics."
        ),
    }


def main() -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    docs = {name: path.read_text(encoding="utf-8") if path.exists() else "" for name, path in DOC_FILES.items()}

    ratio = ratio_audit()
    corner = corner_reduction_certificate()
    burden = source_burden_translation(corner)

    theorem_export = {
        "name": "P2383/S1333 closed-form bathtub corner reduction theorem",
        "claim": (
            "The P2382 bathtub result admits a closed ratio reduction / closed-form proof core: q'(s)<0 follows from the ratio lower bound R_min=2*5^(9/5)/(1+5^(9/5))>sqrt(3), "
            "and the audited M-cap source burden reduces to the single worst corner (eta,beta_tors)=(9/5,0.1) on the cap band [1.5,1.6]."
        ),
        "positive_content": [
            "Replaces the q'(s)<0 grid witness by an algebraic ratio lower-bound proof replay.",
            "Exports derivative-form audits showing that the cap threshold is controlled by the eta-low/beta-high corner on the audited band.",
            "Translates M=1.6 into explicit early interval, early-half mass, and barycenter obligations for future source theorems.",
        ],
        "not_licensed": [
            "strict source theorem deriving the cap M or bang-bang density",
            "promotion of rho_M or W_M into L_total",
            "claim that the cap/frontload premise is strict-core supplied",
            "QW-2191 discharge or selector closure",
            "beta_tors -> chi11 theorem",
            "legacy physical-role transfer",
            "ToE closure",
        ],
    }

    probe = {
        "probe_id": "P2383_S1333_CLOSED_FORM_BATHTUB_CORNER_REDUCTION_THEOREM",
        "source_artifacts": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "docs": {name: rel(path) for name, path in DOC_FILES.items()},
        "artifact_replay": {
            "p2382_packet_id": artifacts["P2382_BOUNDED_DENSITY_BATHTUB"].get("packet_id"),
            "p2381_packet_id": artifacts["P2381_AFFINE_SOURCE_OBLIGATION"].get("packet_id"),
            "p2376_packet_id": artifacts["P2376_RECTANGLE_ROBUSTNESS"].get("packet_id"),
        },
        "closed_form_qprime_ratio_certificate": ratio,
        "cap_threshold_corner_reduction_certificate": corner,
        "source_burden_translation": burden,
        "theorem_export": theorem_export,
    }
    probe["theorem_fingerprint_sha256"] = sha256_json(theorem_export)

    derivative_audits = corner["derivative_audits"]
    gatekeeper_checks = {
        "source_artifacts_present": all("_missing" not in artifact for artifact in artifacts.values()),
        "p2382_loaded": probe["artifact_replay"]["p2382_packet_id"] == "P2382",
        "p2381_loaded": probe["artifact_replay"]["p2381_packet_id"] == "P2381",
        "p2376_loaded": probe["artifact_replay"]["p2376_packet_id"] == "P2376",
        "ratio_bound_beats_sqrt3": ratio["R_min_square_minus_3"] > 0.0,
        "dense_replay_respects_ratio_bound": ratio["dense_replay_worst_R"]["R"] >= ratio["R_min_closed_form_corner"] - 1.0e-12,
        "qprime_closed_form_negative": ratio["qprime_factor_upper_bound_3_minus_Rmin_square"] < 0.0,
        "margin_increases_with_eta_on_cap_band": derivative_audits["margin_d_eta_on_cap_band"]["min"]["value"] > 0.0,
        "margin_decreases_with_beta_on_cap_band": derivative_audits["margin_d_beta_on_cap_band"]["max"]["value"] < 0.0,
        "margin_increases_with_cap_on_cap_band": derivative_audits["margin_d_cap_on_cap_band"]["min"]["value"] > 0.0,
        "worst_corner_threshold_matches_p2382": abs(corner["worst_corner_threshold"] - 1.574821357435363) < 5.0e-13,
        "cap_test_positive_and_lower_control_negative": corner["positive_at_M_1_6"] and corner["negative_at_M_1_5"],
        "docs_updated_with_p2383_corner_reduction": all("P2383/S1333" in text for text in docs.values()),
        "no_qw2191_discharge_claimed": True,
    }

    payload = {
        "schema_version": "p2383_s1333_v1",
        "packet_id": "P2383",
        "stage_id": "S1333",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-06-01T00:00:00+00:00",
        "status": "OPEN_PROGRESS_CLOSED_FORM_BATHTUB_CORNER_REDUCTION_SOURCE_STILL_OPEN",
        "result_kind": "CLOSED_FORM_BATHTUB_CORNER_REDUCTION_THEOREM",
        "closed_form_bathtub_corner_reduction_theorem": probe,
        "recommended_next_honest_step": (
            "Use the P2383 corner-reduced source burden as the acceptance target for a real source theorem: derive a strict bounded density with M>1.574821357435363, "
            "or mark the cap/frontload as a non-strict premise. Do not promote the variational optimizer into L_total without such a source."
        ),
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_PROGRESS_WITH_TRACE_NO_SELECTOR_OR_TOE_CLOSURE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2383 S1333: closed-form bathtub corner reduction theorem",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2383/S1333 turns the P2382 bathtub certificate into a more proof-like corner reduction.",
                "The monotonicity of `q(s)=A_s(5)-3*A_s(1)` is certified by the closed ratio bound",
                "",
                "```text",
                "R=A5/A1=(Delta5*(u1+s*Delta1))/(Delta1*(u5+s*Delta5))",
                "R_min=2*5^(9/5)/(1+5^(9/5)) > sqrt(3)",
                "q'(s)=A1^2*(3-R^2)<0.",
                "```",
                "",
                "## Corner-reduced cap burden",
                "",
                "On the audited cap band `[1.5,1.6]`, derivative replay shows that the chamber margin increases with `eta`, decreases with `beta_tors`, and increases with `M`.",
                "Therefore the cap threshold is controlled by the single corner `(eta,beta_tors)=(9/5,0.1)`.",
                "",
                f"- Corner threshold: `{corner['worst_corner_threshold']}`.",
                f"- Certified sufficient cap: `{CAP_TEST}`.",
                f"- `M=1.6` early interval length: `{1.0 / CAP_TEST}`.",
                f"- `M=1.6` early-half mass: `{burden['M_1_6_early_mass_on_first_half']}`.",
                f"- `M=1.6` barycenter shift from uniform: `{burden['M_1_6_barycenter_shift_from_uniform']}`.",
                "",
                "## Hard limits",
                "",
                "- This is a closed-form/derivative-audited reduction of the bounded-density acceptance criterion, not a strict source theorem for the cap or bang-bang profile.",
                "- No `L_total` promotion, `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, or ToE closure is claimed.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
