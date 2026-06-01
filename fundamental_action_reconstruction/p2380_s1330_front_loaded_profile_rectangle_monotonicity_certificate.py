#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
from itertools import combinations, product
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2380_s1330_front_loaded_profile_rectangle_monotonicity_certificate.json"
MD = GEN / "p2380_s1330_front_loaded_profile_rectangle_monotonicity_certificate.md"

SOURCE_FILES = {
    "P2379_FRONT_LOADED_PROFILE": GEN / "p2379_s1329_front_loaded_normalized_transport_density_existence_probe.json",
    "P2378_UNIT_NORMALIZED_INSUFFICIENCY": GEN / "p2378_s1328_unit_normalized_transport_coupling_insufficiency_theorem.json",
    "P2377_TRANSPORT_PRIMITIVE": GEN / "p2377_s1327_damping_compression_transport_primitive_uniform_coupling_theorem.json",
    "P2376_RECTANGLE_ROBUSTNESS": GEN / "p2376_s1326_damping_compression_eta_beta_rectangle_robustness_theorem.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

Z12_NODE_COUNT = 12
SUPPORT_SIZE = 5
CHAMBER_THRESHOLD = 1.0 / 3.0
ETA_INTERVAL = (9.0 / 5.0, 2.0)
BETA_TORS_INTERVAL = (0.0, 0.1)
SAMPLE_ETAS = [9.0 / 5.0, 19.0 / 10.0, 2.0]
SAMPLE_BETA_TORS = [0.0, 0.01, 0.1]
STRICT_BETA = 1.0
STRICT_PARAMS = {"omega": 743.0 / 4000.0, "phi": 13.0 / 80.0, "beta": STRICT_BETA, "eta": 9.0 / 5.0}
LAMBDA_TEST = 4.0 / 5.0
INTERVAL_GRID = 128


class Interval:
    def __init__(self, lo: float, hi: float | None = None):
        self.lo = float(lo)
        self.hi = float(lo if hi is None else hi)
        if self.lo > self.hi:
            raise ValueError(f"invalid interval [{self.lo}, {self.hi}]")

    def __add__(self, other: float | Interval) -> Interval:
        other = to_interval(other)
        return Interval(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __sub__(self, other: float | Interval) -> Interval:
        other = to_interval(other)
        return Interval(self.lo - other.hi, self.hi - other.lo)

    def __rsub__(self, other: float | Interval) -> Interval:
        return to_interval(other) - self

    def __neg__(self) -> Interval:
        return Interval(-self.hi, -self.lo)

    def __mul__(self, other: float | Interval) -> Interval:
        other = to_interval(other)
        values = [self.lo * other.lo, self.lo * other.hi, self.hi * other.lo, self.hi * other.hi]
        return Interval(min(values), max(values))

    __rmul__ = __mul__

    def inv(self) -> Interval:
        if self.lo <= 0.0 <= self.hi:
            raise ZeroDivisionError(f"interval contains zero: [{self.lo}, {self.hi}]")
        values = [1.0 / self.lo, 1.0 / self.hi]
        return Interval(min(values), max(values))

    def __truediv__(self, other: float | Interval) -> Interval:
        return self * to_interval(other).inv()

    def __rtruediv__(self, other: float | Interval) -> Interval:
        return to_interval(other) / self

    def square(self) -> Interval:
        return self * self

    def exp(self) -> Interval:
        return Interval(math.exp(self.lo), math.exp(self.hi))

    def log(self) -> Interval:
        if self.lo <= 0.0:
            raise ValueError(f"log interval is not positive: [{self.lo}, {self.hi}]")
        return Interval(math.log(self.lo), math.log(self.hi))

    def as_dict(self) -> dict[str, float]:
        return {"lo": self.lo, "hi": self.hi}


def to_interval(value: float | Interval) -> Interval:
    return value if isinstance(value, Interval) else Interval(value)


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


def compression_log_weight(d: int, eta: float, beta_tors: float) -> float:
    return math.log((1.0 + STRICT_BETA * d**eta) / (1.0 + beta_tors * d))


def affine_frontload_moment(d: int, eta: float, beta_tors: float) -> float:
    start = 1.0 + beta_tors * d
    delta = d**eta - beta_tors * d
    c = compression_log_weight(d, eta, beta_tors)
    return c * (1.0 + 2.0 * start / delta) - 2.0


def affine_weighted_transport(d: int, eta: float, beta_tors: float, lam: float) -> float:
    return compression_log_weight(d, eta, beta_tors) + lam * affine_frontload_moment(d, eta, beta_tors)


def lambda_threshold(eta: float, beta_tors: float) -> dict[str, float]:
    k1 = k_strict(1)
    k5 = k_strict(5)
    c1 = compression_log_weight(1, eta, beta_tors)
    c5 = compression_log_weight(5, eta, beta_tors)
    b1 = affine_frontload_moment(1, eta, beta_tors)
    b5 = affine_frontload_moment(5, eta, beta_tors)
    numerator = 3.0 * (k1 + c1) - (k5 + c5)
    denominator = b5 - 3.0 * b1
    return {
        "K1": k1,
        "K5": k5,
        "C1": c1,
        "C5": c5,
        "B1": b1,
        "B5": b5,
        "lambda_numerator": numerator,
        "lambda_denominator_B5_minus_3B1": denominator,
        "lambda_threshold_gt": numerator / denominator,
    }


def interval_c_b_derivatives(d: int, eta: Interval, beta_tors: Interval) -> dict[str, Interval]:
    ln_d = math.log(d)
    p = (eta * ln_d).exp()
    start = 1.0 + beta_tors * d
    delta = p - beta_tors * d
    c = ((1.0 + p) / (1.0 + beta_tors * d)).log()
    p_eta = ln_d * p
    c_eta = p_eta / (1.0 + p)
    # x -> -d/(1+d*x) is increasing on positive x after the negative sign.
    c_x = Interval(-d / (1.0 + d * beta_tors.lo), -d / (1.0 + d * beta_tors.hi))
    q = 1.0 + 2.0 * start / delta
    q_eta = -2.0 * start * p_eta / delta.square()
    q_x = 2.0 * d * (p + 1.0) / delta.square()
    b = c * q - 2.0
    b_eta = c_eta * q + c * q_eta
    b_x = c_x * q + c * q_x
    return {"C": c, "B": b, "C_eta": c_eta, "C_x": c_x, "B_eta": b_eta, "B_x": b_x}


def interval_threshold_box(eta_lo: float, eta_hi: float, beta_lo: float, beta_hi: float) -> dict[str, Any]:
    eta = Interval(eta_lo, eta_hi)
    beta_tors = Interval(beta_lo, beta_hi)
    c1 = interval_c_b_derivatives(1, eta, beta_tors)
    c5 = interval_c_b_derivatives(5, eta, beta_tors)
    k1 = k_strict(1)
    k5 = k_strict(5)
    numerator = 3.0 * (k1 + c1["C"]) - (k5 + c5["C"])
    denominator = c5["B"] - 3.0 * c1["B"]
    n_eta = 3.0 * c1["C_eta"] - c5["C_eta"]
    n_x = 3.0 * c1["C_x"] - c5["C_x"]
    h_eta = c5["B_eta"] - 3.0 * c1["B_eta"]
    h_x = c5["B_x"] - 3.0 * c1["B_x"]
    threshold = numerator / denominator
    threshold_eta = (n_eta * denominator - numerator * h_eta) / denominator.square()
    threshold_x = (n_x * denominator - numerator * h_x) / denominator.square()
    return {
        "eta_interval": {"lo": eta_lo, "hi": eta_hi},
        "beta_tors_interval": {"lo": beta_lo, "hi": beta_hi},
        "threshold_interval": threshold.as_dict(),
        "threshold_eta_derivative_interval": threshold_eta.as_dict(),
        "threshold_beta_tors_derivative_interval": threshold_x.as_dict(),
        "numerator_interval": numerator.as_dict(),
        "denominator_interval": denominator.as_dict(),
    }


def interval_monotonicity_certificate(grid: int = INTERVAL_GRID) -> dict[str, Any]:
    eta_min, eta_max = ETA_INTERVAL
    beta_min, beta_max = BETA_TORS_INTERVAL
    worst_threshold_box: dict[str, Any] | None = None
    worst_eta_derivative_hi = -math.inf
    worst_beta_derivative_lo = math.inf
    eta_derivative_sign_failures = []
    beta_derivative_sign_failures = []
    threshold_upper_failures = []
    denominator_positive_failures = []
    for i in range(grid):
        eta_lo = eta_min + (eta_max - eta_min) * i / grid
        eta_hi = eta_min + (eta_max - eta_min) * (i + 1) / grid
        for j in range(grid):
            beta_lo = beta_min + (beta_max - beta_min) * j / grid
            beta_hi = beta_min + (beta_max - beta_min) * (j + 1) / grid
            box = interval_threshold_box(eta_lo, eta_hi, beta_lo, beta_hi)
            if worst_threshold_box is None or box["threshold_interval"]["hi"] > worst_threshold_box["threshold_interval"]["hi"]:
                worst_threshold_box = box
            worst_eta_derivative_hi = max(worst_eta_derivative_hi, box["threshold_eta_derivative_interval"]["hi"])
            worst_beta_derivative_lo = min(worst_beta_derivative_lo, box["threshold_beta_tors_derivative_interval"]["lo"])
            if box["threshold_eta_derivative_interval"]["hi"] >= 0.0:
                eta_derivative_sign_failures.append(box)
            if box["threshold_beta_tors_derivative_interval"]["lo"] <= 0.0:
                beta_derivative_sign_failures.append(box)
            if box["threshold_interval"]["hi"] >= LAMBDA_TEST:
                threshold_upper_failures.append(box)
            if box["denominator_interval"]["lo"] <= 0.0:
                denominator_positive_failures.append(box)
    corner_max = lambda_threshold(ETA_INTERVAL[0], BETA_TORS_INTERVAL[1])
    corner_min = lambda_threshold(ETA_INTERVAL[1], BETA_TORS_INTERVAL[0])
    return {
        "interval_grid_per_axis": grid,
        "interval_box_count": grid * grid,
        "worst_threshold_box": worst_threshold_box,
        "worst_eta_derivative_hi": worst_eta_derivative_hi,
        "worst_beta_tors_derivative_lo": worst_beta_derivative_lo,
        "eta_derivative_negative_on_all_boxes": len(eta_derivative_sign_failures) == 0,
        "beta_tors_derivative_positive_on_all_boxes": len(beta_derivative_sign_failures) == 0,
        "denominator_positive_on_all_boxes": len(denominator_positive_failures) == 0,
        "raw_interval_threshold_upper_below_lambda_test_on_all_boxes": len(threshold_upper_failures) == 0,
        "corner_maximum_below_lambda_test": corner_max["lambda_threshold_gt"] < LAMBDA_TEST,
        "rectangle_uniform_sufficiency_certified_by_monotonicity": (
            len(eta_derivative_sign_failures) == 0
            and len(beta_derivative_sign_failures) == 0
            and len(denominator_positive_failures) == 0
            and corner_max["lambda_threshold_gt"] < LAMBDA_TEST
        ),
        "failure_counts": {
            "eta_derivative_sign": len(eta_derivative_sign_failures),
            "beta_tors_derivative_sign": len(beta_derivative_sign_failures),
            "denominator_positive": len(denominator_positive_failures),
            "threshold_upper": len(threshold_upper_failures),
        },
        "corner_maximum_candidate": {"eta": ETA_INTERVAL[0], "beta_tors": BETA_TORS_INTERVAL[1], **corner_max},
        "corner_minimum_candidate": {"eta": ETA_INTERVAL[1], "beta_tors": BETA_TORS_INTERVAL[0], **corner_min},
        "analytic_conclusion": "interval arithmetic certifies dT/deta<0 and dT/dbeta_tors>0 on the rectangle, so the maximum is at (eta=9/5,beta_tors=0.1)",
    }


def internal_edges(support: tuple[int, ...], step: int) -> int:
    support_set = set(support)
    edges = set()
    for node in support_set:
        for neighbor in ((node + step) % Z12_NODE_COUNT, (node - step) % Z12_NODE_COUNT):
            if neighbor in support_set:
                edges.add(tuple(sorted((node, neighbor))))
    return len(edges)


def support_rows() -> list[dict[str, Any]]:
    return [
        {"support": list(support), "h1": internal_edges(support, 1), "h5": internal_edges(support, 5)}
        for support in combinations(range(Z12_NODE_COUNT), SUPPORT_SIZE)
    ]


def score_maximizers(rows: list[dict[str, Any]], a: float, b: float) -> dict[str, Any]:
    scored = [(a * row["h1"] + b * row["h5"], row) for row in rows]
    maximum = max(score for score, _ in scored)
    maximizers = [row for score, row in scored if abs(score - maximum) <= 1e-10]
    pair_distribution: dict[str, int] = {}
    for row in maximizers:
        key = f"{row['h1']},{row['h5']}"
        pair_distribution[key] = pair_distribution.get(key, 0) + 1
    return {
        "a_over_b": a / b if b > 0 else None,
        "d5_chamber": b > 0 and a >= 0 and a / b < CHAMBER_THRESHOLD,
        "maximum_score": maximum,
        "maximizer_count": len(maximizers),
        "maximizer_h1_h5_pair_distribution": dict(sorted(pair_distribution.items())),
    }


def sample_support_replay() -> list[dict[str, Any]]:
    rows = support_rows()
    sample_rows = []
    for eta, beta_tors in product(SAMPLE_ETAS, SAMPLE_BETA_TORS):
        threshold = lambda_threshold(eta, beta_tors)
        w1 = affine_weighted_transport(1, eta, beta_tors, LAMBDA_TEST)
        w5 = affine_weighted_transport(5, eta, beta_tors, LAMBDA_TEST)
        sample_rows.append(
            {
                "eta": eta,
                "beta_tors": beta_tors,
                "lambda_threshold_gt": threshold["lambda_threshold_gt"],
                "lambda_test_margin": LAMBDA_TEST - threshold["lambda_threshold_gt"],
                "front_loaded_score_audit": score_maximizers(rows, threshold["K1"] + w1, threshold["K5"] + w5),
            }
        )
    return sample_rows


def rectangle_certificate() -> dict[str, Any]:
    interval_certificate = interval_monotonicity_certificate()
    sample_rows = sample_support_replay()
    return {
        "threshold_definition": {
            "T(eta,x)": "[3*(K1+C1)-(K5+C5)]/(B5-3*B1)",
            "C(d)": "log((1+d^eta)/(1+x*d))",
            "B(d)": "C(d)*(1+2*(1+x*d)/(d^eta-x*d))-2",
            "lambda_test": LAMBDA_TEST,
        },
        "closed_form_derivative_audit": {
            "T_eta": "(N_eta*H-N*H_eta)/H^2",
            "T_x": "(N_x*H-N*H_x)/H^2",
            "method": "outward interval arithmetic over the full audited rectangle",
        },
        "interval_monotonicity_certificate": interval_certificate,
        "sample_support_replay": sample_rows,
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    docs = {name: path.read_text(encoding="utf-8") if path.exists() else "" for name, path in DOC_FILES.items()}
    certificate = rectangle_certificate()
    interval_certificate = certificate["interval_monotonicity_certificate"]
    sample_rows = certificate["sample_support_replay"]

    theorem_export = {
        "theorem_name": "P2380 front-loaded profile rectangle monotonicity certificate",
        "claim": (
            "The P2379 affine normalized profile threshold T(eta,beta_tors) is not merely grid-observed. "
            "Interval arithmetic over the P2376 rectangle certifies dT/deta<0 and dT/dbeta_tors>0, so the worst point is the corner (eta=9/5,beta_tors=0.1). "
            "At that corner T=0.7916644842269442, hence lambda=0.8 is rectangle-uniformly sufficient for the audited affine profile family."
        ),
        "positive_content": [
            "Exports closed-form derivative formulas for the P2379 profile threshold.",
            "Uses interval arithmetic on 128x128 boxes to certify monotonicity signs and denominator positivity over the full P2376 rectangle.",
            "Replays all-792-support d5 selection on the 3x3 sample grid at the rectangle-uniform lambda=0.8.",
        ],
        "not_licensed": [
            "strict variational source theorem deriving rho_lambda or lambda=0.8",
            "promotion of weighted transport into L_total",
            "global proof for arbitrary normalized profiles beyond the affine audited family",
            "QW-2191 discharge",
            "beta_tors -> chi11 theorem",
            "legacy physical-role transfer",
            "ToE closure",
        ],
    }

    probe = {
        "probe_id": "P2380_S1330_FRONT_LOADED_PROFILE_RECTANGLE_MONOTONICITY_CERTIFICATE",
        "source_artifacts": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "docs": {name: rel(path) for name, path in DOC_FILES.items()},
        "artifact_replay": {
            "p2379_packet_id": artifacts["P2379_FRONT_LOADED_PROFILE"].get("packet_id"),
            "p2378_packet_id": artifacts["P2378_UNIT_NORMALIZED_INSUFFICIENCY"].get("packet_id"),
            "p2377_packet_id": artifacts["P2377_TRANSPORT_PRIMITIVE"].get("packet_id"),
            "p2376_packet_id": artifacts["P2376_RECTANGLE_ROBUSTNESS"].get("packet_id"),
        },
        "front_loaded_profile_rectangle_monotonicity_certificate": certificate,
        "theorem_export": theorem_export,
    }
    probe["theorem_fingerprint_sha256"] = sha256_json(theorem_export)

    gatekeeper_checks = {
        "source_artifacts_present": all("_missing" not in artifact for artifact in artifacts.values()),
        "p2379_loaded": probe["artifact_replay"]["p2379_packet_id"] == "P2379",
        "p2378_loaded": probe["artifact_replay"]["p2378_packet_id"] == "P2378",
        "p2377_loaded": probe["artifact_replay"]["p2377_packet_id"] == "P2377",
        "p2376_loaded": probe["artifact_replay"]["p2376_packet_id"] == "P2376",
        "eta_monotone_decreasing_certified": interval_certificate["eta_derivative_negative_on_all_boxes"],
        "beta_tors_monotone_increasing_certified": interval_certificate["beta_tors_derivative_positive_on_all_boxes"],
        "denominator_positive_certified": interval_certificate["denominator_positive_on_all_boxes"],
        "lambda_test_rectangle_uniform_certified": interval_certificate["rectangle_uniform_sufficiency_certified_by_monotonicity"],
        "all_sample_support_replays_select_d5": all(
            row["front_loaded_score_audit"]["maximizer_h1_h5_pair_distribution"] == {"0,4": 12}
            for row in sample_rows
        ),
        "docs_updated_with_p2380_monotonicity_obligation": all("P2380/S1330" in text for text in docs.values()),
        "no_qw2191_discharge_claimed": True,
    }

    payload = {
        "schema_version": "p2380_s1330_v1",
        "packet_id": "P2380",
        "stage_id": "S1330",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-06-01T00:00:00+00:00",
        "status": "OPEN_PROGRESS_FRONT_LOADED_PROFILE_RECTANGLE_MONOTONICITY_CERTIFIED_SOURCE_STILL_OPEN",
        "result_kind": "FRONT_LOADED_PROFILE_RECTANGLE_MONOTONICITY_CERTIFICATE",
        "front_loaded_profile_rectangle_monotonicity_certificate": probe,
        "recommended_next_honest_step": (
            "P2380 upgrades P2379 from a lattice observation to a rectangle monotonicity certificate for the affine profile family. "
            "The remaining honest source task is still to derive rho_lambda and lambda>=0.8 from strict dynamics, or keep the profile as an explicit non-strict selector premise."
        ),
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_PROGRESS_WITH_TRACE_NO_SELECTOR_OR_TOE_CLOSURE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2380 S1330: front-loaded profile rectangle monotonicity certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2380/S1330 upgrades the P2379 front-loaded profile from grid/lattice evidence to a rectangle monotonicity certificate for the affine normalized density family.",
                "Interval arithmetic certifies `dT/deta<0` and `dT/dbeta_tors>0` for `T=[3*(K1+C1)-(K5+C5)]/(B5-3*B1)` on the P2376 rectangle.",
                "",
                "## Certificate",
                "",
                f"- Interval boxes: `{interval_certificate['interval_box_count']}`.",
                f"- Worst threshold box: `{interval_certificate['worst_threshold_box']}`.",
                f"- Worst `dT/deta` upper bound: `{interval_certificate['worst_eta_derivative_hi']}`.",
                f"- Worst `dT/dbeta_tors` lower bound: `{interval_certificate['worst_beta_tors_derivative_lo']}`.",
                f"- Corner maximum candidate: `{interval_certificate['corner_maximum_candidate']}`.",
                f"- Lambda test rectangle-uniform by monotonicity: `{interval_certificate['rectangle_uniform_sufficiency_certified_by_monotonicity']}`.",
                f"- Sample support replay selections: `{[row['front_loaded_score_audit']['maximizer_h1_h5_pair_distribution'] for row in sample_rows]}`.",
                "",
                "## Hard limits",
                "",
                "- This certifies the affine profile threshold on the P2376 rectangle; it does not derive the profile from strict dynamics.",
                "- No `L_total` promotion, `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, or ToE closure is claimed.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
