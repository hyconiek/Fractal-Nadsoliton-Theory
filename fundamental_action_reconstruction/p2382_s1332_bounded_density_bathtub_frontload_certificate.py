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

OUT = GEN / "p2382_s1332_bounded_density_bathtub_frontload_certificate.json"
MD = GEN / "p2382_s1332_bounded_density_bathtub_frontload_certificate.md"

SOURCE_FILES = {
    "P2381_AFFINE_SOURCE_OBLIGATION": GEN / "p2381_s1331_affine_frontload_source_obligation_theorem.json",
    "P2380_RECTANGLE_MONOTONICITY": GEN / "p2380_s1330_front_loaded_profile_rectangle_monotonicity_certificate.json",
    "P2379_FRONT_LOADED_PROFILE": GEN / "p2379_s1329_front_loaded_normalized_transport_density_existence_probe.json",
    "P2377_TRANSPORT_PRIMITIVE": GEN / "p2377_s1327_damping_compression_transport_primitive_uniform_coupling_theorem.json",
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
CAP_TEST = 8.0 / 5.0
CAP_NEGATIVE_EPSILON = 1.0e-9
INTERVAL_GRID = 128
QPRIME_GRID = 32


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


def u0(d: int, beta_tors: float) -> float:
    return 1.0 + beta_tors * d


def delta(d: int, eta: float, beta_tors: float) -> float:
    return d**eta - beta_tors * d


def transport_speed(d: int, eta: float, beta_tors: float, s: float) -> float:
    return delta(d, eta, beta_tors) / (u0(d, beta_tors) + s * delta(d, eta, beta_tors))


def q_contrast(eta: float, beta_tors: float, s: float) -> float:
    return transport_speed(5, eta, beta_tors, s) - 3.0 * transport_speed(1, eta, beta_tors, s)


def q_contrast_prime(eta: float, beta_tors: float, s: float) -> float:
    d1 = delta(1, eta, beta_tors)
    d5 = delta(5, eta, beta_tors)
    u1s = u0(1, beta_tors) + s * d1
    u5s = u0(5, beta_tors) + s * d5
    return 3.0 * d1 * d1 / (u1s * u1s) - d5 * d5 / (u5s * u5s)


def capped_bangbang_weight(d: int, eta: float, beta_tors: float, cap: float) -> float:
    if cap < 1.0:
        raise ValueError("cap must be at least one for a unit-mass density on [0,1]")
    t = 1.0 / cap
    start = u0(d, beta_tors)
    return cap * math.log((start + t * delta(d, eta, beta_tors)) / start)


def capped_chamber_margin(eta: float, beta_tors: float, cap: float) -> dict[str, Any]:
    w1 = capped_bangbang_weight(1, eta, beta_tors, cap)
    w5 = capped_bangbang_weight(5, eta, beta_tors, cap)
    a = k_strict(1) + w1
    b = k_strict(5) + w5
    margin = b - 3.0 * a
    return {
        "eta": eta,
        "beta_tors": beta_tors,
        "density_cap_M": cap,
        "early_interval_length_1_over_M": 1.0 / cap,
        "W1": w1,
        "W5": w5,
        "a": a,
        "b": b,
        "d5_margin_b_minus_3a": margin,
        "a_over_b": a / b,
        "d5_chamber": b > 0.0 and a >= 0.0 and a / b < CHAMBER_THRESHOLD,
    }


def target_contrast() -> float:
    return 3.0 * k_strict(1) - k_strict(5)


def cap_threshold(eta: float, beta_tors: float) -> dict[str, float]:
    target = target_contrast()
    lo = 1.0
    hi = 2.0
    while capped_bangbang_weight(5, eta, beta_tors, hi) - 3.0 * capped_bangbang_weight(1, eta, beta_tors, hi) <= target:
        hi *= 2.0
    for _ in range(96):
        mid = 0.5 * (lo + hi)
        contrast = capped_bangbang_weight(5, eta, beta_tors, mid) - 3.0 * capped_bangbang_weight(1, eta, beta_tors, mid)
        if contrast > target:
            hi = mid
        else:
            lo = mid
    return {
        "eta": eta,
        "beta_tors": beta_tors,
        "cap_threshold_gt": hi,
        "early_interval_length_lt": 1.0 / hi,
        "early_support_fraction_lt": 1.0 / hi,
        "endpoint_cap_over_uniform": hi,
    }


def interval_weight(d: int, eta: Interval, beta_tors: Interval, cap: float) -> Interval:
    p = (eta * math.log(d)).exp()
    start = 1.0 + beta_tors * d
    dlt = p - beta_tors * d
    t = 1.0 / cap
    return cap * ((start + t * dlt) / start).log()


def interval_margin_box(eta_lo: float, eta_hi: float, beta_lo: float, beta_hi: float, cap: float) -> dict[str, Any]:
    eta = Interval(eta_lo, eta_hi)
    beta_tors = Interval(beta_lo, beta_hi)
    w1 = interval_weight(1, eta, beta_tors, cap)
    w5 = interval_weight(5, eta, beta_tors, cap)
    margin = (k_strict(5) + w5) - 3.0 * (k_strict(1) + w1)
    return {
        "eta_interval": {"lo": eta_lo, "hi": eta_hi},
        "beta_tors_interval": {"lo": beta_lo, "hi": beta_hi},
        "margin_interval": margin.as_dict(),
        "W1_interval": w1.as_dict(),
        "W5_interval": w5.as_dict(),
    }


def interval_qprime_box(eta_lo: float, eta_hi: float, beta_lo: float, beta_hi: float, s_lo: float, s_hi: float) -> dict[str, Any]:
    eta = Interval(eta_lo, eta_hi)
    beta_tors = Interval(beta_lo, beta_hi)
    s = Interval(s_lo, s_hi)
    p5 = (eta * math.log(5)).exp()
    d1 = 1.0 - beta_tors
    d5 = p5 - 5.0 * beta_tors
    u1s = 1.0 + beta_tors + s * d1
    u5s = 1.0 + 5.0 * beta_tors + s * d5
    qprime = 3.0 * (d1 / u1s).square() - (d5 / u5s).square()
    return {
        "eta_interval": {"lo": eta_lo, "hi": eta_hi},
        "beta_tors_interval": {"lo": beta_lo, "hi": beta_hi},
        "s_interval": {"lo": s_lo, "hi": s_hi},
        "qprime_interval": qprime.as_dict(),
    }


def rectangle_cap_certificate(grid: int = INTERVAL_GRID) -> dict[str, Any]:
    eta_min, eta_max = ETA_INTERVAL
    beta_min, beta_max = BETA_TORS_INTERVAL
    worst_margin_box: dict[str, Any] | None = None
    margin_failures = []
    for i in range(grid):
        eta_lo = eta_min + (eta_max - eta_min) * i / grid
        eta_hi = eta_min + (eta_max - eta_min) * (i + 1) / grid
        for j in range(grid):
            beta_lo = beta_min + (beta_max - beta_min) * j / grid
            beta_hi = beta_min + (beta_max - beta_min) * (j + 1) / grid
            box = interval_margin_box(eta_lo, eta_hi, beta_lo, beta_hi, CAP_TEST)
            if worst_margin_box is None or box["margin_interval"]["lo"] < worst_margin_box["margin_interval"]["lo"]:
                worst_margin_box = box
            if box["margin_interval"]["lo"] <= 0.0:
                margin_failures.append(box)
    corner = capped_chamber_margin(ETA_INTERVAL[0], BETA_TORS_INTERVAL[1], CAP_TEST)
    return {
        "interval_grid_per_axis": grid,
        "interval_box_count": grid * grid,
        "density_cap_test_M": CAP_TEST,
        "worst_margin_box": worst_margin_box,
        "strictly_positive_margin_on_all_boxes": len(margin_failures) == 0,
        "failure_count": len(margin_failures),
        "worst_corner_replay": corner,
    }


def qprime_monotonicity_certificate(grid: int = QPRIME_GRID) -> dict[str, Any]:
    eta_min, eta_max = ETA_INTERVAL
    beta_min, beta_max = BETA_TORS_INTERVAL
    worst_qprime_box: dict[str, Any] | None = None
    failures = []
    for i in range(grid):
        eta_lo = eta_min + (eta_max - eta_min) * i / grid
        eta_hi = eta_min + (eta_max - eta_min) * (i + 1) / grid
        for j in range(grid):
            beta_lo = beta_min + (beta_max - beta_min) * j / grid
            beta_hi = beta_min + (beta_max - beta_min) * (j + 1) / grid
            for k in range(grid):
                s_lo = k / grid
                s_hi = (k + 1) / grid
                box = interval_qprime_box(eta_lo, eta_hi, beta_lo, beta_hi, s_lo, s_hi)
                if worst_qprime_box is None or box["qprime_interval"]["hi"] > worst_qprime_box["qprime_interval"]["hi"]:
                    worst_qprime_box = box
                if box["qprime_interval"]["hi"] >= 0.0:
                    failures.append(box)
    return {
        "interval_grid_per_axis": grid,
        "interval_box_count": grid * grid * grid,
        "worst_qprime_box": worst_qprime_box,
        "qprime_strictly_negative_on_all_boxes": len(failures) == 0,
        "failure_count": len(failures),
        "consequence": "q(s)=A5(s)-3*A1(s) is strictly decreasing on the audited rectangle, so the bathtub maximizer under 0<=rho<=M puts all permitted density at earliest transport time.",
    }


def support_rows() -> list[dict[str, Any]]:
    def internal_edges(support: tuple[int, ...], step: int) -> int:
        support_set = set(support)
        edges = set()
        for node in support_set:
            for neighbor in ((node + step) % Z12_NODE_COUNT, (node - step) % Z12_NODE_COUNT):
                if neighbor in support_set:
                    edges.add(tuple(sorted((node, neighbor))))
        return len(edges)

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


def sample_table() -> list[dict[str, Any]]:
    rows = support_rows()
    table = []
    for eta, beta_tors in product(SAMPLE_ETAS, SAMPLE_BETA_TORS):
        threshold = cap_threshold(eta, beta_tors)
        margin = capped_chamber_margin(eta, beta_tors, CAP_TEST)
        score = score_maximizers(rows, margin["a"], margin["b"])
        table.append(
            {
                **threshold,
                "cap_test_margin": CAP_TEST - threshold["cap_threshold_gt"],
                "cap_test_chamber_margin": margin["d5_margin_b_minus_3a"],
                "cap_test_support_score": score,
            }
        )
    return table


def dense_cap_threshold_audit(points: int = 81) -> dict[str, Any]:
    eta_min, eta_max = ETA_INTERVAL
    beta_min, beta_max = BETA_TORS_INTERVAL
    max_row: dict[str, float] | None = None
    for i in range(points):
        eta = eta_min + (eta_max - eta_min) * i / (points - 1)
        for j in range(points):
            beta_tors = beta_min + (beta_max - beta_min) * j / (points - 1)
            row = cap_threshold(eta, beta_tors)
            if max_row is None or row["cap_threshold_gt"] > max_row["cap_threshold_gt"]:
                max_row = row
    return {
        "grid_points_per_axis": points,
        "evaluated_points": points * points,
        "max_cap_threshold_row": max_row,
        "cap_test_exceeds_all_grid_thresholds": bool(max_row and CAP_TEST > max_row["cap_threshold_gt"]),
    }


def midpoint_replay(eta: float, beta_tors: float, cap: float, steps: int = 20000) -> dict[str, float]:
    total_q = 0.0
    inv_steps = 1.0 / steps
    cutoff = 1.0 / cap
    for i in range(steps):
        s = (i + 0.5) * inv_steps
        rho = cap if s <= cutoff else 0.0
        total_q += rho * q_contrast(eta, beta_tors, s)
    closed_q = capped_bangbang_weight(5, eta, beta_tors, cap) - 3.0 * capped_bangbang_weight(1, eta, beta_tors, cap)
    return {
        "steps": steps,
        "midpoint_q_integral": total_q * inv_steps,
        "closed_form_q_integral": closed_q,
        "absolute_error": abs(total_q * inv_steps - closed_q),
    }


def build_certificate() -> dict[str, Any]:
    corner_eta = ETA_INTERVAL[0]
    corner_beta = BETA_TORS_INTERVAL[1]
    corner_threshold = cap_threshold(corner_eta, corner_beta)
    below_threshold = corner_threshold["cap_threshold_gt"] - CAP_NEGATIVE_EPSILON
    below_margin = capped_chamber_margin(corner_eta, corner_beta, below_threshold)
    positive_margin = capped_chamber_margin(corner_eta, corner_beta, CAP_TEST)
    rows = support_rows()
    return {
        "bathtub_variational_reduction": {
            "contrast_function": "q(s)=A_s(5)-3*A_s(1)",
            "qprime_formula": "q'(s)=3*Delta1^2/u1(s)^2-Delta5^2/u5(s)^2",
            "density_class": "0<=rho(s)<=M, integral_0^1 rho(s) ds=1",
            "maximizer_when_q_decreases": "rho_M(s)=M on 0<=s<=1/M and 0 otherwise",
            "closed_form_weight": "W_M(d)=M*log((1+beta_tors*d+(d^eta-beta_tors*d)/M)/(1+beta_tors*d))",
        },
        "qprime_monotonicity_certificate": qprime_monotonicity_certificate(),
        "rectangle_cap_acceptance_certificate": rectangle_cap_certificate(),
        "rectangle_worst_cap_source_obligation": {
            **corner_threshold,
            "below_threshold_negative_control": below_margin,
            "cap_test_positive_replay": positive_margin,
            "below_threshold_support_score": score_maximizers(rows, below_margin["a"], below_margin["b"]),
            "cap_test_support_score": score_maximizers(rows, positive_margin["a"], positive_margin["b"]),
        },
        "dense_cap_threshold_audit": dense_cap_threshold_audit(),
        "sample_local_cap_table": sample_table(),
        "midpoint_replay_at_worst_corner": midpoint_replay(corner_eta, corner_beta, CAP_TEST),
    }


def main() -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    docs = {name: path.read_text(encoding="utf-8") if path.exists() else "" for name, path in DOC_FILES.items()}
    certificate = build_certificate()
    obligation = certificate["rectangle_worst_cap_source_obligation"]

    theorem_export = {
        "theorem_name": "P2382 bounded-density bathtub frontload certificate",
        "claim": (
            "For the audited P2377 homotopy, the d5-vs-h1 decision for arbitrary normalized bounded densities reduces to a one-dimensional bathtub problem. "
            "Interval arithmetic certifies q(s)=A5(s)-3*A1(s) is decreasing on the P2376 rectangle, so the optimal M-capped unit density is the early bang-bang profile. "
            "At the worst corner (eta=9/5,beta_tors=0.1), the cap threshold is about 1.5748213574; the explicit cap M=1.6 is rectangle-uniformly sufficient."
        ),
        "positive_content": [
            "Replaces the affine-only profile with a variational reduction over all normalized densities satisfying 0<=rho<=M.",
            "Exports the closed-form capped bang-bang transport weight W_M(d).",
            "Provides interval acceptance for M=1.6 and a below-threshold negative control at the worst corner.",
        ],
        "not_licensed": [
            "strict variational source theorem deriving an M-capped bang-bang density from nadsoliton dynamics",
            "promotion of rho_M or W_M into L_total",
            "claim that the cap M source is strict-core supplied",
            "QW-2191 discharge",
            "beta_tors -> chi11 theorem",
            "legacy physical-role transfer",
            "ToE closure",
        ],
    }

    probe = {
        "probe_id": "P2382_S1332_BOUNDED_DENSITY_BATHTUB_FRONTLOAD_CERTIFICATE",
        "source_artifacts": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "docs": {name: rel(path) for name, path in DOC_FILES.items()},
        "artifact_replay": {
            "p2381_packet_id": artifacts["P2381_AFFINE_SOURCE_OBLIGATION"].get("packet_id"),
            "p2380_packet_id": artifacts["P2380_RECTANGLE_MONOTONICITY"].get("packet_id"),
            "p2379_packet_id": artifacts["P2379_FRONT_LOADED_PROFILE"].get("packet_id"),
            "p2377_packet_id": artifacts["P2377_TRANSPORT_PRIMITIVE"].get("packet_id"),
        },
        "bounded_density_bathtub_frontload_certificate": certificate,
        "theorem_export": theorem_export,
    }
    probe["theorem_fingerprint_sha256"] = sha256_json(theorem_export)

    gatekeeper_checks = {
        "source_artifacts_present": all("_missing" not in artifact for artifact in artifacts.values()),
        "p2381_loaded": probe["artifact_replay"]["p2381_packet_id"] == "P2381",
        "p2380_loaded": probe["artifact_replay"]["p2380_packet_id"] == "P2380",
        "p2379_loaded": probe["artifact_replay"]["p2379_packet_id"] == "P2379",
        "p2377_loaded": probe["artifact_replay"]["p2377_packet_id"] == "P2377",
        "qprime_negative_certified": certificate["qprime_monotonicity_certificate"]["qprime_strictly_negative_on_all_boxes"],
        "cap_test_rectangle_margin_positive": certificate["rectangle_cap_acceptance_certificate"]["strictly_positive_margin_on_all_boxes"],
        "worst_cap_threshold_below_cap_test": obligation["cap_threshold_gt"] < CAP_TEST,
        "below_threshold_negative_control_fails_d5": not obligation["below_threshold_negative_control"]["d5_chamber"],
        "cap_test_selects_d5_supports": obligation["cap_test_support_score"]["maximizer_h1_h5_pair_distribution"] == {"0,4": 12},
        "docs_updated_with_p2382_bathtub_certificate": all("P2382/S1332" in text for text in docs.values()),
        "no_qw2191_discharge_claimed": True,
    }

    payload = {
        "schema_version": "p2382_s1332_v1",
        "packet_id": "P2382",
        "stage_id": "S1332",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-06-01T00:00:00+00:00",
        "status": "OPEN_PROGRESS_BOUNDED_DENSITY_BATHTUB_FRONTLOAD_CERTIFICATE_SOURCE_STILL_OPEN",
        "result_kind": "BOUNDED_DENSITY_BATHTUB_FRONTLOAD_CERTIFICATE",
        "bounded_density_bathtub_frontload_theorem": probe,
        "recommended_next_honest_step": (
            "P2382 reduces arbitrary normalized M-bounded profile selection to an early bang-bang source obligation. "
            "A future strict theorem must derive a bounded source with M exceeding the certified cap threshold, or keep that cap/frontload as an explicit non-strict premise."
        ),
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_PROGRESS_WITH_TRACE_NO_SELECTOR_OR_TOE_CLOSURE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2382 S1332: bounded-density bathtub frontload certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2382/S1332 generalizes the P2379-P2381 affine profile work to all normalized densities with a pointwise cap `0<=rho<=M`.",
                "Because the audited contrast `q(s)=A_s(5)-3*A_s(1)` is strictly decreasing, the bathtub maximizer is the earliest bang-bang profile.",
                "",
                "## Certificate",
                "",
                f"- Worst-corner cap threshold: `{obligation['cap_threshold_gt']}`.",
                f"- Certified sufficient cap: `{CAP_TEST}`.",
                f"- Worst-corner early interval length for `M=1.6`: `{1.0 / CAP_TEST}`.",
                f"- Below-threshold negative control selects: `{obligation['below_threshold_support_score']['maximizer_h1_h5_pair_distribution']}`.",
                f"- Cap `M=1.6` replay selects: `{obligation['cap_test_support_score']['maximizer_h1_h5_pair_distribution']}`.",
                "",
                "## Hard limits",
                "",
                "- This is a bounded-density variational reduction/acceptance certificate, not a strict source theorem deriving the density cap or bang-bang profile.",
                "- No `L_total` promotion, `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, or ToE closure is claimed.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
