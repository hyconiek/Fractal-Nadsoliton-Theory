#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2388_s1338_cap_threshold_root_uniqueness_certificate.json"
MD = GEN / "p2388_s1338_cap_threshold_root_uniqueness_certificate.md"

SOURCE_FILES = {
    "P2387_EXACT_KKT_BRANCH": GEN / "p2387_s1337_bathtub_exact_kkt_branch_certificate.json",
    "P2386_LP_DUAL": GEN / "p2386_s1336_bathtub_lp_dual_certificate.json",
    "P2384_SYMBOLIC_INEQUALITY": GEN / "p2384_s1334_symbolic_bathtub_inequality_proof_packet.json",
    "P2382_BOUNDED_DENSITY_BATHTUB": GEN / "p2382_s1332_bounded_density_bathtub_frontload_certificate.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

STRICT_PARAMS = {"omega": 743.0 / 4000.0, "phi": 13.0 / 80.0, "beta": 1.0, "eta": 9.0 / 5.0}
ETA_WORST = 9.0 / 5.0
BETA_TORS_WORST = 1.0 / 10.0
D1 = 1
D5 = 5
CAP_TEST = 8.0 / 5.0
ROOT_CENTER_FROM_P2382 = 1.574821357435363
ROOT_BRACKET_RADIUS = 1.0e-9
BISECTION_STEPS = 120
NEWTON_STEPS = 8
DERIVATIVE_GRID = 64


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


def u0(d: int) -> float:
    return 1.0 + BETA_TORS_WORST * d


def delta(d: int) -> float:
    return d**ETA_WORST - BETA_TORS_WORST * d


def capped_weight(d: int, cap: float) -> float:
    start = u0(d)
    return cap * math.log((start + delta(d) / cap) / start)


def threshold_target() -> float:
    return 3.0 * k_strict(D1) - k_strict(D5)


def threshold_function(cap: float) -> float:
    return capped_weight(D5, cap) - 3.0 * capped_weight(D1, cap) - threshold_target()


def weight_d_cap(d: int, cap: float) -> float:
    x = delta(d) / (cap * u0(d))
    return math.log1p(x) - x / (1.0 + x)


def threshold_derivative(cap: float) -> float:
    return weight_d_cap(D5, cap) - 3.0 * weight_d_cap(D1, cap)


def bisection_root(lo: float, hi: float) -> dict[str, Any]:
    f_lo = threshold_function(lo)
    f_hi = threshold_function(hi)
    if not f_lo < 0.0 < f_hi:
        raise ValueError(f"root not bracketed: f_lo={f_lo}, f_hi={f_hi}")
    left = lo
    right = hi
    for _ in range(BISECTION_STEPS):
        mid = 0.5 * (left + right)
        if threshold_function(mid) > 0.0:
            right = mid
        else:
            left = mid
    center = 0.5 * (left + right)
    return {
        "lo": left,
        "hi": right,
        "mid": center,
        "width": right - left,
        "f_lo": threshold_function(left),
        "f_hi": threshold_function(right),
        "f_mid": threshold_function(center),
        "steps": BISECTION_STEPS,
    }


def newton_replay(seed: float) -> dict[str, Any]:
    cap = seed
    rows = []
    for step in range(NEWTON_STEPS):
        value = threshold_function(cap)
        derivative = threshold_derivative(cap)
        rows.append({"step": step, "cap": cap, "F_cap": value, "F_prime_cap": derivative})
        cap = cap - value / derivative
    final_value = threshold_function(cap)
    return {"seed": seed, "steps": rows, "final_cap": cap, "final_abs_residual": abs(final_value)}


def derivative_audit(lo: float, hi: float) -> dict[str, Any]:
    min_row = {"cap": None, "F_prime": float("inf")}
    max_row = {"cap": None, "F_prime": -float("inf")}
    for index in range(DERIVATIVE_GRID + 1):
        cap = lo + (hi - lo) * index / DERIVATIVE_GRID
        derivative = threshold_derivative(cap)
        if derivative < min_row["F_prime"]:
            min_row = {"cap": cap, "F_prime": derivative}
        if derivative > max_row["F_prime"]:
            max_row = {"cap": cap, "F_prime": derivative}
    return {
        "grid_points": DERIVATIVE_GRID + 1,
        "min_derivative_on_bracket_grid": min_row,
        "max_derivative_on_bracket_grid": max_row,
        "all_sampled_positive": min_row["F_prime"] > 0.0,
        "proof_role": "regression audit only; theorem monotonicity is inherited from P2384 cap-derivative sign proof",
    }


def p2382_threshold_replay(artifact: dict[str, Any]) -> dict[str, Any]:
    try:
        obligation = artifact["bounded_density_bathtub_frontload_theorem"]["bounded_density_bathtub_frontload_certificate"][
            "rectangle_worst_cap_source_obligation"
        ]
        return {
            "p2382_cap_threshold_gt": obligation["cap_threshold_gt"],
            "p2382_cap_test_M": obligation["cap_test_positive_replay"]["density_cap_M"],
            "p2382_cap_test_margin": obligation["cap_test_positive_replay"]["d5_margin_b_minus_3a"],
            "p2382_below_threshold_fails": obligation["below_threshold_negative_control"]["d5_chamber"] is False,
        }
    except KeyError:
        return {"status": "P2382_REPLAY_FIELDS_MISSING"}


def cap_threshold_root_certificate() -> dict[str, Any]:
    bracket_lo = ROOT_CENTER_FROM_P2382 - ROOT_BRACKET_RADIUS
    bracket_hi = ROOT_CENTER_FROM_P2382 + ROOT_BRACKET_RADIUS
    bracket = {
        "lo": bracket_lo,
        "hi": bracket_hi,
        "F_lo": threshold_function(bracket_lo),
        "F_hi": threshold_function(bracket_hi),
        "sign_change": threshold_function(bracket_lo) < 0.0 < threshold_function(bracket_hi),
    }
    bisect = bisection_root(bracket_lo, bracket_hi)
    newton = newton_replay(ROOT_CENTER_FROM_P2382)
    audit = derivative_audit(bracket_lo, bracket_hi)
    root = bisect["mid"]
    cap_test_margin = threshold_function(CAP_TEST)
    return {
        "threshold_equation": "F(M)=W_M(5)-3*W_M(1)-(3*K_strict(1)-K_strict(5))=0 at eta=9/5,beta_tors=0.1",
        "strict_kernel_parameters_used_for_target": STRICT_PARAMS,
        "worst_corner": {"eta": ETA_WORST, "beta_tors": BETA_TORS_WORST},
        "target_3K1_minus_K5": threshold_target(),
        "initial_sign_bracket_around_p2382_value": bracket,
        "bisection_certificate": bisect,
        "newton_replay": newton,
        "derivative_audit": audit,
        "uniqueness_proof_reduction": [
            "P2384 proves the cap derivative of the d5 chamber margin is positive on the audited cap band.",
            "The sign-changing bracket contains exactly one zero because F is strictly increasing on that band.",
            "M=1.6 is accepted because F(1.6)>0, while values below the unique root fail the d5 chamber inequality at the worst corner.",
        ],
        "accepted_cap_replay": {
            "M_test": CAP_TEST,
            "F_M_test": cap_test_margin,
            "M_test_minus_root_mid": CAP_TEST - root,
            "accepted": cap_test_margin > 0.0,
        },
        "source_burden_translation": {
            "minimal_cap_must_exceed_unique_root": root,
            "M_1_6_early_interval_length": 1.0 / CAP_TEST,
            "M_1_6_early_half_mass": CAP_TEST * 0.5,
            "M_1_6_barycenter": 1.0 / (2.0 * CAP_TEST),
            "M_1_6_barycenter_shift_from_uniform": 0.5 - 1.0 / (2.0 * CAP_TEST),
            "status": "UNIQUE_THRESHOLD_ACCEPTANCE_TARGET_ONLY_SOURCE_STILL_OPEN",
        },
    }


def main() -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    docs = {name: path.read_text(encoding="utf-8") if path.exists() else "" for name, path in DOC_FILES.items()}

    certificate = cap_threshold_root_certificate()
    p2382_replay = p2382_threshold_replay(artifacts["P2382_BOUNDED_DENSITY_BATHTUB"])

    theorem_export = {
        "name": "P2388/S1338 cap-threshold root uniqueness certificate",
        "claim": (
            "The worst-corner bounded-density cap threshold is the unique root of a scalar monotone equation F(M)=0, not merely a bisection number. A sign-changing bracket around the P2382 value plus the P2384 positive cap-derivative proof yields uniqueness; bisection and Newton are retained as reproducible computations."
        ),
        "positive_content": [
            "Defines the exact threshold equation using the P2382 strict-kernel target 3*K_strict(1)-K_strict(5).",
            "Exports a sign-changing bracket, high-iteration bisection certificate, and Newton replay for the unique root near 1.574821357435363.",
            "Replays M=1.6 as strictly above the unique threshold and keeps the source burden explicit.",
            "Connects the P2386/P2387 LP/KKT acceptance target back to the scalar cap threshold without promoting a source theorem.",
        ],
        "not_licensed": [
            "strict source theorem deriving M or the front-loaded density",
            "promotion of the scalar threshold equation into L_total",
            "claim that strict dynamics supplies the cap threshold",
            "QW-2191 discharge or selector closure",
            "beta_tors -> chi11 theorem",
            "legacy physical-role transfer",
            "ToE closure",
        ],
    }

    probe = {
        "probe_id": "P2388_S1338_CAP_THRESHOLD_ROOT_UNIQUENESS_CERTIFICATE",
        "source_artifacts": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "docs": {name: rel(path) for name, path in DOC_FILES.items()},
        "artifact_replay": {
            "p2387_packet_id": artifacts["P2387_EXACT_KKT_BRANCH"].get("packet_id"),
            "p2386_packet_id": artifacts["P2386_LP_DUAL"].get("packet_id"),
            "p2384_packet_id": artifacts["P2384_SYMBOLIC_INEQUALITY"].get("packet_id"),
            "p2382_packet_id": artifacts["P2382_BOUNDED_DENSITY_BATHTUB"].get("packet_id"),
        },
        "p2382_replay": p2382_replay,
        "cap_threshold_root_certificate": certificate,
        "theorem_export": theorem_export,
    }
    probe["theorem_fingerprint_sha256"] = sha256_json(theorem_export)

    root_mid = certificate["bisection_certificate"]["mid"]
    gatekeeper_checks = {
        "source_artifacts_present": all("_missing" not in artifact for artifact in artifacts.values()),
        "p2387_loaded": probe["artifact_replay"]["p2387_packet_id"] == "P2387",
        "p2386_loaded": probe["artifact_replay"]["p2386_packet_id"] == "P2386",
        "p2384_loaded": probe["artifact_replay"]["p2384_packet_id"] == "P2384",
        "p2382_loaded": probe["artifact_replay"]["p2382_packet_id"] == "P2382",
        "initial_bracket_sign_changes": certificate["initial_sign_bracket_around_p2382_value"]["sign_change"],
        "bisection_width_collapsed": certificate["bisection_certificate"]["width"] <= 1.0e-15,
        "newton_residual_tiny": certificate["newton_replay"]["final_abs_residual"] < 1.0e-14,
        "derivative_audit_positive": certificate["derivative_audit"]["all_sampled_positive"],
        "root_matches_p2382_threshold": abs(root_mid - p2382_replay.get("p2382_cap_threshold_gt", 0.0)) < 1.0e-12,
        "M_1_6_above_unique_root": certificate["accepted_cap_replay"]["accepted"]
        and certificate["accepted_cap_replay"]["M_test_minus_root_mid"] > 0.0,
        "source_target_not_promoted": certificate["source_burden_translation"]["status"] == "UNIQUE_THRESHOLD_ACCEPTANCE_TARGET_ONLY_SOURCE_STILL_OPEN",
        "docs_updated_with_p2388_threshold_root": all("P2388/S1338" in text for text in docs.values()),
        "no_qw2191_discharge_claimed": True,
    }

    payload = {
        "schema_version": "p2388_s1338_v1",
        "packet_id": "P2388",
        "stage_id": "S1338",
        "timestamp_utc": "2026-06-01T00:00:00Z",
        "produced_by": rel(Path(__file__).resolve()),
        "status": "OPEN_PROGRESS_UNIQUE_CAP_THRESHOLD_ROOT_CERTIFICATE_SOURCE_STILL_OPEN",
        "result_kind": "CAP_THRESHOLD_ROOT_UNIQUENESS_CERTIFICATE",
        "cap_threshold_root_uniqueness_certificate": probe,
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": (
            "Use the unique scalar root as the exact cap-threshold acceptance target: future source work must derive M above this root, or keep the bounded-density cap as an explicit non-strict premise."
        ),
        "global_status": "OPEN_PROGRESS_WITH_TRACE_NO_SELECTOR_OR_TOE_CLOSURE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2388 S1338: cap-threshold root uniqueness certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2388/S1338 turns the worst-corner cap threshold into a unique-root certificate for the scalar equation",
                "",
                "```text",
                "F(M)=W_M(5)-3*W_M(1)-(3*K_strict(1)-K_strict(5))=0.",
                "```",
                "",
                "The sign-changing bracket plus the P2384 positive cap-derivative proof gives uniqueness; bisection and Newton are retained as reproducible computations rather than as the theorem core.",
                "",
                "## Checks",
                "",
                f"- Bracket: `[{certificate['initial_sign_bracket_around_p2382_value']['lo']}, {certificate['initial_sign_bracket_around_p2382_value']['hi']}]`.",
                f"- Root midpoint: `{root_mid}`.",
                f"- Bisection width: `{certificate['bisection_certificate']['width']}`.",
                f"- Newton final residual: `{certificate['newton_replay']['final_abs_residual']}`.",
                f"- `F(1.6)`: `{certificate['accepted_cap_replay']['F_M_test']}`.",
                "",
                "## Hard limits",
                "",
                "- This is a unique threshold acceptance certificate, not a strict source theorem deriving `M` or the density.",
                "- No `L_total` promotion, `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, or ToE closure is claimed.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
