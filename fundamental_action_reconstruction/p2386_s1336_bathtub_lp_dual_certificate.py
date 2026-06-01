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

OUT = GEN / "p2386_s1336_bathtub_lp_dual_certificate.json"
MD = GEN / "p2386_s1336_bathtub_lp_dual_certificate.md"

SOURCE_FILES = {
    "P2385_EXACT_SUPPORT_CHAMBER": GEN / "p2385_s1335_exact_z12_support_chamber_theorem.json",
    "P2384_SYMBOLIC_INEQUALITY": GEN / "p2384_s1334_symbolic_bathtub_inequality_proof_packet.json",
    "P2382_BOUNDED_DENSITY_BATHTUB": GEN / "p2382_s1332_bounded_density_bathtub_frontload_certificate.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

STRICT_BETA = 1.0
OMEGA = math.pi / 4.0
PHI = math.pi / 6.0
D1 = 1
D5 = 5
CAP_M = 8.0 / 5.0
ETA_WORST = 9.0 / 5.0
BETA_TORS_WORST = 1.0 / 10.0
SAMPLE_GRID = 64
COMPLEMENTARITY_GRID = 41


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
    return math.cos(OMEGA * d + PHI) / (1.0 + STRICT_BETA * d**ETA_WORST)


def u0(d: int, beta_tors: float = BETA_TORS_WORST) -> float:
    return 1.0 + beta_tors * d


def delta(d: int, eta: float = ETA_WORST, beta_tors: float = BETA_TORS_WORST) -> float:
    return d**eta - beta_tors * d


def transport_speed(d: int, s: float, eta: float = ETA_WORST, beta_tors: float = BETA_TORS_WORST) -> float:
    return delta(d, eta, beta_tors) / (u0(d, beta_tors) + s * delta(d, eta, beta_tors))


def contrast_q(s: float) -> float:
    return transport_speed(D5, s) - 3.0 * transport_speed(D1, s)


def contrast_q_prime(s: float) -> float:
    d1 = delta(D1)
    d5 = delta(D5)
    u1 = u0(D1) + s * d1
    u5 = u0(D5) + s * d5
    return 3.0 * d1 * d1 / (u1 * u1) - d5 * d5 / (u5 * u5)


def capped_bangbang_weight(d: int, cap: float = CAP_M) -> float:
    t = 1.0 / cap
    start = u0(d)
    dlt = delta(d)
    return cap * math.log((start + t * dlt) / start)


def primitive_q_integral(lo: float, hi: float) -> float:
    return math.log((u0(D5) + hi * delta(D5)) / (u0(D5) + lo * delta(D5))) - 3.0 * math.log(
        (u0(D1) + hi * delta(D1)) / (u0(D1) + lo * delta(D1))
    )


def p2382_replay_summary(artifact: dict[str, Any]) -> dict[str, Any]:
    try:
        cert = artifact["bounded_density_bathtub_frontload_theorem"]["bounded_density_bathtub_frontload_certificate"]
        obligation = cert["rectangle_worst_cap_source_obligation"]
        replay = obligation["cap_test_positive_replay"]
        return {
            "cap_M": replay["density_cap_M"],
            "eta": replay["eta"],
            "beta_tors": replay["beta_tors"],
            "W5_minus_3W1": replay["W5"] - 3.0 * replay["W1"],
            "d5_margin_b_minus_3a": replay["d5_margin_b_minus_3a"],
            "a_over_b": replay["a_over_b"],
            "d5_chamber": replay["d5_chamber"],
        }
    except KeyError:
        return {"status": "P2382_REPLAY_FIELDS_MISSING"}


def lp_dual_certificate(cap: float = CAP_M) -> dict[str, Any]:
    t = 1.0 / cap
    lam = contrast_q(t)
    primal_value = cap * primitive_q_integral(0.0, t)
    dual_slack_integral = primitive_q_integral(0.0, t) - t * lam
    dual_value = lam + cap * dual_slack_integral
    w_contrast = capped_bangbang_weight(D5, cap) - 3.0 * capped_bangbang_weight(D1, cap)

    samples = []
    max_dual_gap_violation = 0.0
    max_monotonicity_prime = -1.0e100
    for index in range(SAMPLE_GRID + 1):
        s = index / SAMPLE_GRID
        q = contrast_q(s)
        mu = max(q - lam, 0.0)
        dual_gap = lam + mu - q
        qprime = contrast_q_prime(s)
        max_dual_gap_violation = max(max_dual_gap_violation, -dual_gap)
        max_monotonicity_prime = max(max_monotonicity_prime, qprime)
        if index in {0, int(round(t * SAMPLE_GRID)), SAMPLE_GRID}:
            samples.append(
                {
                    "s": s,
                    "q_s": q,
                    "lambda": lam,
                    "mu_s": mu,
                    "dual_gap_lambda_plus_mu_minus_q": dual_gap,
                    "q_prime_s": qprime,
                }
            )

    complementarity_rows = []
    max_complementarity_error = 0.0
    for index in range(COMPLEMENTARITY_GRID + 1):
        s = index / COMPLEMENTARITY_GRID
        rho = cap if s < t else 0.0
        q = contrast_q(s)
        mu = max(q - lam, 0.0)
        dual_gap = lam + mu - q
        density_gap = cap - rho
        comp_mu_density = mu * density_gap
        comp_dual_density = dual_gap * rho
        max_complementarity_error = max(max_complementarity_error, abs(comp_mu_density), abs(comp_dual_density))
        if index in {0, int(t * COMPLEMENTARITY_GRID), int(t * COMPLEMENTARITY_GRID) + 1, COMPLEMENTARITY_GRID}:
            complementarity_rows.append(
                {
                    "s": s,
                    "rho_star_s": rho,
                    "mu_s": mu,
                    "dual_gap": dual_gap,
                    "mu_times_cap_minus_rho": comp_mu_density,
                    "dual_gap_times_rho": comp_dual_density,
                }
            )

    return {
        "linear_program": {
            "primal": "maximize integral_0^1 q(s)*rho(s) ds subject to 0<=rho<=M and integral rho=1",
            "dual": "minimize lambda + M*integral_0^1 mu(s) ds subject to lambda+mu(s)>=q(s), mu(s)>=0",
            "q_s": "A5(s)-3*A1(s)",
        },
        "cap_M": cap,
        "early_cut_t_1_over_M": t,
        "lambda_q_t": lam,
        "dual_mu_rule": "mu(s)=max(q(s)-q(1/M),0)",
        "rho_star_rule": "rho*(s)=M on 0<=s<1/M and rho*(s)=0 on 1/M<s<=1",
        "primal_value": primal_value,
        "dual_slack_integral": dual_slack_integral,
        "dual_value": dual_value,
        "closed_form_W5_minus_3W1": w_contrast,
        "absolute_primal_dual_gap": abs(primal_value - dual_value),
        "absolute_primal_weight_gap": abs(primal_value - w_contrast),
        "sampled_dual_feasibility": {
            "grid_points": SAMPLE_GRID + 1,
            "max_negative_dual_gap_violation": max_dual_gap_violation,
            "max_q_prime_on_grid": max_monotonicity_prime,
            "samples": samples,
        },
        "sampled_complementarity": {
            "grid_points": COMPLEMENTARITY_GRID + 1,
            "max_complementarity_error": max_complementarity_error,
            "rows": complementarity_rows,
        },
        "uniqueness_statement": (
            "Since P2383/P2384 certify q'(s)<0 on the audited rectangle, equality in the LP dual can only saturate the cap on the earliest interval of length 1/M, up to the endpoint of measure zero."
        ),
        "source_burden_translation_for_M_1_6": {
            "early_interval_length": t,
            "early_half_mass": min(1.0, cap * 0.5),
            "barycenter": t / 2.0,
            "shift_from_uniform_barycenter": 0.5 - t / 2.0,
            "status": "SOURCE_STILL_OPEN_NON_STRICT_PREMISE_UNTIL_DERIVED",
        },
    }


def main() -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    docs = {name: path.read_text(encoding="utf-8") if path.exists() else "" for name, path in DOC_FILES.items()}

    certificate = lp_dual_certificate(CAP_M)
    p2382_summary = p2382_replay_summary(artifacts["P2382_BOUNDED_DENSITY_BATHTUB"])

    theorem_export = {
        "name": "P2386/S1336 bathtub LP dual certificate",
        "claim": (
            "The bounded-density bathtub step can be written as an infinite-dimensional linear program with an explicit dual certificate. "
            "For M=1.6 at the P2382/P2384 worst corner, lambda=q(1/M) and mu(s)=max(q(s)-lambda,0) match the early bang-bang primal value exactly, so the cap/frontload optimization is certified by strong-duality/KKT bookkeeping rather than by an informal rearrangement argument."
        ),
        "positive_content": [
            "Exports the primal bounded-density LP and its nonnegative-slack dual.",
            "Computes the exact closed-form primal value M*int_0^(1/M) q(s) ds and matches it to W_M(5)-3*W_M(1).",
            "Checks sampled dual feasibility, sampled q'(s)<0 consistency, and KKT complementarity for the early bang-bang optimizer.",
            "Keeps the finite support theorem P2385 separate from the analytic LP certificate.",
        ],
        "not_licensed": [
            "strict variational source theorem deriving the cap M or the density rho*",
            "promotion of rho*, q, mu, or the LP objective into L_total",
            "claim that the bounded-density source is supplied by the strict core",
            "QW-2191 discharge or selector closure",
            "beta_tors -> chi11 theorem",
            "legacy physical-role transfer",
            "ToE closure",
        ],
    }

    probe = {
        "probe_id": "P2386_S1336_BATHTUB_LP_DUAL_CERTIFICATE",
        "source_artifacts": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "docs": {name: rel(path) for name, path in DOC_FILES.items()},
        "artifact_replay": {
            "p2385_packet_id": artifacts["P2385_EXACT_SUPPORT_CHAMBER"].get("packet_id"),
            "p2384_packet_id": artifacts["P2384_SYMBOLIC_INEQUALITY"].get("packet_id"),
            "p2382_packet_id": artifacts["P2382_BOUNDED_DENSITY_BATHTUB"].get("packet_id"),
        },
        "p2382_cap_replay_summary": p2382_summary,
        "lp_dual_certificate": certificate,
        "theorem_export": theorem_export,
    }
    probe["theorem_fingerprint_sha256"] = sha256_json(theorem_export)

    gatekeeper_checks = {
        "source_artifacts_present": all("_missing" not in artifact for artifact in artifacts.values()),
        "p2385_loaded": probe["artifact_replay"]["p2385_packet_id"] == "P2385",
        "p2384_loaded": probe["artifact_replay"]["p2384_packet_id"] == "P2384",
        "p2382_loaded": probe["artifact_replay"]["p2382_packet_id"] == "P2382",
        "p2382_cap_is_m_1_6": abs(p2382_summary.get("cap_M", 0.0) - CAP_M) < 1.0e-12,
        "p2382_cap_chamber_still_true": p2382_summary.get("d5_chamber") is True,
        "lp_primal_dual_gap_closed": certificate["absolute_primal_dual_gap"] < 1.0e-12,
        "lp_matches_closed_form_weight_contrast": certificate["absolute_primal_weight_gap"] < 1.0e-12,
        "sampled_dual_feasibility_ok": certificate["sampled_dual_feasibility"]["max_negative_dual_gap_violation"] < 1.0e-12,
        "sampled_monotonicity_consistent": certificate["sampled_dual_feasibility"]["max_q_prime_on_grid"] < 0.0,
        "sampled_kkt_complementarity_ok": certificate["sampled_complementarity"]["max_complementarity_error"] < 1.0e-12,
        "source_target_not_promoted": certificate["source_burden_translation_for_M_1_6"]["status"] == "SOURCE_STILL_OPEN_NON_STRICT_PREMISE_UNTIL_DERIVED",
        "docs_updated_with_p2386_lp_dual": all("P2386/S1336" in text for text in docs.values()),
        "no_qw2191_discharge_claimed": True,
    }

    payload = {
        "schema_version": "p2386_s1336_v1",
        "packet_id": "P2386",
        "stage_id": "S1336",
        "timestamp_utc": "2026-06-01T00:00:00Z",
        "produced_by": rel(Path(__file__).resolve()),
        "status": "OPEN_PROGRESS_PROOF_SIDE_LP_DUAL_CERTIFICATE_SOURCE_STILL_OPEN",
        "result_kind": "BATHTUB_LP_DUAL_CERTIFICATE",
        "bathtub_lp_dual_certificate": probe,
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": (
            "Use the LP-dual/KKT certificate as the acceptance target for any future strict source theorem: the next source-side result must derive a density whose support/cap data meet this dual-saturated early-frontload profile, or explicitly mark the cap as a non-strict premise."
        ),
        "global_status": "OPEN_PROGRESS_WITH_TRACE_NO_SELECTOR_OR_TOE_CLOSURE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2386 S1336: bathtub LP dual certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2386/S1336 rewrites the bounded-density bathtub step as a primal/dual LP certificate.",
                "The primal maximizes `int q(s)*rho(s) ds` under `0<=rho<=M` and `int rho=1`; the dual minimizes `lambda + M*int mu(s) ds` under `lambda+mu(s)>=q(s)` and `mu>=0`.",
                "For `M=1.6`, the certificate uses `lambda=q(1/M)` and `mu(s)=max(q(s)-lambda,0)`.",
                "",
                "## Numerical/closed-form checks",
                "",
                f"- Early cut `1/M`: `{certificate['early_cut_t_1_over_M']}`.",
                f"- Primal value: `{certificate['primal_value']}`.",
                f"- Dual value: `{certificate['dual_value']}`.",
                f"- Closed-form `W_M(5)-3*W_M(1)`: `{certificate['closed_form_W5_minus_3W1']}`.",
                f"- Absolute primal-dual gap: `{certificate['absolute_primal_dual_gap']}`.",
                f"- Sampled KKT complementarity max error: `{certificate['sampled_complementarity']['max_complementarity_error']}`.",
                "",
                "## Hard limits",
                "",
                "- This is a proof-side LP/KKT certificate for the bounded-density acceptance criterion, not a strict source theorem deriving the cap or density.",
                "- No `L_total` promotion, `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, or ToE closure is claimed.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
