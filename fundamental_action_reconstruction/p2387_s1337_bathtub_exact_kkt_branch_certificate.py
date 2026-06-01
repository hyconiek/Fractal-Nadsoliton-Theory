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

OUT = GEN / "p2387_s1337_bathtub_exact_kkt_branch_certificate.json"
MD = GEN / "p2387_s1337_bathtub_exact_kkt_branch_certificate.md"

SOURCE_FILES = {
    "P2386_LP_DUAL": GEN / "p2386_s1336_bathtub_lp_dual_certificate.json",
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
BRANCH_AUDIT_GRID = 32


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


def u0(d: int) -> float:
    return 1.0 + BETA_TORS_WORST * d


def delta(d: int) -> float:
    return d**ETA_WORST - BETA_TORS_WORST * d


def transport_speed(d: int, s: float) -> float:
    return delta(d) / (u0(d) + s * delta(d))


def q_contrast(s: float) -> float:
    return transport_speed(D5, s) - 3.0 * transport_speed(D1, s)


def q_prime(s: float) -> float:
    d1 = delta(D1)
    d5 = delta(D5)
    u1s = u0(D1) + s * d1
    u5s = u0(D5) + s * d5
    return 3.0 * d1 * d1 / (u1s * u1s) - d5 * d5 / (u5s * u5s)


def q_primitive(lo: float, hi: float) -> float:
    return math.log((u0(D5) + hi * delta(D5)) / (u0(D5) + lo * delta(D5))) - 3.0 * math.log(
        (u0(D1) + hi * delta(D1)) / (u0(D1) + lo * delta(D1))
    )


def exact_branch_certificate(cap: float = CAP_M) -> dict[str, Any]:
    t = 1.0 / cap
    lam = q_contrast(t)
    q0 = q_contrast(0.0)
    q1 = q_contrast(1.0)
    qt = q_contrast(t)
    left_margin = q0 - lam
    right_margin = lam - q1
    primal_value = cap * q_primitive(0.0, t)
    dual_value = lam + cap * (q_primitive(0.0, t) - t * lam)

    branch_cells = [
        {
            "branch": "left_s_less_than_t",
            "domain": "0 <= s < t",
            "monotone_consequence": "q(s) > lambda",
            "rho_star": "M",
            "mu": "q(s)-lambda",
            "dual_gap_lambda_plus_mu_minus_q": "0",
            "cap_slack_M_minus_rho": "0",
            "kkt_products": {"mu_times_cap_slack": "0", "dual_gap_times_rho": "0"},
        },
        {
            "branch": "cut_point_measure_zero",
            "domain": "s = t",
            "monotone_consequence": "q(s) = lambda",
            "rho_star": "arbitrary endpoint convention; measure zero",
            "mu": "0",
            "dual_gap_lambda_plus_mu_minus_q": "0",
            "cap_slack_M_minus_rho": "endpoint convention",
            "kkt_products": {"mu_times_cap_slack": "0", "dual_gap_times_rho": "0"},
        },
        {
            "branch": "right_s_greater_than_t",
            "domain": "t < s <= 1",
            "monotone_consequence": "q(s) < lambda",
            "rho_star": "0",
            "mu": "0",
            "dual_gap_lambda_plus_mu_minus_q": "lambda-q(s)",
            "cap_slack_M_minus_rho": "M",
            "kkt_products": {"mu_times_cap_slack": "0", "dual_gap_times_rho": "0"},
        },
    ]

    branch_audit_rows = []
    min_left_q_minus_lambda = float("inf")
    min_right_lambda_minus_q = float("inf")
    max_q_prime = -float("inf")
    for index in range(BRANCH_AUDIT_GRID + 1):
        s = index / BRANCH_AUDIT_GRID
        qp = q_prime(s)
        max_q_prime = max(max_q_prime, qp)
        if s < t:
            min_left_q_minus_lambda = min(min_left_q_minus_lambda, q_contrast(s) - lam)
        if s > t:
            min_right_lambda_minus_q = min(min_right_lambda_minus_q, lam - q_contrast(s))
        if index in {0, int(t * BRANCH_AUDIT_GRID), int(t * BRANCH_AUDIT_GRID) + 1, BRANCH_AUDIT_GRID}:
            branch_audit_rows.append(
                {
                    "s": s,
                    "q_s": q_contrast(s),
                    "q_prime_s": qp,
                    "q_minus_lambda": q_contrast(s) - lam,
                    "lambda_minus_q": lam - q_contrast(s),
                    "branch": "left" if s < t else "right_or_cut",
                }
            )

    return {
        "cap_M": cap,
        "cut_t_1_over_M": t,
        "lambda_q_t": lam,
        "endpoint_values": {"q_0": q0, "q_t": qt, "q_1": q1},
        "strict_endpoint_margins_from_monotone_q": {
            "q_0_minus_lambda": left_margin,
            "lambda_minus_q_1": right_margin,
            "both_positive": left_margin > 0.0 and right_margin > 0.0,
        },
        "closed_branch_kkt_certificate": branch_cells,
        "proof_reduction": [
            "P2384/P2386 provide q'(s)<0 for the audited rectangle and cap target.",
            "Therefore q(s)>q(t)=lambda on 0<=s<t and q(s)<q(t)=lambda on t<s<=1.",
            "The branch definitions mu=max(q-lambda,0) and rho*=M*1_[0,t) make dual feasibility and all KKT products identities on each branch.",
            "The primal-dual values match by the closed primitive integral, so no grid sample is part of the theorem statement.",
        ],
        "closed_value_identity": {
            "primal_value_M_int_0_t_q": primal_value,
            "dual_value_lambda_plus_M_int_mu": dual_value,
            "absolute_gap": abs(primal_value - dual_value),
        },
        "computational_audit_not_theorem_core": {
            "grid_points": BRANCH_AUDIT_GRID + 1,
            "max_q_prime_on_grid": max_q_prime,
            "min_left_q_minus_lambda_on_grid_excluding_cut": min_left_q_minus_lambda,
            "min_right_lambda_minus_q_on_grid_excluding_cut": min_right_lambda_minus_q,
            "rows": branch_audit_rows,
        },
        "source_burden_acceptance_target": {
            "cap_M": cap,
            "early_interval_length": t,
            "required_branch_saturation": "rho=M on [0,1/M), rho=0 on (1/M,1] up to null-set conventions",
            "status": "EXACT_KKT_ACCEPTANCE_TARGET_ONLY_SOURCE_STILL_OPEN",
        },
    }


def extract_p2386_gap(artifact: dict[str, Any]) -> dict[str, Any]:
    try:
        cert = artifact["bathtub_lp_dual_certificate"]["lp_dual_certificate"]
        return {
            "packet_value_gap": cert["absolute_primal_dual_gap"],
            "packet_weight_gap": cert["absolute_primal_weight_gap"],
            "packet_cap_M": cert["cap_M"],
        }
    except KeyError:
        return {"status": "P2386_REPLAY_FIELDS_MISSING"}


def main() -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    docs = {name: path.read_text(encoding="utf-8") if path.exists() else "" for name, path in DOC_FILES.items()}

    certificate = exact_branch_certificate(CAP_M)
    p2386_replay = extract_p2386_gap(artifacts["P2386_LP_DUAL"])

    theorem_export = {
        "name": "P2387/S1337 bathtub exact KKT branch certificate",
        "claim": (
            "The P2386 LP-dual certificate can be upgraded from sampled KKT evidence to an exact branch certificate: once P2384/P2386 provide q'(s)<0, the cut t=1/M partitions the domain into q>lambda, q=lambda, and q<lambda branches, making dual feasibility and complementarity algebraic identities."
        ),
        "positive_content": [
            "Exports the three KKT branches for s<t, s=t, and s>t.",
            "Checks strict endpoint margins q(0)-lambda and lambda-q(1), then uses monotonicity to certify branch signs without sampling as theorem core.",
            "Replays the P2386 primal-dual value gap and recomputes the closed value identity.",
            "Leaves the source question as an exact KKT acceptance target rather than a sourced density theorem.",
        ],
        "not_licensed": [
            "strict source theorem deriving rho* or cap M",
            "promotion of the LP/KKT branch certificate into L_total",
            "claim that strict nadsoliton dynamics supplies the front-loaded density",
            "QW-2191 discharge or selector closure",
            "beta_tors -> chi11 theorem",
            "legacy physical-role transfer",
            "ToE closure",
        ],
    }

    probe = {
        "probe_id": "P2387_S1337_BATHTUB_EXACT_KKT_BRANCH_CERTIFICATE",
        "source_artifacts": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "docs": {name: rel(path) for name, path in DOC_FILES.items()},
        "artifact_replay": {
            "p2386_packet_id": artifacts["P2386_LP_DUAL"].get("packet_id"),
            "p2385_packet_id": artifacts["P2385_EXACT_SUPPORT_CHAMBER"].get("packet_id"),
            "p2384_packet_id": artifacts["P2384_SYMBOLIC_INEQUALITY"].get("packet_id"),
            "p2382_packet_id": artifacts["P2382_BOUNDED_DENSITY_BATHTUB"].get("packet_id"),
        },
        "p2386_replay": p2386_replay,
        "exact_kkt_branch_certificate": certificate,
        "theorem_export": theorem_export,
    }
    probe["theorem_fingerprint_sha256"] = sha256_json(theorem_export)

    gatekeeper_checks = {
        "source_artifacts_present": all("_missing" not in artifact for artifact in artifacts.values()),
        "p2386_loaded": probe["artifact_replay"]["p2386_packet_id"] == "P2386",
        "p2385_loaded": probe["artifact_replay"]["p2385_packet_id"] == "P2385",
        "p2384_loaded": probe["artifact_replay"]["p2384_packet_id"] == "P2384",
        "p2382_loaded": probe["artifact_replay"]["p2382_packet_id"] == "P2382",
        "strict_endpoint_branch_margins_positive": certificate["strict_endpoint_margins_from_monotone_q"]["both_positive"],
        "closed_value_identity_gap_zero": certificate["closed_value_identity"]["absolute_gap"] < 1.0e-12,
        "p2386_replay_gap_zero": p2386_replay.get("packet_value_gap", 1.0) < 1.0e-12,
        "audit_grid_qprime_negative": certificate["computational_audit_not_theorem_core"]["max_q_prime_on_grid"] < 0.0,
        "audit_grid_branch_margins_positive": certificate["computational_audit_not_theorem_core"]["min_left_q_minus_lambda_on_grid_excluding_cut"] > 0.0
        and certificate["computational_audit_not_theorem_core"]["min_right_lambda_minus_q_on_grid_excluding_cut"] > 0.0,
        "source_target_not_promoted": certificate["source_burden_acceptance_target"]["status"] == "EXACT_KKT_ACCEPTANCE_TARGET_ONLY_SOURCE_STILL_OPEN",
        "docs_updated_with_p2387_exact_kkt": all("P2387/S1337" in text for text in docs.values()),
        "no_qw2191_discharge_claimed": True,
    }

    payload = {
        "schema_version": "p2387_s1337_v1",
        "packet_id": "P2387",
        "stage_id": "S1337",
        "timestamp_utc": "2026-06-01T00:00:00Z",
        "produced_by": rel(Path(__file__).resolve()),
        "status": "OPEN_PROGRESS_EXACT_KKT_BRANCH_CERTIFICATE_SOURCE_STILL_OPEN",
        "result_kind": "BATHTUB_EXACT_KKT_BRANCH_CERTIFICATE",
        "bathtub_exact_kkt_branch_certificate": probe,
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": (
            "Treat the P2387 branch certificate as the exact source-side acceptance target: future work must derive the front-loaded density/cap from strict dynamics or explicitly keep it as a non-strict premise."
        ),
        "global_status": "OPEN_PROGRESS_WITH_TRACE_NO_SELECTOR_OR_TOE_CLOSURE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2387 S1337: bathtub exact KKT branch certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2387/S1337 upgrades the P2386 LP-dual audit from sampled KKT evidence to a branchwise KKT certificate.",
                "Given the already-certified monotonicity `q'(s)<0`, the cut `t=1/M=0.625` gives `q(s)>lambda` for `s<t`, `q(s)=lambda` at the cut, and `q(s)<lambda` for `s>t`.",
                "Therefore `mu=max(q-lambda,0)` and `rho*=M*1_[0,t)` make dual feasibility and complementarity algebraic branch identities.",
                "",
                "## Checks",
                "",
                f"- `q(0)-lambda`: `{certificate['strict_endpoint_margins_from_monotone_q']['q_0_minus_lambda']}`.",
                f"- `lambda-q(1)`: `{certificate['strict_endpoint_margins_from_monotone_q']['lambda_minus_q_1']}`.",
                f"- Closed value gap: `{certificate['closed_value_identity']['absolute_gap']}`.",
                f"- Audit-grid max `q_prime`: `{certificate['computational_audit_not_theorem_core']['max_q_prime_on_grid']}`.",
                "",
                "## Hard limits",
                "",
                "- This is an exact KKT acceptance certificate, not a strict source theorem deriving `rho*` or `M`.",
                "- No `L_total` promotion, `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, or ToE closure is claimed.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
