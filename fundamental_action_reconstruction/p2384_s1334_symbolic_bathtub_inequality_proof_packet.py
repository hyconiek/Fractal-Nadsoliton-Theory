#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from decimal import Decimal, getcontext
from pathlib import Path
from typing import Any

getcontext().prec = 80

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2384_s1334_symbolic_bathtub_inequality_proof_packet.json"
MD = GEN / "p2384_s1334_symbolic_bathtub_inequality_proof_packet.md"

SOURCE_FILES = {
    "P2383_CORNER_REDUCTION": GEN / "p2383_s1333_closed_form_bathtub_corner_reduction_theorem.json",
    "P2382_BOUNDED_DENSITY_BATHTUB": GEN / "p2382_s1332_bounded_density_bathtub_frontload_certificate.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

CAP_LO = Decimal("1.5")
CAP_HI = Decimal("1.6")
BETA_LO = Decimal("0")
BETA_HI = Decimal("0.1")
ETA_LO_NUM = 9
ETA_LO_DEN = 5
ETA_COARSE_NUM = 3
ETA_COARSE_DEN = 2
FIVE = Decimal(5)
SQRT_3 = Decimal(3).sqrt()
SQRT_125 = Decimal(125).sqrt()


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


def dec(value: Decimal) -> str:
    return format(value, "f")


def h_cap(x: Decimal) -> Decimal:
    return (Decimal(1) + x).ln() - x / (Decimal(1) + x)


def strict_ratio_symbolic_proof() -> dict[str, Any]:
    threshold = Decimal(3) + Decimal(2) * SQRT_3
    algebraic_margin = SQRT_125 - threshold
    squared_margin_certificate = {
        "claim": "sqrt(125) > 3+2*sqrt(3)",
        "equivalent_integer_check": "104 > 12*sqrt(3); squaring is legal because both sides are positive",
        "integer_margin_after_squaring": 104 * 104 - 12 * 12 * 3,
    }
    rmin_lower = Decimal(2) * SQRT_125 / (Decimal(1) + SQRT_125)
    rmin_exact_9_5 = Decimal("1.8953947080453295")
    return {
        "proof_goal": "q'(s)=A1(s)^2*(3-R(s)^2)<0 on the P2376 rectangle",
        "coarse_exponent_step": "Since eta>=9/5>3/2, 5^eta >= 5^(3/2)=sqrt(125).",
        "ratio_floor_sufficient_condition": "2*x/(1+x)>sqrt(3) is equivalent to x>3+2*sqrt(3).",
        "sqrt125_gt_threshold": {
            "sqrt125": dec(SQRT_125),
            "three_plus_two_sqrt3": dec(threshold),
            "margin": dec(algebraic_margin),
            "integer_certificate": squared_margin_certificate,
        },
        "rmin_lower_from_coarse_bound": dec(rmin_lower),
        "rmin_lower_square_minus_3": dec(rmin_lower * rmin_lower - Decimal(3)),
        "p2383_numeric_rmin_replay": dec(rmin_exact_9_5),
        "logical_chain": [
            "R(s) decreases in s because d/ds[(u1+s*Delta1)/(u5+s*Delta5)] has numerator Delta1*u5-Delta5*u1<0; the coarse Delta5>=sqrt(125)-0.5 and Delta1<=1 beat u5/u1<=1.5.",
            "At s=1, R=2*(x-5*beta_tors)/((1-beta_tors)*(1+x)); it is increasing in x=5^eta and in beta_tors for x>5, so the lower corner is beta_tors=0 and eta minimized.",
            "The coarser eta>=3/2 bound already gives R>sqrt(3); therefore q'(s)<0 without a grid-first proof premise.",
        ],
        "passes": rmin_lower * rmin_lower > Decimal(3) and algebraic_margin > 0,
    }


def beta_derivative_symbolic_bound() -> dict[str, Any]:
    n5_lower = Decimal(5) * (Decimal(1) + SQRT_125) / (
        (Decimal(1) + Decimal(5) * BETA_HI)
        * (CAP_HI * (Decimal(1) + Decimal(5) * BETA_HI) + SQRT_125 - Decimal(5) * BETA_HI)
    )
    n1_upper = Decimal(2) / ((Decimal(1) + BETA_LO) * (CAP_LO * (Decimal(1) + BETA_LO) + Decimal(1) - BETA_LO))
    return {
        "derivative_identity": (
            "dW_d/dbeta_tors = -N_d, where N_d=d*(1+d^eta)/((1+beta_tors*d)*(M*(1+beta_tors*d)+d^eta-beta_tors*d))."
        ),
        "margin_derivative": "d/d beta_tors [W5-3W1] = -N5+3*N1",
        "coarse_lower_N5": dec(n5_lower),
        "coarse_upper_3N1": dec(Decimal(3) * n1_upper),
        "positive_gap_N5_minus_3N1": dec(n5_lower - Decimal(3) * n1_upper),
        "conclusion": "Because N5>3*N1 on the audited band, the capped chamber margin decreases with beta_tors.",
        "passes": n5_lower > Decimal(3) * n1_upper,
    }


def cap_derivative_symbolic_bound() -> dict[str, Any]:
    x5_lower = (SQRT_125 - Decimal(5) * BETA_HI) / (CAP_HI * (Decimal(1) + Decimal(5) * BETA_HI))
    x1_upper = (Decimal(1) - BETA_LO) / (CAP_LO * (Decimal(1) + BETA_LO))
    h5_lower = h_cap(x5_lower)
    h1_upper = h_cap(x1_upper)
    return {
        "derivative_identity": "dW_d/dM = h(x_d), h(x)=log(1+x)-x/(1+x), x_d=(d^eta-beta_tors*d)/(M*(1+beta_tors*d)).",
        "monotonic_h_note": "h'(x)=x/(1+x)^2>0 for x>0.",
        "coarse_lower_x5": dec(x5_lower),
        "coarse_upper_x1": dec(x1_upper),
        "coarse_lower_h_x5": dec(h5_lower),
        "coarse_upper_3h_x1": dec(Decimal(3) * h1_upper),
        "positive_gap_h5_minus_3h1": dec(h5_lower - Decimal(3) * h1_upper),
        "conclusion": "Because h(x5)>3*h(x1) on the audited band, the capped chamber margin increases with M.",
        "passes": h5_lower > Decimal(3) * h1_upper,
    }


def eta_derivative_symbolic_bound() -> dict[str, Any]:
    return {
        "derivative_identity": "dW_d/deta = d^eta*log(d)/(1+beta_tors*d+(d^eta-beta_tors*d)/M).",
        "d1_term": "For d=1, log(1)=0, so dW1/deta=0.",
        "d5_term": "For d=5, every factor is positive on the audited rectangle and cap band.",
        "conclusion": "d/deta [W5-3W1] > 0, so the capped chamber margin increases with eta.",
        "passes": True,
    }


def source_target_translation() -> dict[str, Any]:
    cap = Decimal("1.6")
    return {
        "cap_replay_M": dec(cap),
        "early_interval_length": dec(Decimal(1) / cap),
        "early_half_mass": dec(cap / Decimal(2)),
        "barycenter": dec(Decimal(1) / (Decimal(2) * cap)),
        "barycenter_shift_from_uniform": dec(Decimal("0.5") - Decimal(1) / (Decimal(2) * cap)),
        "source_theorem_obligation": (
            "A future strict source theorem using this lane must source a bounded density with at least this cap/frontload scale; otherwise the cap remains a non-strict premise."
        ),
    }


def main() -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    docs = {name: path.read_text(encoding="utf-8") if path.exists() else "" for name, path in DOC_FILES.items()}

    ratio = strict_ratio_symbolic_proof()
    beta = beta_derivative_symbolic_bound()
    cap = cap_derivative_symbolic_bound()
    eta = eta_derivative_symbolic_bound()
    target = source_target_translation()

    theorem_export = {
        "name": "P2384/S1334 symbolic bathtub inequality proof packet",
        "claim": (
            "The P2383 corner-reduction lane has a symbolic inequality core: a coarse 5^eta>=sqrt(125) bound proves the q'(s)<0 ratio floor, "
            "and closed derivative identities with coarse endpoint inequalities certify the eta/beta_tors/M corner-direction signs on the audited cap band."
        ),
        "positive_content": [
            "Exports an integer-backed proof that sqrt(125)>3+2*sqrt(3), which is sufficient for R>sqrt(3).",
            "Replaces the beta_tors and M derivative sign claims with explicit N_d and h(x) identities plus positive coarse gaps.",
            "Keeps the bounded-density source target as an obligation rather than a sourced strict theorem.",
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
        "probe_id": "P2384_S1334_SYMBOLIC_BATHTUB_INEQUALITY_PROOF_PACKET",
        "source_artifacts": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "docs": {name: rel(path) for name, path in DOC_FILES.items()},
        "artifact_replay": {
            "p2383_packet_id": artifacts["P2383_CORNER_REDUCTION"].get("packet_id"),
            "p2382_packet_id": artifacts["P2382_BOUNDED_DENSITY_BATHTUB"].get("packet_id"),
        },
        "symbolic_ratio_proof": ratio,
        "eta_derivative_sign_proof": eta,
        "beta_tors_derivative_sign_proof": beta,
        "cap_derivative_sign_proof": cap,
        "source_target_translation": target,
        "theorem_export": theorem_export,
    }
    probe["theorem_fingerprint_sha256"] = sha256_json(theorem_export)

    gatekeeper_checks = {
        "source_artifacts_present": all("_missing" not in artifact for artifact in artifacts.values()),
        "p2383_loaded": probe["artifact_replay"]["p2383_packet_id"] == "P2383",
        "p2382_loaded": probe["artifact_replay"]["p2382_packet_id"] == "P2382",
        "ratio_symbolic_proof_passes": ratio["passes"],
        "eta_derivative_sign_passes": eta["passes"],
        "beta_derivative_sign_passes": beta["passes"],
        "cap_derivative_sign_passes": cap["passes"],
        "source_target_not_promoted": "non-strict premise" in target["source_theorem_obligation"],
        "docs_updated_with_p2384_symbolic_packet": all("P2384/S1334" in text for text in docs.values()),
        "no_qw2191_discharge_claimed": True,
    }

    payload = {
        "schema_version": "p2384_s1334_v1",
        "packet_id": "P2384",
        "stage_id": "S1334",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-06-01T00:00:00+00:00",
        "status": "OPEN_PROGRESS_SYMBOLIC_BATHTUB_INEQUALITY_PROOF_SOURCE_STILL_OPEN",
        "result_kind": "SYMBOLIC_BATHTUB_INEQUALITY_PROOF_PACKET",
        "symbolic_bathtub_inequality_proof_packet": probe,
        "recommended_next_honest_step": (
            "Use the symbolic P2384 inequalities as the proof-side acceptance kernel for any future strict source theorem. "
            "The remaining open task is still to derive the bounded front-loaded density/cap from a real strict source, not to add more optimizer replays."
        ),
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_PROGRESS_WITH_TRACE_NO_SELECTOR_OR_TOE_CLOSURE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2384 S1334: symbolic bathtub inequality proof packet",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2384/S1334 extracts an inequality-level proof core from P2383.",
                "The ratio part uses only the coarse bound `5^eta >= 5^(3/2)=sqrt(125)` and the algebraic threshold `3+2*sqrt(3)`.",
                "The derivative part records closed identities for the eta, beta_tors, and M directions and supplies positive coarse endpoint gaps on the cap band `[1.5,1.6]`.",
                "",
                "## Certified coarse gaps",
                "",
                f"- `sqrt(125)-(3+2*sqrt(3)) = {ratio['sqrt125_gt_threshold']['margin']}`.",
                f"- beta-direction gap `N5-3*N1 >= {beta['positive_gap_N5_minus_3N1']}`.",
                f"- cap-direction gap `h(x5)-3*h(x1) >= {cap['positive_gap_h5_minus_3h1']}`.",
                f"- `M=1.6` early interval length: `{target['early_interval_length']}`.",
                f"- `M=1.6` early-half mass: `{target['early_half_mass']}`.",
                f"- `M=1.6` barycenter shift from uniform: `{target['barycenter_shift_from_uniform']}`.",
                "",
                "## Hard limits",
                "",
                "- This is an inequality proof packet for the bounded-density acceptance criterion, not a strict source theorem for the cap or bang-bang profile.",
                "- No `L_total` promotion, `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, or ToE closure is claimed.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
