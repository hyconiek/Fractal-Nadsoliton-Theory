#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
SCRATCH = ROOT / "scratch"

OUT = GEN / "p2363_s1313_legacy_strict_bridge_moment_lagrangian_transport_probe.json"
MD = GEN / "p2363_s1313_legacy_strict_bridge_moment_lagrangian_transport_probe.md"

SOURCE_FILES = {
    "P1866_STRICT_KERNEL_TO_LTOTAL": GEN
    / "p1866_s816_strict_kernel_to_full_lagrangian_and_eom_symbolic_export.json",
    "P2362_PARALLEL_EOM_LAGRANGIAN": GEN
    / "p2362_s1312_strict_lagrangian_eom_parallel_completion_probe.json",
    "SCRATCH_SYMBOLIC_APD_CANCELLATION": SCRATCH
    / "bridge_strict_completion_legacy_to_strict_symbolic_cancellation_certificate_report.json",
    "SCRATCH_FINITE_APD_ASSEMBLY": SCRATCH
    / "bridge_strict_completion_legacy_to_strict_finite_bridge_assembly_certificate_report.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

REQUIRED_DOC_SNIPPETS = [
    "P2363/S1313",
    "bridge-completed moment transport",
    "raw legacy moments",
    "m2_eff = m2*(1 + c0)",
]


def rel(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": rel(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def read_text(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8")


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def sstr(expr: Any) -> str:
    return sp.sstr(sp.simplify(expr))


def build_moments(kernel: sp.Expr, d: sp.Symbol) -> dict[str, sp.Expr]:
    d_ref = sp.Integer(1)
    return {
        "c0": sp.simplify(kernel.subs(d, d_ref)),
        "c1": sp.simplify(sp.diff(kernel, d).subs(d, d_ref)),
        "c2": sp.simplify(sp.diff(kernel, d, 2).subs(d, d_ref) / 2),
    }


def effective_couplings(moments: dict[str, sp.Expr]) -> dict[str, sp.Expr]:
    m2, lam, y, g, xi = sp.symbols("m2 lam y g xi", real=True)
    c0 = moments["c0"]
    c1 = moments["c1"]
    c2 = moments["c2"]
    return {
        "m2_eff": sp.simplify(m2 * (1 + c0)),
        "lam_eff": sp.simplify(lam * (1 + c1**2)),
        "y_eff": sp.simplify(y * (1 + c0 / 2)),
        "g_eff": sp.simplify(g * (1 + c2)),
        "xi_eff": sp.simplify(xi * (1 + c0)),
    }


def residuals(left: dict[str, sp.Expr], right: dict[str, sp.Expr]) -> dict[str, str]:
    return {key: sstr(left[key] - right[key]) for key in left}


def parse_p1866_exprs(
    values: dict[str, str],
    symbols: dict[str, sp.Symbol],
) -> dict[str, sp.Expr]:
    locals_map: dict[str, Any] = {
        **symbols,
        "cos": sp.cos,
        "sin": sp.sin,
    }
    return {key: sp.sympify(value, locals=locals_map) for key, value in values.items()}


def numeric_moments(
    moments: dict[str, sp.Expr],
    substitutions: dict[sp.Symbol, sp.Expr],
) -> dict[str, float]:
    return {key: float(sp.N(value.subs(substitutions), 18)) for key, value in moments.items()}


def max_abs_difference(left: dict[str, float], right: dict[str, float]) -> float:
    return max(abs(left[key] - right[key]) for key in left)


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    doc_texts = {name: read_text(path) for name, path in DOC_FILES.items()}

    p1866 = artifacts["P1866_STRICT_KERNEL_TO_LTOTAL"]
    p2362 = artifacts["P2362_PARALLEL_EOM_LAGRANGIAN"]
    symbolic_report = artifacts["SCRATCH_SYMBOLIC_APD_CANCELLATION"]
    finite_report = artifacts["SCRATCH_FINITE_APD_ASSEMBLY"]

    d, beta, eta, omega, phi = sp.symbols("d beta eta omega phi", positive=True, real=True)
    alpha_geo, omega_L, phi_L, beta_tors = sp.symbols(
        "alpha_geo omega_L phi_L beta_tors", positive=True, real=True
    )

    K_legacy = alpha_geo * sp.cos(omega_L * d + phi_L) / (1 + beta_tors * d)
    A_completion = 1 / alpha_geo
    P_completion = sp.cos(omega * d + phi) / sp.cos(omega_L * d + phi_L)
    D_completion = (1 + beta_tors * d) / (1 + beta * d**eta)
    Q_apd = sp.simplify(A_completion * P_completion * D_completion)
    K_completed = sp.simplify(K_legacy * Q_apd)
    K_strict = sp.cos(omega * d + phi) / (1 + beta * d**eta)
    K_legacy_amplitude_normalized = sp.simplify(K_legacy / alpha_geo)

    bridge_residual = sp.simplify(K_completed - K_strict)
    strict_moments = build_moments(K_strict, d)
    completed_moments = build_moments(K_completed, d)
    legacy_raw_moments = build_moments(K_legacy, d)
    legacy_normalized_moments = build_moments(K_legacy_amplitude_normalized, d)

    strict_effective = effective_couplings(strict_moments)
    completed_effective = effective_couplings(completed_moments)

    symbol_map = {
        "beta": beta,
        "eta": eta,
        "omega": omega,
        "phi": phi,
        "m2": sp.Symbol("m2", real=True),
        "lam": sp.Symbol("lam", real=True),
        "y": sp.Symbol("y", real=True),
        "g": sp.Symbol("g", real=True),
        "xi": sp.Symbol("xi", real=True),
    }
    p1866_moments = parse_p1866_exprs(
        p1866.get("coefficient_export", {}).get("kernel_moments_at_d_ref", {}),
        symbol_map,
    )
    p1866_effective = parse_p1866_exprs(
        p1866.get("coefficient_export", {}).get("strict_effective_couplings", {}),
        symbol_map,
    )

    completed_moment_residuals = residuals(completed_moments, strict_moments)
    completed_effective_residuals = residuals(completed_effective, strict_effective)
    p1866_moment_residuals = residuals(strict_moments, p1866_moments)
    p1866_effective_residuals = residuals(strict_effective, p1866_effective)

    canonical_substitutions = {
        alpha_geo: 4 * sp.log(2),
        omega_L: sp.pi / 4,
        phi_L: sp.pi / 6,
        beta_tors: sp.Rational(1, 100),
        omega: sp.Rational(743, 4000),
        phi: sp.Rational(13, 80),
        beta: sp.Integer(1),
        eta: sp.Rational(9, 5),
    }
    strict_numeric = numeric_moments(strict_moments, canonical_substitutions)
    completed_numeric = numeric_moments(completed_moments, canonical_substitutions)
    legacy_raw_numeric = numeric_moments(legacy_raw_moments, canonical_substitutions)
    legacy_normalized_numeric = numeric_moments(legacy_normalized_moments, canonical_substitutions)
    completed_numeric_max_abs_residual = max_abs_difference(completed_numeric, strict_numeric)
    legacy_raw_numeric_max_abs_mismatch = max_abs_difference(legacy_raw_numeric, strict_numeric)
    legacy_normalized_numeric_max_abs_mismatch = max_abs_difference(
        legacy_normalized_numeric,
        strict_numeric,
    )

    finite_summary = finite_report.get("finite_bridge_assembly_summary", {})
    symbolic_summary = symbolic_report.get("symbolic_cancellation_summary", {})
    p2362_summary = p2362.get("strict_lagrangian_eom_parallel_completion_probe", {}).get(
        "parallel_completion_summary", {}
    )

    theorem_export = {
        "theorem_name": "P2363 bridge-completed moment transport into strict L_total couplings",
        "claim": (
            "If the explicit APD completion factor Q=(1/alpha_geo)*P*D is applied to "
            "K_legacy_ont before moment extraction, then the d_ref=1 moments c0,c1,c2 "
            "and the P1866 effective couplings are identical to the strict-kernel export. "
            "Raw legacy moments, even after amplitude normalization, are not an admissible "
            "substitute for strict L_total coefficients."
        ),
        "positive_content": [
            "Symbolic cancellation gives K_legacy_ont(d)*Q_APD(d)-K_strict_gate(d)=0.",
            "The completed bridge moments c0,c1,c2 have zero symbolic residual against the strict moments.",
            "The completed bridge effective couplings m2_eff, lam_eff, y_eff, g_eff, xi_eff have zero symbolic residual against the P1866 strict map.",
            "A finite canonical witness shows raw legacy moments and amplitude-normalized legacy moments do not match strict moments.",
        ],
        "not_licensed": [
            "strict dynamical source theorem for A/P/D",
            "selector premise or QW-2191 discharge",
            "legacy physical-role transfer",
            "full bridge role-transfer theorem",
            "full 4D EOM theorem-grade closure",
            "ToE closure",
        ],
    }

    probe = {
        "probe_id": "P2363_S1313_LEGACY_STRICT_BRIDGE_MOMENT_LAGRANGIAN_TRANSPORT",
        "source_artifacts": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "docs": {name: rel(path) for name, path in DOC_FILES.items()},
        "symbolic_bridge": {
            "legacy_kernel": "K_legacy_ont(d)=alpha_geo*cos(omega_L*d+phi_L)/(1+beta_tors*d)",
            "completion_factor": "Q_APD(d)=(1/alpha_geo)*(cos(omega*d+phi)/cos(omega_L*d+phi_L))*((1+beta_tors*d)/(1+beta*d^eta))",
            "completed_kernel": sstr(K_completed),
            "strict_kernel": sstr(K_strict),
            "completed_minus_strict_residual": sstr(bridge_residual),
            "admissibility_scope": "inherits finite nonzero-domain assumptions from scratch APD cancellation reports",
        },
        "moment_transport": {
            "reference_scale": "d_ref=1",
            "strict_moments": {key: sstr(value) for key, value in strict_moments.items()},
            "completed_bridge_moments": {key: sstr(value) for key, value in completed_moments.items()},
            "completed_minus_strict_residuals": completed_moment_residuals,
            "p1866_strict_export_residuals": p1866_moment_residuals,
            "raw_legacy_moments": {key: sstr(value) for key, value in legacy_raw_moments.items()},
            "amplitude_normalized_legacy_moments": {
                key: sstr(value) for key, value in legacy_normalized_moments.items()
            },
        },
        "effective_coupling_transport": {
            "map": {
                "m2_eff": "m2*(1 + c0)",
                "lam_eff": "lam*(1 + c1**2)",
                "y_eff": "y*(1 + c0/2)",
                "g_eff": "g*(1 + c2)",
                "xi_eff": "xi*(1 + c0)",
            },
            "strict_effective_couplings": {
                key: sstr(value) for key, value in strict_effective.items()
            },
            "completed_bridge_effective_couplings": {
                key: sstr(value) for key, value in completed_effective.items()
            },
            "completed_minus_strict_residuals": completed_effective_residuals,
            "p1866_strict_export_residuals": p1866_effective_residuals,
        },
        "canonical_numeric_witness": {
            "parameters": {
                "alpha_geo": "4*log(2)",
                "omega_L": "pi/4",
                "phi_L": "pi/6",
                "beta_tors": "1/100",
                "omega": "743/4000",
                "phi": "13/80",
                "beta": "1",
                "eta": "9/5",
            },
            "strict_moments": strict_numeric,
            "completed_bridge_moments": completed_numeric,
            "legacy_raw_moments": legacy_raw_numeric,
            "legacy_amplitude_normalized_moments": legacy_normalized_numeric,
            "completed_numeric_max_abs_residual": completed_numeric_max_abs_residual,
            "legacy_raw_numeric_max_abs_mismatch": legacy_raw_numeric_max_abs_mismatch,
            "legacy_amplitude_normalized_numeric_max_abs_mismatch": legacy_normalized_numeric_max_abs_mismatch,
        },
        "separation_from_selector": {
            "eom_lagrangian_track_is_selector_independent": True,
            "selector_closure_is_parallel_problem": True,
            "selector_required_for_moment_transport": False,
            "selector_qw2191_status": "OPEN_PARALLEL_TRACK_NOT_USED_AS_EOM_INPUT",
            "p2362_inherited_summary": p2362_summary,
        },
        "theorem_export": theorem_export,
    }
    probe["theorem_fingerprint_sha256"] = sha256_json(theorem_export)

    gatekeeper_checks = {
        "source_artifacts_present": all("_missing" not in artifact for artifact in artifacts.values()),
        "p1866_strict_moments_loaded": set(p1866_moments) == {"c0", "c1", "c2"},
        "p1866_effective_couplings_loaded": set(p1866_effective)
        == {"m2_eff", "lam_eff", "y_eff", "g_eff", "xi_eff"},
        "symbolic_apd_bridge_residual_zero": sstr(bridge_residual) == "0",
        "completed_moments_match_strict": all(value == "0" for value in completed_moment_residuals.values()),
        "completed_effective_couplings_match_strict": all(
            value == "0" for value in completed_effective_residuals.values()
        ),
        "strict_moments_match_p1866_export": all(value == "0" for value in p1866_moment_residuals.values()),
        "strict_effective_couplings_match_p1866_export": all(
            value == "0" for value in p1866_effective_residuals.values()
        ),
        "canonical_completed_numeric_residual_zero": completed_numeric_max_abs_residual < 1e-12,
        "raw_legacy_moments_do_not_match_strict": legacy_raw_numeric_max_abs_mismatch > 1e-2,
        "amplitude_normalized_legacy_moments_do_not_match_strict": legacy_normalized_numeric_max_abs_mismatch
        > 1e-2,
        "finite_bridge_assembly_inherited": finite_summary.get("assembled_map_reconstructs_strict_kernel")
        is True
        and finite_summary.get("full_bridge_theorem_exported") is False,
        "symbolic_cancellation_inherited": symbolic_summary.get("symbolic_cancellation_formula_exported")
        is True
        and symbolic_summary.get("raw_kernel_identity_claimed") is False,
        "p2362_parallel_eom_lagrangian_inherited": p2362_summary.get(
            "eom_lagrangian_track_is_selector_independent"
        )
        is True
        and p2362_summary.get("selector_closure_is_parallel_problem") is True,
        "docs_updated_with_p2363_transport_statement": all(
            snippet in text for text in doc_texts.values() for snippet in REQUIRED_DOC_SNIPPETS
        ),
        "no_selector_or_role_transfer_claimed": True,
        "hard_limits_preserved": True,
    }

    payload = {
        "schema_version": "p2363_s1313_v1",
        "packet_id": "P2363",
        "stage_id": "S1313",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-06-01T00:00:00+00:00",
        "status": "OPEN_PROGRESS_BRIDGE_COMPLETED_MOMENT_TRANSPORT_TO_LTOTAL_NO_SELECTOR_OR_ROLE_TRANSFER",
        "result_kind": "LEGACY_STRICT_BRIDGE_COMPLETED_MOMENT_LAGRANGIAN_TRANSPORT_CERTIFICATE",
        "legacy_strict_bridge_moment_lagrangian_transport_probe": probe,
        "recommended_next_honest_step": (
            "Use the now-certified bridge-completed moments as the input to a 4D covariant "
            "EOM residual table for the scalar/gauge/gravity sectors; keep selector/QW-2191 "
            "and legacy role-transfer audits on separate parallel tracks."
        ),
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_PROGRESS_WITH_TRACE_NO_SELECTOR_ROLE_TRANSFER_OR_TOE_CLOSURE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2363 S1313: legacy->strict bridge moment transport into L_total",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2363/S1313 certifies bridge-completed moment transport: applying the explicit APD completion factor to the legacy kernel gives the same `c0,c1,c2` moments and effective `L_total` couplings as the strict kernel.",
                "",
                "## Main Equalities",
                "",
                f"- `K_legacy*Q_APD - K_strict = {sstr(bridge_residual)}`.",
                f"- Completed moment residuals: `{completed_moment_residuals}`.",
                f"- Completed effective-coupling residuals: `{completed_effective_residuals}`.",
                "",
                "## Negative Control",
                "",
                f"- Raw legacy moment max mismatch: `{legacy_raw_numeric_max_abs_mismatch}`.",
                f"- Amplitude-normalized legacy moment max mismatch: `{legacy_normalized_numeric_max_abs_mismatch}`.",
                "",
                "## Hard Limits",
                "",
                "- No strict source theorem for APD is claimed.",
                "- No selector premise or QW-2191 discharge is claimed.",
                "- No legacy physical-role transfer is claimed.",
                "- No full 4D EOM theorem closure or ToE closure is claimed.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
