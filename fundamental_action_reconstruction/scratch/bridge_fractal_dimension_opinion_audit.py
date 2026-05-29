#!/usr/bin/env python3
"""Scratch audit of the external AI opinion about a fractal-dimension bridge.

The opinion under review suggests that D_F≈1.73 might be the missing link
between eta_eff≈1.74 and strict eta=1.8, and proposes several anomalous-
diffusion / Hausdorff / Z12 formulas.  This audit separates:

1. repo-supported facts,
2. numerically plausible but unsupported formulas,
3. tautological parameter choices,
4. theorem obligations that remain open.

No bridge theorem or physical-role transfer is exported.
"""
from __future__ import annotations

import json
import math
import subprocess
from pathlib import Path
from typing import Any

import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FAR = ROOT / "fundamental_action_reconstruction"
F151 = FAR / "F151_FIRST_LEGACY_TO_STRICT_KERNEL_BRIDGE_TARGET_PACKET.md"
DF_AUDIT = HERE / "bridge_df_eta_eff_report.json"
TRANSPORT_AUDIT = HERE / "bridge_df_transport_exponent_report.json"
OUT_JSON = HERE / "bridge_fractal_dimension_opinion_audit_report.json"
OUT_MD = HERE / "bridge_fractal_dimension_opinion_audit_report.md"

STRICT_ETA = 1.8
ALPHA_GEO = 4.0 * math.log(2.0)
DF_MINUS_ONE = ALPHA_GEO - 1.0
Z12_N = 12.0


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rg_hits(pattern: str, *, include_scratch: bool = False, limit: int = 80) -> list[str]:
    cmd = [
        "rg",
        "-n",
        "-S",
        pattern,
        "fundamental_action_reconstruction",
        "-g",
        "*.md",
        "-g",
        "*.py",
        "-g",
        "*.json",
    ]
    if not include_scratch:
        cmd.extend(["-g", "!fundamental_action_reconstruction/scratch/**"])
    proc = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, check=False)
    return proc.stdout.splitlines()[:limit]


def rel_gap(a: float, b: float) -> float:
    return abs(a - b) / max(abs(b), 1e-15)


def formula_audit(eta_eff: float) -> dict[str, Any]:
    # Opinion formulas evaluated both on a hypothetical D≈eta_eff and on the
    # repo-supported D_f-1 candidate.  The goal is not to endorse them but to
    # make their assumptions explicit.
    hypothetical_df = eta_eff
    embed_dim = 2.0
    spectral_dim_assumed = 2.0
    z12_lambda_required_for_eta = Z12_N ** (1.0 / STRICT_ETA)
    z12_lambda_required_for_eta_eff = Z12_N ** (1.0 / eta_eff)
    z12_lambda_required_for_df_minus_one = Z12_N ** (1.0 / DF_MINUS_ONE)

    return {
        "direct_dimension_identification": {
            "formula": "eta = D_candidate",
            "eta_eff_as_dimension_gap_to_strict_eta": STRICT_ETA - eta_eff,
            "D_f_minus_1_gap_to_strict_eta": STRICT_ETA - DF_MINUS_ONE,
            "verdict": "numerically plausible for D_f-1, but still requires a theorem identifying the damping exponent with a fractal/transport dimension",
        },
        "spectral_dimension_relation": {
            "formula": "eta = 2*D_candidate/d_s",
            "with_D_candidate_eta_eff_and_d_s_2": 2.0 * hypothetical_df / spectral_dim_assumed,
            "with_D_candidate_D_f_minus_1_and_d_s_2": 2.0 * DF_MINUS_ONE / spectral_dim_assumed,
            "required_d_s_for_strict_eta_from_D_f_minus_1": 2.0 * DF_MINUS_ONE / STRICT_ETA,
            "verdict": "unsupported here: repo grep does not export a spectral dimension d_s bridge for this kernel exponent",
        },
        "walk_dimension_relation_from_opinion": {
            "formula": "eta = 2 - D_candidate/D_embed",
            "with_D_candidate_eta_eff_and_D_embed_2": 2.0 - hypothetical_df / embed_dim,
            "with_D_candidate_D_f_minus_1_and_D_embed_2": 2.0 - DF_MINUS_ONE / embed_dim,
            "verdict": "numerically too small for eta≈1.8 under D_embed=2; not the live bridge route",
        },
        "floor_fraction_formula": {
            "formula": "eta = D_candidate - floor(D_candidate) + 1",
            "with_repo_D_f_4ln2": ALPHA_GEO - math.floor(ALPHA_GEO) + 1.0,
            "with_D_candidate_eta_eff": hypothetical_df - math.floor(hypothetical_df) + 1.0,
            "verdict": "for repo D_f=4ln2 this equals D_f-1 exactly; this is consistent with the D_f-1 route but does not independently prove it",
        },
        "z12_hausdorff_scaling_claim": {
            "formula": "D = ln(12)/ln(lambda)",
            "lambda_required_for_D_equals_strict_eta": z12_lambda_required_for_eta,
            "lambda_required_for_D_equals_eta_eff": z12_lambda_required_for_eta_eff,
            "lambda_required_for_D_equals_D_f_minus_1": z12_lambda_required_for_df_minus_one,
            "independent_lambda_export_found": False,
            "verdict": "choosing lambda=12^(1/eta) is tautological unless repo exports an independent scaling lambda",
        },
    }


def symbolic_checks() -> dict[str, str]:
    d, beta, Df, lam = sp.symbols("d beta D_f lambda", positive=True, real=True)
    gamma = Df - 1
    denominator = 1 + beta * d**gamma
    log_slope = sp.simplify(d * sp.diff(sp.log(beta * d**gamma), d))
    z12_dim = sp.log(12) / sp.log(lam)
    lambda_sol = sp.solve(sp.Eq(z12_dim, sp.Rational(9, 5)), lam)[0]
    return {
        "transport_denominator_from_gamma_D_f_minus_1": sp.sstr(denominator),
        "coupling_log_slope": sp.sstr(log_slope),
        "z12_dimension_formula": sp.sstr(z12_dim),
        "lambda_solving_z12_dimension_equals_9_over_5": sp.sstr(lambda_sol),
    }


def main() -> None:
    df_audit = load_json(DF_AUDIT)
    transport = load_json(TRANSPORT_AUDIT)
    eta_eff = float(df_audit["candidate"]["numerical_constants"]["eta_eff_from_bridge_missing_coupling"])
    f151_text = F151.read_text(encoding="utf-8") if F151.exists() else ""

    grep_audit = {
        "non_scratch_D_f_semantic_hits": rg_hits(r"\bD_F\b|\bD_f\b|D_f\s*≡|alpha_geo\s*≡\s*4 ln 2", limit=120),
        "non_scratch_D_f_minus_one_transport_hits": rg_hits(r"D_f\s*-\s*1|D_F\s*-\s*1|p\^\{D_f-1\}|d\^\{D_f - 1\}|N\(d\).*D_f", limit=120),
        "non_scratch_1p73_1p74_eta_eff_hits": rg_hits(r"eta_eff|1\.73[0-9]*|1\.74[0-9]*", limit=120),
        "non_scratch_hausdorff_walk_spectral_hits": rg_hits(r"Hausdorff|hausdorff|walk.?dim|spectral dimension|d_s\b|d_w\b|subdiffusion|anomalous.diffusion|fractal_dim", limit=120),
        "scratch_bridge_hits": rg_hits(r"eta_eff|D_f-1|transport exponent|bridge_df", include_scratch=True, limit=120),
    }

    f151_audit = {
        "F151_file_exists": F151.exists(),
        "bridge_target_present": "B_legacy_strict_bridge_target_v1" in f151_text,
        "component_packet_present": "Gamma_bridge_components_target_v1" in f151_text,
        "damping_component_present": "R_damp_renorm_target_v1" in f151_text,
        "future_only_no_discharge_scope_present": "actual bridge discharge" in f151_text and "no false pass" in f151_text.lower(),
    }

    formulas = formula_audit(eta_eff)
    symbolic = symbolic_checks()
    opinion_verdict = {
        "claim_DF_1p73_absent_from_repo": {
            "verdict": "partly_true_only_if_read_as_no_exported_symbol_D_F_equals_1p73_in_non_scratch_FAR; false as a dismissal of the D_f link because non-scratch FAR exports D_f and D_f-1 measure/path-growth structures, and scratch now contains eta_eff evidence",
            "non_scratch_1p73_1p74_hit_count": len(grep_audit["non_scratch_1p73_1p74_eta_eff_hits"]),
        },
        "claim_fractal_dimension_is_missing_idea": {
            "verdict": "too_strong; F151 exports the bridge target and non-scratch FAR contains D_f/D_f-1 structures, but no discharged theorem yet ties them to eta",
        },
        "claim_component_2_damping_is_the_live_target": {
            "verdict": "supported_by_F151_target_scope; R_damp_renorm_target_v1 is present, but only future-only and below bridge discharge",
        },
        "best_supported_next_formula": {
            "verdict": "D_f_minus_1_transport_denominator remains the best supported route: it has repo provenance, a symbolic denominator insertion ansatz, and the prior transport probe captures >=99% of linear-to-free exponent improvement across tested windows",
            "transport_probe_min_capture_fraction": transport["aggregate_summary"]["min_capture_fraction"],
            "transport_probe_min_delta_bic_linear_minus_df": transport["aggregate_summary"]["min_delta_bic_linear_minus_df"],
        },
    }

    report = {
        "status": "OPEN_OPINION_AUDIT_WITH_THEOREM_PREP_TRACE_NO_BRIDGE_THEOREM",
        "result_kind": "SCRATCH_FRACTAL_DIMENSION_OPINION_AUDIT__NOT_A_THEOREM",
        "source_reports": {
            "df_eta_eff_audit": str(DF_AUDIT.relative_to(ROOT)),
            "df_transport_exponent_audit": str(TRANSPORT_AUDIT.relative_to(ROOT)),
            "F151_bridge_target_packet": str(F151.relative_to(ROOT)),
        },
        "constants": {
            "eta_eff": eta_eff,
            "strict_eta": STRICT_ETA,
            "D_f_4ln2": ALPHA_GEO,
            "D_f_minus_1": DF_MINUS_ONE,
            "gap_eta_eff_to_D_f_minus_1": DF_MINUS_ONE - eta_eff,
            "gap_D_f_minus_1_to_strict_eta": STRICT_ETA - DF_MINUS_ONE,
            "rel_gap_midpoint_to_D_f_minus_1": df_audit["candidate"]["deltas"]["rel_midpoint_eta_eff_strict_minus_D_f_minus_1"],
        },
        "f151_audit": f151_audit,
        "grep_audit": grep_audit,
        "formula_audit": formulas,
        "symbolic_checks": symbolic,
        "opinion_verdict": opinion_verdict,
        "hard_limits": [
            "No Hausdorff, spectral-dimension, walk-dimension, or Z12 scaling formula is promoted to theorem status here.",
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role formula is transferred to the strict kernel.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Stay on the D_f-1 damping-renormalization component: prove or falsify that the FAR/A4 measure factor or legacy path-counting factor enters the kernel damping denominator, and reject Z12/Hausdorff formulas unless an independent scaling lambda or spectral/walk dimension is exported.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch audit: external fractal-dimension bridge opinion\n\n"
        "Status: opinion checked; theorem-prep trace strengthened; no bridge theorem.\n\n"
        f"- F151 bridge target exists `{f151_audit['F151_file_exists']}` and damping component `R_damp_renorm_target_v1` exists `{f151_audit['damping_component_present']}`.\n"
        f"- Constants: `eta_eff={eta_eff:.15f}`, `D_f-1={DF_MINUS_ONE:.15f}`, strict `eta={STRICT_ETA:.15f}`.\n"
        f"- Best-supported route remains `D_f-1` transport denominator; prior min capture `{transport['aggregate_summary']['min_capture_fraction']:.12f}`, min ΔBIC `{transport['aggregate_summary']['min_delta_bic_linear_minus_df']:.6f}`.\n"
        "- Hausdorff/spectral/walk/Z12 formulas are not promoted: no independent `lambda`, `d_s`, or `d_w` export was found by this audit.\n"
        "- No false pass: no kernel identity, no physical-role transfer, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
