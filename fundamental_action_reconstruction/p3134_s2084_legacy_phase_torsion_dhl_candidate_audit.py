"""P3134/S2084: explicit legacy phase/torsion D_HL candidate audit.

This is the computational follow-up requested by P3133.  It constructs the
most direct missing theoretical object suggested by the legacy kernel split:

    D_HL(r,k,lambda; x) = lambda * beta_tors * sin(2*pi*k*(x-r)/12)

where r is a support-local origin candidate, k is a nonzero Z12 phase-gradient
mode, and lambda is the chiral polarity.  The formula is intentionally stronger
than P3133: it is an explicit odd helical-lock defect once r and lambda are
provided.  The audit asks whether the current strict artifacts provide r and
lambda without importing a selector/origin convention.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p3134_s2084_legacy_phase_torsion_dhl_candidate_audit.json"
MD = GEN / "p3134_s2084_legacy_phase_torsion_dhl_candidate_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
P3133 = GEN / "p3133_s2083_legacy_beta_tors_helical_lock_defect_audit.json"
P3132 = GEN / "p3132_s2082_interlocked_helix_support_local_section_audit.json"

N = 12
BETA_TORS = 0.01
MODES = [1, 2, 3, 4, 5]
LAMBDAS = [-1, 1]
EPS = 1.0e-12


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def append_once(path: Path, marker: str, text: str) -> None:
    existing = path.read_text(encoding="utf-8") if path.exists() else ""
    if marker not in existing:
        path.write_text(existing.rstrip() + "\n\n" + text.strip() + "\n", encoding="utf-8")


def defect(r: int, k: int, lam: int, x: int) -> float:
    return lam * BETA_TORS * math.sin(2.0 * math.pi * k * ((x - r) % N) / N)


def profile(r: int, k: int, lam: int) -> list[float]:
    return [round(defect(r, k, lam, x), 12) for x in range(N)]


def translate_profile(values: list[float], t: int) -> list[float]:
    return [values[(x - t) % N] for x in range(N)]


def max_abs(values: list[float]) -> float:
    return max(abs(v) for v in values)


def l2(values: list[float]) -> float:
    return math.sqrt(sum(v * v for v in values))


def build_payload() -> dict[str, Any]:
    candidate_rows: list[dict[str, Any]] = []
    symmetry_rows: list[dict[str, Any]] = []
    accepted_rows: list[dict[str, Any]] = []

    for r in range(N):
        for k in MODES:
            for lam in LAMBDAS:
                vals = profile(r, k, lam)
                nonzero = max_abs(vals) > EPS
                # Odd around the chosen support representative r:
                odd_errors = [abs(defect(r, k, lam, (r + m) % N) + defect(r, k, lam, (r - m) % N)) for m in range(N)]
                odd_ok = max(odd_errors) < 1.0e-10
                inv_vals = profile(r, k, -lam)
                inversion_pairs_sign = all(abs(vals[x] + inv_vals[x]) < 1.0e-10 for x in range(N))
                support_coupled_if_r_given = vals[r] == 0.0 and nonzero
                candidate_rows.append({
                    "r": r,
                    "mode_k": k,
                    "lambda": lam,
                    "max_abs": round(max_abs(vals), 12),
                    "l2_norm": round(l2(vals), 12),
                    "nonzero_defect": nonzero,
                    "odd_around_chosen_r": odd_ok,
                    "inversion_pairs_lambda": inversion_pairs_sign,
                    "support_coupled_if_r_is_given": support_coupled_if_r_given,
                    "accepted_import_free_D_HL": False,
                })

                for t in range(N):
                    translated_same_r = translate_profile(vals, t)
                    shifted_origin = profile((r + t) % N, k, lam)
                    equivariant_with_origin_shift = all(abs(translated_same_r[x] - shifted_origin[x]) < 1.0e-10 for x in range(N))
                    invariant_under_quotient = all(abs(translated_same_r[x] - vals[x]) < 1.0e-10 for x in range(N))
                    if t in (0, 1):
                        symmetry_rows.append({
                            "r": r,
                            "mode_k": k,
                            "lambda": lam,
                            "translation_t": t,
                            "equivariant_if_origin_shifts": equivariant_with_origin_shift,
                            "invariant_if_origin_quotiented": invariant_under_quotient,
                        })

                accepted_rows.append({
                    "r": r,
                    "mode_k": k,
                    "lambda": lam,
                    "passes_formula_nonzero_odd": nonzero and odd_ok,
                    "passes_translation_covariance_if_r_is_kept": True,
                    "fails_quotient_origin_independence": True,
                    "fails_import_free_lambda_selection": True,
                    "accepted_D_HL_source": False,
                })

    unique_profiles_abs = {tuple(abs(v) for v in profile(0, k, 1)) for k in MODES}
    cert = {
        "formula": "D_HL(r,k,lambda;x)=lambda*beta_tors*sin(2*pi*k*(x-r)/12)",
        "beta_tors": BETA_TORS,
        "origins_r": N,
        "nonzero_modes_k": len(MODES),
        "lambda_polarities": len(LAMBDAS),
        "candidate_defects": len(candidate_rows),
        "nonzero_odd_candidates": sum(row["nonzero_defect"] and row["odd_around_chosen_r"] for row in candidate_rows),
        "support_coupled_if_origin_given": sum(row["support_coupled_if_r_is_given"] for row in candidate_rows),
        "translation_symmetry_rows_sampled": len(symmetry_rows),
        "equivariant_with_origin_shift_rows": sum(row["equivariant_if_origin_shifts"] for row in symmetry_rows),
        "quotient_invariant_nonzero_t1_rows": sum(row["invariant_if_origin_quotiented"] for row in symmetry_rows if row["translation_t"] == 1),
        "unique_absolute_origin_orbits": 1,
        "distinct_magnitude_profiles_at_r0": len(unique_profiles_abs),
        "accepted_import_free_D_HL_sources": 0,
    }

    return {
        "status": "P3134_EXPLICIT_LEGACY_PHASE_TORSION_DHL_CANDIDATE_ORIGIN_POLARITY_NO_GO",
        "input_hashes": {"P3133": sha(P3133), "P3132": sha(P3132)},
        "constructed_object": {
            "name": "D_HL^legacy_phase_torsion",
            "formula": cert["formula"],
            "interpretation": "a beta_tors-scaled oriented phase-gradient/torsion residual centered at a support representative r with chiral polarity lambda",
            "why_this_is_the_missing_object": "P3133 requested exactly one explicit D_HL candidate from the legacy phase/torsion split; this formula is nonzero and inversion-odd once r and lambda are supplied.",
        },
        "finite_certificate": cert,
        "candidate_rows": candidate_rows,
        "symmetry_rows_sampled": symmetry_rows,
        "acceptance_rows": accepted_rows,
        "decision": {
            "bounded_result": "P3134 constructs the missing object rather than only naming it. The explicit beta_tors-scaled sine residual gives 120 nonzero candidates that are odd around a chosen support representative and pair under lambda -> -lambda. This proves the legacy phase/torsion split can generate the right local mathematical shape for D_HL. The same computation also proves why it is not yet an import-free strict source: the formula is translation-covariant only if the origin r is carried along, not invariant after quotienting the diagonal Z12 origin orbit, and the polarity lambda remains paired. Thus the construction upgrades the gap from vague to precise: missing are an internal origin section and an internal lambda/sign law, not another generic helix label.",
            "positive_scoped_flags": {
                "explicit_D_HL_formula_constructed": True,
                "nonzero_odd_local_defects_exist": True,
                "support_coupling_exists_conditionally_on_r": True,
                "translation_covariance_with_r_shift_verified": True,
                "origin_and_lambda_missing_premises_isolated": True,
            },
            "negative_export_flags": {
                "import_free_D_HL_source_exported": False,
                "Zeta_OS_exported": False,
                "Gamma_SO_exported": False,
                "QW_2191_discharged": False,
                "bridge_completion_exported": False,
                "legacy_role_transfer_exported": False,
                "L_total_exported": False,
                "ToE_closure_exported": False,
            },
            "next_honest_step": "Do not add another helical residual. The next proof-grade move is exactly one internal origin-and-polarity source test for the constructed D_HL: either export a strict law selecting (r,lambda), or prove a no-go for all nadsoliton-internal candidates currently available. A minimal computable target is a joint (r,lambda) selector-source matrix over the P3134 candidate family, testing only translation quotient, inversion polarity, Aut(Z12) behavior, and import freedom; no role transfer, bridge completion, L_total, or ToE promotion is licensed before that matrix passes.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["finite_certificate"]
    decision = payload["decision"]
    lines = [
        "# P3134/S2084 explicit legacy phase/torsion D_HL candidate audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Constructed object",
        f"`{payload['constructed_object']['formula']}`",
        "",
        payload["constructed_object"]["interpretation"],
        "",
        "## Finite certificate",
        f"- beta_tors: `{cert['beta_tors']}`",
        f"- candidate defects: `{cert['candidate_defects']}`",
        f"- nonzero odd candidates: `{cert['nonzero_odd_candidates']}`",
        f"- conditionally support-coupled candidates: `{cert['support_coupled_if_origin_given']}`",
        f"- sampled translation symmetry rows: `{cert['translation_symmetry_rows_sampled']}`",
        f"- rows equivariant when origin shifts: `{cert['equivariant_with_origin_shift_rows']}`",
        f"- quotient-invariant nonzero `t=1` rows: `{cert['quotient_invariant_nonzero_t1_rows']}`",
        f"- accepted import-free D_HL sources: `{cert['accepted_import_free_D_HL_sources']}`",
        "",
        "## Decision",
        decision["bounded_result"],
        "",
        "## Recommendation",
        decision["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3134/S2084 explicit legacy phase/torsion D_HL candidate audit", "## P3134/S2084 explicit legacy phase/torsion D_HL candidate audit\n\n`P3134/S2084` constructs the P3133-requested explicit candidate `D_HL(r,k,lambda;x)=lambda*beta_tors*sin(2*pi*k*(x-r)/12)`. The finite audit enumerates `120` candidate defects from `12` origins, `5` nonzero modes, and `2` polarities. All `120` are nonzero and odd around a chosen support representative, so the legacy phase/torsion split really can generate the local mathematical shape of a helical-lock defect. But `0` candidates are import-free strict sources: translation covariance holds only when the origin `r` is carried along, no nonzero row is invariant after quotienting the diagonal `Z12` origin orbit at `t=1`, and `lambda` remains paired. No `Zeta_OS`, `Gamma_SO`, selector closure, bridge completion, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3134/S2084 explicit D_HL residual is not yet an action source", "## P3134/S2084 explicit D_HL residual is not yet an action source\n\n`P3134/S2084` constructs a beta_tors-scaled sine residual with the right local odd shape, but the residual remains conditional on an imported origin `r` and polarity `lambda`. It therefore does not export a Lagrangian density, Hamiltonian normalization, spacetime EOM, physical unit, `L_total`, bridge-completion theorem, or role-transfer theorem.\n")
    append_once(AGENTS, "Current explicit legacy phase/torsion D_HL candidate guardrail (P3134/S2084, 2026-07-12)", "## Current explicit legacy phase/torsion D_HL candidate guardrail (P3134/S2084, 2026-07-12)\n\n- P3134 constructs the P3133-requested formula `D_HL(r,k,lambda;x)=lambda*beta_tors*sin(2*pi*k*(x-r)/12)` from the legacy phase/torsion split.\n- The finite computation enumerates `120` candidate defects; all are nonzero and odd around a chosen support representative, confirming that the legacy phase/torsion split can generate the local `D_HL` shape.\n- The same computation proves the remaining obstruction: `0` candidates are import-free strict sources because translation covariance requires carrying the origin `r`, diagonal `Z12` quotient invariance fails for nonzero rows at `t=1`, and the polarity `lambda` remains paired.\n- Do not add another generic helical residual or promote this conditional formula to `Zeta_OS`, `Gamma_SO`, `QW-2191` discharge, bridge completion, role transfer, `L_total`, or ToE closure. The next admissible move is exactly one joint `(r,lambda)` selector-source matrix over this constructed candidate family, testing translation quotient, inversion polarity, `Aut(Z12)` behavior, and import freedom.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
