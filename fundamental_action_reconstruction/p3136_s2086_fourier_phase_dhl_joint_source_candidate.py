"""P3136/S2086: Fourier-phase J_DHL candidate audit.

P3135 demanded one formula-level joint source law J_DHL(support, field)->(r,lambda).
This audit constructs the sharpest natural candidate for the P3134 sine residual:
read the phase of its discrete Fourier coefficient.  The computation proves what
this object can and cannot source.
"""

from __future__ import annotations

import cmath
import hashlib
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p3136_s2086_fourier_phase_dhl_joint_source_candidate.json"
MD = GEN / "p3136_s2086_fourier_phase_dhl_joint_source_candidate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
P3135 = GEN / "p3135_s2085_dhl_origin_polarity_selector_source_matrix.json"
P3134 = GEN / "p3134_s2084_legacy_phase_torsion_dhl_candidate_audit.json"

N = 12
BETA_TORS = 0.01
MODES = [1, 2, 3, 4, 5]
LAMBDAS = [-1, 1]
EPS = 1e-9


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def append_once(path: Path, marker: str, text: str) -> None:
    existing = path.read_text(encoding="utf-8") if path.exists() else ""
    if marker not in existing:
        path.write_text(existing.rstrip() + "\n\n" + text.strip() + "\n", encoding="utf-8")


def defect_profile(r: int, k: int, lam: int) -> list[float]:
    return [lam * BETA_TORS * math.sin(2 * math.pi * k * ((x - r) % N) / N) for x in range(N)]


def fourier_coeff(vals: list[float], k: int) -> complex:
    return sum(vals[x] * cmath.exp(-2j * math.pi * k * x / N) for x in range(N))


def same_profile(a: list[float], b: list[float]) -> bool:
    return all(abs(x - y) < EPS for x, y in zip(a, b))


def matching_pairs(vals: list[float], k: int) -> list[tuple[int, int]]:
    return [(r, lam) for r in range(N) for lam in LAMBDAS if same_profile(vals, defect_profile(r, k, lam))]


def translate(vals: list[float], t: int) -> list[float]:
    return [vals[(x - t) % N] for x in range(N)]


def build_payload() -> dict[str, Any]:
    rows = []
    representative_rows = []
    for k in MODES:
        for r in range(N):
            for lam in LAMBDAS:
                vals = defect_profile(r, k, lam)
                coeff = fourier_coeff(vals, k)
                matches = matching_pairs(vals, k)
                alias_size = len(matches)
                full_pair = alias_size == 1 and matches[0] == (r, lam)
                # Coefficient is nonzero and its phase changes under translation; this is the chart-relative extraction witness.
                shifted_coeff = fourier_coeff(translate(vals, 1), k)
                phase_rotates = abs(shifted_coeff - coeff * cmath.exp(-2j * math.pi * k / N)) < 1e-9
                row = {
                    "mode_k": k,
                    "r": r,
                    "lambda": lam,
                    "gcd_k_12": math.gcd(k, N),
                    "coefficient_abs": round(abs(coeff), 12),
                    "nonzero_coefficient": abs(coeff) > EPS,
                    "matching_pair_count": alias_size,
                    "chart_relative_full_pair_recovery": full_pair,
                    "translation_phase_rotates_at_t1": phase_rotates,
                    "strict_import_free_recovery": False,
                }
                rows.append(row)
                if r == 0 and lam == 1:
                    representative_rows.append({**row, "matching_pairs": [list(pair) for pair in matches]})

    unit_mode_rows = [row for row in rows if row["gcd_k_12"] == 1]
    nonunit_mode_rows = [row for row in rows if row["gcd_k_12"] != 1]
    return {
        "status": "P3136_FOURIER_PHASE_J_DHL_CANDIDATE_CHART_RELATIVE_POSITIVE_STRICT_SOURCE_NO_GO",
        "input_hashes": {"P3135": sha(P3135), "P3134": sha(P3134)},
        "repo_backscan_summary": [
            "P2992/P2994 Fourier-character source-atom lanes already warn that exact character receivers do not provide strict provenance or atom-specific source coupling.",
            "P3039 already tested chi_i=sin(2*pi*i/12): the torsor is real and inversion-odd, but translations move the phase origin and current receivers choose chart-relative representatives.",
            "P2869 warns that Aut-character Fourier idempotents can represent endpoints only by importing projector/polarity coefficients.",
        ],
        "constructed_object": {
            "name": "J_DHL^FourierPhase",
            "formula": "J reads C_k(D)=sum_x D(x)*exp(-2*pi*i*k*x/12); phase(C_k) reconstructs r modulo gcd(k,12), and sign convention reconstructs lambda when a Fourier frame is fixed.",
            "scope": "formula-level extractor for the P3134 D_HL sine residual family",
        },
        "finite_certificate": {
            "candidate_profiles_tested": len(rows),
            "nonzero_fourier_coefficients": sum(row["nonzero_coefficient"] for row in rows),
            "unit_mode_profiles": len(unit_mode_rows),
            "nonunit_mode_profiles": len(nonunit_mode_rows),
            "chart_relative_full_pair_recoveries": sum(row["chart_relative_full_pair_recovery"] for row in rows),
            "aliased_nonunit_recoveries": sum(row["matching_pair_count"] > 1 for row in rows),
            "translation_phase_rotation_witnesses": sum(row["translation_phase_rotates_at_t1"] for row in rows),
            "accepted_import_free_J_DHL_sources": 0,
            "mode_alias_counts": {str(k): sorted({row["matching_pair_count"] for row in rows if row["mode_k"] == k}) for k in MODES},
        },
        "representative_rows": representative_rows,
        "extraction_rows": rows,
        "decision": {
            "bounded_result": "P3136 constructs an actual formula-level J_DHL candidate rather than another matrix label. The Fourier phase extractor is mathematically strong: all 120 P3134 profiles have nonzero Fourier coefficient at their generating mode. It also exposes a sharper obstruction than expected: even the primitive unit modes k=1 and k=5 have a half-period ambiguity (r,lambda) ~ (r+6,-lambda), while nonunit modes have larger aliases. Thus no row recovers an import-free full pair. The extraction assumes a labelled Fourier frame, a chosen mode/character, and a polarity convention; translation rotates the coefficient phase in all 120 rows. P3136 therefore proves a positive receiver theorem and a stricter joint-source no-go boundary at the same time.",
            "positive_scoped_flags": {
                "formula_level_J_DHL_candidate_constructed": True,
                "nonzero_fourier_receiver_verified": True,
                "unit_mode_half_period_aliasing_verified": True,
                "nonunit_mode_aliasing_computed": True,
                "translation_phase_rotation_obstruction_verified": True,
            },
            "negative_export_flags": {
                "import_free_J_DHL_source_exported": False,
                "Fourier_frame_source_exported": False,
                "mode_character_source_exported": False,
                "polarity_convention_source_exported": False,
                "D_HL_source_exported": False,
                "Zeta_OS_exported": False,
                "Gamma_SO_exported": False,
                "QW_2191_discharged": False,
                "bridge_completion_exported": False,
                "legacy_role_transfer_exported": False,
                "L_total_exported": False,
                "ToE_closure_exported": False,
            },
            "next_honest_step": "Do not repeat Fourier receiver extraction. The next admissible proof-grade object is one strict Fourier-frame/source law F_DHL that selects a primitive mode/character and phase-zero reference without importing chart labels. Its audit should be a finite frame-source obstruction/witness test against prior Fourier-character no-go lanes P2992/P2994, the chi_i localizer boundary P3039, and the Aut-character idempotent warning P2869. Without such F_DHL, preserve the P3135-P3136 conditional-positive/no-strict-source certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["finite_certificate"]
    decision = payload["decision"]
    lines = [
        "# P3136/S2086 Fourier-phase J_DHL candidate audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Constructed object",
        f"`{payload['constructed_object']['formula']}`",
        "",
        "## Repo backscan",
    ]
    for hit in payload["repo_backscan_summary"]:
        lines.append(f"- {hit}")
    lines.extend([
        "",
        "## Finite certificate",
        f"- candidate profiles tested: `{cert['candidate_profiles_tested']}`",
        f"- nonzero Fourier coefficients: `{cert['nonzero_fourier_coefficients']}`",
        f"- unit-mode profiles: `{cert['unit_mode_profiles']}`",
        f"- nonunit-mode profiles: `{cert['nonunit_mode_profiles']}`",
        f"- chart-relative full pair recoveries: `{cert['chart_relative_full_pair_recoveries']}`",
        f"- aliased nonunit recoveries: `{cert['aliased_nonunit_recoveries']}`",
        f"- translation phase-rotation witnesses: `{cert['translation_phase_rotation_witnesses']}`",
        f"- accepted import-free J_DHL sources: `{cert['accepted_import_free_J_DHL_sources']}`",
        f"- mode alias counts: `{cert['mode_alias_counts']}`",
        "",
        "## Decision",
        decision["bounded_result"],
        "",
        "## Recommendation",
        decision["next_honest_step"],
    ])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3136/S2086 Fourier-phase J_DHL candidate audit", "## P3136/S2086 Fourier-phase J_DHL candidate audit\n\n`P3136/S2086` constructs the formula-level joint source candidate requested after P3135: the Fourier phase extractor `C_k(D)=sum_x D(x)*exp(-2*pi*i*k*x/12)` for the P3134 sine residual.  The finite audit tests all `120` profiles.  Every generating-mode coefficient is nonzero; even the `48` primitive/unit-mode rows (`k=1,5`) retain the half-period ambiguity `(r,lambda) ~ (r+6,-lambda)`, while nonunit modes have larger aliases.  Translation rotates the coefficient phase in all `120` rows, so the extractor assumes a labelled Fourier frame, selected mode/character, and polarity convention.  No import-free `J_DHL`, `D_HL` source, `Zeta_OS`, `Gamma_SO`, selector closure, bridge completion, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3136/S2086 Fourier-phase extractor is not yet a variational source", "## P3136/S2086 Fourier-phase extractor is not yet a variational source\n\n`P3136/S2086` gives a conditional positive reconstruction theorem for `D_HL` profiles through Fourier phase, but the reconstruction depends on an imported Fourier frame, mode/character, and polarity convention.  Therefore it does not export a Lagrangian density, Hamiltonian normalization, spacetime EOM, physical unit, `L_total`, bridge-completion theorem, or role-transfer theorem.\n")
    append_once(AGENTS, "Current Fourier-phase J_DHL candidate guardrail (P3136/S2086, 2026-07-12)", "## Current Fourier-phase J_DHL candidate guardrail (P3136/S2086, 2026-07-12)\n\n- P3136 constructs a formula-level `J_DHL` candidate using the Fourier coefficient phase of the P3134 sine residual.\n- The finite audit verifies a strong conditional result: all `120` profiles have nonzero generating-mode coefficients, but the `48` primitive mode rows `k=1,5` still retain the half-period ambiguity `(r,lambda) ~ (r+6,-lambda)`; `72` nonunit rows have larger aliasing.\n- The same audit proves the strict-source boundary: translation rotates coefficient phase in all rows, and the extraction imports a labelled Fourier frame, selected mode/character, and polarity convention.\n- Do not replay Fourier receiver extraction as `D_HL`, `Zeta_OS`, `Gamma_SO`, `QW-2191` discharge, bridge completion, role transfer, `L_total`, or ToE closure. The next admissible object is exactly one strict Fourier-frame/source law `F_DHL` selecting primitive character and phase-zero reference, audited against P2992/P2994, P3039, and P2869.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
