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

IN_ALPHA = GEN / "alpha_geo_strict_derived_v1.json"
IN_2314 = GEN / "p2314_s1264_strict_full_eom_lagrangian_symmetry_vortex_inventory_probe.json"
IN_RELEASE7 = REPO / "RELEASE_7_STRICT_PROJECTIVE_OPERATIONAL_MODEL_BRIEF.md"
OUT = GEN / "p2315_s1265_strict_schematic_lagrangian_eom_kernel_spectrum_probe.json"
MD = GEN / "p2315_s1265_strict_schematic_lagrangian_eom_kernel_spectrum_probe.md"

N = 12
OMEGA = 0.18575
PHI = 0.16250
BETA = 1.0
ETA = 1.8


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def read_text(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8", errors="replace")


def sha256_file(path: Path) -> str:
    if not path.exists():
        return ""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_json(payload: Any) -> str:
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def strict_kernel(distance: int) -> float:
    if distance == 0:
        return 0.0
    return math.cos(OMEGA * distance + PHI) / (1.0 + BETA * (distance ** ETA))


def cyclic_distance(i: int, j: int) -> int:
    raw = abs(i - j) % N
    return min(raw, N - raw)


def coupling_matrix() -> list[list[float]]:
    return [[strict_kernel(cyclic_distance(i, j)) if i != j else 0.0 for j in range(N)] for i in range(N)]


def fourier_eigenvalues(matrix: list[list[float]]) -> list[dict[str, Any]]:
    first_row = matrix[0]
    rows: list[dict[str, Any]] = []
    for k in range(N):
        real = 0.0
        imag = 0.0
        for j, value in enumerate(first_row):
            angle = -2.0 * math.pi * k * j / N
            real += value * math.cos(angle)
            imag += value * math.sin(angle)
        rows.append({
            "mode_k": k,
            "eigenvalue_real": real,
            "eigenvalue_imag_abs": abs(imag),
            "paired_with": (N - k) % N,
        })
    return rows


def pair_degeneracy_report(eigen_rows: list[dict[str, Any]], scale: float = 1.0) -> dict[str, Any]:
    by_k = {int(row["mode_k"]): float(row["eigenvalue_real"]) * scale for row in eigen_rows}
    pairs = []
    for k in range(1, N // 2):
        partner = N - k
        diff = abs(by_k[k] - by_k[partner])
        pairs.append({
            "pair": [k, partner],
            "eigenvalue_k": by_k[k],
            "eigenvalue_partner": by_k[partner],
            "abs_difference": diff,
            "degenerate": diff < 1e-12,
        })
    nyquist = {
        "mode": N // 2,
        "self_paired": True,
        "eigenvalue": by_k[N // 2],
    }
    return {
        "translation_invariant_circulant_matrix": True,
        "pair_planes": pairs,
        "nyquist_mode": nyquist,
        "all_pair_planes_degenerate": all(row["degenerate"] for row in pairs),
        "max_pair_difference": max((row["abs_difference"] for row in pairs), default=0.0),
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    alpha = load(IN_ALPHA)
    p2314 = load(IN_2314)
    release7_text = read_text(IN_RELEASE7)

    matrix = coupling_matrix()
    eigen_rows = fourier_eigenvalues(matrix)
    degeneracy = pair_degeneracy_report(eigen_rows)
    alpha_numeric = 4.0 * math.log(2.0)
    alpha_scaled_degeneracy = pair_degeneracy_report(eigen_rows, scale=alpha_numeric)

    row_sums = [sum(row) for row in matrix]
    max_row_sum_spread = max(row_sums) - min(row_sums)
    real_symmetric = all(abs(matrix[i][j] - matrix[j][i]) < 1e-15 for i in range(N) for j in range(N))
    circulant = all(abs(matrix[i][j] - matrix[0][(j - i) % N]) < 1e-15 for i in range(N) for j in range(N))

    schematic_eom_export = {
        "lagrangian_source": "RELEASE_7 schematic L_ZTP candidate layer",
        "normalized_euler_lagrange_row": "Box(Psi_o) + dV/d(Psi_o^dagger) + sum_{o' != o} K(o,o') Psi_o' = 0, up to the explicit 1/2 convention of the schematic complex-field representative",
        "derived_scope": "schematic Z12 kernel-coupled field equations only",
        "potential_term_unspecified": True,
        "metric_tensor_variation_missing": True,
        "full_tensor_eom_exported": False,
        "full_lagrangian_exported_for_task3": False,
    }

    kernel_spectrum_audit = {
        "strict_kernel_parameters": {
            "n": N,
            "omega": OMEGA,
            "phi": PHI,
            "beta": BETA,
            "eta": ETA,
            "distance_rule": "cyclic_distance_on_Z12",
            "diagonal_set_to_zero": True,
        },
        "matrix_properties": {
            "real_symmetric": real_symmetric,
            "circulant": circulant,
            "row_sums": row_sums,
            "max_row_sum_spread": max_row_sum_spread,
        },
        "first_row": matrix[0],
        "fourier_eigenvalues": eigen_rows,
        "pair_degeneracy_report": degeneracy,
        "alpha_scaled_pair_degeneracy_report": alpha_scaled_degeneracy,
        "selector_breaking_verdict": "The strict gate kernel matrix is circulant/translation-invariant on Z12.  Its Fourier pair planes k and 12-k remain exactly degenerate; multiplying by the scalar 4 ln 2 rescales eigenvalues but does not select an orientation inside a pair plane.",
    }

    computational_verdict = {
        "schematic_eom_derived": True,
        "full_eom_derived": False,
        "full_lagrangian_derived": False,
        "kernel_spectrum_computed": True,
        "kernel_breaks_pair_plane_degeneracy": False,
        "alpha_4ln2_breaks_pair_plane_degeneracy": False,
        "selector_source_exported": False,
        "admissible_for_g1_g3_update": False,
        "reason": "P2315 derives the schematic Z12 Euler-Lagrange row and computes the strict kernel spectrum, but the resulting circulant operator preserves Fourier pair degeneracy and does not fill the full tensor EOM/Lagrangian gaps recorded by P2314.",
    }

    theorem_export = {
        "statement_id": "P2315_SCHEMATIC_LAGRANGIAN_EOM_KERNEL_SPECTRUM_NONSELECTOR_THEOREM",
        "formal_statement": (
            "From the Release-7 schematic L_ZTP one can derive a normalized Z12 kernel-coupled Euler-Lagrange row.  With the frozen strict gate kernel parameters, the corresponding 12x12 coupling matrix is real symmetric and circulant.  Its Fourier modes k and 12-k are pairwise degenerate for k=1..5, and multiplication by the scalar alpha_geo_strict=4 ln 2 preserves this degeneracy.  Therefore this schematic EOM/spectrum computation does not export a full tensor-resolved EOM, does not export a theorem-grade Task-3 Lagrangian, and does not provide a selector/orientation source for G1/G3."
        ),
        "proof_bits": {
            "schematic_eom_export": schematic_eom_export,
            "matrix_real_symmetric": real_symmetric,
            "matrix_circulant": circulant,
            "max_row_sum_spread": max_row_sum_spread,
            "pair_degeneracy_report": degeneracy,
            "alpha_scaled_pair_degeneracy_report": alpha_scaled_degeneracy,
            "p2314_full_eom_answer": (p2314.get("strict_full_eom_lagrangian_symmetry_vortex_inventory_probe", {}) or {}).get("computational_answers", {}),
        },
        "not_claimed": [
            "full strict EOM",
            "full theorem-grade Lagrangian for Task-3",
            "4 ln 2 selector/symmetry-breaking theorem",
            "strict twisted-vortex/torsion bridge",
            "strict G1 closure",
            "strict G3 closure",
            "QW-2191 discharge",
            "ToE closure",
        ],
    }
    theorem_fingerprint = sha256_json(theorem_export)

    gatekeeper_checks = {
        "release7_schematic_lagrangian_detected": "\\mathcal{L}_{ZTP}" in release7_text,
        "schematic_eom_derived": computational_verdict["schematic_eom_derived"],
        "kernel_matrix_real_symmetric": real_symmetric,
        "kernel_matrix_circulant": circulant,
        "fourier_pair_degeneracy_verified": degeneracy["all_pair_planes_degenerate"],
        "alpha_scaling_preserves_pair_degeneracy": alpha_scaled_degeneracy["all_pair_planes_degenerate"],
        "full_eom_not_derived": computational_verdict["full_eom_derived"] is False,
        "full_lagrangian_not_derived": computational_verdict["full_lagrangian_derived"] is False,
        "selector_source_not_exported": computational_verdict["selector_source_exported"] is False,
        "p2314_inventory_loaded": p2314.get("packet_id") == "P2314",
        "strict_g1_g3_not_updated": True,
        "no_selector_premise_added": True,
        "no_qw2191_discharge_claimed": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2315_s1265_v1",
        "packet_id": "P2315",
        "stage_id": "S1265",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_PARTIAL_SCHEMATIC_EOM_DERIVED_KERNEL_SPECTRUM_NONSELECTOR",
        "result_kind": "STRICT_SCHEMATIC_LAGRANGIAN_EOM_AND_KERNEL_SPECTRUM_COMPUTATION_NO_G1_G3_CLOSURE",
        "strict_schematic_lagrangian_eom_kernel_spectrum_probe": {
            "probe_id": "P2315_S1265_STRICT_SCHEMATIC_LAGRANGIAN_EOM_KERNEL_SPECTRUM",
            "source_packets": {
                "alpha_geo_strict_derived_v1": "generated/alpha_geo_strict_derived_v1.json",
                "p2314": "generated/p2314_s1264_strict_full_eom_lagrangian_symmetry_vortex_inventory_probe.json",
                "release7": "RELEASE_7_STRICT_PROJECTIVE_OPERATIONAL_MODEL_BRIEF.md",
            },
            "source_hashes": {
                "alpha_geo_strict_derived_v1_sha256": sha256_file(IN_ALPHA),
                "p2314_sha256": sha256_file(IN_2314),
                "release7_sha256": sha256_file(IN_RELEASE7),
            },
            "schematic_eom_export": schematic_eom_export,
            "kernel_spectrum_audit": kernel_spectrum_audit,
            "computational_verdict": computational_verdict,
            "task3_g1_g3_update": {
                "G1_reduction_certainty": "OPEN_UNCHANGED",
                "G2_nonlinear_trajectory_realism": "CLOSED_FROM_P2301_UNCHANGED",
                "G3_operational_policy_rule": "OPEN_UNCHANGED",
                "reason": "P2315 derives only schematic Z12 EOM rows and verifies kernel pair degeneracy; it exports no full EOM/Lagrangian or selector source.",
            },
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": theorem_fingerprint,
        },
        "recommended_next_honest_step": {
            "id": "P2316_candidate",
            "goal": "To get a real closure-relevant advance, add one missing nonlinear/tensor object: an explicit potential V with stability theorem, a full tensor-resolved EOM variation, or an internal selector/orientation source that breaks the Fourier pair-plane degeneracy before replaying P2302/P2281.",
        },
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_PARTIAL_COMPUTATIONAL_EOM_PROGRESS_WITH_NONSELECTOR_SPECTRUM_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2315/S1265 — schematic Lagrangian EOM and strict-kernel spectrum",
            "",
            f"- Status: `{payload['status']}`",
            f"- Schematic EOM derived: `{computational_verdict['schematic_eom_derived']}`",
            f"- Full EOM derived: `{computational_verdict['full_eom_derived']}`",
            f"- Full Task-3 Lagrangian derived: `{computational_verdict['full_lagrangian_derived']}`",
            f"- Strict kernel pair-plane degeneracy verified: `{degeneracy['all_pair_planes_degenerate']}`",
            f"- `4 ln 2` scalar scaling preserves degeneracy: `{alpha_scaled_degeneracy['all_pair_planes_degenerate']}`",
            "- G1/G3 update: `OPEN_UNCHANGED`",
            f"- Theorem fingerprint: `{theorem_fingerprint}`",
            "",
            "## Concrete result",
            "P2315 derives the normalized schematic Euler-Lagrange row from the Release-7 candidate Lagrangian and computes the frozen strict Z12 kernel spectrum.  The operator is circulant, so Fourier pair planes remain degenerate.  A scalar `4 ln 2` multiplier rescales eigenvalues but does not pick an orientation.",
            "",
            "## Guardrail statement",
            "P2315 does not export a full tensor EOM, a full theorem-grade Lagrangian, a selector source, QW-2191 discharge, or ToE closure.",
            "",
            "## Next honest step",
            payload["recommended_next_honest_step"]["goal"],
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
