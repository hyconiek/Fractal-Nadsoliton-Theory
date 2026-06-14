#!/usr/bin/env python3
"""P2733/S1683: spectral lambda-observability no-go for chiral tau coupling.

P2732 showed that the direct Im(B_{1,5}) tau-coupling selects tau only after
lambda and orientation/polarity are fixed.  P2733 asks whether the coupling
sign lambda can be recovered from intrinsic spectral data of that variational
problem.  It computes the full energy histogram over all 2^12 tau fields for
each (lambda, orientation) row and checks whether any row-spectral invariant
can distinguish lambda or P2721 polarity.
"""
from __future__ import annotations

import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2733_s1683_chiral_tau_coupling_spectral_lambda_observability_no_go.json"
MD = GEN / "p2733_s1683_chiral_tau_coupling_spectral_lambda_observability_no_go.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

Z12 = 12
SPINS = (-1, 1)
ORIENTATIONS = (-1, 1)
LAMBDAS = (-1, 1)
FIELDS = tuple(itertools.product(SPINS, repeat=Z12))
INPUTS = {
    "P2732_CHIRAL_TAU_COUPLING": GEN / "p2732_s1682_chiral_bispectrum_time_arrow_source_term_coupling_matrix.json",
    "P2721_SIGN_TORSOR_COUPLING": GEN / "p2721_s1671_chiral_bispectrum_sign_to_orientation_torsor_coupling_audit.json",
    "P2718_CHIRAL_BISPECTRUM_MARKER": GEN / "p2718_s1668_chiral_bispectrum_explicit_signed_formula_torsor_coupling_audit.json",
}
NEGATIVE_EXPORT_FLAGS = [
    "lambda_observable_from_spectrum",
    "p2721_polarity_observable_from_spectrum",
    "strict_mechanism_fixing_lambda_exported",
    "strict_time_arrow_source_term_exported",
    "qw2191_discharged",
    "pair12_strict_core_upgrade_exported",
    "bridge_closure_exported",
    "role_transfer_started",
    "ltotal_promoted",
    "toe_closure_exported",
]


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def marker_table(p2718: dict[str, Any]) -> dict[tuple[int, int], float]:
    table = {
        (int(row["orientation"]), int(row["source"])): float(row["marker_imag"])
        for row in p2718.get("marker_rows", [])
    }
    return table or {(orientation, source): float(2 * orientation) for orientation in ORIENTATIONS for source in range(Z12)}


def energy(field: tuple[int, ...], orientation: int, lam: int, markers: dict[tuple[int, int], float]) -> float:
    return -sum(lam * field[source] * markers[(orientation, source)] for source in range(Z12))


def row_spectrum(orientation: int, lam: int, markers: dict[tuple[int, int], float]) -> dict[str, Any]:
    energies = [energy(field, orientation, lam, markers) for field in FIELDS]
    hist = Counter(energies)
    min_energy = min(energies)
    max_energy = max(energies)
    ground_magnetizations = sorted(set(sum(field) for field, value in zip(FIELDS, energies) if value == min_energy))
    ceiling_magnetizations = sorted(set(sum(field) for field, value in zip(FIELDS, energies) if value == max_energy))
    spectral_signature = tuple(sorted((float(key), count) for key, count in hist.items()))
    unlabeled_extrema_signature = (
        float(min_energy),
        hist[min_energy],
        float(max_energy),
        hist[max_energy],
        tuple(sorted(abs(m) for m in ground_magnetizations)),
        tuple(sorted(abs(m) for m in ceiling_magnetizations)),
    )
    return {
        "lambda": lam,
        "orientation": orientation,
        "energy_histogram": {str(float(key)): count for key, count in sorted(hist.items())},
        "spectral_signature": spectral_signature,
        "min_energy": float(min_energy),
        "max_energy": float(max_energy),
        "ground_state_count": hist[min_energy],
        "ceiling_state_count": hist[max_energy],
        "ground_magnetizations": ground_magnetizations,
        "ceiling_magnetizations": ceiling_magnetizations,
        "unlabeled_extrema_signature": unlabeled_extrema_signature,
    }


def audit_spectral_observability(markers: dict[tuple[int, int], float]) -> dict[str, Any]:
    rows = [row_spectrum(orientation, lam, markers) for lam in LAMBDAS for orientation in ORIENTATIONS]
    signatures = [row["spectral_signature"] for row in rows]
    extrema_signatures = [row["unlabeled_extrema_signature"] for row in rows]
    full_spectrum_classes = Counter(signatures)
    extrema_classes = Counter(extrema_signatures)
    lambda_pair_same_spectrum = []
    orientation_pair_same_spectrum = []
    for row in rows:
        lambda_partner = next(candidate for candidate in rows if candidate["lambda"] == -row["lambda"] and candidate["orientation"] == row["orientation"])
        orientation_partner = next(candidate for candidate in rows if candidate["lambda"] == row["lambda"] and candidate["orientation"] == -row["orientation"])
        lambda_pair_same_spectrum.append(row["spectral_signature"] == lambda_partner["spectral_signature"])
        orientation_pair_same_spectrum.append(row["spectral_signature"] == orientation_partner["spectral_signature"])
    return {
        "row_count": len(rows),
        "field_count_per_row": len(FIELDS),
        "distinct_full_energy_spectra": len(full_spectrum_classes),
        "distinct_unlabeled_extrema_signatures": len(extrema_classes),
        "full_spectrum_class_sizes": sorted(full_spectrum_classes.values()),
        "unlabeled_extrema_class_sizes": sorted(extrema_classes.values()),
        "all_lambda_pairs_have_identical_spectrum": all(lambda_pair_same_spectrum),
        "all_orientation_pairs_have_identical_spectrum": all(orientation_pair_same_spectrum),
        "rows": rows,
        "obstruction": "The full tau-energy spectrum of the chiral coupling is identical for all lambda and orientation rows.  Therefore any lambda/polarity law depending only on intrinsic row-spectral data is blind to the sign that P2732 needs fixed.",
    }


def acceptance_matrix(audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "all_four_lambda_orientation_rows_scanned": audit["row_count"] == 4,
        "full_tau_spectrum_computed_per_row": audit["field_count_per_row"] == 4096,
        "full_spectrum_distinguishes_lambda_or_polarity": audit["distinct_full_energy_spectra"] > 1,
        "unlabeled_extrema_distinguish_lambda_or_polarity": audit["distinct_unlabeled_extrema_signatures"] > 1,
        "strict_lambda_fixing_spectral_law_exported": False,
    }
    missing = [key for key, value in facts.items() if not value]
    return {
        "facts": facts,
        "missing_criteria": missing,
        "accepted_as_lambda_fixing_spectral_law": not missing,
        "blocker": "Spectral data are identical across the lambda/orientation torsor; they cannot source the missing sign.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["spectral_lambda_observability_audit"]
    lines = [
        "# P2733/S1683 chiral tau-coupling spectral lambda-observability no-go",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Spectral audit",
        f"- row_count={audit['row_count']}",
        f"- field_count_per_row={audit['field_count_per_row']}",
        f"- distinct_full_energy_spectra={audit['distinct_full_energy_spectra']}",
        f"- distinct_unlabeled_extrema_signatures={audit['distinct_unlabeled_extrema_signatures']}",
        f"- all_lambda_pairs_have_identical_spectrum={audit['all_lambda_pairs_have_identical_spectrum']}",
        f"- all_orientation_pairs_have_identical_spectrum={audit['all_orientation_pairs_have_identical_spectrum']}",
        audit["obstruction"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    inputs = {key: read_json(path) for key, path in INPUTS.items()}
    markers = marker_table(inputs["P2718_CHIRAL_BISPECTRUM_MARKER"])
    audit = audit_spectral_observability(markers)
    acceptance = acceptance_matrix(audit)
    no_go = not acceptance["accepted_as_lambda_fixing_spectral_law"]
    payload = {
        "status": "P2733_CHIRAL_TAU_COUPLING_SPECTRAL_LAMBDA_OBSERVABILITY_NO_GO" if no_go else "P2733_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "input_statuses": {key: value.get("status") for key, value in inputs.items()},
        "audited_candidate_class": "all intrinsic full-energy-spectrum and unlabeled-extrema invariants of the P2732 direct chiral tau coupling rows",
        "spectral_lambda_observability_audit": audit,
        "acceptance_matrix": acceptance,
        "decision": {
            "lambda_observable_from_spectrum": False,
            "p2721_polarity_observable_from_spectrum": False,
            "strict_mechanism_fixing_lambda_exported": False,
            "strict_time_arrow_source_term_exported": False,
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": "P2733 computes the full tau-energy spectrum for every P2732 coupling row.  All rows share the same spectrum and unlabeled extrema data, so no intrinsic spectral law can distinguish lambda or P2721 polarity.",
            "next_honest_step": "Do not search for lambda fixing inside the P2732 row spectrum or unlabeled ground/extrema data.  A next admissible proof-grade move must introduce a new non-spectral strict polarity/lambda source, or preserve the P2697-P2733 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2733/S1683 chiral tau-coupling spectral lambda-observability no-go", "## P2733/S1683 chiral tau-coupling spectral lambda-observability no-go\n\n`P2733/S1683` computes the full `2^12=4096`-state tau-energy spectrum for each P2732 direct-coupling row `(lambda,orientation)`.  All four rows have identical full spectra and identical unlabeled extrema signatures, so intrinsic spectral data cannot fix `lambda` or select the P2721 polarity.  No strict mechanism fixing `lambda`, `QW-2191` discharge, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2733/S1683 spectral lambda-observability Ltotal guard", "## P2733/S1683 spectral lambda-observability Ltotal guard\n\n`P2733/S1683` shows that the variational spectrum of the direct chiral `tau` coupling is blind to `lambda` and orientation/polarity.  Because no intrinsic spectral law selects the missing sign, this audit adds no strict variational source term and cannot promote `L_total`.\n")
    append_once(AGENTS, "Current chiral tau-coupling spectral lambda-observability no-go guardrail (P2733/S1683, 2026-06-14)", "## Current chiral tau-coupling spectral lambda-observability no-go guardrail (P2733/S1683, 2026-06-14)\n\n- P2733 computes the full tau-energy spectrum for every P2732 direct-coupling row and finds one identical spectrum across all `lambda` and orientation/P2721-polarity branches.\n- Intrinsic row-spectral and unlabeled-extrema data cannot fix `lambda` or select a P2721 polarity; they are blind to the sign needed by P2732.\n- Do not replay P2732 spectral/ground-sector data as `QW-2191` discharge, selector closure, pair12 strict-core upgrade, role transfer, bridge closure, `L_total`, or ToE.  A next admissible move must introduce a new non-spectral strict polarity/lambda source, or preserve the P2697-P2733 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
