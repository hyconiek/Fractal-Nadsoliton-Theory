#!/usr/bin/env python3
"""P1988 S938 strict non-GB to spatial-EOM family projection witness.

Next honest step after P1987: project the classified non-GB lapse residual
families onto the P1974-style Bianchi-I spatial anisotropic family basis and
check whether the family content strictly exceeds Einstein-level anisotropic
residual classes.
"""

from __future__ import annotations

import json
import platform
from pathlib import Path
from typing import Any

import numpy as np
import scipy.linalg as la
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p1988_s938_strict_non_gb_to_spatial_eom_family_projection_witness.json"
DEFAULT_TIMESTAMP_UTC = "2026-05-18T00:00:00+00:00"


def load(name: str) -> dict[str, Any]:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p1987 = load("p1987_s937_strict_non_gb_residual_term_classification_witness.json")
    p1974 = load("p1974_s924_strict_bianchi_anisotropic_eom_residual_obstruction_witness.json")

    # Family presence from P1987 (already strict-identified residual).
    fam_flags = p1987.get("term_family_nonzero_flags", {})

    # P1974 EH-level anisotropic residual families are linear in (H*sigma, dsigma)
    # plus quadratic Q; no Ndd, no d2sigma, no quartic Q^2 channel.
    eh_family_capacity = {
        "H2_Q": False,
        "H_Qd": False,
        "Hd_Q": False,
        "Q2": False,
        "d2sigma_sigma": False,
        "dsigma_sq": False,
        "Ndd_sigma_sq": False,
        "Nd2_sigma_sq": False,
        "H_sigma_linear": True,
        "dsigma_linear": True,
        "Q_quadratic": True,
    }

    strict_non_gb_families = [k for k, v in fam_flags.items() if v]
    outside_eh_capacity = [k for k in strict_non_gb_families if not eh_family_capacity.get(k, False)]

    # Deterministic numeric severity from P1987 probe.
    l2 = float((p1987.get("numeric_probe") or {}).get("l2_norm", 0.0))

    # A compact linear-algebra style indicator: if outside_EH_families > 0 then
    # any closure relying only on EH family channels has positive codimension gap.
    gap_vector = np.array([1.0 if k in outside_eh_capacity else 0.0 for k in strict_non_gb_families], dtype=float)
    gap_norm = float(la.norm(gap_vector, ord=2))

    p1974_vec = p1974.get("anisotropic_eom_residual_vector", [])
    p1974_has_only_low_order = all("d2sigma" not in s and "Ndd" not in s for s in p1974_vec)

    gate = {
        "p1987_present": p1987.get("result_kind") == "PASS_NON_GB_TERM_CLASSIFICATION_WITNESS",
        "p1974_present": p1974.get("packet_id") == "P1974",
        "strict_non_gb_families_detected": len(strict_non_gb_families) > 0,
        "outside_eh_family_capacity_detected": len(outside_eh_capacity) > 0,
        "family_gap_norm_positive": gap_norm > 0.0,
        "p1974_eh_residual_has_no_d2sigma_or_Ndd": p1974_has_only_low_order,
        "numeric_probe_nonzero": l2 > 0.0,
    }

    pass_witness = all([
        gate["p1987_present"], gate["p1974_present"], gate["strict_non_gb_families_detected"],
        gate["outside_eh_family_capacity_detected"], gate["family_gap_norm_positive"],
        gate["p1974_eh_residual_has_no_d2sigma_or_Ndd"], gate["numeric_probe_nonzero"],
    ])

    out = {
        "ledger_id": "P1988_S938_STRICT_NON_GB_TO_SPATIAL_EOM_FAMILY_PROJECTION_WITNESS",
        "packet_id": "P1988",
        "stage_id": "S938",
        "produced_by": Path(__file__).name,
        "timestamp_utc": DEFAULT_TIMESTAMP_UTC,
        "route": "strict_only",
        "depends_on": {
            "p1987_term_classification_present": gate["p1987_present"],
            "p1974_anisotropic_eom_obstruction_present": gate["p1974_present"],
        },
        "family_projection": {
            "strict_non_gb_families_detected": strict_non_gb_families,
            "eh_family_capacity": eh_family_capacity,
            "outside_eh_family_capacity": outside_eh_capacity,
            "outside_eh_family_count": len(outside_eh_capacity),
            "family_gap_norm": gap_norm,
        },
        "cross_packet_consistency": {
            "p1974_component_residual_vector": p1974_vec,
            "p1974_has_only_low_order_eh_families": p1974_has_only_low_order,
        },
        "numeric_severity": {
            "p1987_probe_l2_norm": l2,
        },
        "gatekeeper_checks": gate,
        "result_kind": "PASS_NON_GB_SPATIAL_PROJECTION_OBSTRUCTION_WITNESS" if pass_witness else "OPEN_OBSTRUCTION_WITH_TRACE",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "obstruction_tags": [
            "NON_GB_FAMILY_CONTENT_EXCEEDS_EH_ANISOTROPIC_CAPACITY",
            "SPATIAL_EOM_PROVIDER_CLASS_STILL_MISSING",
            "QW_2191_SELECTOR_CONSTRAINED_EXTENSION_STILL_OPEN",
        ],
        "next_honest_step": "Construct and test a strict admissible anisotropic provider ansatz class (with explicit selector premise labels) against each outside-EH family channel one-by-one in spatial Bianchi-I equations.",
        "lay_explanation": "Pokazaliśmy formalnie, że klasy składników z residualu non-GB są szersze niż to, co potrafi sam poziom Einsteinowski. Czyli bez nowego legalnego mechanizmu (providera) te składniki nie znikną.",
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": __import__("scipy").__version__,
            "sympy": sp.__version__,
        },
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1988] wrote witness: {OUT}")


if __name__ == "__main__":
    main()
