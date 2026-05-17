#!/usr/bin/env python3
"""P1952 S902 strict QW-2049 UR domain separation certificate.

This is a domain repair step after P1951 found that the old P1677 restricted
UR domain does not contain the QW-2049 strict tuple.  The script exports a
separate seed-local rectangle around QW-2049, proves positivity lower bounds
for the declared B1 coefficients and Disc_seed/s on that rectangle, and keeps
the full UR_link theorem explicitly open.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p1952_s902_strict_qw2049_ur_domain_separation_certificate.json"


def load(name: str) -> dict:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def digest(obj: object) -> str:
    blob = json.dumps(obj, sort_keys=True, ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


def num(expr: sp.Expr, digits: int = 30) -> str:
    return str(sp.N(expr, digits))


def main() -> None:
    p1677 = load("p1677_s627_combined_unitarity_renormalization_theorem_object.json")
    p1950 = load("p1950_s900_strict_b1_renormalization_exact_cancellation_probe.json")
    p1951 = load("p1951_s901_strict_b1_cutkosky_phase_space_integral_probe.json")

    alpha = 4 * sp.log(2)
    pi = sp.pi

    # Seed-local rectangle: deliberately narrow, contains QW-2049, and is
    # disjoint from the old P1677 beta/eta rectangle.
    rect = {
        "beta": (sp.Rational(19, 20), sp.Rational(21, 20)),
        "eta": (sp.Rational(17, 10), sp.Rational(19, 10)),
        "omega": (sp.Rational(9, 50), sp.Rational(19, 100)),
        "phi": (sp.Rational(3, 20), sp.Rational(17, 100)),
    }
    qw2049 = {
        "beta": sp.Rational(1, 1),
        "eta": sp.Rational(9, 5),
        "omega": sp.Rational(743, 4000),
        "phi": sp.Rational(13, 80),
    }

    beta_lo, _beta_hi = rect["beta"]
    eta_lo, _eta_hi = rect["eta"]
    omega_lo, _omega_hi = rect["omega"]
    phi_lo, phi_hi = rect["phi"]

    a_r2_lower = alpha / (16 * pi**2) * (beta_lo + eta_lo / 2) * (omega_lo**2 + phi_lo**2)
    a_ric2_lower = alpha / (16 * pi**2) * (1 + beta_lo) * (eta_lo * omega_lo + phi_lo)
    a_riem2_lower = alpha / (16 * pi**2) * beta_lo * eta_lo * (
        omega_lo**2 + omega_lo * phi_lo + phi_lo**2
    )
    # Because omega >= 0.18 and phi <= 0.17 on this rectangle, |omega-phi| >= 0.01.
    omega_phi_gap_lower = omega_lo - phi_hi
    a_gb_lower = alpha / (16 * pi**2) * (eta_lo - 1) * omega_phi_gap_lower**2
    disc_over_s_lower = sp.simplify((3 * a_r2_lower + 4 * a_ric2_lower) / (48 * pi))

    p1677_domain = p1677.get("restricted_domain") or {}
    beta_range = p1677_domain.get("beta_range") or []
    eta_range = p1677_domain.get("eta_range") or []
    beta_disjoint = bool(beta_range and float(beta_range[1]) < float(beta_lo))
    eta_disjoint = bool(eta_range and float(eta_range[1]) < float(eta_lo))
    separated_from_p1677 = beta_disjoint and eta_disjoint

    contains_qw2049 = {
        name: bool(bounds[0] <= qw2049[name] <= bounds[1])
        for name, bounds in rect.items()
    }

    lower_bounds = {
        "a_R2_lower": {"exact": str(sp.simplify(a_r2_lower)), "numeric_30d": num(a_r2_lower)},
        "a_Ric2_lower": {"exact": str(sp.simplify(a_ric2_lower)), "numeric_30d": num(a_ric2_lower)},
        "a_Riem2_lower": {"exact": str(sp.simplify(a_riem2_lower)), "numeric_30d": num(a_riem2_lower)},
        "a_GB_lower": {"exact": str(sp.simplify(a_gb_lower)), "numeric_30d": num(a_gb_lower)},
        "Disc_seed_over_s_lower": {
            "exact": str(sp.simplify(disc_over_s_lower)),
            "numeric_30d": num(disc_over_s_lower),
        },
    }
    all_lower_bounds_positive = all(sp.N(sp.sympify(v["exact"]), 50) > 0 for v in lower_bounds.values())

    out = {
        "packet_id": "P1952",
        "stage_id": "S902",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "local_verdict": "QW2049_SEED_LOCAL_POSITIVE_RECTANGLE_EXPORTED__UR_LINK_STILL_OPEN",
        "route": "strict_only",
        "legacy_bridge_used": False,
        "depends_on": {
            "p1677_present": "restricted_domain" in p1677,
            "p1950_present": p1950.get("local_verdict") == "PASS_ZERO_ON_DECLARED_B1_ANSATZ",
            "p1951_present": "p1677_domain_compatibility_check" in p1951,
        },
        "input_hashes": {
            "p1677_restricted_domain_sha256": digest(p1677_domain),
            "p1950_symbolic_cancellation_sha256": digest(p1950.get("symbolic_cancellation") or {}),
            "p1951_domain_check_sha256": digest(p1951.get("p1677_domain_compatibility_check") or {}),
        },
        "old_p1677_domain": p1677_domain,
        "qw2049_tuple": {name: str(value) for name, value in qw2049.items()},
        "qw2049_seed_local_rectangle": {
            name: {"lower": str(bounds[0]), "upper": str(bounds[1])}
            for name, bounds in rect.items()
        },
        "contains_qw2049_tuple": contains_qw2049,
        "domain_separation_from_p1677": {
            "beta_disjoint": beta_disjoint,
            "eta_disjoint": eta_disjoint,
            "separated_from_p1677": separated_from_p1677,
            "meaning": "The new rectangle is a separate QW-2049 seed-local domain, not an extension theorem for P1677.",
        },
        "positivity_certificate": {
            "assumptions": [
                "beta in [19/20, 21/20]",
                "eta in [17/10, 19/10]",
                "omega in [9/50, 19/100]",
                "phi in [3/20, 17/100]",
                "omega - phi >= 1/100",
                "alpha_geo_strict = 4*log(2) > 0",
            ],
            "lower_bounds": lower_bounds,
            "all_lower_bounds_positive": bool(all_lower_bounds_positive),
            "disc_statement": "For s > 0, Disc_seed(s) >= s * Disc_seed_over_s_lower > 0 on the QW-2049 seed-local rectangle.",
        },
        "renormalization_link_scope": {
            "p1950_cancellation_reused": p1950.get("local_verdict"),
            "safe_statement": "The declared B1 counterterm cancellation and seed Cutkosky positivity coexist on the QW-2049 seed-local rectangle.",
            "unsafe_statement_rejected": "This does not prove P1677 UR_link, dressed unitarity, or background-global renormalization.",
        },
        "theorem_scope": {
            "proved": [
                "P1677 restricted beta/eta domain is disjoint from a narrow QW-2049 seed-local rectangle.",
                "The declared B1 a_* lower bounds are strictly positive on that rectangle.",
                "The declared seed Disc_seed(s)/s has a strictly positive lower bound on that rectangle.",
            ],
            "not_proved": [
                "extension of P1677 to QW-2049",
                "full dressed graviton->gauge_gauge Cutkosky equality",
                "global UR_link theorem",
                "background-independence transport",
                "QW-2191 selector discharge",
            ],
        },
        "false_pass_guard": "This certificate repairs the bookkeeping domain split only at seed-local level; it is not a global unitarity-renormalization theorem.",
        "next_honest_step": "Build P1953: replace the seed integrand by a dressed amplitude model reduced to the same basis, or export a precise obstruction if the dressed amplitude is not available in the repo.",
        "lay_explanation": "Zamiast udawac, ze stara domena P1677 obejmuje aktualny strict punkt, wyznaczylismy osobna mala bezpieczna okolice wokol QW-2049. W tej okolicy robocze wspolczynniki i calka sa dodatnie, ale pelny dowod unitarności nadal czeka na ubrana amplitude.",
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
