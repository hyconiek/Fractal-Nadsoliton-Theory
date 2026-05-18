#!/usr/bin/env python3
"""P1990 S940 strict extended provider-family probe witness.

Next honest step after P1989: test one extended anisotropic-provider ansatz
family against representative outside-EH channels, with explicit non-strict
selector-premise labeling per guardrails.

This packet does NOT claim strict closure; it reports solvability structure for
one extended class only.
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
OUT = GEN / "p1990_s940_strict_extended_provider_family_probe_witness.json"
DEFAULT_TIMESTAMP_UTC = "2026-05-18T00:00:00+00:00"


def load(name: str) -> dict[str, Any]:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p1987 = load("p1987_s937_strict_non_gb_residual_term_classification_witness.json")
    p1989 = load("p1989_s939_strict_minimal_provider_family_no_go_witness.json")

    s1, s2, sd1, sd2, sdd1, sdd2 = sp.symbols("sigma1 sigma2 dsigma1 dsigma2 d2sigma1 d2sigma2", real=True)
    N, Nd, Ndd, V, H, Hd = sp.symbols("N Nd Ndd V H Hd", real=True)

    locals_map = {
        "N": N, "Nd": Nd, "Ndd": Ndd, "V": V, "H": H, "Hd": Hd,
        "sigma1": s1, "sigma2": s2, "dsigma1": sd1, "dsigma2": sd2, "d2sigma1": sdd1, "d2sigma2": sdd2,
        "pi": sp.pi, "log": sp.log, "ln": sp.log,
    }
    residual = sp.sympify(p1987.get("scaled_expression_N6_over_V", "0"), locals=locals_map)

    # Extended provider family with explicit selector-premise marker:
    # P_ext = minimal EH-like + eta1*(d2sigma*sigma) + eta2*(Ndd*Q) + eta3*(Q^2)
    # NOTE: eta_i are marked as selector-premise-augmented (non-strict) unless a
    # strict bridge is exported later.
    c1, c2, c3, eta1, eta2, eta3 = sp.symbols("c1 c2 c3 eta1 eta2 eta3", real=True)
    Q = s1**2 + s1 * s2 + s2**2
    Q2 = Q**2
    d2sig_sig = 2 * sdd1 * s1 + sdd1 * s2 + sdd2 * s1 + 2 * sdd2 * s2
    NddQ = Ndd * Q
    provider_ext = sp.expand(c1 * (3 * H * (s1 + s2)) + c2 * (sd1 + sd2) + c3 * Q + eta1 * d2sig_sig + eta2 * NddQ + eta3 * Q2)

    target = sp.expand(residual - provider_ext)

    # Representative outside-EH channels to cancel.
    ch = {
        "d2sigma1_sigma1": sp.expand(target).coeff(sdd1 * s1),
        "Ndd_sigma1_sq": sp.expand(target).coeff(Ndd * s1**2),
        "Qquartic_s1_4": sp.expand(target).coeff(s1**4),
    }

    # Solve only for eta1/eta2/eta3 to check if augmented channels can match these
    # representative outside-EH coefficients.
    eqs = [sp.Eq(ch["d2sigma1_sigma1"], 0), sp.Eq(ch["Ndd_sigma1_sq"], 0), sp.Eq(ch["Qquartic_s1_4"], 0)]
    sol_aug = sp.solve(eqs, (eta1, eta2, eta3), dict=True)
    has_aug_solution = len(sol_aug) > 0

    # Substitute first solution if exists; check residual remains in other channels.
    residual_after = target
    if has_aug_solution:
        residual_after = sp.expand(target.subs(sol_aug[0]))

    # Detect that remaining terms involving H*Qd/Hd*Q-like channels persist.
    mixed_flags = {
        "H_sigma1_dsigma1": bool(residual_after.coeff(H * s1 * sd1) != 0),
        "Hd_sigma1_sq": bool(residual_after.coeff(Hd * s1**2) != 0),
        "Nd2_sigma1_sq": bool(residual_after.coeff(Nd**2 * s1**2) != 0),
    }

    # Numeric replay.
    sample = {
        N: sp.Rational(1, 1), Nd: sp.Rational(1, 20), Ndd: sp.Rational(-1, 200),
        V: sp.Rational(1, 1), H: sp.Rational(1, 1), Hd: sp.Rational(1, 10),
        s1: sp.Rational(1, 10), s2: sp.Rational(-1, 20), sd1: sp.Rational(1, 100), sd2: sp.Rational(-1, 200),
        sdd1: sp.Rational(1, 1000), sdd2: sp.Rational(-1, 2000),
        c1: 0, c2: 0, c3: 0,
    }
    if has_aug_solution:
        sample.update({eta1: sp.nsimplify(sol_aug[0][eta1]), eta2: sp.nsimplify(sol_aug[0][eta2]), eta3: sp.nsimplify(sol_aug[0][eta3])})
    else:
        sample.update({eta1: 0, eta2: 0, eta3: 0})

    val = sp.simplify(residual_after.subs(sample))
    l2 = float(la.norm(np.array([float(sp.N(val, 40))], dtype=float), ord=2))

    gate = {
        "p1987_present": p1987.get("result_kind") == "PASS_NON_GB_TERM_CLASSIFICATION_WITNESS",
        "p1989_minimal_no_go_present": p1989.get("result_kind") == "PASS_MINIMAL_PROVIDER_NO_GO_WITNESS",
        "augmented_channel_solution_exists": has_aug_solution,
        "mixed_channels_still_persist": any(mixed_flags.values()),
        "numeric_residual_after_augmented_fit_nonzero": l2 > 0.0,
    }

    sol_aug_json = [{str(k): str(v) for k, v in d.items()} for d in sol_aug]

    out = {
        "ledger_id": "P1990_S940_STRICT_EXTENDED_PROVIDER_FAMILY_PROBE_WITNESS",
        "packet_id": "P1990",
        "stage_id": "S940",
        "produced_by": Path(__file__).name,
        "timestamp_utc": DEFAULT_TIMESTAMP_UTC,
        "route": "strict_only_with_explicit_non_strict_selector_labeling",
        "depends_on": {
            "p1987_present": gate["p1987_present"],
            "p1989_present": gate["p1989_minimal_no_go_present"],
        },
        "selector_premise_label": {
            "status": "NON_STRICT_AUGMENTED_CLASS",
            "reason": "eta1/eta2/eta3 extended channels are selector-premise-augmented and must not be claimed strict-core closure without a bridge.",
        },
        "extended_provider_family": {
            "definition": "P_ext = P_min + eta1*(d2sigma*sigma) + eta2*(Ndd*Q) + eta3*(Q^2)",
            "P_min": "c1*(3H(sigma1+sigma2)) + c2*(dsigma1+dsigma2) + c3*Q",
            "unknowns": ["c1", "c2", "c3", "eta1", "eta2", "eta3"],
        },
        "representative_outside_channel_equations": {k: str(v) for k, v in ch.items()},
        "augmented_channel_solver": {
            "equations": [str(e) for e in eqs],
            "solutions_for_eta": sol_aug_json,
            "has_solution": has_aug_solution,
        },
        "post_fit_residual_channel_flags": mixed_flags,
        "numeric_probe": {
            "sample": {str(k): str(v) for k, v in sample.items()},
            "residual_value_after_augmented_fit": str(val),
            "l2_norm": l2,
        },
        "gatekeeper_checks": gate,
        "result_kind": "PASS_EXTENDED_PROVIDER_PROBE_WITH_NONSTRICT_LABEL" if all(gate.values()) else "OPEN_OBSTRUCTION_WITH_TRACE",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "obstruction_tags": [
            "AUGMENTED_CLASS_CAN_FIT_SOME_OUTSIDE_CHANNELS",
            "MIXED_CHANNELS_STILL_PERSIST_AFTER_FIT",
            "STRICT_CORE_CLOSURE_STILL_BLOCKED_BY_SELECTOR_GUARDRAIL",
        ],
        "false_pass_guard": "P1990 does not prove strict-core closure. It probes one explicitly non-strict augmented class and still finds persistent mixed channels.",
        "next_honest_step": "Derive full spatial Bianchi-I component equations for the same augmented class and test componentwise cancellation map; keep verdict non-strict unless a strict selector bridge is exported.",
        "lay_explanation": "Nawet gdy dopuścimy bogatszy, jawnie oznaczony (niestety jeszcze nierygorystycznie strict) mechanizm, da się skasować tylko część problemu. Inne mieszane składniki nadal zostają.",
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": __import__("scipy").__version__,
            "sympy": sp.__version__,
        },
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1990] wrote witness: {OUT}")


if __name__ == "__main__":
    main()
