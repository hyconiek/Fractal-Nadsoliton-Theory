#!/usr/bin/env python3
"""P1891 S841 strict first-gate promotion contract probe."""
from __future__ import annotations
import json
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def load(name: str) -> dict:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def main() -> None:
    p1890 = load("p1890_s840_strict_qg_gate_closure_scoreboard_probe.json")

    beta, eta, omega, phi = sp.symbols("beta eta omega phi", positive=True, real=True)
    lam, y, kappa2 = sp.symbols("lam y kappa2", positive=True, real=True)
    eps, mu2, s = sp.symbols("eps mu2 s", positive=True, real=True)

    c0 = sp.cos(omega + phi) / (beta + 1)
    c1 = (-beta * eta * sp.cos(omega + phi) - omega * (beta + 1) * sp.sin(omega + phi)) / (beta + 1) ** 2
    lam_eff = lam * (1 + c1**2)
    y_eff = y * (1 + c0 / 2)

    G_R = lam_eff / (16 * sp.pi**2) * (1 / eps + sp.log(mu2 / s))
    G_U = y_eff**2 / (16 * sp.pi**2) * (1 / eps + sp.log(mu2 / s) + sp.Rational(1, 2))
    G_BI = (y_eff**2 + kappa2) / (16 * sp.pi**2) * (1 / eps + sp.log(mu2 / s) + sp.Rational(1, 2))

    promotion_contract = {
        "gate_candidate": "renormalization_gate",
        "required_evidence": [
            "diagram_resolved_coefficients_for_G_R",
            "finite_part_scheme_lock",
            "cross-check_against_same_EOM_registry",
        ],
        "promotion_condition": "G_R witness closed while preserving OPEN status of U/BI/selector gates",
        "no_false_pass_rule": "single-gate promotion must not be misread as strict-core closure",
    }

    qw2049 = {beta: 1.0, eta: 1.8, omega: 0.18575, phi: 0.16250}

    out = {
        "packet_id": "P1891",
        "stage_id": "S841",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {"p1890_present": "qg_gate_scoreboard" in p1890},
        "strict_chain_step": "K_strict -> coefficients -> full L_SM+L_GR -> EOM -> gate scoreboard -> first-gate promotion contract",
        "effective_coefficients": {
            "lam_eff": str(lam_eff),
            "y_eff": str(y_eff),
        },
        "gate_symbols": {
            "G_R": str(sp.expand(G_R)),
            "G_U": str(sp.expand(G_U)),
            "G_BI": str(sp.expand(G_BI)),
        },
        "first_gate_promotion_contract": promotion_contract,
        "qw2049_trace": {
            "lam_eff_over_lam": str(sp.N((lam_eff / lam).subs(qw2049), 12)),
            "y_eff_over_y": str(sp.N((y_eff / y).subs(qw2049), 12)),
        },
        "strict_core_closure_missing_items": {
            "renormalization": "Gate-promotion requires explicit solved G_R witness package.",
            "unitarity": "G_U still OPEN; needs Cutkosky solved inputs.",
            "background_independence": "G_BI still OPEN; needs transport theorem/witness.",
            "selector_obstruction": "QW-2191 remains open.",
        },
        "false_pass_guard": "First-gate promotion contract does not imply full strict-core closure.",
        "next_honest_step": "Deliver explicit G_R witness packet and update scoreboard with one promoted gate while preserving open blockers elsewhere.",
        "lay_explanation": "To strategia krokowa: najpierw domknąć jedną bramę (renormalizację) bardzo rygorystycznie, zamiast udawać od razu pełne domknięcie całej teorii.",
    }

    path = GEN / "p1891_s841_strict_first_gate_promotion_contract_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
