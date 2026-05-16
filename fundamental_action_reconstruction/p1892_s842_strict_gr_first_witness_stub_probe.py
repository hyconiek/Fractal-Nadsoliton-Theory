#!/usr/bin/env python3
"""P1892 S842 strict renormalization-gate first witness stub probe."""
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
    p1891 = load("p1891_s841_strict_first_gate_promotion_contract_probe.json")

    beta, eta, omega, phi = sp.symbols("beta eta omega phi", positive=True, real=True)
    lam = sp.symbols("lam", positive=True, real=True)
    eps, mu2, s = sp.symbols("eps mu2 s", positive=True, real=True)

    c1 = (-beta * eta * sp.cos(omega + phi) - omega * (beta + 1) * sp.sin(omega + phi)) / (beta + 1) ** 2
    lam_eff = lam * (1 + c1**2)

    # First explicit renormalization witness stub with pole + finite part.
    A_div, A_fin = sp.symbols("A_div A_fin", real=True)
    loop_scalar4 = lam_eff**2 / (16 * sp.pi**2) * (A_div / eps + A_fin + sp.log(mu2 / s))
    ct_scalar4 = -lam_eff**2 / (16 * sp.pi**2) * (A_div / eps + A_fin)
    ren_scalar4 = sp.expand(loop_scalar4 + ct_scalar4)

    pass_metric = sp.simplify(ren_scalar4 - lam_eff**2 * sp.log(mu2 / s) / (16 * sp.pi**2))

    qw2049 = {beta: 1.0, eta: 1.8, omega: 0.18575, phi: 0.16250}

    out = {
        "packet_id": "P1892",
        "stage_id": "S842",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {"p1891_present": "first_gate_promotion_contract" in p1891},
        "strict_chain_step": "K_strict -> coefficients -> full L_total quartic sector -> first renormalization witness stub",
        "effective_coefficients": {
            "lam_eff": str(lam_eff),
        },
        "first_gate_witness_stub": {
            "loop_scalar4": str(loop_scalar4),
            "ct_scalar4": str(ct_scalar4),
            "ren_scalar4": str(ren_scalar4),
            "pass_metric_symbol": str(pass_metric),
            "interpretation": "Pole+finite cancellation contract for one representative term of full L_total.",
        },
        "promotion_gate_target": "G_R",
        "promotion_readiness_condition": "pass_metric_symbol == 0 and same-scheme consistency across remaining sectors",
        "qw2049_trace": {
            "lam_eff_over_lam": str(sp.N((lam_eff / lam).subs(qw2049), 12)),
        },
        "strict_core_closure_missing_items": {
            "renormalization": "Need full-sector extension beyond quartic term and fixed (A_div, A_fin) from explicit diagrams.",
            "unitarity": "G_U remains open.",
            "background_independence": "G_BI remains open.",
            "selector_obstruction": "QW-2191 remains open.",
        },
        "false_pass_guard": "Single-term witness stub is not full renormalization gate closure.",
        "next_honest_step": "Extend this stub to gauge/fermion/nonminimal/gravity sectors and then attempt true G_R promotion update.",
        "lay_explanation": "To pierwszy bardzo konkretny test renormalizacji: pokazujemy na jednym składniku, jak mają znosić się nieskończoności i co zostaje fizycznie.",
    }

    path = GEN / "p1892_s842_strict_gr_first_witness_stub_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
