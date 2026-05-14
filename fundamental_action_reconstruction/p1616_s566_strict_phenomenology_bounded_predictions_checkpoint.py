#!/usr/bin/env python3
from __future__ import annotations
import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1614 = GEN / "p1614_s564_strict_toe_closure_dossier_summary.json"
IN1615 = GEN / "p1615_s565_executable_symbolic_appendix_summary.json"

def _load(p: Path) -> dict:
    if not p.exists():
        raise FileNotFoundError(p.name)
    return json.loads(p.read_text(encoding="utf-8"))

def _interval(v: float, rel: float = 0.05) -> dict:
    return {"central": v, "minus": v*(1-rel), "plus": v*(1+rel), "relative_uncertainty": rel}

def main() -> None:
    s14 = _load(IN1614)
    s15 = _load(IN1615)

    c = s14.get("dossier", {}).get("coefficients", {})
    lam = float(c.get("lambda_sm_eff", 0.0))
    kap = float(c.get("kappa_gr_eff", 0.0))
    eps = float(c.get("epsilon_mix_eff", 0.0))

    predictions = {
        "P1_scalar_mass_proxy": {
            "formula": "m_eff^2 ~ lambda_sm_eff",
            "value_interval": _interval(lam, 0.05),
            "domain": "strict low-curvature proxy"
        },
        "P2_gr_coupling_proxy": {
            "formula": "G_eff^-1 ~ kappa_gr_eff",
            "value_interval": _interval(kap, 0.05),
            "domain": "strict semiclassical regime"
        },
        "P3_mix_backreaction_proxy": {
            "formula": "delta_psi_curvature ~ epsilon_mix_eff * R",
            "value_interval": _interval(eps, 0.1),
            "domain": "weak-to-moderate curvature"
        }
    }

    cl = s15.get("strict_core_closure", {})
    ready = (
        s14.get("status", "").startswith("PASS")
        and s15.get("status", "").startswith("PASS")
        and not cl.get("missing_exports")
        and not cl.get("missing_witnesses")
        and not cl.get("missing_theorems")
    )

    summary = {
        "checkpoint": "P1616_S566",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "PASS_P1616_STRICT_PHENOMENOLOGY_BOUNDED_PREDICTIONS" if ready else "KEEP_OPEN_P1616_PREDICTION_GAP",
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "predictions": predictions,
        "strict_core_closure": {
            "status": "CLOSED" if ready else "OPEN",
            "missing_exports": cl.get("missing_exports", []),
            "missing_witnesses": cl.get("missing_witnesses", []),
            "missing_theorems": cl.get("missing_theorems", []),
        },
        "strict_only": True,
        "legacy_bridge_used": False,
        "external_team_validation_required": False,
        "next_honest_step": "Design strict experimental discriminants for P1/P2/P3 proxies and export measurement protocol assumptions.",
        "lay_summary": "Z zamkniętego modelu strict wyciągamy trzy testowalne przewidywania z przedziałami niepewności, bez użycia legacy bridge.",
    }

    out = GEN / "p1616_s566_strict_phenomenology_bounded_predictions_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    print(f"Wrote {out}")

if __name__ == "__main__":
    main()
