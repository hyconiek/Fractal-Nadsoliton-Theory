#!/usr/bin/env python3
from __future__ import annotations
import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1617 = GEN / "p1617_s567_strict_experimental_discriminants_summary.json"
IN1616 = GEN / "p1616_s566_strict_phenomenology_bounded_predictions_summary.json"

def _load(p: Path) -> dict:
    if not p.exists():
        raise FileNotFoundError(p.name)
    return json.loads(p.read_text(encoding="utf-8"))

def main() -> None:
    s17 = _load(IN1617)
    s16 = _load(IN1616)

    protocol = {
        "measurement_model": {
            "D1": "D1_obs = f1(omega, phi, beta, eta) + eps1",
            "D2": "D2_obs = f2(omega, phi, beta, eta) + eps2",
            "D3": "D3_obs = f3(omega, phi, beta, eta) + eps3",
        },
        "error_model": {
            "eps1": "Gaussian(0, sigma1^2)",
            "eps2": "Gaussian(0, sigma2^2)",
            "eps3": "Gaussian(0, sigma3^2)",
            "systematics_included": ["calibration", "curvature proxy", "model truncation"],
        },
        "identifiability": {
            "omega": "high from D1+D3",
            "phi": "medium from D1 phase-sensitive fit",
            "beta": "high from D2",
            "eta": "medium from D3 curvature scaling",
        },
        "posterior_target": "p(omega,phi,beta,eta | D1,D2,D3, strict-model)",
    }

    cl = s17.get("strict_core_closure", {})
    ready = (
        s17.get("status", "").startswith("PASS")
        and s16.get("status", "").startswith("PASS")
        and not cl.get("missing_exports")
        and not cl.get("missing_witnesses")
        and not cl.get("missing_theorems")
    )

    summary = {
        "checkpoint": "P1618_S568",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "PASS_P1618_STRICT_DATA_ASSIMILATION_PROTOCOL" if ready else "KEEP_OPEN_P1618_ASSIMILATION_GAP",
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "protocol": protocol,
        "strict_core_closure": {
            "status": "CLOSED" if ready else "OPEN",
            "missing_exports": cl.get("missing_exports", []),
            "missing_witnesses": cl.get("missing_witnesses", []),
            "missing_theorems": cl.get("missing_theorems", []),
        },
        "strict_only": True,
        "legacy_bridge_used": False,
        "external_team_validation_required": False,
        "next_honest_step": "Implement Bayesian posterior sampler checkpoint with synthetic strict datasets for D1/D2/D3.",
        "lay_summary": "To instrukcja od danych do parametrów modelu: jak z pomiarów D1-D3 odtworzyć liczby kernela strict.",
    }

    out = GEN / "p1618_s568_strict_data_assimilation_protocol_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    print(f"Wrote {out}")

if __name__ == "__main__":
    main()
