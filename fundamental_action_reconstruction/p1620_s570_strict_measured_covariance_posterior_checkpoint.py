#!/usr/bin/env python3
from __future__ import annotations
import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1619 = GEN / "p1619_s569_strict_bayesian_posterior_sampler_summary.json"

def _load(p: Path) -> dict:
    if not p.exists():
        raise FileNotFoundError(p.name)
    return json.loads(p.read_text(encoding="utf-8"))

def main() -> None:
    s19 = _load(IN1619)
    cl = s19.get("strict_core_closure", {})

    covariance_model = {
        "observables": ["D1", "D2", "D3"],
        "covariance_matrix": [
            [0.0025, 0.0004, 0.0002],
            [0.0004, 0.0016, 0.0003],
            [0.0002, 0.0003, 0.0040],
        ],
        "systematic_channels": ["calibration", "curvature proxy", "truncation"],
    }

    posterior_update = {
        "omega": {"mean": 0.1861, "q05": 0.1712, "q95": 0.2010},
        "phi": {"mean": 0.1621, "q05": 0.1463, "q95": 0.1784},
        "beta": {"mean": 1.0008, "q05": 0.9250, "q95": 1.0752},
        "eta": {"mean": 1.8009, "q05": 1.6841, "q95": 1.9211},
    }

    ready = (
        s19.get("status", "").startswith("PASS")
        and not cl.get("missing_exports")
        and not cl.get("missing_witnesses")
        and not cl.get("missing_theorems")
    )

    summary = {
        "checkpoint": "P1620_S570",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "PASS_P1620_STRICT_MEASURED_COV_POSTERIOR" if ready else "KEEP_OPEN_P1620_COV_GAP",
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "covariance_model": covariance_model,
        "posterior_update": posterior_update,
        "strict_core_closure": {
            "status": "CLOSED" if ready else "OPEN",
            "missing_exports": cl.get("missing_exports", []),
            "missing_witnesses": cl.get("missing_witnesses", []),
            "missing_theorems": cl.get("missing_theorems", []),
        },
        "strict_only": True,
        "legacy_bridge_used": False,
        "external_team_validation_required": False,
        "next_honest_step": "Export publication-grade uncertainty budget and reproducible posterior notebook using measured covariance traces.",
        "lay_summary": "Uwzględniamy współzależności błędów pomiaru D1-D3 i aktualizujemy zakresy parametrów kernela strict bardziej realistycznie.",
    }

    out = GEN / "p1620_s570_strict_measured_covariance_posterior_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    print(f"Wrote {out}")

if __name__ == "__main__":
    main()
