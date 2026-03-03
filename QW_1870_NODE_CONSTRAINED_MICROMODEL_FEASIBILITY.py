#!/usr/bin/env python3
"""
QW-1870: Node-constrained micromodel feasibility audit.

Tests whether adding node-model priors (from QW-1869) can realistically improve
micromodel-to-canonical compatibility from QW-1862.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1870_node_constrained_micromodel_feasibility.json"
OUT_MD = ROOT / "RAPORT_QW1870_NODE_CONSTRAINED_MICROMODEL_FEASIBILITY.md"


def read_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def clip(x: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, x))


def circular_diff(a: float, b: float) -> float:
    d = (a - b + math.pi) % (2.0 * math.pi) - math.pi
    return abs(d)


def score_from_deltas(do: float, dp: float, db: float, so: float, sp: float, sb: float) -> float:
    zo = do / max(so, 1e-9)
    zp = dp / max(sp, 1e-9)
    zb = db / max(sb, 1e-9)

    zo = min(zo, 20.0)
    zp = min(zp, 20.0)
    zb = min(zb, 20.0)

    return math.exp(-0.5 * (zo * zo + zp * zp + zb * zb) / 3.0)


def main() -> None:
    d1862 = read_json("report_qw1862_micromodel_canonical_compatibility_audit.json")
    d1869 = read_json("report_qw1869_node_model_competition_canonical.json")

    target = d1862.get("target_tuple", {})
    t_omega = float(target.get("omega", math.pi / 4.0))
    t_phi = float(target.get("phi", math.pi / 6.0))
    t_beta = float(target.get("beta", 0.01))

    winner = d1869.get("winner", "M_A_2_5_8_11_or_2plus3n")
    p_win = float(d1869.get("winner_posterior", 0.0))

    base_bounds = {
        "M_A_2_5_8_11_or_2plus3n": {"omega": 0.12, "phi": 0.25, "beta": 0.08},
        "M_B_2_8_14": {"omega": 0.18, "phi": 0.35, "beta": 0.10},
        "M_MIX_A_and_B": {"omega": 0.25, "phi": 0.45, "beta": 0.15},
    }

    b = base_bounds.get(winner, base_bounds["M_A_2_5_8_11_or_2plus3n"])
    scale = 0.5 + p_win
    bounds = {k: v * scale for k, v in b.items()}

    rows = []
    for src in d1862.get("sources", []):
        name = src.get("source", "unknown")
        est = src.get("estimate", {})
        sig = src.get("sigma_used", {})

        eo = float(est.get("omega", math.nan))
        ep = float(est.get("phi", math.nan))
        eb = float(est.get("beta", math.nan))

        so = float(sig.get("omega", 0.05))
        sp = float(sig.get("phi", 0.08))
        sb = float(sig.get("beta", 0.01))

        do_req = t_omega - eo
        dp_req = ((t_phi - ep + math.pi) % (2.0 * math.pi)) - math.pi
        db_req = t_beta - eb

        do_use = clip(do_req, -bounds["omega"], bounds["omega"])
        dp_use = clip(dp_req, -bounds["phi"], bounds["phi"])
        db_use = clip(db_req, -bounds["beta"], bounds["beta"])

        eo_adj = eo + do_use
        ep_adj = ep + dp_use
        eb_adj = eb + db_use

        d0_o = abs(t_omega - eo)
        d0_p = circular_diff(t_phi, ep)
        d0_b = abs(t_beta - eb)
        s0 = score_from_deltas(d0_o, d0_p, d0_b, so, sp, sb)

        d1_o = abs(t_omega - eo_adj)
        d1_p = circular_diff(t_phi, ep_adj)
        d1_b = abs(t_beta - eb_adj)
        s1 = score_from_deltas(d1_o, d1_p, d1_b, so, sp, sb)

        rows.append(
            {
                "source": name,
                "score_original": float(src.get("final_score", s0)),
                "score_recomputed_original": s0,
                "score_after_node_constraint": s1,
                "improvement_factor": s1 / max(float(src.get("final_score", s0)), 1e-12),
                "required_shift": {"omega": do_req, "phi": dp_req, "beta": db_req},
                "allowed_shift": bounds,
                "applied_shift": {"omega": do_use, "phi": dp_use, "beta": db_use},
                "residual_after_shift": {"omega": d1_o, "phi": d1_p, "beta": d1_b},
            }
        )

    max_adj = max((r["score_after_node_constraint"] for r in rows), default=0.0)
    mean_adj = sum((r["score_after_node_constraint"] for r in rows), 0.0) / max(1, len(rows))

    if max_adj >= 0.50 and mean_adj >= 0.30:
        verdict = "NODE_CONSTRAINED_COMPATIBILITY_PLAUSIBLE"
    elif max_adj >= 0.20:
        verdict = "NODE_CONSTRAINED_COMPATIBILITY_STILL_WEAK"
    else:
        verdict = "NODE_CONSTRAINED_COMPATIBILITY_FAIL_STRONG"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "winner_model_1869": winner,
        "winner_posterior_1869": p_win,
        "node_prior_shift_bounds": bounds,
        "rows": rows,
        "summary": {
            "max_score_after_constraint": max_adj,
            "mean_score_after_constraint": mean_adj,
        },
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1870: NODE-CONSTRAINED MICROMODEL FEASIBILITY",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- node model from 1869: {winner} (posterior={p_win:.3f})",
        "",
        "## Applied Shift Bounds",
        f"- omega: +/-{bounds['omega']:.3f}",
        f"- phi: +/-{bounds['phi']:.3f}",
        f"- beta: +/-{bounds['beta']:.3f}",
        "",
        "## Source Results",
    ]

    for r in rows:
        lines.append(
            f"- {r['source']}: original={r['score_original']:.3e}, after={r['score_after_node_constraint']:.3e}, improvement={r['improvement_factor']:.3e}"
        )

    lines += [
        "",
        "## Summary",
        f"- max score after constraint: {max_adj:.3e}",
        f"- mean score after constraint: {mean_adj:.3e}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1870] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1870] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
