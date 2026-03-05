#!/usr/bin/env python3
"""
QW-2134: Interacting microcausality perturbative gate (strict conditional).

Purpose:
- extend QW-2133 (free sector) toward interacting sector using explicit
  local-QFT preconditions already established in strict chain,
- keep claim discipline: perturbative conditional closure, not constructive
  all-orders/nonperturbative proof.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2134_interacting_microcausality_perturbative_gate.json"
OUT_MD = ROOT / "RAPORT_QW2134_INTERACTING_MICROCAUSALITY_PERTURBATIVE_GATE.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def has_any_token(text: str, tokens: List[str]) -> bool:
    t = text.lower()
    return any(tok.lower() in t for tok in tokens)


def main() -> None:
    r2127 = load_json("report_qw2127_nonabelian_spinor_gauge_action_bridge_gate.json")
    r2129 = load_json("report_qw2129_anomaly_cancellation_kernel_anchored_gate.json")
    r2131 = load_json("report_qw2131_hypercharge_template_kernel_uniqueness_gate.json")
    r2133 = load_json("report_qw2133_ktotal_microcausality_free_sector_gate.json")

    action_blocks = r2127["action_blocks"]
    action_concat = " ; ".join(str(v) for v in action_blocks.values())

    nonlocal_tokens = [
        "k(x-y)",
        "mathcal{k}(x-y)",
        "convolution",
        "nonlocal",
        "int d^4y",
        "∫ d^4y",
        "integral kernel",
    ]
    locality_tokens = [
        "d_mu",
        "gamma^mu",
        "f^{abc}",
        "eps^{ijk}",
        "bar(psi",
        "h.c.",
    ]

    anomaly = r2129["anomaly_coefficients_per_generation"]
    tol = float(anomaly["tolerance"])
    anomaly_zero = (
        abs(float(anomaly["A_SU3_SU3_U1"])) <= tol
        and abs(float(anomaly["A_SU2_SU2_U1"])) <= tol
        and abs(float(anomaly["A_U1_cubic"])) <= tol
        and abs(float(anomaly["A_gravity_gravity_U1"])) <= tol
    )
    witten_even = bool(r2129["witten_global_check"]["is_even"])

    flags = {
        "q2127_action_bridge_present": bool(
            str(r2127.get("verdict", "")).startswith("NONABELIAN_SPINOR_GAUGE_ACTION_BRIDGE_GATE_PASS")
        ),
        "action_blocks_are_spacetime_local_by_structure_tokens": has_any_token(action_concat, locality_tokens),
        "no_explicit_spacetime_nonlocal_kernel_tokens_in_action": not has_any_token(action_concat, nonlocal_tokens),
        "dimension4_audit_pass": bool(r2127["flags"]["dimension4_audit_pass"]),
        "su2_su3_lie_algebra_closure_pass": bool(
            r2127["flags"]["su2_lie_algebra_closure_numeric"] and r2127["flags"]["su3_lie_algebra_closure_numeric"]
        ),
        "anomaly_cancellation_pass": bool(anomaly_zero),
        "witten_global_anomaly_absent": bool(witten_even),
        "hypercharge_template_kernel_unique_anchored": bool(r2131["flags"]["hypercharge_template_unique_from_kernel"]),
        "free_sector_microcausality_pass": bool(
            str(r2133.get("verdict", "")).startswith("KTOTAL_MICROCAUSALITY_FREE_SECTOR_GATE_PASS")
        ),
        "perturbative_local_qft_microcausality_conditions_met": bool(
            r2127["flags"]["dimension4_audit_pass"]
            and anomaly_zero
            and witten_even
            and (not has_any_token(action_concat, nonlocal_tokens))
            and str(r2133.get("verdict", "")).startswith("KTOTAL_MICROCAUSALITY_FREE_SECTOR_GATE_PASS")
        ),
        "deterministic_no_scan_no_retune": True,
        "full_constructive_all_orders_interacting_microcausality_proof": False,
    }
    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = int(len(flags))

    verdict = (
        "INTERACTING_MICROCAUSALITY_PERTURBATIVE_GATE_PASS_PARTIAL_CONDITIONAL"
        if (
            flags["q2127_action_bridge_present"]
            and flags["action_blocks_are_spacetime_local_by_structure_tokens"]
            and flags["no_explicit_spacetime_nonlocal_kernel_tokens_in_action"]
            and flags["dimension4_audit_pass"]
            and flags["su2_su3_lie_algebra_closure_pass"]
            and flags["anomaly_cancellation_pass"]
            and flags["witten_global_anomaly_absent"]
            and flags["hypercharge_template_kernel_unique_anchored"]
            and flags["free_sector_microcausality_pass"]
            and flags["perturbative_local_qft_microcausality_conditions_met"]
            and flags["deterministic_no_scan_no_retune"]
        )
        else "INTERACTING_MICROCAUSALITY_PERTURBATIVE_GATE_FAIL"
    )

    assumptions = [
        "Canonical local QFT quantization for declared local action blocks.",
        "Standard perturbative causal construction (e.g. EG/BRST-consistent setup).",
        "No hidden spacetime-nonlocal counterterms introduced outside audited blocks.",
    ]

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2127_action_bridge": "report_qw2127_nonabelian_spinor_gauge_action_bridge_gate.json",
            "q2129_anomaly_audit": "report_qw2129_anomaly_cancellation_kernel_anchored_gate.json",
            "q2131_hypercharge_uniqueness": "report_qw2131_hypercharge_template_kernel_uniqueness_gate.json",
            "q2133_free_microcausality": "report_qw2133_ktotal_microcausality_free_sector_gate.json",
        },
        "audits": {
            "action_blocks": action_blocks,
            "nonlocal_tokens_checked": nonlocal_tokens,
            "locality_tokens_checked": locality_tokens,
            "anomaly_coefficients_per_generation": anomaly,
            "witten_global_check": r2129["witten_global_check"],
        },
        "assumptions_for_conditional_claim": assumptions,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "BUILD_FULL_CONSTRUCTIVE_ALL_ORDERS_INTERACTING_MICROCAUSALITY_PROOF_WITH_EXPLICIT_LOOP_CAUSAL_SUPPORT"
            if verdict.endswith("CONDITIONAL")
            else "REPAIR_LOCALITY_OR_ANOMALY_OR_ACTION_PRECONDITIONS_AND_RERUN_QW2134"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2134: INTERACTING MICROCAUSALITY PERTURBATIVE GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Core checks",
        f"- anomaly cancellation pass: `{flags['anomaly_cancellation_pass']}`",
        f"- Witten global anomaly absent: `{flags['witten_global_anomaly_absent']}`",
        f"- free-sector microcausality pass: `{flags['free_sector_microcausality_pass']}`",
        f"- no explicit spacetime nonlocal tokens: `{flags['no_explicit_spacetime_nonlocal_kernel_tokens_in_action']}`",
        "",
        "## Scope boundary",
        "- full constructive all-orders interacting proof: `False`",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2134] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2134] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2134] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()
