#!/usr/bin/env python3
"""QW-2447: strict anti-false-pass integrity gate.

Cross-checks QW-2440..QW-2446 for:
- no false theorem/full-closure pass claims,
- explicit blocker preservation,
- runtime/provisioning status coherence,
- preserved boundary: all_strict_obligations_fully_closed == False.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent

REPORT_FILES = {
    "q2440": "report_qw2440_grep_frontier_single_foundation_audit_gate.json",
    "q2442": "report_qw2442_dual_nadsoliton_single_foundation_provider_execution_status_gate.json",
    "q2443": "report_qw2443_dual_nadsoliton_single_foundation_minimal_blocker_cut_gate.json",
    "q2444": "report_qw2444_lean_runtime_discovery_gate.json",
    "q2445": "report_qw2445_dual_nadsoliton_single_foundation_provider_execution_status_v2_gate.json",
    "q2446": "report_qw2446_lean_runtime_provisioning_attempt_gate.json",
}

SPEC_2443 = "spec_qw2443_dual_nadsoliton_single_foundation_minimal_blocker_cut.json"
READINESS_MD = "RAPORT_STAN_TEORII_FIN_V5_1_READINESS_2026-03-05.md"

EXPECTED_L12 = "RG_CanonicalAction_to_WellPosedness_EXPORT"
EXPECTED_L5 = "QFT_CanonicalAction_to_Positivity_EXPORT"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def load_json(name: str) -> dict[str, Any] | None:
    path = ROOT / name
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    loaded = {k: load_json(v) for k, v in REPORT_FILES.items()}
    missing_reports = [name for name, obj in loaded.items() if obj is None]

    spec_2443 = load_json(SPEC_2443)
    readiness_path = ROOT / READINESS_MD
    readiness_text = readiness_path.read_text(encoding="utf-8") if readiness_path.exists() else ""

    verdicts = {k: (v or {}).get("verdict", "") for k, v in loaded.items()}

    overclaim_tokens = [
        "FULL_CLOSURE_READY",
        "THEORY_FULLY_CLOSED",
        "PASS_TERMINAL_THEORY_CLOSED",
        "THEOREM_LEVEL_FULL_PASS",
    ]
    no_overclaim_tokens_in_verdicts = all(
        all(tok not in verdict for tok in overclaim_tokens) for verdict in verdicts.values()
    )

    l12_unknown = set((loaded["q2442"] or {}).get("unknown_identifiers", {}).get("l12", []))
    l5_unknown = set((loaded["q2442"] or {}).get("unknown_identifiers", {}).get("l5", []))
    cut_symbols = set()
    if spec_2443 is not None:
        for row in spec_2443.get("minimal_blocker_cut", []):
            sym = row.get("symbol")
            if isinstance(sym, str):
                cut_symbols.add(sym)

    q2444_verdict = verdicts["q2444"]
    q2444_runtime_available = q2444_verdict == "LEAN_RUNTIME_DISCOVERY_GATE_PASS_RUNTIME_AVAILABLE"
    q2444_runtime_unavailable = q2444_verdict == "LEAN_RUNTIME_DISCOVERY_GATE_PASS_PARTIAL_BLOCKED_BY_NO_RUNTIME"

    q2445_flags = (loaded["q2445"] or {}).get("flags", {})
    q2446_flags = (loaded["q2446"] or {}).get("flags", {})

    q2445_coherent_with_q2444 = (
        (q2444_runtime_available and q2445_flags.get("q2444_runtime_available") is True and q2445_flags.get("execution_attempt_performed") is True)
        or (q2444_runtime_unavailable and verdicts["q2445"] == "DUAL_NADSOLITON_SINGLE_FOUNDATION_PROVIDER_EXECUTION_STATUS_V2_GATE_PASS_PARTIAL_BLOCKED_BY_NO_RUNTIME")
    )
    q2446_coherent_with_q2444 = (
        (q2444_runtime_available and q2446_flags.get("q2444_runtime_available") is True and q2446_flags.get("q2444_runtime_unavailable_confirmed") is False)
        or (q2444_runtime_unavailable and q2446_flags.get("q2444_runtime_unavailable_confirmed") is True)
    )

    all_not_fully_closed_false = True
    for obj in loaded.values():
        if obj is None:
            continue
        flags = obj.get("flags", {})
        if flags.get("all_strict_obligations_fully_closed") is not False:
            all_not_fully_closed_false = False
            break

    allowed_q2446_verdicts = {
        "LEAN_RUNTIME_PROVISIONING_ATTEMPT_GATE_PASS_SKIPPED_RUNTIME_ALREADY_AVAILABLE",
        "LEAN_RUNTIME_PROVISIONING_ATTEMPT_GATE_PASS_RUNTIME_INSTALLED_LOCALLY",
        "LEAN_RUNTIME_PROVISIONING_ATTEMPT_GATE_PASS_PARTIAL_BLOCKED_BY_NETWORK_DNS",
        "LEAN_RUNTIME_PROVISIONING_ATTEMPT_GATE_PASS_PARTIAL_BLOCKED_BY_ENVIRONMENT",
    }

    flags = {
        "all_required_reports_present": len(missing_reports) == 0,
        "q2443_spec_present": spec_2443 is not None,
        "readiness_declares_not_ready": "RELEASE_5_1_FULL_CLOSURE_NOT_READY" in readiness_text,
        "readiness_no_ready_overclaim": "RELEASE_5_1_FULL_CLOSURE_READY" not in readiness_text,
        "no_overclaim_tokens_in_verdicts": no_overclaim_tokens_in_verdicts,
        "q2442_expected_unknown_symbols_isolated": EXPECTED_L12 in l12_unknown and EXPECTED_L5 in l5_unknown,
        "q2443_minimal_cut_matches_expected_dual_symbols": cut_symbols == {EXPECTED_L12, EXPECTED_L5},
        "q2445_coherent_with_q2444": q2445_coherent_with_q2444,
        "q2446_coherent_with_q2444": q2446_coherent_with_q2444,
        "q2446_verdict_allowed_under_strict_runtime_semantics": verdicts["q2446"] in allowed_q2446_verdicts,
        "all_strict_obligations_fully_closed_false_everywhere": all_not_fully_closed_false,
        "core_blockers_still_explicit": (
            "BLOCKED_BY_MISSING_CANONICAL_EXPORT_SYMBOLS" in verdicts["q2442"]
            and "BLOCKED_BY_MISSING_CANONICAL_EXPORT_SYMBOLS" in verdicts["q2445"]
        ),
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    hard_requirements = [
        "all_required_reports_present",
        "q2443_spec_present",
        "readiness_declares_not_ready",
        "readiness_no_ready_overclaim",
        "no_overclaim_tokens_in_verdicts",
        "q2442_expected_unknown_symbols_isolated",
        "q2443_minimal_cut_matches_expected_dual_symbols",
        "q2445_coherent_with_q2444",
        "q2446_coherent_with_q2444",
        "q2446_verdict_allowed_under_strict_runtime_semantics",
        "all_strict_obligations_fully_closed_false_everywhere",
        "core_blockers_still_explicit",
    ]

    hard_ok = all(flags[k] for k in hard_requirements)
    if hard_ok:
        verdict = "STRICT_ANTI_FALSE_PASS_INTEGRITY_GATE_PASS_WITH_BLOCKERS_EXPLICIT"
        required_next_step = "DISCHARGE_DUAL_CANONICAL_EXPORT_PROVIDER_SYMBOLS_NON_AXIOMATICALLY"
    else:
        verdict = "STRICT_ANTI_FALSE_PASS_INTEGRITY_GATE_FAIL_INCONSISTENT_CHAIN_OR_OVERCLAIM_RISK"
        required_next_step = "REPAIR_INCONSISTENCIES_AND_RERUN_QW2440_TO_QW2447"

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "reports": REPORT_FILES,
            "spec_2443": SPEC_2443,
            "readiness": READINESS_MD,
        },
        "missing_reports": missing_reports,
        "verdicts": verdicts,
        "unknown_identifiers": {
            "q2442_l12": sorted(l12_unknown),
            "q2442_l5": sorted(l5_unknown),
        },
        "cut_symbols_q2443": sorted(cut_symbols),
        "hard_requirements": hard_requirements,
        "flags": flags,
        "scope_boundary": {
            "integrity_audit_only": True,
            "theorem_discharge_completed": False,
            "all_strict_obligations_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }
    proof_path = ROOT / "proof_object_qw2447_strict_anti_false_pass_integrity_gate.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "reports": REPORT_FILES,
            "spec_2443": SPEC_2443,
            "readiness": READINESS_MD,
        },
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": required_next_step,
    }

    out_json = ROOT / "report_qw2447_strict_anti_false_pass_integrity_gate.json"
    out_md = ROOT / "RAPORT_QW2447_STRICT_ANTI_FALSE_PASS_INTEGRITY_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2447: STRICT ANTI-FALSE-PASS INTEGRITY GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- missing_reports: `{missing_reports}`",
                f"- q2442_unknown_symbols: `{sorted(l12_unknown | l5_unknown)}`",
                f"- q2443_cut_symbols: `{sorted(cut_symbols)}`",
                "",
                "## Wniosek rygorystyczny",
                "- Zachowano granice anty-overclaim: brak podstaw do full-closure PASS.",
                "- Blockery canonical-export (RG/QFT) pozostaja jawne i aktywne.",
                "- Flaga `all_strict_obligations_fully_closed` pozostaje `False`.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
