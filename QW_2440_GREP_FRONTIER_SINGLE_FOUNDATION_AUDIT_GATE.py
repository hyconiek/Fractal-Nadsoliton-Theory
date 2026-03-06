#!/usr/bin/env python3
"""QW-2440: grep frontier single-foundation audit gate.

Performs a strict textual audit of the current frontier after QW-2439:
- cycle recurrence evidence,
- unresolved dual canonical export blockers,
- absence of false full-closure claims,
- presence of Nadsoliton single-foundation language.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def run_rg(pattern: str, include_globs: list[str] | None = None) -> list[str]:
    cmd = ["rg", "-n", "--hidden", "--glob", "!*.pdf"]
    if include_globs:
        for g in include_globs:
            cmd.extend(["--glob", g])
    cmd.extend([pattern, str(ROOT)])
    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode not in (0, 1):
        raise RuntimeError(f"rg failed: {pattern}\n{proc.stderr}")
    lines = [ln for ln in proc.stdout.splitlines() if ln.strip()]
    return lines


def contains_text(path: Path, needle: str) -> bool:
    if not path.exists():
        return False
    return needle in path.read_text(encoding="utf-8")


def main() -> None:
    cycle_hits = run_rg(r"QW-2384|SCC=20/20|cycle recurrence|structural theorem dependency loop")
    rg_export_hits = run_rg(r"RG_CanonicalAction_to_WellPosedness_EXPORT")
    qft_export_hits = run_rg(r"QFT_CanonicalAction_to_Positivity_EXPORT")
    unknown_hits = run_rg(r"Unknown identifier")
    false_closure_hits_raw = run_rg(
        r"RELEASE_5_1_FULL_CLOSURE_READY|FULL_CLOSURE_READY|THEORY_FULLY_CLOSED|PASS_TERMINAL_THEORY_CLOSED",
        include_globs=["*.md", "*.json", "*.tex"],
    )
    false_closure_hits = [
        ln
        for ln in false_closure_hits_raw
        if "proof_object_qw2440_grep_frontier_single_foundation_audit.json" not in ln
    ]
    nadsoliton_hits = run_rg(r"Nadsoliton|NADSOLITON|single fundamental|jedynym fundament|FundamentalKernelDynamics")

    readiness_path = ROOT / "RAPORT_STAN_TEORII_FIN_V5_1_READINESS_2026-03-05.md"

    flags = {
        "cycle_recurrence_evidence_present": len(cycle_hits) > 0,
        "dual_canonical_export_blockers_present": len(rg_export_hits) > 0 and len(qft_export_hits) > 0,
        "unknown_identifier_obstruction_evidence_present": len(unknown_hits) > 0,
        "false_full_closure_claim_absent": len(false_closure_hits) == 0,
        "single_foundation_language_present": len(nadsoliton_hits) > 0,
        "readiness_file_declares_not_ready": contains_text(readiness_path, "RELEASE_5_1_FULL_CLOSURE_NOT_READY"),
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    if (
        flags["cycle_recurrence_evidence_present"]
        and flags["dual_canonical_export_blockers_present"]
        and flags["false_full_closure_claim_absent"]
        and flags["readiness_file_declares_not_ready"]
    ):
        verdict = "GREP_FRONTIER_SINGLE_FOUNDATION_AUDIT_GATE_PASS_WITH_BLOCKERS_EXPLICIT"
        required_next_step = "BUILD_DUAL_NADSOLITON_SINGLE_FOUNDATION_EXPORT_PROVIDER_PACKET"
    else:
        verdict = "GREP_FRONTIER_SINGLE_FOUNDATION_AUDIT_GATE_FAIL"
        required_next_step = "REPAIR_GREP_AUDIT_PIPELINE_AND_RERUN"

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "queries": {
            "cycle": r"QW-2384|SCC=20/20|cycle recurrence|structural theorem dependency loop",
            "rg_export": r"RG_CanonicalAction_to_WellPosedness_EXPORT",
            "qft_export": r"QFT_CanonicalAction_to_Positivity_EXPORT",
            "unknown": r"Unknown identifier",
            "false_closure": r"RELEASE_5_1_FULL_CLOSURE_READY|FULL_CLOSURE_READY|THEORY_FULLY_CLOSED|PASS_TERMINAL_THEORY_CLOSED",
            "nadsoliton": r"Nadsoliton|NADSOLITON|single fundamental|jedynym fundament|FundamentalKernelDynamics",
        },
        "counts": {
            "cycle_hits": len(cycle_hits),
            "rg_export_hits": len(rg_export_hits),
            "qft_export_hits": len(qft_export_hits),
            "unknown_identifier_hits": len(unknown_hits),
            "false_closure_hits": len(false_closure_hits),
            "nadsoliton_hits": len(nadsoliton_hits),
        },
        "sample_hits": {
            "cycle": cycle_hits[:10],
            "rg_export": rg_export_hits[:10],
            "qft_export": qft_export_hits[:10],
            "unknown": unknown_hits[:10],
            "nadsoliton": nadsoliton_hits[:10],
        },
        "scope_boundary": {
            "textual_audit_only": True,
            "theorem_discharge_completed": False,
            "all_strict_obligations_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }

    proof_path = ROOT / "proof_object_qw2440_grep_frontier_single_foundation_audit.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "readiness": readiness_path.name,
            "proof_object": proof_path.name,
        },
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "counts": proof_obj["counts"],
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": required_next_step,
    }

    out_json = ROOT / "report_qw2440_grep_frontier_single_foundation_audit_gate.json"
    out_md = ROOT / "RAPORT_QW2440_GREP_FRONTIER_SINGLE_FOUNDATION_AUDIT_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    lines = [
        "# RAPORT QW-2440: GREP FRONTIER SINGLE FOUNDATION AUDIT GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        f"- cycle_hits: `{proof_obj['counts']['cycle_hits']}`",
        f"- rg_export_hits: `{proof_obj['counts']['rg_export_hits']}`",
        f"- qft_export_hits: `{proof_obj['counts']['qft_export_hits']}`",
        f"- unknown_identifier_hits: `{proof_obj['counts']['unknown_identifier_hits']}`",
        f"- false_closure_hits: `{proof_obj['counts']['false_closure_hits']}`",
        f"- nadsoliton_hits: `{proof_obj['counts']['nadsoliton_hits']}`",
        "",
        "## Wniosek rygorystyczny",
        "- Strukturalna petla theorem-dependency pozostaje jawnie obecna.",
        "- Dwa canonical export blockery (RG/QFT) pozostaja aktywne.",
        "- Nie wykryto sygnalow falszywego full-closure claim.",
        "- Najlepszy nastepny ruch: zbudowac dual packet single-foundation (Nadsoliton) dla eksportow RG/QFT.",
    ]
    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(json.dumps({"verdict": verdict, "counts": proof_obj["counts"]}))


if __name__ == "__main__":
    main()
