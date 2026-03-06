#!/usr/bin/env python3
"""QW-2273: RG residual non-axiomatic provider evidence gate.

Performs strict lexical evidence audit for residual symbol:
- theorem definition presence,
- anti-axiomatic filters,
- no *_DerivedOrPending token dependency in candidate file.
"""

from __future__ import annotations

import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent


def load(path: str) -> dict[str, Any]:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def theorem_defined_in_file(text: str, symbol: str) -> bool:
    return bool(re.search(rf"^theorem\s+{re.escape(symbol)}\s*:", text, flags=re.M))


def file_has_axioms(text: str) -> bool:
    return bool(re.search(r"^\s*axiom\s+", text, flags=re.M))


def file_has_derived_or_pending_tokens(text: str) -> bool:
    return "_DerivedOrPending" in text


def main() -> None:
    spec = load("spec_qw2269_rg_residual_core_blocker_discharge_packet.json")
    obligations = spec.get("obligations", [])
    if not obligations:
        raise SystemExit("No obligations in spec_qw2269")

    derived_symbol = obligations[0].get("required_outcome", {}).get("introduce_symbol", "")
    target_symbol = obligations[0].get("target_symbol", "")

    rows: list[dict[str, Any]] = []
    candidate_files = []
    strict_non_axiomatic_candidates = []

    for p in sorted(ROOT.glob("*.lean")):
        txt = p.read_text(encoding="utf-8", errors="ignore")
        has_theorem = theorem_defined_in_file(txt, derived_symbol)
        if not has_theorem:
            continue

        has_axiom = file_has_axioms(txt)
        has_dop = file_has_derived_or_pending_tokens(txt)
        passes = (not has_axiom) and (not has_dop)

        candidate_files.append(p.name)
        if passes:
            strict_non_axiomatic_candidates.append(p.name)

        rows.append(
            {
                "file": p.name,
                "defines_symbol": has_theorem,
                "has_axiom_tokens": has_axiom,
                "has_derived_or_pending_tokens": has_dop,
                "passes_strict_non_axiomatic_filter": passes,
            }
        )

    flags = {
        "residual_spec_present": True,
        "derived_symbol_target_present": bool(derived_symbol),
        "target_symbol_present": bool(target_symbol),
        "candidate_files_found": len(candidate_files) > 0,
        "strict_non_axiomatic_candidate_found": len(strict_non_axiomatic_candidates) > 0,
        "o1c_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    verdict = (
        "RG_RESIDUAL_NON_AXIOMATIC_PROVIDER_EVIDENCE_GATE_PASS_PARTIAL_NO_STRICT_CANDIDATE"
        if flags["residual_spec_present"] and flags["derived_symbol_target_present"]
        else "RG_RESIDUAL_NON_AXIOMATIC_PROVIDER_EVIDENCE_GATE_FAIL"
    )

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": "spec_qw2269_rg_residual_core_blocker_discharge_packet.json",
        "target_symbol": target_symbol,
        "derived_symbol": derived_symbol,
        "theorem_definition_evidence": rows,
        "candidate_files": candidate_files,
        "strict_non_axiomatic_candidates": strict_non_axiomatic_candidates,
        "scope_boundary": {
            "strict_non_axiomatic_provider_found": len(strict_non_axiomatic_candidates) > 0,
            "non_axiomatic_discharge_completed": False,
            "o1c_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }

    proof_path = ROOT / "proof_object_qw2273_rg_residual_non_axiomatic_provider_evidence.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": proof_obj["source"],
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "target_symbol": target_symbol,
        "derived_symbol": derived_symbol,
        "n_candidate_files": len(candidate_files),
        "n_strict_non_axiomatic_candidates": len(strict_non_axiomatic_candidates),
        "candidate_files": candidate_files,
        "strict_non_axiomatic_candidates": strict_non_axiomatic_candidates,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "PRODUCE_STRICT_NON_AXIOMATIC_DERIVED_PROVIDER_FOR_RG_RESIDUAL_SYMBOL",
    }

    out_json = ROOT / "report_qw2273_rg_residual_non_axiomatic_provider_evidence_gate.json"
    out_md = ROOT / "RAPORT_QW2273_RG_RESIDUAL_NON_AXIOMATIC_PROVIDER_EVIDENCE_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2273: RG RESIDUAL NON-AXIOMATIC PROVIDER EVIDENCE GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- target_symbol: `{target_symbol}`",
                f"- derived_symbol: `{derived_symbol}`",
                f"- n_candidate_files: `{len(candidate_files)}`",
                f"- n_strict_non_axiomatic_candidates: `{len(strict_non_axiomatic_candidates)}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "verdict": verdict,
                "n_candidate_files": len(candidate_files),
                "n_strict_non_axiomatic_candidates": len(strict_non_axiomatic_candidates),
            }
        )
    )


if __name__ == "__main__":
    main()
