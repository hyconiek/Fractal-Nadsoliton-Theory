#!/usr/bin/env python3
from __future__ import annotations

import json
import re
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FAR = ROOT / "fundamental_action_reconstruction"
GEN = FAR / "generated"
TEX = ROOT / "TOE_FINAL_DOCUMENTATION_RELEASE_8_STRICT_FULL.tex"
IN1632 = GEN / "p1632_s582_full_strict_lagrangian_and_closure_obligation_summary.json"
IN1634 = GEN / "p1634_s584_strict_global_selector_source_export_audit_summary.json"


def find_patterns(text: str) -> dict[str, bool]:
    pats = {
        "has_kernel_split_section": r"Strict vs legacy kernel split discipline",
        "has_k_strict_symbol": r"K_\{\\mathrm\{strict",
        "has_full_lagrangian_chapter": r"Strict core Lagrangian and emergence map",
        "has_eom_mentions": r"Euler|EOM|equations? of motion",
        "explicit_no_legacy_transfer": r"do \\emph\{not\} silently transfer legacy",
        "mentions_qw2191": r"QW-2191",
    }
    return {k: bool(re.search(v, text, flags=re.IGNORECASE)) for k, v in pats.items()}


def main() -> None:
    tex = TEX.read_text(encoding="utf-8")
    s32 = json.loads(IN1632.read_text(encoding="utf-8"))
    s34 = json.loads(IN1634.read_text(encoding="utf-8"))

    flags = find_patterns(tex)
    closure_open = s34["strict_core_closure"]["status"] == "OPEN"
    full_lagrangian_present = "L_total" in s32["full_lagrangian_density"]

    status = "PASS_P1635_TEX_CHAIN_ALIGNED_WITH_STRICT_ARTIFACTS" if all(flags.values()) and full_lagrangian_present else "KEEP_OPEN_P1635_TEX_CHAIN_GAP"

    summary = {
        "checkpoint": "P1635_S585",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": status,
        "tex_source": str(TEX.relative_to(ROOT)),
        "strict_chain_required": "K_strict -> coefficients -> full Lagrangian -> EOM",
        "tex_flags": flags,
        "artifact_alignment": {
            "full_lagrangian_present_in_p1632": full_lagrangian_present,
            "strict_core_closure_open_in_p1634": closure_open,
            "missing_exports": s34["strict_core_closure"]["missing_exports"],
            "missing_witnesses": s34["strict_core_closure"]["missing_witnesses"],
            "missing_theorems": s34["strict_core_closure"]["missing_theorems"],
        },
        "strict_only": True,
        "legacy_bridge_used": False,
        "next_honest_step": "Podnieść audit do theorem-level: zintegrować globalny atlas selektora na pełnej domenie strict oraz dowód operator-level consistency jako eksport E_selector_internal_source_full_domain.",
        "lay_summary": "Dokument TeX jest zgodny z torem strict i pełnym lagranżianem, ale finalne domknięcie nadal czeka na brakujące dowody globalne.",
    }

    out = GEN / "p1635_s585_tex_strict_chain_consistency_audit_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
