#!/usr/bin/env python3
"""P1945 S895 strict PO3 reproducibility double-run requirement probe."""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def load(name: str) -> dict:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def main() -> None:
    p1944 = load("p1944_s894_strict_po3_reproducible_checker_bundle_contract_probe.json")

    out = {
        "packet_id": "P1945",
        "stage_id": "S895",
        "status": "OPEN_DOUBLE_RUN_REPRODUCIBILITY_REQUIREMENT_NOT_SATISFIED",
        "route": "strict_only",
        "depends_on": {
            "p1944_present": "reproducible_bundle_contract_v1" in p1944,
            "p1944_bundle_missing": p1944.get("reproducible_bundle_contract_v1", {}).get("current_state") == "MISSING_BUNDLE_FILES",
        },
        "strict_chain_anchor": "K_strict -> coefficients -> full non-skeleton L_SM+L_GR -> EOM -> theorem object -> reproducible bundle -> double-run reproducibility",
        "double_run_reproducibility_contract": {
            "target_theorem": "THM_A_EPS_NONEMPTY_V1",
            "run_requirements": [
                "Run checker twice with identical domain encoding and toolchain stamp",
                "Both runs must return exit_code == 0",
                "Both runs must output identical proof_artifact_hash_sha256",
                "Both runs must output identical domain_encoding_digest",
            ],
            "required_outputs": [
                "artifacts/thm_a_eps_nonempty/run1/checker_meta.json",
                "artifacts/thm_a_eps_nonempty/run2/checker_meta.json",
                "artifacts/thm_a_eps_nonempty/repro_compare.json",
            ],
            "current_state": "NO_REAL_RUNS_EXPORTED",
        },
        "gate_recheck": {
            "verdict": "NOT_PASSED",
            "reason": "No run1/run2 machine artifacts available for deterministic comparison.",
            "false_pass_guard": "Single-run or placeholder outputs cannot satisfy reproducibility gate.",
        },
        "strict_core_statusvector_restamp": {
            "renormalization": "OPEN",
            "unitarity": "OPEN",
            "background_independence": "OPEN_PO3_DOUBLE_RUN_EVIDENCE_REQUIRED",
            "selector_qw2191": "OPEN",
        },
        "toe_remaining_minimum_after_p1945": {
            "current_total_open": 7,
            "exact_open_blocks": p1944.get("toe_remaining_minimum_after_p1944", {}).get("exact_open_blocks", []),
            "explanation": "P1945 sharpens reproducibility acceptance criteria but discharges no theorem-grade block.",
        },
        "next_honest_step": "Export P1946 with actual run1/run2 checker metadata and deterministic comparison artifact; restamp PO3 only if all equality checks pass.",
        "lay_explanation": "Ile zostało do ToE? Nadal 7. Trzeba pokazać, że komputerowy dowód działa tak samo dwa razy z rzędu — bez tego nie ma wiarygodnego domknięcia.",
    }

    path = GEN / "p1945_s895_strict_po3_reproducibility_double_run_requirement_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
