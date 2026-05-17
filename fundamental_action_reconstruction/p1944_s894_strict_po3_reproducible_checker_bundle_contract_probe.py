#!/usr/bin/env python3
"""P1944 S894 strict PO3 reproducible checker-bundle contract probe."""
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
    p1943 = load("p1943_s893_strict_po3_real_checker_run_gate_request_probe.json")

    out = {
        "packet_id": "P1944",
        "stage_id": "S894",
        "status": "OPEN_REPRODUCIBLE_CHECKER_BUNDLE_CONTRACT_NOT_FULFILLED",
        "route": "strict_only",
        "depends_on": {
            "p1943_present": "real_checker_run_contract" in p1943,
            "p1943_no_real_run": p1943.get("real_checker_run_contract", {}).get("current_state") == "NO_REAL_RUN_ARTIFACTS",
        },
        "strict_chain_anchor": "K_strict -> coefficients -> full non-skeleton L_SM+L_GR -> EOM -> theorem object -> reproducible checker bundle",
        "reproducible_bundle_contract_v1": {
            "target_theorem": "THM_A_EPS_NONEMPTY_V1",
            "required_files": [
                "artifacts/thm_a_eps_nonempty/checker_stdout.log",
                "artifacts/thm_a_eps_nonempty/checker_stderr.log",
                "artifacts/thm_a_eps_nonempty/checker_meta.json",
                "artifacts/thm_a_eps_nonempty/proof_artifact.sha256",
                "scripts/reproduce_thm_a_eps_nonempty.sh",
            ],
            "checker_meta_required_keys": [
                "exit_code",
                "toolchain_version",
                "domain_encoding_digest",
                "proof_artifact_hash_sha256",
                "run_timestamp_utc",
            ],
            "acceptance_rule": "Bundle is admissible only if files exist, meta keys are complete, and exit_code==0 with hash consistency.",
            "current_state": "MISSING_BUNDLE_FILES",
        },
        "gate_recheck": {
            "verdict": "NOT_PASSED",
            "reason": "No reproducible artifact bundle exported yet.",
            "false_pass_guard": "Do not infer machine verification from schema-only declarations.",
        },
        "strict_core_statusvector_restamp": {
            "renormalization": "OPEN",
            "unitarity": "OPEN",
            "background_independence": "OPEN_PO3_REPRODUCIBLE_BUNDLE_REQUIRED",
            "selector_qw2191": "OPEN",
        },
        "toe_remaining_minimum_after_p1944": {
            "current_total_open": 7,
            "exact_open_blocks": p1943.get("toe_remaining_minimum_after_p1943", {}).get("exact_open_blocks", []),
            "explanation": "P1944 adds reproducibility contract details but no theorem-grade block is discharged.",
        },
        "next_honest_step": "Export P1945 with actual bundle files and run reproducibility script twice to confirm deterministic hash and checker exit code.",
        "lay_explanation": "Ile zostało do ToE? Nadal 7 bloków. Teraz dokładnie wiadomo jakie pliki i metadane muszą zostać dostarczone, żeby komputerowy dowód był powtarzalny i wiarygodny.",
    }

    path = GEN / "p1944_s894_strict_po3_reproducible_checker_bundle_contract_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
