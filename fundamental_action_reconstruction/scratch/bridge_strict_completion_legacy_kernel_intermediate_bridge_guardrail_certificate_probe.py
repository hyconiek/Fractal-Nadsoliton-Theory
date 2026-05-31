#!/usr/bin/env python3
"""Scratch probe: legacy-kernel intermediate bridge guardrail certificate.

Question addressed: what is the honest current status of the bridge from
K_legacy_ont to K_strict_gate after the guardrail correction?

This probe does not claim a finished identity theorem.  It audits the updated
repo guardrails as a finite text/logic certificate: the legacy kernel is
restored as an intermediate/incomplete bridge kernel, the strict kernel is the
completed/enriched strict continuation, and any finished bridge must be followed
by a separate legacy-role transfer audit because the legacy kernel did not carry
all strict nadsoliton characteristics, especially strict compression.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FAR = ROOT / "fundamental_action_reconstruction"
OUT_JSON = HERE / "bridge_strict_completion_legacy_kernel_intermediate_bridge_guardrail_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_legacy_kernel_intermediate_bridge_guardrail_certificate_report.md"

AGENTS = ROOT / "AGENTS.md"
K1 = FAR / "K1_LEGACY_ONTOLOGICAL_KERNEL_VS_STRICT_GATE_KERNEL_SPLIT_NOTE.md"
K2 = FAR / "K2_STRICT_GATE_KERNEL_DERIVATION_CHAIN_NOTE.md"
F2 = FAR / "F2_STRICT_GATE_KERNEL_PROVENANCE_AND_FAR_INPUT_CLASSIFICATION_PACKET.md"
F3 = FAR / "F3_CURRENT_FAR_FRONTIER_KERNEL_ARTIFACT_SENSITIVITY_CLASSIFICATION_PACKET.md"
S2 = FAR / "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md"
CHAIN = HERE / "bridge_strict_completion_certificate_chain_integrity_report.json"
LEGACY_TORSION_AUDIT = HERE / "bridge_legacy_torsion_chi11_opinion_audit_report.md"


def read_text(path: Path) -> str:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite text: {path}")
    return path.read_text(encoding="utf-8")


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def contains_all(text: str, needles: list[str]) -> bool:
    return all(needle in text for needle in needles)


def build_payload() -> dict[str, Any]:
    agents = read_text(AGENTS)
    k1 = read_text(K1)
    k2 = read_text(K2)
    f2 = read_text(F2)
    f3 = read_text(F3)
    s2 = read_text(S2)
    torsion_audit = read_text(LEGACY_TORSION_AUDIT)
    chain = load_json(CHAIN)

    evidence_rows = [
        {
            "source": "AGENTS",
            "claim_checked": "automation guardrails restore legacy as intermediate bridge kernel and require role-transfer audit",
            "required_markers": [
                "restored as an **intermediate bridge kernel**",
                "completed/enriched strict working kernel",
                "omits strict-side characteristics of the nadsoliton",
                "audit legacy role transfer separately",
            ],
            "passes": contains_all(agents, [
                "restored as an **intermediate bridge kernel**",
                "completed/enriched strict working kernel",
                "omits strict-side characteristics of the nadsoliton",
                "audit legacy role transfer separately",
            ]),
        },
        {
            "source": "K1/K2",
            "claim_checked": "raw identity is still not silently exported",
            "required_markers": [
                "This note does **not** claim a bridge theorem.",
                "K_legacy_ont(d) == K_strict_gate(d)",
                "not directly derived from the old 4.4",
            ],
            "passes": contains_all(k1, [
                "This note does **not** claim a bridge theorem.",
                "K_legacy_ont(d) == K_strict_gate(d)",
            ]) and "not directly derived from the old 4.4" in k2,
        },
        {
            "source": "S2",
            "claim_checked": "legacy kernel is restored as bridge-intermediate and strict is completed/enriched continuation",
            "required_markers": [
                "the legacy kernel is **not** a discarded dead end",
                "intermediate kernel on the path toward identifying the strict kernel",
                "completed / enriched strict kernel",
                "legacy kernel + missing strict-side characteristics -> strict kernel",
                "role-transfer audit after bridge completion",
            ],
            "passes": contains_all(s2, [
                "the legacy kernel is **not** a discarded dead end",
                "intermediate kernel on the path toward identifying the strict kernel",
                "completed / enriched strict kernel",
                "legacy kernel + missing strict-side characteristics -> strict kernel",
                "role-transfer audit after bridge completion",
            ]),
        },
        {
            "source": "S2_compression_clause",
            "claim_checked": "legacy kernel incompleteness includes strict compression/damping characteristics",
            "required_markers": [
                "nonlinear `d^eta`",
                "legacy linear `beta_tors*d` denominator",
                "strict-side additions must remain visible",
            ],
            "passes": contains_all(s2, [
                "nonlinear `d^eta`",
                "legacy linear `beta_tors*d` denominator",
                "strict-side additions must remain visible",
            ]),
        },
        {
            "source": "F2/F3",
            "claim_checked": "routing remains disciplined: no silent substitution, prefer robust/bridge-completion routes",
            "required_markers": [
                "must not silently substitute", "kernel-split-robust",
            ],
            "passes": "must not silently substitute" in f2 and "kernel-split-robust" in f3,
        },
        {
            "source": "legacy_torsion_chi11_opinion_audit",
            "claim_checked": "beta_tors -> chi_11 remains a candidate bridge/source theorem target, not a proven transfer",
            "required_markers": [
                "CANDIDATE_BRIDGE_HYPOTHESIS_NOT_THEOREM",
                "role_transfer_control",
                "No beta_tors -> chi_11 quantization/collapse theorem is asserted.",
            ],
            "passes": contains_all(torsion_audit, [
                "CANDIDATE_BRIDGE_HYPOTHESIS_NOT_THEOREM",
                "role_transfer_control",
                "No beta_tors -> chi_11 quantization/collapse theorem is asserted.",
            ]),
        },
        {
            "source": "strict_completion_chain_integrity",
            "claim_checked": "current chain still blocks unqualified identity, role transfer, QW-2191 discharge, and ToE closure",
            "required_markers": [
                "No unqualified identity K_legacy_ont == K_strict_gate is claimed.",
                "No legacy physical-role transfer to K_strict_gate is claimed.",
                "No QW-2191 selector discharge is claimed.",
                "No ToE closure is claimed.",
            ],
            "passes": all(marker in chain["hard_limits"] for marker in [
                "No unqualified identity K_legacy_ont == K_strict_gate is claimed.",
                "No legacy physical-role transfer to K_strict_gate is claimed.",
                "No QW-2191 selector discharge is claimed.",
                "No ToE closure is claimed.",
            ]),
        },
    ]

    summary = {
        "legacy_kernel_restored_as_intermediate": evidence_rows[0]["passes"] and evidence_rows[2]["passes"],
        "strict_kernel_treated_as_completed_legacy_continuation": evidence_rows[0]["passes"] and evidence_rows[2]["passes"],
        "raw_identity_bridge_still_not_silent": evidence_rows[1]["passes"] and evidence_rows[6]["passes"],
        "legacy_kernel_incomplete_for_strict_characteristics": evidence_rows[3]["passes"],
        "strict_compression_missing_from_legacy_recorded": evidence_rows[3]["passes"],
        "role_transfer_audit_required_after_full_bridge": evidence_rows[0]["passes"] and evidence_rows[2]["passes"],
        "beta_tors_to_chi11_remains_candidate_not_theorem": evidence_rows[5]["passes"],
        "q2191_remains_open": "No QW-2191 selector discharge is claimed." in chain["hard_limits"],
        "toe_closure_not_claimed": "No ToE closure is claimed." in chain["hard_limits"],
    }
    summary["intermediate_bridge_guardrail_certificate_passes"] = all(summary.values())

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_LEGACY_KERNEL_INTERMEDIATE_BRIDGE_GUARDRAIL_CERTIFICATE__TEXT_LOGIC_AUDIT",
        "status": "legacy-kernel-restored-as-intermediate-bridge-kernel-role-transfer-audit-required-no-silent-identity",
        "source_reports": {
            "AGENTS_guardrails": str(AGENTS.relative_to(ROOT)),
            "K1_kernel_split_note": str(K1.relative_to(ROOT)),
            "K2_strict_gate_derivation_chain_note": str(K2.relative_to(ROOT)),
            "F2_provenance_packet": str(F2.relative_to(ROOT)),
            "F3_frontier_classification_packet": str(F3.relative_to(ROOT)),
            "S2_strategic_priority_reorientation_packet": str(S2.relative_to(ROOT)),
            "legacy_torsion_chi11_opinion_audit": str(LEGACY_TORSION_AUDIT.relative_to(ROOT)),
            "strict_completion_chain_integrity_report": str(CHAIN.relative_to(ROOT)),
        },
        "grep_disambiguation": {
            "searched_terms": [
                "legacy strict kernel bridge",
                "K_legacy_ont K_strict_gate completed legacy",
                "intermediate bridge kernel",
                "role-transfer theorem",
                "strict compression legacy kernel",
                "beta_tors chi_11 candidate bridge",
            ],
            "finding": "Repo grep finds older retirement/non-bridge and bridge-intuition artifacts; the current guardrail update supersedes retirement-only routing by restoring K_legacy_ont as an intermediate incomplete bridge kernel while preserving no-silent-identity and no-silent-role-transfer controls.",
        },
        "evidence_rows": evidence_rows,
        "legacy_kernel_intermediate_bridge_summary": summary,
        "closure_plan_implications": [
            "Do spend bridge work on the explicit completion map K_legacy_ont -> K_strict_gate rather than treating the kernels as unrelated.",
            "Do not treat the legacy kernel as already containing all strict nadsoliton characteristics; strict compression and certified phase/topology remain completion data.",
            "After a full bridge is specified, run a separate role-transfer audit for every legacy physical role.",
            "Treat beta_tors -> chi_11 as a concrete bridge/source theorem target, not as an already exported theorem.",
            "QW-2191 and ToE closure remain open until a real strict selector/source theorem is exported.",
        ],
        "proof_certificate": {
            "grep_step": "Repo grep separates older non-bridge/retirement language from the current restored intermediate-kernel guardrail.",
            "guardrail_step": "AGENTS and S2 now restore K_legacy_ont as an intermediate bridge kernel and treat K_strict_gate as the completed/enriched strict continuation only through explicit completion evidence.",
            "compression_step": "S2 records that the legacy kernel omits strict-side characteristics, including nonlinear d^eta compression absent from the legacy beta_tors*d denominator.",
            "role_transfer_step": "A full legacy->strict bridge must be followed by a separate role-transfer audit before legacy physical-role claims can move onto K_strict_gate.",
            "selector_step": "beta_tors->chi_11 remains a candidate theorem target and QW-2191 remains open; no selector discharge is claimed.",
            "verdict_step": "The active route is bridge-completion with guarded role-transfer, not retirement-only non-bridge and not raw identity.",
        },
        "hard_limits": [
            "No raw identity K_legacy_ont == K_strict_gate is claimed.",
            "No claim that the legacy kernel already contains all strict nadsoliton characteristics is made.",
            "No legacy physical-role transfer onto K_strict_gate is claimed before the role-transfer audit.",
            "No beta_tors -> chi_11 quantization/collapse theorem is asserted.",
            "No QW-2191 selector discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    summary = payload["legacy_kernel_intermediate_bridge_summary"]
    lines = [
        "# Legacy kernel intermediate bridge guardrail certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "The legacy kernel is restored as an intermediate incomplete bridge kernel.",
        "The strict kernel is treated as a completed/enriched continuation only",
        "through explicit completion evidence.  Role transfer remains a separate",
        "post-bridge audit obligation.",
        "",
        "## Summary",
        "",
    ]
    for key, value in summary.items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Closure-plan implications", ""])
    for implication in payload["closure_plan_implications"]:
        lines.append(f"- {implication}")
    lines.extend(["", "## Hard limits", ""])
    for limit in payload["hard_limits"]:
        lines.append(f"- {limit}")
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    payload = build_payload()
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_markdown(payload)
    print(json.dumps(payload["legacy_kernel_intermediate_bridge_summary"], indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
