#!/usr/bin/env python3
"""Scratch probe: audit the proposed beta_tors -> chi_11 torsion-remnant opinion.

The user asked whether an AI opinion is justified: that legacy continuous torsion
``beta_tors`` is physically the source/topological remnant of the strict-side
``chi_11`` unit-character bit, and that this supplies the missing bridge into
``eta=9/5`` / ``q^5=256/243`` / the balanced ledger.

This probe does not try to prove that bridge.  It performs a finite, repo-local
claim audit against the guardrails and the existing chi_11 reports.  The honest
answer is intentionally split:

* the intuition is a plausible candidate direction to formalize;
* the current repo does not export the bridge theorem claimed by the opinion;
* S2 explicitly retires the legacy kernel as archival, so beta_tors cannot be
  silently promoted into a strict selector source;
* chi_11 remains an extra unit-axis/character premise unless a new strict source
  theorem is added.

No false pass: no legacy-to-strict kernel identity, no beta_tors -> chi_11
collapse theorem, no selector discharge, and no ToE closure are claimed.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FAR = ROOT / "fundamental_action_reconstruction"
OUT_JSON = HERE / "bridge_legacy_torsion_chi11_opinion_audit_report.json"
OUT_MD = HERE / "bridge_legacy_torsion_chi11_opinion_audit_report.md"

K1 = FAR / "K1_LEGACY_ONTOLOGICAL_KERNEL_VS_STRICT_GATE_KERNEL_SPLIT_NOTE.md"
S2 = FAR / "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md"
COHOMOLOGY_REPORT = HERE / "bridge_strict_alpha_unit_character_cohomology_obstruction_audit_report.json"
META_CHARACTER_REPORT = HERE / "bridge_strict_alpha_unit_character_meta_symmetry_uniqueness_obstruction_report.json"
UNORIENTED_BRIDGE_REPORT = HERE / "bridge_legacy_strict_unoriented_equation_without_selector_bit_audit_report.json"
REYNOLDS_REPORT = HERE / "bridge_strict_alpha_reynolds_annihilator_chi11_matrix_certificate_report.json"
PREMISE_REPORT = HERE / "bridge_strict_alpha_chi11_premise_dependency_lattice_audit_report.json"
CONVEX_SELECTOR_REPORT = HERE / "bridge_strict_alpha_symmetric_convex_selector_majorization_certificate_report.json"
FRACTAL_DIMENSION_REPORT = HERE / "bridge_strict_alpha_fractal_dimension_report.json"

TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
BALANCED_LEDGER = [2, 2, 2, 1, 1]


def read_text(path: Path) -> str:
    if not path.exists():
        raise FileNotFoundError(f"missing required text file: {path}")
    return path.read_text(encoding="utf-8")


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing required report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def contains_all(text: str, needles: list[str]) -> bool:
    return all(needle in text for needle in needles)


def claim_rows(k1_text: str, s2_text: str, reports: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    cohomology = reports["cohomology"]
    meta = reports["meta_character"]
    bridge = reports["unoriented_bridge"]
    reynolds = reports["reynolds"]
    premise = reports["premise"]
    convex = reports["convex_selector"]
    fractal = reports["fractal_dimension"]

    return [
        {
            "claim_id": "legacy_beta_tors_exists_as_torsion_damping",
            "opinion_claim": "The legacy kernel contains beta_tors and treats it as torsion/damping data.",
            "verdict": "SUPPORTED_AS_LEGACY_ARCHIVAL_FACT",
            "evidence": [
                contains_all(k1_text, ["K_legacy_ont(d) := alpha_geo * cos", "`beta_tors = 0.01` as hierarchy/torsion damping datum"]),
                "K1 records beta_tors as a legacy hierarchy/torsion damping datum, not as a current strict selector source.",
            ],
            "allowed_use": "May motivate an archival/candidate bridge question if all role-transfer hypotheses are explicit.",
        },
        {
            "claim_id": "strict_kernel_has_different_beta_eta_form",
            "opinion_claim": "The strict kernel has beta=1.0 and eta=1.8 rather than the old beta_tors*d denominator.",
            "verdict": "SUPPORTED_STRUCTURAL_SPLIT_ONLY",
            "evidence": [
                contains_all(k1_text, ["K_strict_gate(d) := cos", "beta = 1.0", "eta = 1.8"]),
                "The split is documented, but the split alone does not identify beta_tors with any strict character bit.",
            ],
            "allowed_use": "May compare forms; may not silently substitute or transfer legacy roles.",
        },
        {
            "claim_id": "torsion_implies_parity_chirality_breaking",
            "opinion_claim": "Because the legacy word torsion appears, beta_tors physically breaks parity/chirality.",
            "verdict": "NOT_EXPORTED_BY_CURRENT_REPO",
            "evidence": [
                "K1 says hierarchy/torsion damping datum; this audit found no current strict theorem equating beta_tors with a parity-breaking unit character.",
                not contains_all(k1_text, ["beta_tors", "chi_11"]),
            ],
            "allowed_use": "At most a conjectural physical analogy until a formal orientation/chirality map is supplied.",
        },
        {
            "claim_id": "chi11_is_missing_unit_character_bit",
            "opinion_claim": "The strict branch problem contains a missing one-bit chi_11/unit-axis datum.",
            "verdict": "SUPPORTED_IN_REDUCED_ALGEBRAIC_AUDITS",
            "evidence": [
                cohomology["interpretation"]["direct_result"],
                meta["interpretation"]["direct_result"],
                reynolds["matrix_certificate"]["chi11_numerator_rank"] == 13,
            ],
            "allowed_use": "Can state as a reduced-support/character-premise result, not as a strict source theorem.",
        },
        {
            "claim_id": "beta_tors_collapses_to_chi11",
            "opinion_claim": "Continuous legacy beta_tors quantizes/collapses to the strict chi_11 bit.",
            "verdict": "CANDIDATE_BRIDGE_HYPOTHESIS_NOT_THEOREM",
            "evidence": [
                bridge["interpretation"]["what_cannot_be_done_now"],
                premise["interpretation"]["honest_negative"],
                contains_all(s2_text, ["legacy kernel is **not** a discarded dead end", "intermediate kernel on the path toward identifying the strict kernel"]),
            ],
            "allowed_use": "Only as an explicitly labelled new theorem target, not as a repo conclusion.",
        },
        {
            "claim_id": "legacy_kernel_bridge_is_now_the_main_bridge",
            "opinion_claim": "This observation supplies the legacy->strict bridge the repo has been seeking.",
            "verdict": "REACTIVATED_AS_BRIDGE_TARGET_BUT_NOT_PROVEN_AS_FACT",
            "evidence": [
                contains_all(s2_text, ["legacy -> strict` completion bridge", "role-transfer audit after bridge completion"]),
                "S2 restores bridge-completion priority and keeps legacy role-transfer forbidden until an explicit post-bridge role-transfer audit.",
            ],
            "allowed_use": "Can motivate the restored legacy->strict completion bridge, but the bridge still requires an explicit map, compression account, selector/source theorem, and role-transfer audit.",
        },
        {
            "claim_id": "chi11_triggers_eta_and_balanced_ledger",
            "opinion_claim": "Supplying chi_11 directly triggers eta=9/5, q^5=256/243, and the [2,2,2,1,1] allocation.",
            "verdict": "OVERSTATED_CONDITIONAL_DOWNSTREAM_ONLY",
            "evidence": [
                convex["target_eta_9_5_certificate"]["balanced_ledger"] == BALANCED_LEDGER,
                convex["target_eta_9_5_certificate"]["q_power_product"] == TARGET_Q_POWER,
                fractal["hard_limits"],
            ],
            "allowed_use": "May cite conditional selector/ledger certificates; may not claim a strict physical trigger theorem.",
        },
        {
            "claim_id": "full_aut_invariant_source_exports_chi11",
            "opinion_claim": "The missing bit was already present physically, so full-Aut strict data can export chi_11 once read correctly.",
            "verdict": "REFUTED_FOR_CURRENT_FULL_AUT_SUPPORT_DATA",
            "evidence": [
                reynolds["matrix_certificate"]["reynolds_times_chi11_numerator_is_zero_matrix"],
                reynolds["interpretation"]["honest_negative"],
                premise["minimal_premise_antichains"]["strict_full_aut_internal_chi11_polarity_source"] == [],
            ],
            "allowed_use": "Requires a new non-full-Aut or strict internal orientation source; current full-Aut support data annihilates it.",
        },
    ]


def build_payload() -> dict[str, Any]:
    k1_text = read_text(K1)
    s2_text = read_text(S2)
    reports = {
        "cohomology": load_json(COHOMOLOGY_REPORT),
        "meta_character": load_json(META_CHARACTER_REPORT),
        "unoriented_bridge": load_json(UNORIENTED_BRIDGE_REPORT),
        "reynolds": load_json(REYNOLDS_REPORT),
        "premise": load_json(PREMISE_REPORT),
        "convex_selector": load_json(CONVEX_SELECTOR_REPORT),
        "fractal_dimension": load_json(FRACTAL_DIMENSION_REPORT),
    }
    rows = claim_rows(k1_text, s2_text, reports)
    verdict_counts: dict[str, int] = {}
    for row in rows:
        verdict_counts[row["verdict"]] = verdict_counts.get(row["verdict"], 0) + 1

    bridge_hypothesis_requirements = [
        {
            "requirement": "orientation_map",
            "needed_theorem": "Define a strict or explicitly bridge-level map from legacy torsion orientation to a U12 character, including why the kernel is {1,11} rather than chi_5 or chi_7.",
            "current_status": reports["meta_character"]["interpretation"]["negative_result"],
        },
        {
            "requirement": "role_transfer_control",
            "needed_theorem": "Prove that using beta_tors as a selector source is not a silent transfer of legacy physical roles onto K_strict_gate.",
            "current_status": "S2 restores K_legacy_ont as an intermediate bridge kernel and requires a post-bridge role-transfer audit before any legacy physical role can move onto K_strict_gate.",
        },
        {
            "requirement": "eta_pipeline_link",
            "needed_theorem": "Connect the character bit to eta=9/5 and q^5=256/243 through already certified selector premises without treating conditional certificates as strict closure.",
            "current_status": "Existing eta/ledger reports are conditional and keep hard limits against exact strict eta theorem/selector closure.",
        },
        {
            "requirement": "full_aut_obstruction_escape",
            "needed_theorem": "Explain how the proposed source avoids the Reynolds annihilator and premise-lattice no-source obstruction for full-Aut invariant data.",
            "current_status": reports["reynolds"]["interpretation"]["honest_negative"],
        },
    ]

    return {
        "result_kind": "SCRATCH_LEGACY_TORSION_CHI11_OPINION_AUDIT__CANDIDATE_NOT_THEOREM",
        "status": "ai-opinion-classified-as-interesting-but-overstated-bridge-hypothesis",
        "finite_model_cross_checks": {
            "ring": "Z_12",
            "support_count": reports["cohomology"]["finite_model"]["support_count"],
            "automorphism_units": reports["cohomology"]["finite_model"]["automorphism_units"],
            "chi11_character_count": reports["cohomology"]["character_table_audit"]["character_count"],
            "chi11_rank_from_reynolds_probe": reports["reynolds"]["matrix_certificate"]["chi11_numerator_rank"],
            "reynolds_rank_from_reynolds_probe": reports["reynolds"]["matrix_certificate"]["reynolds_numerator_rank"],
            "strict_full_aut_internal_chi11_antichain": reports["premise"]["minimal_premise_antichains"]["strict_full_aut_internal_chi11_polarity_source"],
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
            "balanced_ledger": BALANCED_LEDGER,
        },
        "opinion_audit": {
            "overall_verdict": "PARTLY_USEFUL_HEURISTIC_BUT_NOT_A_CURRENT_REPO_THEOREM",
            "verdict_counts": verdict_counts,
            "claims": rows,
        },
        "bridge_hypothesis_requirements": bridge_hypothesis_requirements,
        "recommended_next_honest_step": {
            "if_pursuing_this_intuition": "Write a deliberately non-strict theorem-target spec: beta_tors_orientation_source -> U12 character, then try to prove or refute it against the full-Aut Reynolds obstruction and meta-character ambiguity.",
            "if_following_current_S2_priority": "Treat beta_tors as a candidate bridge/source datum and continue bridge-completion work without importing legacy torsion roles before a role-transfer audit.",
        },
        "exact_proof_certificate": {
            "supported_part": "K1 supports beta_tors as a legacy torsion/damping datum; existing chi_11 audits support a missing nontrivial unit-character premise in reduced models.",
            "unsupported_part": "No loaded report exports beta_tors -> chi_11, torsion -> parity-breaking unit character, or chi_11 -> eta=9/5 as an unconditional strict theorem.",
            "current_obstruction": "Full-Aut Reynolds averaging annihilates chi_11 and the premise lattice has no minimal set for a strict full-Aut internal chi_11 polarity source.",
            "s2_boundary": "S2 restores the legacy kernel as an intermediate bridge kernel; any beta_tors bridge is now an explicit candidate bridge/source theorem target, not an asserted current closure.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself remains the primordial information in a solitonic state.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted.",
            "No beta_tors -> chi_11 quantization/collapse theorem is asserted.",
            "No legacy physical-role transfer onto K_strict_gate is used.",
            "No theorem derives chi_11, shell-labels, unit-axis bit, exact-cover clauses, or cardinality 5 from strict geometry.",
            "No theorem derives eta=9/5 or q^5=256/243 from legacy torsion.",
            "No QW-2191 discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> str:
    checks = payload["finite_model_cross_checks"]
    audit = payload["opinion_audit"]
    lines = [
        "# Legacy torsion -> chi_11 opinion audit",
        "",
        f"Status: `{payload['status']}`",
        f"Overall verdict: `{audit['overall_verdict']}`",
        "",
        "## Cross-checks",
        "",
        f"- Ring/supports: `{checks['ring']}` / `{checks['support_count']}`",
        f"- Automorphism units: `{checks['automorphism_units']}`",
        f"- Character count: `{checks['chi11_character_count']}`",
        f"- Reynolds rank / chi_11 rank: `{checks['reynolds_rank_from_reynolds_probe']}` / `{checks['chi11_rank_from_reynolds_probe']}`",
        f"- Strict full-Aut internal chi_11 antichain: `{checks['strict_full_aut_internal_chi11_antichain']}`",
        f"- Target replay kept conditional: `q^5={checks['target_q_power']}`, eta `{checks['target_eta']}`, ledger `{checks['balanced_ledger']}`",
        "",
        "## Claim verdicts",
        "",
    ]
    for row in audit["claims"]:
        lines.extend(
            [
                f"### {row['claim_id']}",
                f"- Opinion claim: {row['opinion_claim']}",
                f"- Verdict: `{row['verdict']}`",
                f"- Allowed use: {row['allowed_use']}",
                "",
            ]
        )
    lines.extend(["## Bridge-hypothesis requirements", ""])
    for row in payload["bridge_hypothesis_requirements"]:
        lines.extend(
            [
                f"- `{row['requirement']}`: {row['needed_theorem']}",
                f"  - Current status: {row['current_status']}",
            ]
        )
    lines.extend(["", "## Proof certificate", ""])
    for key, value in payload["exact_proof_certificate"].items():
        lines.append(f"- `{key}`: {value}")
    lines.extend(["", "## Hard limits", ""])
    lines.extend(f"- {item}" for item in payload["hard_limits"])
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    payload = build_payload()
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    OUT_MD.write_text(write_markdown(payload), encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
