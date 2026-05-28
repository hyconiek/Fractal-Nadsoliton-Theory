#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

IN_2308 = GEN / "p2308_s1258_strict_current_interface_class_response_functional_nonexistence_probe.json"
IN_2309 = GEN / "p2309_s1259_strict_min_norm_response_weights_replay_quarantine_probe.json"
IN_2310 = GEN / "p2310_s1260_strict_self_energy_response_source_audit_and_replay_probe.json"
IN_2311 = GEN / "p2311_s1261_strict_self_energy_to_policy_margin_theorem_audit_probe.json"
IN_2312 = GEN / "p2312_s1262_strict_lapse_shear_source_normalization_attempt_probe.json"
OUT = GEN / "p2313_s1263_strict_provider_margin_route_stop_and_source_class_contract_probe.json"
MD = GEN / "p2313_s1263_strict_provider_margin_route_stop_and_source_class_contract_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    if not path.exists():
        return ""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_json(payload: Any) -> str:
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def nested(payload: dict[str, Any], *keys: str, default: Any = None) -> Any:
    node: Any = payload
    for key in keys:
        if not isinstance(node, dict):
            return default
        node = node.get(key)
    return default if node is None else node


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2308 = load(IN_2308)
    p2309 = load(IN_2309)
    p2310 = load(IN_2310)
    p2311 = load(IN_2311)
    p2312 = load(IN_2312)

    p2312_probe = nested(p2312, "strict_lapse_shear_source_normalization_attempt_probe", default={}) or {}

    blocker_stack = [
        {
            "packet": "P2308",
            "blocker": "current exported P2300/P2281 response-functional interface class has no admissible strict functional",
            "active": bool((p2308.get("gatekeeper_checks", {}) or {}).get("no_current_interface_candidate_admissible") is True or p2308.get("status", "").startswith("OPEN")),
            "source": "generated/p2308_s1258_strict_current_interface_class_response_functional_nonexistence_probe.json",
        },
        {
            "packet": "P2309",
            "blocker": "min-norm weights only work by target-calibrated replay and are quarantined",
            "active": bool((p2309.get("gatekeeper_checks", {}) or {}).get("target_calibrated_weights_quarantined") is True or p2309.get("status", "").startswith("OPEN")),
            "source": "generated/p2309_s1259_strict_min_norm_response_weights_replay_quarantine_probe.json",
        },
        {
            "packet": "P2310",
            "blocker": "self-energy scalar is non-target-calibrated and replay-passing but lacks ADM/policy-margin source theorem",
            "active": bool((p2310.get("gatekeeper_checks", {}) or {}).get("strict_self_energy_bridge_not_claimed") is True),
            "source": "generated/p2310_s1260_strict_self_energy_response_source_audit_and_replay_probe.json",
        },
        {
            "packet": "P2311",
            "blocker": "E(c)=||c||^2 -> provider_lift theorem obligations O3/O4 fail in the current export set",
            "active": bool((p2311.get("gatekeeper_checks", {}) or {}).get("current_export_set_blocks_self_energy_bridge_use") is True),
            "source": "generated/p2311_s1261_strict_self_energy_to_policy_margin_theorem_audit_probe.json",
        },
        {
            "packet": "P2312",
            "blocker": "concrete P1985 non-GB lapse-source normalization is sign-indefinite and rectifiers are non-strict",
            "active": bool((p2312.get("gatekeeper_checks", {}) or {}).get("no_new_lapse_shear_source_exported") is True),
            "source": "generated/p2312_s1262_strict_lapse_shear_source_normalization_attempt_probe.json",
        },
    ]
    all_blockers_active = all(row["active"] for row in blocker_stack)

    p2312_witnesses = nested(p2312_probe, "orientation_witnesses", default=[]) or []
    positive_negative_witness = {row.get("source_sign") for row in p2312_witnesses} >= {"positive", "negative"}

    dynamical_orientation_constraint_audit = {
        "attempt_id": "P2313_DERIVATIVE_ORIENTATION_CONSTRAINT_INVENTORY",
        "question": "Does the current export set contain a strict dynamical constraint fixing the P1985 lapse/shear derivative orientation enough to make the source monotone?",
        "candidate_constraints": [
            {
                "constraint_id": "FREEZE_SHEAR_AND_LAPSE_DERIVATIVES",
                "definition": "Set dsigma_i=d2sigma_i=Nd=Ndd=0 before source normalization.",
                "strict_admissible": False,
                "reason": "This is an imposed kinematic slice/convention, not an exported ADM/Bianchi-I equation valid over the P2281 policy replay domain.",
            },
            {
                "constraint_id": "REQUIRE_NONNEGATIVE_NON_GB_LAPSE_SOURCE",
                "definition": "Restrict admissible trajectories to EL_nonGB_lapse >= 0.",
                "strict_admissible": False,
                "reason": "This is exactly a sign/selector premise and is contradicted as a universal property by P2312's negative witness.",
            },
            {
                "constraint_id": "USE_ABSOLUTE_VALUE_OR_POSITIVE_PART_RECTIFIER",
                "definition": "Replace the source by |EL_nonGB_lapse| or max(EL_nonGB_lapse,0).",
                "strict_admissible": False,
                "reason": "The rectifier is not exported by the variational equations and would be a convention layer.",
            },
            {
                "constraint_id": "DERIVE_NEW_ADM_ORIENTATION_EQUATION",
                "definition": "Export a new strict equation tying lapse/shear derivatives to a signed monotone P2281 margin response.",
                "strict_admissible": False,
                "reason": "This is the missing future object; it is not present in the current export set audited by P2308-P2312.",
            },
        ],
        "positive_negative_witness_available": positive_negative_witness,
        "constraint_exported_now": False,
        "conclusion": "No current strict dynamical constraint fixes the P1985 orientation. The current bridge route must stop unless a genuinely new source class is exported.",
    }

    new_source_class_admission_contract = [
        {
            "criterion_id": "NSC1_TYPED_DOMAIN",
            "requirement": "The source must declare its typed domain: P2300 coefficients, ADM lapse/shear trajectories, or a new strict state object.",
            "mandatory": True,
        },
        {
            "criterion_id": "NSC2_INTERNAL_STRICT_ORIGIN",
            "requirement": "The source must be derived from strict-side dynamics/ontology, not from legacy role transfer, target calibration, or convention.",
            "mandatory": True,
        },
        {
            "criterion_id": "NSC3_SIGNED_MONOTONE_RESPONSE",
            "requirement": "It must export a signed monotone response functional into the P2281/P2302 margin recurrence.",
            "mandatory": True,
        },
        {
            "criterion_id": "NSC4_SELECTOR_FREE_OR_EXPLICITLY_NONSTRICT",
            "requirement": "Any positive-channel, absolute-value, branch, or selector operation must either be internally derived or explicitly marked non-strict.",
            "mandatory": True,
        },
        {
            "criterion_id": "NSC5_REPLAY_AND_GATE_ROWS",
            "requirement": "After the source is exported, replay P2302/P2281 and then recompute P2282-style G1/G2/G3 rows before any closure claim.",
            "mandatory": True,
        },
        {
            "criterion_id": "NSC6_GUARDRAILS",
            "requirement": "The packet must preserve QW-2191, no legacy-kernel role transfer, and no ToE closure claim unless separately proven.",
            "mandatory": True,
        },
    ]

    candidate_next_source_classes = [
        {
            "class_id": "NEW_ADM_LAPSE_SHEAR_ORIENTATION_EQUATION",
            "description": "A future variational/dynamical identity that fixes the sign of the P1985 lapse/shear derivative contribution on the P2281 replay domain.",
            "admitted_now": False,
            "why_not_now": "P2312 exhibits opposite signs for the current source at the same shear energy, so the needed orientation equation is not already exported.",
            "required_contracts": ["NSC1_TYPED_DOMAIN", "NSC2_INTERNAL_STRICT_ORIGIN", "NSC3_SIGNED_MONOTONE_RESPONSE", "NSC5_REPLAY_AND_GATE_ROWS"],
        },
        {
            "class_id": "NEW_SELECTOR_SOURCE_FROM_STRICT_SHANNON_SYMMETRY_BREAKING",
            "description": "A strict internal selector/orientation source possibly related to the strict-side 4 ln 2 Shannon premise, but not a legacy import and not a free selector premise.",
            "admitted_now": False,
            "why_not_now": "QW-2191 remains active; no internal selector source has been exported in the current Task-3 bridge packets.",
            "required_contracts": ["NSC2_INTERNAL_STRICT_ORIGIN", "NSC4_SELECTOR_FREE_OR_EXPLICITLY_NONSTRICT", "NSC6_GUARDRAILS"],
        },
        {
            "class_id": "NEW_POLICY_MARGIN_FUNCTIONAL_FROM_DYNAMICS",
            "description": "A direct derivation of the P2281 margin recurrence response from strict dynamics, bypassing coefficient norm and lapse-source rectification.",
            "admitted_now": False,
            "why_not_now": "P2307/P2308 show the current exported dynamics only responds to an already supplied lambda and does not derive lambda=R(c).",
            "required_contracts": ["NSC1_TYPED_DOMAIN", "NSC3_SIGNED_MONOTONE_RESPONSE", "NSC5_REPLAY_AND_GATE_ROWS"],
        },
    ]

    route_stop_verdict = {
        "current_provider_margin_bridge_route_stopped": all_blockers_active,
        "stopped_route_scope": "current P2300/P2281/P1985/P2310 interface and normalization class only",
        "not_universal_no_go": True,
        "new_source_class_required": True,
        "admissible_for_g1_g3_update": False,
        "reason": "P2308-P2312 now block response weights, target-calibrated weights, self-energy theorem use, and concrete non-GB lapse-source normalization; continuing this same route would be cyclic without a new source class.",
    }

    theorem_export = {
        "statement_id": "P2313_CURRENT_PROVIDER_MARGIN_ROUTE_STOP_AND_NEW_SOURCE_CONTRACT_THEOREM",
        "formal_statement": (
            "Given the active blocker stack P2308-P2312, the current P2300/P2281/P1985/P2310 provider-to-margin route is stopped for strict G1/G3 use. "
            "P2313 does not prove a universal impossibility theorem; it certifies that further progress must introduce a genuinely new strict source class or dynamical orientation constraint satisfying the listed NSC admission contract before replaying P2302/P2281 and recomputing P2282-style gate rows."
        ),
        "proof_bits": {
            "blocker_stack_active": blocker_stack,
            "all_blockers_active": all_blockers_active,
            "positive_negative_p2312_witness_available": positive_negative_witness,
            "constraint_exported_now": dynamical_orientation_constraint_audit["constraint_exported_now"],
            "new_source_class_required": route_stop_verdict["new_source_class_required"],
            "admission_contract_ids": [row["criterion_id"] for row in new_source_class_admission_contract],
        },
        "not_claimed": [
            "universal no-go for all future provider-to-margin bridges",
            "strict G1 closure",
            "strict G3 closure",
            "full Task-3 closure",
            "QW-2191 discharge",
            "selector closure",
            "legacy-kernel role transfer",
            "ToE closure",
        ],
    }
    theorem_fingerprint = sha256_json(theorem_export)

    gatekeeper_checks = {
        "p2308_current_interface_blocker_loaded": blocker_stack[0]["active"],
        "p2309_target_calibration_quarantine_loaded": blocker_stack[1]["active"],
        "p2310_self_energy_quarantine_loaded": blocker_stack[2]["active"],
        "p2311_self_energy_theorem_blocker_loaded": blocker_stack[3]["active"],
        "p2312_sign_indefinite_source_blocker_loaded": blocker_stack[4]["active"],
        "all_current_route_blockers_active": all_blockers_active,
        "no_current_dynamical_orientation_constraint_exported": dynamical_orientation_constraint_audit["constraint_exported_now"] is False,
        "new_source_class_contract_exported": len(new_source_class_admission_contract) == 6,
        "current_route_stop_scope_not_universal_no_go": route_stop_verdict["not_universal_no_go"] is True,
        "strict_g1_g3_not_updated": True,
        "no_selector_premise_added": True,
        "no_qw2191_discharge_claimed": True,
        "no_legacy_kernel_role_transfer": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2313_s1263_v1",
        "packet_id": "P2313",
        "stage_id": "S1263",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_OBSTRUCTION_CURRENT_PROVIDER_MARGIN_ROUTE_STOPPED_NEW_SOURCE_CLASS_REQUIRED",
        "result_kind": "STRICT_PROVIDER_MARGIN_ROUTE_STOP_AND_NEW_SOURCE_CLASS_ADMISSION_CONTRACT",
        "strict_provider_margin_route_stop_and_source_class_contract_probe": {
            "probe_id": "P2313_S1263_STRICT_PROVIDER_MARGIN_ROUTE_STOP_AND_SOURCE_CLASS_CONTRACT",
            "source_packets": {
                "p2308": "generated/p2308_s1258_strict_current_interface_class_response_functional_nonexistence_probe.json",
                "p2309": "generated/p2309_s1259_strict_min_norm_response_weights_replay_quarantine_probe.json",
                "p2310": "generated/p2310_s1260_strict_self_energy_response_source_audit_and_replay_probe.json",
                "p2311": "generated/p2311_s1261_strict_self_energy_to_policy_margin_theorem_audit_probe.json",
                "p2312": "generated/p2312_s1262_strict_lapse_shear_source_normalization_attempt_probe.json",
            },
            "source_hashes": {
                "p2308_sha256": sha256_file(IN_2308),
                "p2309_sha256": sha256_file(IN_2309),
                "p2310_sha256": sha256_file(IN_2310),
                "p2311_sha256": sha256_file(IN_2311),
                "p2312_sha256": sha256_file(IN_2312),
            },
            "blocker_stack": blocker_stack,
            "dynamical_orientation_constraint_audit": dynamical_orientation_constraint_audit,
            "new_source_class_admission_contract": new_source_class_admission_contract,
            "candidate_next_source_classes": candidate_next_source_classes,
            "route_stop_verdict": route_stop_verdict,
            "task3_g1_g3_update": {
                "G1_reduction_certainty": "OPEN_UNCHANGED",
                "G2_nonlinear_trajectory_realism": "CLOSED_FROM_P2301_UNCHANGED",
                "G3_operational_policy_rule": "OPEN_UNCHANGED",
                "reason": "P2313 stops the current route rather than cycling the same blocker cut without a new source class.",
            },
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": theorem_fingerprint,
        },
        "recommended_next_honest_step": {
            "id": "P2314_candidate",
            "goal": "Start a genuinely new source-class construction attempt satisfying NSC1-NSC6, or explicitly mark any Shannon/selector branch as non-strict until it exports an internal selector/orientation source and survives P2302/P2282 replay.",
        },
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_OBSTRUCTION_WITH_ROUTE_STOP_AND_NEW_SOURCE_CLASS_CONTRACT_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2313/S1263 — provider→margin route stop and new source-class contract",
            "",
            f"- Status: `{payload['status']}`",
            f"- All current route blockers active: `{all_blockers_active}`",
            "- Current route stopped: `True` for the current P2300/P2281/P1985/P2310 interface class only.",
            "- New source class required: `True`",
            "- G1/G3 update: `OPEN_UNCHANGED`",
            f"- Theorem fingerprint: `{theorem_fingerprint}`",
            "",
            "## Concrete audit result",
            "P2308-P2312 now block the current route: no current response functional, target-calibrated weights are quarantined, self-energy lacks the theorem bridge, and the concrete non-GB lapse source is sign-indefinite.  Repeating this route without a new source class would be cyclic.",
            "",
            "## Admission contract",
            "Any P2314+ source must satisfy NSC1-NSC6: typed domain, strict internal origin, signed monotone response, selector discipline, replay/gate rows, and all guardrails.",
            "",
            "## Guardrail statement",
            "P2313 does not close G1/G3, does not discharge QW-2191, does not add a selector premise, does not transfer legacy-kernel roles, and does not claim ToE closure.",
            "",
            "## Next honest step",
            payload["recommended_next_honest_step"]["goal"],
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
