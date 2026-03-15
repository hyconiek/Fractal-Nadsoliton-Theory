#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

DEFAULT_PROVIDER = GENERATED / "f447_current_strict_t169_qw2122_scalar_to_t168_per_site_value_provider_strict_derived_v1.json"

OUT_JSON = GENERATED / "p448_current_strict_f447_t169_provider_provenance_audit_probe.json"
OUT_SUMMARY = GENERATED / "p448_current_strict_f447_t169_provider_provenance_audit_probe_summary.json"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="P448: audit provenance / strict-derived classification readiness of the F447 T169→T168 per-site provider."
    )
    p.add_argument(
        "--provider",
        default=str(DEFAULT_PROVIDER),
        help="Path to the provider JSON artifact to audit (default: generated F447 provider).",
    )
    p.add_argument(
        "--out-json",
        default=str(OUT_JSON),
        help="Output artifact JSON path.",
    )
    p.add_argument(
        "--out-summary",
        default=str(OUT_SUMMARY),
        help="Output summary JSON path.",
    )
    return p.parse_args()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def is_number(x: Any) -> bool:
    return isinstance(x, (int, float))


def is_finite_number(x: Any) -> bool:
    return isinstance(x, (int, float)) and (x == x) and (x not in (float("inf"), float("-inf")))


def list_len_12_finite(xs: Any) -> bool:
    return (
        isinstance(xs, list)
        and len(xs) == 12
        and all(is_finite_number(x) for x in xs)
    )


def require_path(p: str | None) -> tuple[bool, str | None]:
    if not isinstance(p, str) or not p.strip():
        return False, None
    path = Path(p)
    if not path.is_absolute():
        path = REPO / path
    return path.exists(), str(path.relative_to(REPO)) if path.exists() else str(path)


def audit_provider(provider: dict[str, Any]) -> dict[str, Any]:
    schema_ok = True
    schema_errors: list[str] = []

    classification = str(provider.get("classification") or "")
    status = str(provider.get("status") or "")

    # Required strict inputs for the claimed T169 constrained lift.
    strict_inputs = provider.get("strict_inputs")
    strict_input_checks: dict[str, Any] = {}
    if not isinstance(strict_inputs, dict):
        schema_ok = False
        schema_errors.append("missing strict_inputs dict")
    else:
        for key in ("QW-2122", "alpha_geo_strict_derived_v1", "R14", "R15"):
            exists, resolved = require_path(strict_inputs.get(key))
            strict_input_checks[key] = {"declared_path": strict_inputs.get(key), "resolved_exists": bool(exists), "resolved": resolved}
            if not exists:
                schema_ok = False
                schema_errors.append(f"strict_inputs.{key} missing or path does not exist")

    # Required arrays for T168 consumption.
    vpsi_ok = list_len_12_finite(provider.get("vpsi"))
    g4_ok = list_len_12_finite(provider.get("g4"))
    g6_ok = list_len_12_finite(provider.get("g6"))
    if not vpsi_ok:
        schema_ok = False
        schema_errors.append("missing or invalid vpsi (need 12 finite numbers)")
    if not g4_ok:
        schema_ok = False
        schema_errors.append("missing or invalid g4 (need 12 finite numbers)")
    if not g6_ok:
        schema_ok = False
        schema_errors.append("missing or invalid g6 (need 12 finite numbers)")

    # Required sigma-six-sum keys for T167/P434 consumption.
    sigma_keys = [f"Sigma_psi{i}_psi{i+6}" for i in range(6)]
    sigma_ok_by_key = {k: is_finite_number(provider.get(k)) for k in sigma_keys}
    if not all(sigma_ok_by_key.values()):
        schema_ok = False
        schema_errors.append("missing or invalid sigma-six-sum keys (Sigma_psi{k}_psi{k+6}, k=0..5)")

    # Extract the nontrivial selection premises from provider.lift_definition.
    lift_definition = provider.get("lift_definition")
    if not isinstance(lift_definition, dict):
        schema_ok = False
        schema_errors.append("missing lift_definition dict (selection rule must be explicit)")
        lift_definition = {}

    declared_reference = None
    ref = provider.get("reference")
    if isinstance(ref, dict):
        declared_reference = ref.get("r_ordpow")

    ref_exists, ref_resolved = require_path(declared_reference)

    nontrivial_premises: list[str] = []
    if isinstance(declared_reference, str):
        nontrivial_premises.append("choice_of_reference_distribution_shape: r_ordpow(x) ∝ ord_Z12(x)^(-alpha_geo)")
    if "vpsi_signs" in lift_definition:
        nontrivial_premises.append("sign_selection: minimize mixing energy over {±1}^12 and fix global Z2 by s0=+1")
    if "g6" in lift_definition:
        nontrivial_premises.append("sextic_mapping: g6_i := 0")
    if "g4_uniform_minimal_matching" in lift_definition:
        nontrivial_premises.append("quartic_mapping: uniform g4_i matched to scalar lambda_psi along lifted vacuum ray")

    # Compare against P446 “real discharge blueprint” checklist (conservative).
    # We do not parse theorems; we only check whether the provider references any theorem-level uniqueness/existence proof.
    theorem_refs: Any = provider.get("theorems")
    if theorem_refs is None and isinstance(provider.get("provenance"), dict):
        theorem_refs = (provider.get("provenance") or {}).get("theorems")
    theorem_refs_list = theorem_refs if isinstance(theorem_refs, list) else []
    has_theorem_refs = len(theorem_refs_list) > 0

    required_theorems = {
        "N479_CURRENT_FIRST_STRICT_Z12_ELEMENT_ORDER_REFERENCE_AUT_Z12_INVARIANCE_NO_MARKED_DIRECTION_THEOREM",
        "N483_CURRENT_FIRST_STRICT_T169_POWERLAW_ELEMENT_ORDER_CONSTRAINED_LIFT_EXISTENCE_UNIQUENESS_THEOREM",
    }
    required_present = {t: (t in theorem_refs_list) for t in sorted(required_theorems)}
    theorem_level_pass = classification.lower() == "strict_derived" and all(required_present.values())

    # Conservative recommendation logic:
    # - If classification claims strict_derived but there are no theorem-level references and the lift introduces explicit mapping premises,
    #   recommend downgrading to strict_extension_only until a theorem-level admissibility/uniqueness story is exported.
    recommendation = "review_required"
    recommendation_reason = "Provider is computable and explicit; classification requires manual review against T169/P446 criteria."
    if theorem_level_pass:
        recommendation = "keep_strict_derived"
        recommendation_reason = (
            "Provider declares strict_derived and references the required theorem-level existence/uniqueness + no-direction-slot anchors "
            "(N483/N479)."
        )
    elif classification.lower() == "strict_derived" and (not has_theorem_refs) and nontrivial_premises:
        recommendation = "downgrade_to_strict_extension_only"
        recommendation_reason = (
            "Provider declares strict_derived but does not reference a theorem-level admissibility/existence+uniqueness proof for the constrained lift; "
            "mapping choices remain explicit premises (P446 blueprint)."
        )
    elif classification.lower() == "strict_extension_only":
        recommendation = "keep_strict_extension_only"
        recommendation_reason = "Provider is already labeled strict_extension_only; keep separation from strict core."

    return {
        "schema_ok": schema_ok,
        "schema_errors": schema_errors,
        "provider_status": status,
        "provider_classification": classification,
        "strict_input_checks": strict_input_checks,
        "reference": {
            "declared_r_ordpow": declared_reference,
            "resolved_exists": bool(ref_exists),
            "resolved": ref_resolved,
        },
        "t168_arrays_present": {"vpsi": vpsi_ok, "g4": g4_ok, "g6": g6_ok},
        "t167_sigma_present": sigma_ok_by_key,
        "declared_lift_definition_keys": sorted(list(lift_definition.keys())),
        "nontrivial_selection_premises": nontrivial_premises,
        "theorem_references_present": bool(has_theorem_refs),
        "theorem_references": theorem_refs_list,
        "required_theorems_present": required_present,
        "theorem_level_pass": bool(theorem_level_pass),
        "recommendation": recommendation,
        "recommendation_reason": recommendation_reason,
    }


def main() -> None:
    GENERATED.mkdir(exist_ok=True)
    args = parse_args()

    provider_path = Path(args.provider)
    if not provider_path.is_absolute():
        provider_path = REPO / provider_path

    out_json = Path(args.out_json)
    if not out_json.is_absolute():
        out_json = REPO / out_json

    out_summary = Path(args.out_summary)
    if not out_summary.is_absolute():
        out_summary = REPO / out_summary

    if not provider_path.exists():
        summary = {
            "stage": "P448",
            "status": "NOT_COMPUTABLE_MISSING_PROVIDER",
            "provider_path": str(provider_path.relative_to(REPO)) if provider_path.is_absolute() else str(provider_path),
            "no_false_pass": True,
        }
        out_summary.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(out_summary)
        return

    provider = load_json(provider_path)
    audit = audit_provider(provider)

    artifact = {
        "stage": "P448",
        "goal": "audit_F447_T169_provider_provenance_and_strict_derived_readiness_without_false_promotion",
        "provider": str(provider_path.relative_to(REPO)),
        "audit": audit,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P448",
        "status": "PASS_AUDIT_READY" if audit["schema_ok"] else "FAIL_SCHEMA_ERRORS",
        "provider_classification": audit["provider_classification"],
        "recommendation": audit["recommendation"],
        "theorem_level_pass": bool(audit.get("theorem_level_pass")),
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_summary.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    out_summary.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(out_json)


if __name__ == "__main__":
    main()
