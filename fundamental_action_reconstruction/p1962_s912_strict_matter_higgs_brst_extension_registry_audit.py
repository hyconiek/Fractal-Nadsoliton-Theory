#!/usr/bin/env python3
"""P1962 S912 strict matter/Higgs BRST extension registry audit.

P1961 proved local gauge-sector BRST nilpotency for
SU(3) x SU(2) x U(1).  This executor asks the next narrower question:

    can the current strict exports extend those local rules to the
    matter/Higgs/spinor part of the exported L_total bundle?

The result is an audit, not a global BRST theorem.  It records conditional
matter-module BRST rules and checks whether the repository currently exports
the representation matrices, chiral family registry, Higgs representation,
Yukawa gauge-invariance checks, and anomaly checks needed to make those rules
theorem-grade.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

P1961 = GEN / "p1961_s911_strict_local_nonabelian_brst_differential_and_nilpotency_certificate.json"
P1907 = GEN / "p1907_s857_strict_full_lagrangian_to_eom_witness_matrix_probe.json"
P1837 = GEN / "p1837_s787_strict_full_lagrangian_explicit_anchor_registry_checkpoint.json"
P1855 = GEN / "p1855_s805_strict_b1_gauge_fermion_k5_and_cutkosky_integral_stub_checkpoint.json"
P1856 = GEN / "p1856_s806_strict_b1_representation_and_k5_triangle_seed_checkpoint.json"
P1801 = GEN / "p1801_s751_strict_brst_nilpotency_witness_intake_and_gate_link_checkpoint.json"
P1957 = GEN / "p1957_s907_strict_bv_brst_ghost_sector_nonavailability_theorem.json"

OUT = GEN / "p1962_s912_strict_matter_higgs_brst_extension_registry_audit.json"

IGNORED_FILES = [
    "TOE_FINAL_DOCUMENTATION_RELEASE_8_STRICT_FULL.pdf",
    OUT.name,
    Path(__file__).name,
]


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def digest_path(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bool_from_path(path: Path) -> bool:
    return path.exists()


def artifact_status(obj: dict[str, Any]) -> str:
    return str(obj.get("status", "UNKNOWN_STATUS"))


def representation_seed_audit(p1856: dict[str, Any]) -> dict[str, Any]:
    reps = (
        p1856.get("strict_representation_seed", {})
        .get("fermion_representations_seed", [])
    )
    families = sorted({str(row.get("family", "UNKNOWN")) for row in reps})
    su3_labels = sorted({str(row.get("su3", "MISSING")) for row in reps})
    su2_labels = sorted({str(row.get("su2", "MISSING")) for row in reps})
    hypercharges = [str(row.get("Y", "MISSING")) for row in reps]

    has_field_names = all("field" in row or "field_id" in row for row in reps)
    has_chirality = all("chirality" in row or "chiral_sector" in row for row in reps)
    has_generator_matrices = any(
        "generator_matrices" in row
        or "T_SU3" in row
        or "T_SU2" in row
        or "Y_matrix" in row
        for row in reps
    )
    has_commutator_certificate = "commutator_certificate" in p1856 or "representation_commutator_certificate" in p1856
    has_higgs_representation = any(
        str(row.get("field", row.get("field_id", ""))).lower() in {"h", "higgs", "higgs_h"}
        for row in reps
    )

    return {
        "seed_present": bool(reps),
        "seed_row_count": len(reps),
        "families_seen": families,
        "su3_labels_seen": su3_labels,
        "su2_labels_seen": su2_labels,
        "hypercharges_seen": hypercharges,
        "has_field_names": has_field_names,
        "has_chirality_labels": has_chirality,
        "has_generator_matrices": has_generator_matrices,
        "has_commutator_certificate": has_commutator_certificate,
        "has_higgs_representation": has_higgs_representation,
        "classification": (
            "REPRESENTATION_SEED_ONLY_NOT_THEOREM_GRADE"
            if reps
            else "NO_REPRESENTATION_SEED"
        ),
    }


def ltotal_registry_audit(p1907: dict[str, Any], p1837: dict[str, Any]) -> dict[str, Any]:
    registry = p1907.get("full_lagrangian_term_registry_non_skeleton", {})
    eom_matrix = p1907.get("eom_export_matrix", [])
    eom_statuses = {str(row.get("field", "UNKNOWN")): str(row.get("status", "UNKNOWN")) for row in eom_matrix}
    reduced_anchor = p1837.get("explicit_anchor_registry", {})

    sector_keys = [
        "sm_gauge_sector",
        "sm_fermion_kinetic_sector",
        "higgs_sector",
        "yukawa_sector",
        "gravity_sector",
        "nonminimal_mixed_sector",
    ]
    present_sectors = [key for key in sector_keys if key in registry]

    return {
        "p1907_registry_present": bool(registry),
        "present_sector_keys": present_sectors,
        "missing_sector_keys": [key for key in sector_keys if key not in registry],
        "eom_statuses": eom_statuses,
        "all_eom_concrete": bool(eom_statuses) and all(status == "PASS_EXPLICIT" for status in eom_statuses.values()),
        "p1837_reduced_anchor_present": bool(reduced_anchor.get("reduced_non_skeleton_anchor_present", False)),
        "p1837_promotion_policy": p1837.get("promotion_policy", "UNKNOWN"),
        "classification": "SECTOR_LEVEL_AND_REDUCED_ANCHOR_ONLY_NOT_FULL_NONPROXY_FIELD_REGISTRY",
    }


def requirement_table(
    p1961: dict[str, Any],
    rep_audit: dict[str, Any],
    ltotal_audit: dict[str, Any],
) -> list[dict[str, Any]]:
    p1961_ready = bool(
        p1961.get("ReadinessAfterP1961_strict_v1", {}).get("evaluated_local_gauge_sector_ready", False)
    )
    return [
        {
            "id": "LOCAL_GAUGE_BRST_S2_ZERO",
            "required": "P1961 local gauge-sector s^2=0 for A,c,cbar,B",
            "available": p1961_ready,
            "source": "P1961",
        },
        {
            "id": "FULL_NONPROXY_FIELD_REGISTRY",
            "required": "field-by-field matter/Higgs/spinor registry, not only sector-level L_total text",
            "available": False,
            "source": "P1907/P1837",
            "trace": ltotal_audit["classification"],
        },
        {
            "id": "CHIRAL_FAMILY_REPRESENTATION_REGISTRY",
            "required": "all chiral fermion multiplets with family labels and Higgs multiplet",
            "available": bool(
                rep_audit["has_field_names"]
                and rep_audit["has_chirality_labels"]
                and rep_audit["has_higgs_representation"]
            ),
            "source": "P1856",
            "trace": "P1856 gives a seed table, but no field/chirality-complete Higgs-inclusive registry.",
        },
        {
            "id": "REPRESENTATION_MATRICES",
            "required": "explicit T_SU3, T_SU2, and Y matrices/operators for each multiplet",
            "available": bool(rep_audit["has_generator_matrices"]),
            "source": "P1856",
        },
        {
            "id": "REPRESENTATION_COMMUTATOR_CERTIFICATE",
            "required": "[T_a,T_b]=f_ab^c T_c for every exported representation, plus commuting product factors",
            "available": bool(rep_audit["has_commutator_certificate"]),
            "source": "P1856/P1960",
        },
        {
            "id": "YUKAWA_GAUGE_INVARIANCE_CERTIFICATE",
            "required": "explicit BRST/gauge invariance of qbar H u, qbar H d, lbar H e terms",
            "available": False,
            "source": "not exported",
        },
        {
            "id": "ANOMALY_CANCELLATION_CERTIFICATE",
            "required": "full chiral family anomaly sums in the same representation convention",
            "available": False,
            "source": "P1856 marks this open",
        },
        {
            "id": "GLOBAL_BRST_WITNESS_PACK",
            "required": "P1801 TG2 pack: Q, Q^2, cohomology, ghost consistency, shared freeze",
            "available": False,
            "source": "P1801/P1957",
        },
    ]


def main() -> None:
    p1961 = load_json(P1961)
    p1907 = load_json(P1907)
    p1837 = load_json(P1837)
    p1855 = load_json(P1855)
    p1856 = load_json(P1856)
    p1801 = load_json(P1801)
    p1957 = load_json(P1957)

    rep_audit = representation_seed_audit(p1856)
    ltotal_audit = ltotal_registry_audit(p1907, p1837)
    requirements = requirement_table(p1961, rep_audit, ltotal_audit)
    extension_ready = all(row["available"] for row in requirements)

    out = {
        "packet_id": "P1962",
        "stage_id": "S912",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE__MATTER_HIGGS_BRST_EXTENSION_REPRESENTATION_REGISTRY_NOT_THEOREM_GRADE",
        "route": "strict_only",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "legacy_bridge_used": False,
        "higher_reasoning_used": True,
        "ignored_files": IGNORED_FILES,
        "input_paths": {
            "p1961_json": str(P1961.relative_to(ROOT)),
            "p1907_json": str(P1907.relative_to(ROOT)),
            "p1837_json": str(P1837.relative_to(ROOT)),
            "p1855_json": str(P1855.relative_to(ROOT)),
            "p1856_json": str(P1856.relative_to(ROOT)),
            "p1801_json": str(P1801.relative_to(ROOT)),
            "p1957_json": str(P1957.relative_to(ROOT)),
        },
        "input_hashes": {
            "p1961_sha256": digest_path(P1961),
            "p1907_sha256": digest_path(P1907),
            "p1837_sha256": digest_path(P1837),
            "p1855_sha256": digest_path(P1855),
            "p1856_sha256": digest_path(P1856),
            "p1801_sha256": digest_path(P1801),
            "p1957_sha256": digest_path(P1957),
        },
        "depends_on": {
            "p1961_present": bool_from_path(P1961),
            "p1907_present": bool_from_path(P1907),
            "p1837_present": bool_from_path(P1837),
            "p1855_present": bool_from_path(P1855),
            "p1856_present": bool_from_path(P1856),
            "p1801_present": bool_from_path(P1801),
            "p1957_present": bool_from_path(P1957),
            "p1957_status": artifact_status(p1957),
        },
        "ConditionalMatterHiggsBRSTExtension_strict_B1_v1": {
            "scope": "conditional local matter/Higgs module rule only",
            "rules": {
                "s_Phi_R": "s Phi_R = - C_R Phi_R",
                "C_R": "C_R = g3 c3^A T3_R^A + g2 c2^I T2_R^I + g1 c1 Y_R",
                "s_C_R": "s C_R = - C_R^2",
                "s_Phi_bar_R": "s Phi_bar_R = Phi_bar_R C_R",
            },
            "nilpotency_condition": (
                "s^2 Phi_R=0 follows only after exporting representation operators "
                "satisfying [T_a,T_b]=f_ab^c T_c, cross-factor commutation, and central U(1) action."
            ),
            "not_claimed": [
                "unconditional matter/Higgs BRST nilpotency",
                "full strict L_total BRST invariance",
                "BV master action",
                "BRST charge Q",
                "TG2_BRST_GLOBAL_NILPOTENCY PASS",
            ],
        },
        "RepresentationRegistryAudit_strict_B1_v1": rep_audit,
        "LTotalFieldRegistryAudit_strict_B1_v1": ltotal_audit,
        "ExtensionReadinessRequirements_strict_B1_v1": requirements,
        "ReadinessAfterP1962_strict_v1": {
            "formula_local_extension": (
                "P1961_LOCAL_READY & FULL_NONPROXY_FIELD_REGISTRY & CHIRAL_REP_REGISTRY "
                "& REPRESENTATION_MATRICES & COMMUTATOR_CERTIFICATES "
                "& YUKAWA_GAUGE_INVARIANCE & ANOMALY_CANCELLATION"
            ),
            "evaluated_local_matter_higgs_extension_ready": extension_ready,
            "evaluated_global_TG2_ready": False,
            "tg2_lock_source": "P1767/P1801/P1957: no G_BW PASS, no global Q, no cohomology pack, no shared freeze",
            "p1801_requirements_reused": p1801.get("tg2_pass_requirements", []),
        },
        "P1961_to_P1962_delta": {
            "closed_by_P1961": [
                "local gauge-sector BRST rules for SU3/SU2/U1",
                "local s^2=0 on A,c,cbar,B",
            ],
            "new_in_P1962": [
                "conditional matter/Higgs module BRST rule recorded",
                "full-bundle extension readiness audited against current exported artifacts",
                "precise representation/field-registry obstruction exported",
            ],
            "still_open": [
                "field-by-field full nonproxy matter/Higgs/spinor registry",
                "explicit representation matrices for each multiplet",
                "representation commutator certificates",
                "Higgs representation in the same registry",
                "Yukawa gauge-invariance certificate",
                "full anomaly cancellation witness",
                "global BV/BRST Q and Q^2=0 witness pack",
            ],
        },
        "false_pass_guard": (
            "P1962 is a local extension-readiness obstruction and conditional rule packet. "
            "It does not promote TG2 and does not override the BW->BRST->Cutkosky gate order."
        ),
        "next_honest_step": (
            "Build P1963: export a strict representation registry with explicit SU3/SU2/U1 "
            "operators for each chiral fermion multiplet and Higgs field, then machine-check "
            "commutators, Yukawa gauge invariance, and anomaly sums before reattempting local "
            "matter/Higgs BRST nilpotency."
        ),
        "lay_explanation": (
            "P1961 proved the BRST lock works for gauge fields alone. P1962 checks whether the "
            "lock has the exact keys for matter and Higgs. The current repo has a seed list, "
            "but not the explicit representation matrices and full field registry needed for "
            "a theorem-grade extension."
        ),
    }

    GEN.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
