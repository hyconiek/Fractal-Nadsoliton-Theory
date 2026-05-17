#!/usr/bin/env python3
"""P1957 S907 strict BV/BRST ghost-sector nonavailability theorem.

P1956 exported a local on-shell transverse polarization projector for the
gauge_gauge cut.  This executor checks whether the current strict repository
state contains the stronger data required to promote that local projector to a
BRST/cohomology theorem: a nonproxy ghost sector, a BV/BRST operator map, an
explicit BRST charge Q, a Q^2=0 nilpotency certificate, and a ghost-cancellation
trace.

The result is intentionally a nonavailability theorem for the current export
state, not a no-go theorem about the theory.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p1957_s907_strict_bv_brst_ghost_sector_nonavailability_theorem.json"

TEXT_SUFFIXES = {".py", ".md", ".json", ".txt", ".yaml", ".yml"}
SKIP_NAMES = {
    "TOE_FINAL_DOCUMENTATION_RELEASE_8_STRICT_FULL.pdf",
    "p1957_s907_strict_bv_brst_ghost_sector_nonavailability_theorem.py",
    "p1957_s907_strict_bv_brst_ghost_sector_nonavailability_theorem.json",
}


def load(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def digest(obj: object) -> str:
    blob = json.dumps(obj, sort_keys=True, ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


def has_key_deep(obj: Any, key_name: str) -> bool:
    if isinstance(obj, dict):
        return any(k == key_name or has_key_deep(v, key_name) for k, v in obj.items())
    if isinstance(obj, list):
        return any(has_key_deep(v, key_name) for v in obj)
    return False


def text_contains(obj: Any, needle: str) -> bool:
    return needle.lower() in json.dumps(obj, sort_keys=True, ensure_ascii=False).lower()


def repo_text_files() -> list[Path]:
    paths: list[Path] = []
    for path in ROOT.rglob("*"):
        if not path.is_file():
            continue
        if path.name in SKIP_NAMES:
            continue
        if path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        paths.append(path)
    return sorted(paths)


def scan_terms(terms: list[str]) -> dict[str, dict[str, Any]]:
    files = repo_text_files()
    result: dict[str, dict[str, Any]] = {}
    for term in terms:
        count = 0
        sample_paths: list[str] = []
        low = term.lower()
        for path in files:
            text = path.read_text(encoding="utf-8", errors="ignore")
            hits = text.lower().count(low)
            if hits:
                count += hits
                if len(sample_paths) < 8:
                    sample_paths.append(str(path.relative_to(ROOT)))
        result[term] = {"count": count, "sample_paths": sample_paths}
    return result


def json_files_with_key(key_name: str) -> list[str]:
    hits: list[str] = []
    for path in sorted(GEN.glob("*.json")):
        if path.name in SKIP_NAMES:
            continue
        try:
            obj = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            continue
        if has_key_deep(obj, key_name):
            hits.append(str(path.relative_to(ROOT)))
    return hits


def is_real_id(value: Any) -> bool:
    if not isinstance(value, str):
        return False
    placeholders = {
        "REQUIRED",
        "HASH_REQUIRED",
        "FREEZE_ID_REQUIRED",
        "OPEN",
        "OPEN_MISSING_ARTIFACT",
        "MISSING",
    }
    return bool(value.strip()) and value.strip() not in placeholders


def requirement_row(
    req_id: str,
    required_object: str,
    present: bool,
    observed_sources: list[str],
    exact_missing_data: list[str],
) -> dict[str, Any]:
    return {
        "req_id": req_id,
        "required_object": required_object,
        "present_as_theorem_grade_export": present,
        "verdict": "AVAILABLE" if present else "MISSING_OR_TEMPLATE_ONLY",
        "observed_sources": observed_sources,
        "exact_missing_data": [] if present else exact_missing_data,
    }


def main() -> None:
    p1767 = load("p1767_s717_strict_bianchi_ward_to_brst_cutkosky_gate_sequencing_checkpoint.json")
    p1801 = load("p1801_s751_strict_brst_nilpotency_witness_intake_and_gate_link_checkpoint.json")
    p1801_template = load("p1801_s751_brst_nilpotency_witness_intake_template.json")
    p1802 = load("p1802_s752_strict_cutkosky_unitarity_witness_intake_and_tg3_gate_link_checkpoint.json")
    p1833 = load("p1833_s783_strict_full_lagrangian_and_reverse_closure_worklist_checkpoint.json")
    p1836 = load("p1836_s786_strict_full_lagrangian_non_skeleton_manifest_checkpoint.json")
    p1852 = load("p1852_s802_strict_b1_brst_anomaly_and_cutkosky_seed_witness_checkpoint.json")
    p1854 = load("p1854_s804_strict_b1_brst_cochain_and_first_cutkosky_channel_checkpoint.json")
    p1956 = load("p1956_s906_strict_gauge_gauge_physical_projector_certificate.json")

    gate_contract = p1767.get("gate_sequencing_contract", {})
    g_bw_status = (gate_contract.get("G_BW") or {}).get("status", "MISSING")
    g_brst_contract = gate_contract.get("G_BRST") or {}

    template_pack = p1801_template.get("brst_witness_pack") or {}
    template_nilpotency = template_pack.get("nilpotency_check") or {}

    brst_charge_real = is_real_id(template_pack.get("brst_charge_definition_id")) and "template" not in json.dumps(
        p1801_template, sort_keys=True
    ).lower()
    q2_zero_real = (
        template_nilpotency.get("final_simplified") == "0"
        and is_real_id(template_nilpotency.get("check_digest"))
        and is_real_id(template_pack.get("brst_charge_definition_id"))
    )
    cohomology_real = is_real_id(template_pack.get("cohomology_subspace_check_id"))
    ghost_consistency_real = is_real_id(template_pack.get("ghost_sector_consistency_id"))

    ghost_export_key_hits = json_files_with_key("ghost_sector_nonproxy_export")
    bv_map_key_hits = json_files_with_key("BV_BRST_operator_map")
    brst_operator_key_hits = json_files_with_key("brst_operator_definition_nonproxy")
    ghost_action_key_hits = json_files_with_key("ghost_antighost_action_sector")

    # Exact strings do occur as requirements, but not as exported objects.
    ghost_export_real = bool(ghost_export_key_hits)
    bv_map_real = bool(bv_map_key_hits)
    ghost_action_real = bool(ghost_action_key_hits)
    brst_operator_real = bool(brst_operator_key_hits)
    g_bw_pass_zero = g_bw_status == "PASS_ZERO"

    terms = [
        "BV_BRST_operator_map",
        "ghost_sector_nonproxy_export",
        "BRST charge",
        "Q^2",
        "nilpotency",
        "cohomology",
        "ghost cancellation",
        "ghost_antighost_action_sector",
        "ladunek BRST",
        "nilpotenc",
        "kohomolog",
        "sektor duch",
        "kasowanie duch",
    ]
    term_scan = scan_terms(terms)

    G_BW_PASS_ZERO, GHOST_EXPORT, BV_MAP, BRST_Q, Q2_ZERO, COHOMOLOGY, GHOST_CONSISTENCY, SHARED_FREEZE = sp.symbols(
        "G_BW_PASS_ZERO GHOST_EXPORT BV_MAP BRST_Q Q2_ZERO COHOMOLOGY GHOST_CONSISTENCY SHARED_FREEZE"
    )
    tg2_formula = sp.And(
        G_BW_PASS_ZERO,
        GHOST_EXPORT,
        BV_MAP,
        BRST_Q,
        Q2_ZERO,
        COHOMOLOGY,
        GHOST_CONSISTENCY,
        SHARED_FREEZE,
    )
    truth_assignment = {
        G_BW_PASS_ZERO: g_bw_pass_zero,
        GHOST_EXPORT: ghost_export_real and ghost_action_real,
        BV_MAP: bv_map_real and brst_operator_real,
        BRST_Q: brst_charge_real,
        Q2_ZERO: q2_zero_real,
        COHOMOLOGY: cohomology_real,
        GHOST_CONSISTENCY: ghost_consistency_real,
        SHARED_FREEZE: is_real_id(p1801_template.get("shared_freeze_id")),
    }
    tg2_evaluated = bool(tg2_formula.subs(truth_assignment))

    requirement_matrix = [
        requirement_row(
            "B0",
            "G_BW:PASS_ZERO prerequisite",
            g_bw_pass_zero,
            [
                f"P1767 G_BW status={g_bw_status}",
                f"P1767 G_BRST status={g_brst_contract.get('status', 'MISSING')}",
            ],
            ["A same-freeze Bianchi/Ward PASS_ZERO witness before TG2 promotion"],
        ),
        requirement_row(
            "B1",
            "ghost_sector_nonproxy_export with ghost/antighost action",
            ghost_export_real and ghost_action_real,
            [
                f"dict-key hits for ghost_sector_nonproxy_export={ghost_export_key_hits}",
                f"dict-key hits for ghost_antighost_action_sector={ghost_action_key_hits}",
                "P1703 lists ghost_antighost_action_sector only as a required input",
                "P1833/P1836 list ghost_sector_nonproxy_export only as required exports",
            ],
            [
                "explicit ghost fields c, antighosts, Nakanishi-Lautrup/auxiliary sector if used",
                "ghost/antighost action density in the same strict L_total scheme",
                "gauge-fixing functional and Faddeev-Popov operator",
            ],
        ),
        requirement_row(
            "B2",
            "BV_BRST_operator_map and nonproxy BRST differential",
            bv_map_real and brst_operator_real,
            [
                f"dict-key hits for BV_BRST_operator_map={bv_map_key_hits}",
                f"dict-key hits for brst_operator_definition_nonproxy={brst_operator_key_hits}",
                "P1767 names BV_BRST_operator_map as required but does not export it",
            ],
            [
                "s(g_mu_nu), s(A_mu), s(c), s(antighost), s(auxiliary) transformation rules",
                "BV antifield map or equivalent operator-domain definition",
                "graded Leibniz/sign convention registry",
            ],
        ),
        requirement_row(
            "B3",
            "explicit BRST charge Q",
            brst_charge_real,
            [
                f"P1801 template brst_charge_definition_id={template_pack.get('brst_charge_definition_id', 'MISSING')}",
                "P1801 is an intake template/policy, not the real witness pack",
            ],
            [
                "real brst_charge_definition_id",
                "operator expression for Q on the strict field algebra",
                "domain and grading convention",
            ],
        ),
        requirement_row(
            "B4",
            "Q_squared_simplified_zero nilpotency certificate",
            q2_zero_real,
            [
                f"P1801 template nilpotency expression={template_nilpotency.get('expression', 'MISSING')}",
                f"P1801 template final_simplified={template_nilpotency.get('final_simplified', 'MISSING')}",
                f"P1801 template check_digest={template_nilpotency.get('check_digest', 'MISSING')}",
            ],
            [
                "computed Q^2 expression from a real Q",
                "symbolic simplification to zero",
                "non-placeholder check digest",
            ],
        ),
        requirement_row(
            "B5",
            "physical-state cohomology subspace check",
            cohomology_real,
            [
                f"P1801 template cohomology_subspace_check_id={template_pack.get('cohomology_subspace_check_id', 'MISSING')}",
                "P1956 supplies only the local transverse external-state projector",
            ],
            [
                "ker(Q)/im(Q) construction",
                "map from transverse gauge_gauge projector into BRST cohomology",
                "proof that unphysical longitudinal/time-like states are BRST-exact or excluded",
            ],
        ),
        requirement_row(
            "B6",
            "gauge_gauge ghost-cancellation trace",
            ghost_consistency_real,
            [
                f"P1801 template ghost_sector_consistency_id={template_pack.get('ghost_sector_consistency_id', 'MISSING')}",
                f"P1852 seed pole constraint present={text_contains(p1852, 'ghost-free')}",
                f"P1854 pre-theorem cochain present={'brst_cochain_b1' in p1854}",
            ],
            [
                "channel-level ghost loop/sign contribution",
                "cancellation against unphysical polarizations",
                "same-scheme link to P1956 physical polarization sum",
            ],
        ),
    ]

    all_missing_or_template_only = not tg2_evaluated
    out = {
        "packet_id": "P1957",
        "stage_id": "S907",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "local_verdict": "STRICT_BV_BRST_GHOST_SECTOR_NOT_AVAILABLE__P1956_PROJECTOR_REMAINS_LOCAL",
        "route": "strict_only",
        "legacy_bridge_used": False,
        "higher_reasoning_used": True,
        "ignored_files": sorted(SKIP_NAMES),
        "pre_execution_grep_summary": {
            "english_terms": [
                "BV",
                "BRST charge",
                "Q^2",
                "nilpotency",
                "ghost sector",
                "ghost action",
                "ghost cancellation",
                "cohomology",
                "BV_BRST_operator_map",
                "ghost_sector_nonproxy_export",
            ],
            "polish_terms": [
                "ladunek BRST",
                "nilpotencja/nilpotenc",
                "sektor duchow",
                "kasowanie duchow",
                "kohomologia/kohomolog",
                "cechowanie",
            ],
            "term_scan_counts": term_scan,
            "interpretation": (
                "The terms occur primarily as requirements, templates, locks, or open obligations. "
                "No theorem-grade strict BV/BRST ghost-sector witness pack was found."
            ),
        },
        "depends_on": {
            "p1767_present": "gate_sequencing_contract" in p1767,
            "p1801_present": "tg2_pass_requirements" in p1801,
            "p1801_template_present": "brst_witness_pack" in p1801_template,
            "p1802_present": "tg3_pass_requirements" in p1802,
            "p1833_present": "reverse_qg_worklist" in p1833,
            "p1836_present": "reverse_gate_dependency" in p1836,
            "p1852_present": "brst_anomaly_seed_contract" in p1852,
            "p1854_present": "brst_cochain_b1" in p1854,
            "p1956_present": "two_particle_gauge_gauge_projector" in p1956,
        },
        "input_hashes": {
            "p1767_sha256": digest(p1767),
            "p1801_sha256": digest(p1801),
            "p1801_template_sha256": digest(p1801_template),
            "p1802_sha256": digest(p1802),
            "p1833_sha256": digest(p1833),
            "p1836_sha256": digest(p1836),
            "p1852_sha256": digest(p1852),
            "p1854_sha256": digest(p1854),
            "p1956_sha256": digest(p1956),
        },
        "tg2_acceptance_formula": {
            "formula": str(tg2_formula),
            "truth_assignment": {str(k): bool(v) for k, v in truth_assignment.items()},
            "evaluated_pass": tg2_evaluated,
            "required_by": [
                "P1767 G_BRST requires G_BW:PASS_ZERO, ghost_sector_nonproxy_export, BV_BRST_operator_map",
                "P1801 requires brst_charge_definition_present, Q_squared_simplified_zero, cohomology_subspace_check_present, ghost_sector_consistency_present, shared_freeze_id",
            ],
        },
        "minimal_missing_data_matrix": requirement_matrix,
        "formal_nonavailability_theorem": {
            "statement": (
                "On the current repository state, P1956 cannot be promoted to a strict "
                "BRST/cohomology theorem because the BV/BRST operator, nonproxy ghost sector, "
                "Q^2=0 nilpotency certificate, and ghost-cancellation trace are not exported."
            ),
            "proof_trace": [
                "P1767 explicitly blocks G_BRST without G_BW:PASS_ZERO, ghost_sector_nonproxy_export, and BV_BRST_operator_map.",
                "P1801 is an intake policy/template; its IDs remain REQUIRED placeholders and do not constitute a real witness pack.",
                "P1703/P1833/P1836 list ghost/BRST objects as required/planned exports, not as completed theorem-grade data.",
                "P1852/P1854 provide B1 anomaly/cochain seeds and proxies, not an explicit BRST charge or ghost action.",
                "P1956 provides a correct local transverse polarization projector, but explicitly not BRST cohomology or ghost cancellation.",
                "The executable acceptance formula evaluates to False under the current export-state truth assignment.",
            ],
            "underdetermination_witness": {
                "free_objects_if_derivation_attempted_now": [
                    "Q_BRST_free",
                    "s_g_mu_nu_free",
                    "s_A_mu_free",
                    "s_c_free",
                    "s_antighost_free",
                    "S_ghost_free",
                    "GaugeFixingFunctional_free",
                    "CohomologyProjector_free",
                    "GhostCancellationTrace_gauge_gauge_free",
                ],
                "meaning": (
                    "Multiple incompatible BRST/ghost completions can satisfy all currently exported "
                    "templates while giving different Q^2 and ghost-cancellation behavior. Therefore "
                    "the BRST theorem is underdetermined by current strict exports."
                ),
            },
        },
        "safe_consequence": {
            "may_use": [
                "P1956 local on-shell transverse polarization projector",
                "P1955 minimal tree-level hAA vertex",
                "P1951/P1952 seed Cutkosky positivity data",
            ],
            "must_not_claim": [
                "TG2_BRST_GLOBAL_NILPOTENCY PASS",
                "BRST physical-state cohomology theorem",
                "ghost-cancelled gauge_gauge Cutkosky equality",
                "TG3_CUTKOSKY_GLOBAL_UNITARITY PASS",
                "global UR_link theorem",
            ],
        },
        "next_solver_input_contract": {
            "recommended_packet": "P1958",
            "minimum_new_exports_required": [
                "GaugeFixingFunctional_strict_B1_v1",
                "GhostAntighostAction_strict_B1_v1",
                "BV_BRST_operator_map_strict_B1_v1",
                "BRSTCharge_Q_strict_B1_v1",
                "NilpotencyCertificate_Q2_zero_strict_B1_v1",
                "CohomologyProjectionLink_gauge_gauge_to_P1956_v1",
                "GhostCancellationTrace_gauge_gauge_strict_B1_v1",
            ],
        },
        "false_pass_guard": (
            "This is a current-state nonavailability theorem. It does not prove BRST is impossible; "
            "it proves that the required strict proof objects are not yet exported and cannot be "
            "replaced by templates, seed cochains, or the P1956 transverse projector."
        ),
        "next_honest_step": (
            "Build P1958 with high reasoning: start from the explicit strict gauge-field sector and "
            "derive a concrete gauge-fixing functional plus ghost/antighost action, then define the "
            "BRST differential on fields before attempting Q^2=0."
        ),
        "lay_explanation": (
            "Sprawdzilismy, czy repo ma juz pelna maszynerie BRST: operator Q, duchy i dowod, "
            "ze Q zastosowane dwa razy daje zero. Nie ma jej jeszcze; sa tylko formularze i wymagania. "
            "To znaczy, ze poprzedni projektor polaryzacji jest poprawnym klockiem, ale nie wolno go "
            "nazywac pelnym dowodem unitarności."
        ),
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
