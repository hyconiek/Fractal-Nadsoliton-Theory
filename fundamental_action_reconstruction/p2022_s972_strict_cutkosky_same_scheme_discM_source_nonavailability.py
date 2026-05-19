#!/usr/bin/env python3
"""P2022 S972 strict Cutkosky same-scheme DiscM source nonavailability.

Next honest step after P2021: stop manipulating proxy residue factors and ask
whether the current strict-side repository exports the actual source data needed
to compute DiscM_common_basis for the exact P2020 {plus,cross} normalization.
The answer is a narrow current-state nonavailability theorem, not a physical
no-go theorem.
"""
from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2022_s972_strict_cutkosky_same_scheme_discM_source_nonavailability.json"
TS = "2026-05-18T00:00:00+00:00"
CHANNEL = "graviton->gauge_gauge"


def load(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def digest(obj: object) -> str:
    blob = json.dumps(obj, sort_keys=True, ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


def contains_text(obj: Any, needle: str) -> bool:
    return needle.lower() in json.dumps(obj, sort_keys=True, ensure_ascii=True).lower()


def find_generated_text_hits(terms: list[str], limit: int = 12) -> dict[str, list[str]]:
    hits = {term: [] for term in terms}
    for path in sorted(GEN.glob("*.json")):
        text = path.read_text(encoding="utf-8")
        lower = text.lower()
        for term in terms:
            if term.lower() in lower and len(hits[term]) < limit:
                hits[term].append(path.name)
    return hits


def row(req_id: str, required_object: str, verdict: str, evidence: list[str], missing: list[str]) -> dict[str, Any]:
    return {
        "req_id": req_id,
        "required_object": required_object,
        "verdict": verdict,
        "evidence": evidence,
        "missing_or_insufficient": missing,
    }


def symbolic_underdetermination_witness() -> dict[str, Any]:
    alpha, beta, gamma = sp.symbols("alpha_gr beta_gr gamma_gr")
    j_r2, j_ric2, j_eh = sp.symbols("J_R2 J_Rmunu2 J_EH_mix")
    cut = sp.eye(2) / sp.pi
    disc = sp.Matrix([[alpha * j_r2 + beta * j_ric2 + gamma * j_eh, 0], [0, alpha * j_r2 + beta * j_ric2 + gamma * j_eh]])
    defect = sp.simplify(disc - cut)

    assignment_zero = {alpha: 0, beta: 0, gamma: 0}
    assignment_cut_like = {alpha: 1, beta: 0, gamma: 0, j_r2: 1 / sp.pi}
    defect_zero = sp.simplify(defect.subs(assignment_zero))
    defect_cut_like = sp.simplify(defect.subs(assignment_cut_like))

    return {
        "discM_placeholder_matrix": [[sp.sstr(disc[i, j]) for j in range(2)] for i in range(2)],
        "cutSum_P2020_no_symmetry_matrix": [[sp.sstr(cut[i, j]) for j in range(2)] for i in range(2)],
        "defect_matrix": [[sp.sstr(defect[i, j]) for j in range(2)] for i in range(2)],
        "open_coefficients": ["alpha_gr", "beta_gr", "gamma_gr", "J_R2", "J_Rmunu2", "J_EH_mix"],
        "assignment_zero_defect_trace": sp.sstr(sp.trace(defect_zero)),
        "assignment_cut_like_defect_trace": sp.sstr(sp.simplify(sp.trace(defect_cut_like))),
        "meaning": "With the current symbolic placeholders, DiscM_minus_CutSum can be nonzero or zero depending on unexported loop coefficients/integrals.  Therefore P2020/P2021 cannot decide Cutkosky equality.",
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p1662 = load("p1662_s612_strict_full_lagrangian_explicit_density_summary.json")
    p1852 = load("p1852_s802_strict_b1_brst_anomaly_and_cutkosky_seed_witness_checkpoint.json")
    p1857 = load("p1857_s807_strict_b1_triangle_amplitude_seed_and_k5_instantiation_checkpoint.json")
    p1955 = load("p1955_s905_strict_minimal_hAA_vertex_export.json")
    p2019 = load("p2019_s969_strict_cutkosky_first_transverse_tree_amplitude_component.json")
    p1859 = load("p1859_s809_strict_b1_cutkosky_kernel_and_residue_certificate_stub_checkpoint.json")
    p1862 = load("p1862_s812_strict_b1_dressed_pole_residue_and_projected_disc_evaluation_checkpoint.json")
    p1913 = load("p1913_s863_strict_c1_gr_common_basis_unitarity_background_probe.json")
    p1953 = load("p1953_s903_strict_dressed_cutkosky_amplitude_availability_audit.json")
    p1994 = load("p1994_s944_strict_cutkosky_dressed_amplitude_first_import_witness.json")
    p2020 = load("p2020_s970_strict_cutkosky_p2019_tree_phase_space_cut_sum_witness.json")
    p2021 = load("p2021_s971_strict_cutkosky_p2020_linear_polarization_residue_dressing_factor_witness.json")

    p1913_scheme = ((p1913.get("common_basis_definition") or {}).get("scheme") or "MISSING")
    p1953_contract = p1953.get("dressed_amplitude_interface_contract", {})
    p1955_minimal_vertex_ok = (p1955.get("machine_check") or {}).get("all_local_checks_zero") is True
    p2019_tree_component_ok = p2019.get("result_kind") == "PASS_FIRST_TREE_TRANSVERSE_COMPONENT_WITNESS"
    p1859_required = (p1859.get("residue_certificate_stub") or {}).get("required_exports") or []

    requirements = [
        row(
            "D1",
            "Explicit hAA vertices derived from full non-skeleton L_total in the P2020 normalization",
            "PARTIAL_TREE_LEVEL_SOURCE_AVAILABLE_NOT_LOOP_DISCM_READY",
            [
                f"P1662 density present={contains_text(p1662, 'L_total') or contains_text(p1662, 'lagrangian')}",
                f"P1955 minimal hAA vertex local checks zero={p1955_minimal_vertex_ok}.",
                f"P2019 P1955/P1956 transverse tree component available={p2019_tree_component_ok}.",
                "P2020 transports and integrates that tree component in the real {plus,cross} linear-polarization basis.",
            ],
            [
                "D1 is not a hard missing source at minimal tree level; the previous P2022 wording was too strong.",
                "Still need loop-level hAA/self-energy insertion rules and finite-part scheme data before DiscM_common_basis can be evaluated.",
            ],
        ),
        row(
            "D2",
            "Gauge fixing, ghost sector, and BRST physical-state projector for the same channel",
            "FAIL_PROJECTOR_MISSING",
            [
                f"P1852 mentions BRST={contains_text(p1852, 'BRST')}",
                "P1953 already reports no channel-level BRST physical-state projector.",
            ],
            [
                "Need nilpotent BRST charge or cohomology projector acting on the gauge_gauge final state.",
                "Need ghost-sector exclusion trace in the same scheme.",
            ],
        ),
        row(
            "D3",
            "Loop-derived dressed residue/pole data, not seed inheritance or proxy Z(s)",
            "FAIL_PROXY_OR_SEED_ONLY",
            [
                f"P1859 required exports={p1859_required}",
                "P1862 uses seed inherited dressed residues.",
                "P1994 Z(s) is a mirrored proxy and P2021 rejects it as P1953 input.",
            ],
            [
                "Need dressed propagator pole list computed from the loop-corrected operator.",
                "Need residue values per pole and positivity proof after BRST projection.",
            ],
        ),
        row(
            "D4",
            "DiscM_common_basis for the exact P2020 {plus,cross} matrix basis",
            "FAIL_OPEN_SYMBOLIC_ONLY",
            [
                f"P1913 scheme={p1913_scheme}",
                "P1913 grmix row remains symbolic placeholders rather than evaluated DiscM.",
            ],
            [
                "Need evaluated alpha_gr, beta_gr, gamma_gr and loop master integrals in the P2020 basis.",
                "Need matrix-valued DiscM with the same external linear-polarization index order {plus,cross}.",
            ],
        ),
        row(
            "D5",
            "Same renormalization/gauge scheme lock across P1953/P2020/P2021 and DiscM backend",
            "FAIL_SCHEME_NOT_LOCKED",
            [
                "P1953 requires scheme=MSbar_B1_seed.",
                f"P1913 common-basis scheme={p1913_scheme}.",
            ],
            [
                "Need a single scheme tag and finite-part convention shared by vertex, propagator, DiscM, and CutSum.",
            ],
        ),
    ]

    hard_failures = [item for item in requirements if str(item["verdict"]).startswith("FAIL")]
    partial_sources = [item for item in requirements if str(item["verdict"]).startswith("PARTIAL")]
    source_hits = find_generated_text_hits([
        "hAA",
        "M_dressed_common_basis",
        "DiscM_common_basis",
        "BRST physical-state projector",
        "ghost_sector_exclusion_trace",
        "Z(s)",
        "P2021",
    ])

    gate = {
        "p1955_minimal_hAA_vertex_available": p1955_minimal_vertex_ok,
        "p2019_tree_component_available": p2019_tree_component_ok,
        "p2020_matrix_available": p2020.get("result_kind") == "PASS_TREE_PHASE_SPACE_LINEAR_POLARIZATION_CUT_SUM_COMPONENT_WITNESS",
        "p2021_proxy_rejected": p2021.get("result_kind") == "OPEN_PROXY_RESIDUE_TRANSPORT_SANITY_ONLY_NOT_P1953_ADMISSIBLE",
        "p1953_contract_available": bool(p1953_contract.get("required_fields")),
        "one_tree_source_requirement_partial_not_failed": len(partial_sources) == 1 and partial_sources[0]["req_id"] == "D1",
        "remaining_discM_source_requirements_fail_currently": len(hard_failures) == len(requirements) - 1,
        "symbolic_underdetermination_exported": True,
        "no_cutkosky_promotion": True,
    }

    out = {
        "ledger_id": "P2022_S972_STRICT_CUTKOSKY_SAME_SCHEME_DISCM_SOURCE_NONAVAILABILITY",
        "packet_id": "P2022",
        "stage_id": "S972",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TS,
        "route": "strict_only",
        "legacy_bridge_used": False,
        "channel": CHANNEL,
        "depends_on": {
            "p1953": p1953.get("packet_id"),
            "p2020": p2020.get("ledger_id"),
            "p2021": p2021.get("ledger_id"),
        },
        "input_hashes": {
            "p1852_sha256": digest(p1852),
            "p1857_sha256": digest(p1857),
            "p1859_sha256": digest(p1859),
            "p1955_sha256": digest(p1955),
            "p2019_sha256": digest(p2019),
            "p1862_sha256": digest(p1862),
            "p1913_sha256": digest(p1913),
            "p1953_sha256": digest(p1953),
            "p1994_sha256": digest(p1994),
            "p2020_sha256": digest(p2020),
            "p2021_sha256": digest(p2021),
        },
        "source_search_hits": source_hits,
        "discM_source_requirements": requirements,
        "hard_failure_count": len(hard_failures),
        "partial_source_count": len(partial_sources),
        "correction_to_prior_p2022_d1_wording": "D1 is now marked partial, not hard-fail, because P1955/P2019/P2020 already export the minimal tree hAA source chain.  The obstruction is the missing loop/scheme/BRST data needed for DiscM, not the complete absence of a tree hAA vertex.",
        "symbolic_underdetermination_witness": symbolic_underdetermination_witness(),
        "gatekeeper_checks": gate,
        "result_kind": "OPEN_SAME_SCHEME_DISCM_SOURCE_PARTIAL_TREE_VERTEX_AVAILABLE_WITH_TRACE" if all(gate.values()) else "OPEN_SAME_SCHEME_DISCM_SOURCE_AUDIT_INCOMPLETE",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "false_pass_guard": "P2022 proves only current-state nonavailability of the loop-derived same-scheme DiscM source data needed to compare with P2020, while preserving the already available minimal tree hAA source chain. It is not a no-go theorem, not DiscM=CutSum, not BRST closure, not all-state unitarity, not QW-2191 discharge, and not ToE closure.",
        "next_honest_step": "Build P2023 by deriving or explicitly tabulating the BRST physical-state projector and ghost-sector exclusion trace for graviton->gauge_gauge in one fixed scheme; in parallel, specify the loop-level self-energy/vertex insertion rules that upgrade the existing tree hAA chain into a same-scheme DiscM_common_basis calculation.",
        "toe_progress": "Moves the unitarity route from proxy-factor algebra to the real source-data frontier: exact P2020 CutSum and the minimal tree hAA source chain are ready, but same-scheme loop DiscM cannot be computed until loop-level vertex/self-energy data, BRST projector, ghost trace, and scheme lock are exported.",
        "lay_explanation": "Sprawdziliśmy, czy repo ma już prawdziwe dane potrzebne do policzenia lewej strony twierdzenia optycznego dla tej samej macierzy polaryzacji. Nie ma: są seedy, proxy i symbole, ale brakuje pętlowego wyprowadzenia, projektora BRST i wspólnego schematu. To zawęża następny problem do konkretnych brakujących elementów.",
        "environment": {
            "python": platform.python_version(),
            "sympy": sp.__version__,
        },
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P2022] wrote same-scheme DiscM source nonavailability: {OUT}")


if __name__ == "__main__":
    main()
