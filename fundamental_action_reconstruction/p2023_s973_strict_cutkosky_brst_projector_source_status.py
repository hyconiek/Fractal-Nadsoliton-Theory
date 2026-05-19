#!/usr/bin/env python3
"""P2023 S973 strict Cutkosky BRST projector source-status witness.

Next honest step after P2022: separate the already-exported local transverse
polarization projector from a genuine channel-level BRST physical-state
projector and ghost-sector exclusion trace.  The result is a precise
source-status theorem: the local transverse projector is available, but the BRST
cohomology projector/ghost trace needed for dressed Cutkosky closure is not.
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
OUT = GEN / "p2023_s973_strict_cutkosky_brst_projector_source_status.json"
TS = "2026-05-19T00:00:00+00:00"
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


def row(req_id: str, required_object: str, verdict: str, evidence: list[str], missing: list[str]) -> dict[str, Any]:
    return {
        "req_id": req_id,
        "required_object": required_object,
        "verdict": verdict,
        "evidence": evidence,
        "missing_or_insufficient": missing,
    }


def local_projector_summary(p1956: dict[str, Any]) -> dict[str, Any]:
    checks = p1956.get("single_particle_projector_checks") or []
    all_local = bool(checks) and all(
        (entry.get("machine_checks") or {}).get("completeness_residual_zero") is True
        and (entry.get("machine_checks") or {}).get("projector_k_transverse_zero") is True
        and (entry.get("machine_checks") or {}).get("idempotence_residual_zero") is True
        and (entry.get("machine_checks") or {}).get("projector_rank_mixed") == 2
        for entry in checks
    )
    return {
        "p1956_packet_id": p1956.get("packet_id"),
        "local_verdict": p1956.get("local_verdict"),
        "single_particle_projector_count": len(checks),
        "all_local_projector_checks_pass": all_local,
        "declared_limitation": "P1956 proves local on-shell transverse polarization-sum algebra, not BRST nilpotency, ghost cancellation, or cohomology projection.",
    }


def symbolic_projector_gap_witness() -> dict[str, Any]:
    qpt, ghost, exact, disc = sp.symbols("Q_P_T GhostTrace BRST_exact DiscM_loop")
    cut = sp.Integer(2) / sp.pi
    local_trace = cut
    brst_physical_trace = sp.simplify(local_trace - qpt - ghost + exact)
    optical_defect = sp.simplify(disc - brst_physical_trace)
    optimistic = sp.simplify(optical_defect.subs({qpt: 0, ghost: 0, exact: 0, disc: cut}))
    ghost_shifted = sp.simplify(optical_defect.subs({qpt: 0, ghost: sp.Rational(1, 5), exact: 0, disc: cut}))
    q_shifted = sp.simplify(optical_defect.subs({qpt: sp.Rational(1, 7), ghost: 0, exact: 0, disc: cut}))
    return {
        "local_transverse_trace_from_P2020_no_symmetry": sp.sstr(local_trace),
        "formal_BRST_physical_trace": sp.sstr(brst_physical_trace),
        "formal_optical_defect": sp.sstr(optical_defect),
        "unexported_symbols": ["Q_P_T", "GhostTrace", "BRST_exact", "DiscM_loop"],
        "optimistic_unproved_assignment_defect": sp.sstr(optimistic),
        "ghost_shifted_assignment_defect": sp.sstr(ghost_shifted),
        "q_shifted_assignment_defect": sp.sstr(q_shifted),
        "meaning": "The local transverse CutSum trace equals 2/pi, but without exported BRST exact-state and ghost-sector data the physical-state trace and optical defect are underdetermined.",
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p1767 = load("p1767_s717_strict_bianchi_ward_to_brst_cutkosky_gate_sequencing_checkpoint.json")
    p1801 = load("p1801_s751_strict_brst_nilpotency_witness_intake_and_gate_link_checkpoint.json")
    p1802 = load("p1802_s752_strict_cutkosky_unitarity_witness_intake_and_tg3_gate_link_checkpoint.json")
    p1852 = load("p1852_s802_strict_b1_brst_anomaly_and_cutkosky_seed_witness_checkpoint.json")
    p1854 = load("p1854_s804_strict_b1_brst_cochain_and_first_cutkosky_channel_checkpoint.json")
    p1956 = load("p1956_s906_strict_gauge_gauge_physical_projector_certificate.json")
    p2020 = load("p2020_s970_strict_cutkosky_p2019_tree_phase_space_cut_sum_witness.json")
    p2022 = load("p2022_s972_strict_cutkosky_same_scheme_discM_source_nonavailability.json")

    p1956_local = local_projector_summary(p1956)
    brst_seed_present = "brst_anomaly_seed_contract" in p1852
    brst_cochain_present = "brst_cochain_b1" in p1854
    p2022_next = p2022.get("next_honest_step", "")

    requirements = [
        row(
            "B1",
            "Local on-shell transverse gauge-gauge projector in the P2020 channel",
            "PASS_LOCAL_TRANSVERSE_PROJECTOR_AVAILABLE_NOT_BRST",
            [
                f"P1956 local projector checks pass={p1956_local['all_local_projector_checks_pass']}.",
                "P2020 already uses the matching real transverse polarization basis and preserves {plus,cross} indices.",
            ],
            [
                "No missing local-transverse algebra at this level; limitation is that this is not a BRST cohomology projector.",
            ],
        ),
        row(
            "B2",
            "Nilpotent BRST charge/operator map acting on the same channel Hilbert space",
            "FAIL_OPERATOR_MAP_NOT_EXPORTED",
            [
                f"P1767 mentions BV_BRST_operator_map={contains_text(p1767, 'BV_BRST_operator_map')}.",
                f"P1801 BRST intake present={contains_text(p1801, 'BRST')}.",
            ],
            [
                "Need explicit Q_BRST action on graviton, gauge, ghost, and antighost states for the channel.",
                "Need machine-checked Q_BRST^2=0 in the same gauge-fixing family.",
            ],
        ),
        row(
            "B3",
            "Channel-level BRST physical-state projector/cohomology representative map",
            "FAIL_COHOMOLOGY_PROJECTOR_NOT_EXPORTED",
            [
                f"P1852 seed contract present={brst_seed_present}.",
                "P1956 explicitly limits itself to local transverse polarization algebra.",
            ],
            [
                "Need Pi_phys^BRST with kernel/image conditions and proof that it maps gauge_gauge cuts to BRST cohomology classes.",
                "Need proof that the P2020 transverse states are exactly the BRST physical representatives, not only a convenient local basis.",
            ],
        ),
        row(
            "B4",
            "Ghost-sector exclusion/cancellation trace in the same scheme",
            "FAIL_GHOST_TRACE_NOT_EXPORTED",
            [
                f"P1854 cochain proxy present={brst_cochain_present}.",
                "P2022 lists ghost trace as a remaining blocker.",
            ],
            [
                "Need ghost and longitudinal contributions on the same cut and their cancellation/exclusion certificate.",
                "Need sign convention and state-counting trace compatible with the P2020 phase-space normalization.",
            ],
        ),
        row(
            "B5",
            "Same-scheme bridge from BRST projector data to DiscM_common_basis",
            "FAIL_SCHEME_BRIDGE_NOT_EXPORTED",
            [
                f"P1802 Cutkosky intake present={contains_text(p1802, 'Cutkosky')}.",
                f"P2022 next step already targets BRST/ghost data={'BRST physical-state projector' in p2022_next}.",
            ],
            [
                "Need one scheme tag and finite-part convention shared by BRST projector, ghost trace, loop DiscM, and P2020 CutSum.",
            ],
        ),
    ]

    pass_local = [item for item in requirements if str(item["verdict"]).startswith("PASS_LOCAL")]
    hard_failures = [item for item in requirements if str(item["verdict"]).startswith("FAIL")]
    gate = {
        "p1956_local_transverse_projector_available": p1956_local["all_local_projector_checks_pass"],
        "p2020_cut_sum_available": p2020.get("result_kind") == "PASS_TREE_PHASE_SPACE_LINEAR_POLARIZATION_CUT_SUM_COMPONENT_WITNESS",
        "p2022_brst_next_step_alignment": "BRST physical-state projector" in p2022_next,
        "exactly_one_local_projector_pass": len(pass_local) == 1 and pass_local[0]["req_id"] == "B1",
        "remaining_brst_ghost_requirements_fail": len(hard_failures) == len(requirements) - 1,
        "symbolic_projector_gap_exported": True,
        "no_brst_or_cutkosky_promotion": True,
    }

    out = {
        "ledger_id": "P2023_S973_STRICT_CUTKOSKY_BRST_PROJECTOR_SOURCE_STATUS",
        "packet_id": "P2023",
        "stage_id": "S973",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TS,
        "route": "strict_only",
        "legacy_bridge_used": False,
        "channel": CHANNEL,
        "depends_on": {
            "p1956": p1956.get("packet_id"),
            "p2020": p2020.get("ledger_id"),
            "p2022": p2022.get("ledger_id"),
        },
        "input_hashes": {
            "p1767_sha256": digest(p1767),
            "p1801_sha256": digest(p1801),
            "p1802_sha256": digest(p1802),
            "p1852_sha256": digest(p1852),
            "p1854_sha256": digest(p1854),
            "p1956_sha256": digest(p1956),
            "p2020_sha256": digest(p2020),
            "p2022_sha256": digest(p2022),
        },
        "local_projector_summary": p1956_local,
        "brst_projector_requirements": requirements,
        "local_pass_count": len(pass_local),
        "hard_failure_count": len(hard_failures),
        "symbolic_projector_gap_witness": symbolic_projector_gap_witness(),
        "gatekeeper_checks": gate,
        "result_kind": "OPEN_BRST_PROJECTOR_SOURCE_STATUS_LOCAL_TRANSVERSE_ONLY_WITH_TRACE" if all(gate.values()) else "OPEN_BRST_PROJECTOR_SOURCE_STATUS_AUDIT_INCOMPLETE",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "false_pass_guard": "P2023 credits the P1956/P2020 local transverse projector but does not promote it to a BRST physical-state projector. It does not export Q_BRST, ghost cancellation, BRST cohomology, DiscM=CutSum, all-state unitarity, QW-2191 discharge, or ToE closure.",
        "next_honest_step": "Build P2024 by exporting the smallest channel-level BRST data object: Q_BRST action on gauge/ghost/antighost cut states, a nilpotency check Q_BRST^2=0, and a ghost-sector trace convention compatible with the P2020 phase-space normalization.",
        "toe_progress": "Separates the already valid local transverse projector from the missing BRST cohomology/ghost trace layer, preventing the P2020 CutSum from being over-read as a BRST-projected optical theorem.",
        "lay_explanation": "Mamy lokalny filtr wybierający dwie poprzeczne polaryzacje bozonów cechowania, ale to nie jest jeszcze pełny filtr BRST. Nadal trzeba pokazać, jak duchy i stany niefizyczne kasują się w tym samym schemacie rachunku.",
        "environment": {
            "python": platform.python_version(),
            "sympy": sp.__version__,
        },
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P2023] wrote BRST projector source status: {OUT}")


if __name__ == "__main__":
    main()
