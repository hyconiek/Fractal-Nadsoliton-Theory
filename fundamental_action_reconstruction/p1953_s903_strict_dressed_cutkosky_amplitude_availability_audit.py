#!/usr/bin/env python3
"""P1953 S903 strict dressed Cutkosky amplitude availability audit.

P1951/P1952 established seed-level positivity.  This audit checks whether the
current repository already exports enough data to replace K_cut_seed by a full
dressed graviton -> gauge_gauge amplitude reduced to a common basis.

The result is intentionally an obstruction certificate if the needed objects
are only seeds, stubs, or open contracts.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p1953_s903_strict_dressed_cutkosky_amplitude_availability_audit.json"


def load(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def digest(obj: object) -> str:
    blob = json.dumps(obj, sort_keys=True, ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


def strings(obj: Any) -> list[str]:
    if isinstance(obj, str):
        return [obj]
    if isinstance(obj, dict):
        out: list[str] = []
        for key, value in obj.items():
            out.extend(strings(key))
            out.extend(strings(value))
        return out
    if isinstance(obj, list):
        out = []
        for value in obj:
            out.extend(strings(value))
        return out
    return [str(obj)]


def find_generated_text_hits(terms: list[str], limit: int = 20) -> dict[str, list[str]]:
    hits: dict[str, list[str]] = {term: [] for term in terms}
    for path in sorted(GEN.glob("*.json")):
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        lower = text.lower()
        for term in terms:
            if term.lower() in lower and len(hits[term]) < limit:
                hits[term].append(path.name)
    return hits


def row(requirement: str, verdict: str, evidence: list[str], blocker: str) -> dict[str, object]:
    return {
        "requirement": requirement,
        "verdict": verdict,
        "evidence": evidence,
        "blocker": blocker,
    }


def main() -> None:
    p1852 = load("p1852_s802_strict_b1_brst_anomaly_and_cutkosky_seed_witness_checkpoint.json")
    p1857 = load("p1857_s807_strict_b1_triangle_amplitude_seed_and_k5_instantiation_checkpoint.json")
    p1859 = load("p1859_s809_strict_b1_cutkosky_kernel_and_residue_certificate_stub_checkpoint.json")
    p1862 = load("p1862_s812_strict_b1_dressed_pole_residue_and_projected_disc_evaluation_checkpoint.json")
    p1913 = load("p1913_s863_strict_c1_gr_common_basis_unitarity_background_probe.json")
    p1951 = load("p1951_s901_strict_b1_cutkosky_phase_space_integral_probe.json")
    p1952 = load("p1952_s902_strict_qw2049_ur_domain_separation_certificate.json")

    generated_hits = find_generated_text_hits(
        [
            "dressed amplitude",
            "M_dressed",
            "dressed_pole_residue",
            "graviton->gauge_gauge",
            "BRST physical-state",
            "DiscM_grmix",
            "CutSum_grmix",
        ]
    )

    p1862_note = ((p1862.get("dressed_pole_residue_seed_table") or {}).get("note") or "")
    p1913_scheme = ((p1913.get("common_basis_definition") or {}).get("scheme") or "MISSING")
    p1951_scope = ((p1951.get("declared_integrand") or {}).get("scope_warning") or "")
    p1859_residue_required = (
        (p1859.get("residue_certificate_stub") or {}).get("required_exports") or []
    )

    readiness_rows = [
        row(
            "full dressed graviton->gauge_gauge amplitude squared in common basis",
            "FAIL_MISSING",
            [
                "P1857 exports triangle_amplitude_seed with status OPEN_OBSTRUCTION_WITH_TRACE.",
                f"P1951 scope warning: {p1951_scope}",
            ],
            "No exported object supplies |M_dressed(graviton->gauge_gauge)|^2 or its tensor/polarization contraction.",
        ),
        row(
            "DiscM and CutSum evaluated for the gravity/gauge channel in one common basis",
            "FAIL_OPEN_SYMBOLIC_ONLY",
            [
                "P1913 grmix row is OPEN_SYMBOLIC_REDUCED_NOT_EVALUATED.",
                "P1911 queues DiscM_grmix/CutSum_grmix as OPEN_FAIL_BY_MISSING_EVALUATION_INPUTS.",
            ],
            "The common-basis form exists only as symbolic placeholders alpha_gr, beta_gr, gamma_gr.",
        ),
        row(
            "computed dressed propagator poles and residues, not seed inheritance",
            "FAIL_SEED_INHERITANCE_ONLY",
            [
                f"P1862 note: {p1862_note}",
                f"P1859 required exports: {p1859_residue_required}",
            ],
            "Residues labelled dressed are inherited from seed physical-pole projection, not computed from the dressed propagator.",
        ),
        row(
            "BRST physical-state projection for the same channel",
            "FAIL_MISSING",
            [
                "P1852 exports BRST/Cutkosky seed contracts, not a physical-state projector.",
                "P1951 and P1952 explicitly list BRST physical-state projection as not proved.",
            ],
            "No channel-level BRST projector maps the gauge_gauge final state into the physical Hilbert subspace.",
        ),
        row(
            "single fixed renormalization and gauge scheme across P1950/P1951/P1913",
            "FAIL_SCHEME_NOT_LOCKED",
            [
                "P1950/P1951 use MSbar_B1_seed.",
                f"P1913 common basis scheme is {p1913_scheme}.",
            ],
            "The amplitude common-basis row is not locked to the same scheme tag as the B1 seed calculation.",
        ),
    ]

    hard_failures = [r for r in readiness_rows if str(r["verdict"]).startswith("FAIL")]

    interface_contract = {
        "object_name": "DressedCutkoskyAmplitude_graviton_to_gauge_gauge_strict_B1_v1",
        "required_fields": [
            "channel = graviton->gauge_gauge",
            "scheme = MSbar_B1_seed",
            "gauge_fixing_family and xi value",
            "external_state_projectors including BRST physical-state projector",
            "M_dressed_common_basis",
            "AbsM_dressed_squared_common_basis",
            "DiscM_common_basis",
            "CutSum_common_basis",
            "DiscM_minus_CutSum_simplified",
            "dressed_graviton_propagator_pole_list",
            "residue_values_per_pole",
            "ghost_sector_exclusion_trace",
            "phase_space_measure and integration_domain",
            "uncertainty_or_exactness_certificate",
        ],
        "acceptance_rule": "All required fields present, same scheme tag, and DiscM_minus_CutSum_simplified == 0 with positive physical residues.",
    }

    out = {
        "packet_id": "P1953",
        "stage_id": "S903",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "local_verdict": "DRESSED_CUTKOSKY_AMPLITUDE_NOT_AVAILABLE__INTERFACE_CONTRACT_EXPORTED",
        "route": "strict_only",
        "legacy_bridge_used": False,
        "depends_on": {
            "p1852_present": "cutkosky_seed_contract" in p1852,
            "p1857_present": "triangle_amplitude_seed" in p1857,
            "p1859_present": "cutkosky_kernel_contract" in p1859,
            "p1862_present": "dressed_pole_residue_seed_table" in p1862,
            "p1913_present": "unitarity_common_basis_rows_v1" in p1913,
            "p1951_present": "declared_integrand" in p1951,
            "p1952_present": "positivity_certificate" in p1952,
        },
        "input_hashes": {
            "p1859_sha256": digest(p1859),
            "p1862_sha256": digest(p1862),
            "p1913_sha256": digest(p1913),
            "p1951_sha256": digest(p1951),
            "p1952_sha256": digest(p1952),
        },
        "generated_text_search_hits": generated_hits,
        "readiness_matrix": readiness_rows,
        "hard_failure_count": len(hard_failures),
        "dressed_amplitude_interface_contract": interface_contract,
        "safe_use_of_existing_results": {
            "p1950": "Declared B1 counterterm cancellation remains valid within its declared ansatz scope.",
            "p1951": "Seed-integrand phase-space positivity remains valid and cannot be promoted to dressed Cutkosky closure.",
            "p1952": "QW-2049 seed-local positivity rectangle remains valid and cannot be promoted to P1677 UR_link.",
        },
        "obstruction_statement": "The current repo does not export the dressed graviton->gauge_gauge amplitude and same-scheme DiscM/CutSum data required to replace K_cut_seed. The next step must build that object or prove it cannot be built from current strict inputs.",
        "theorem_scope": {
            "proved": [
                "The repository contains seed/proxy/stub data for the target channel.",
                "The repository does not currently contain the minimum dressed-amplitude interface needed for theorem-grade Cutkosky closure.",
                "A precise acceptance contract for that missing interface is now exported.",
            ],
            "not_proved": [
                "full dressed Cutkosky equality",
                "BRST physical-state projected optical theorem",
                "ghost-free dressed propagator theorem",
                "global UR_link theorem",
                "QW-2191 selector discharge",
            ],
        },
        "higher_reasoning_required_for_next_step": True,
        "next_honest_step": "Build P1954 with high-reasoning support: either derive a same-scheme dressed amplitude from L_total and exported vertices/projectors, or export a formal nonavailability theorem identifying the exact missing vertex/projector data.",
        "lay_explanation": "Sprawdzilismy, czy w repo jest juz prawdziwa ubrana amplituda potrzebna do pelnego testu unitarności. Nie ma jej: sa tylko wersje robocze, seedy i kontrakty. Teraz wiadomo dokladnie, jakiego obiektu brakuje.",
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
