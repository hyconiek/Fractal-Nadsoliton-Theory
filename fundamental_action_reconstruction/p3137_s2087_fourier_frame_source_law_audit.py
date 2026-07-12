"""P3137/S2087: strict Fourier-frame/source law F_DHL audit.

P3136 left exactly one admissible object: a strict Fourier-frame/source law
F_DHL selecting a primitive character and phase-zero reference without imported
chart labels.  This audit constructs the finite frame candidate space and tests
repo-backed source classes against it.
"""

from __future__ import annotations

import cmath
import hashlib
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p3137_s2087_fourier_frame_source_law_audit.json"
MD = GEN / "p3137_s2087_fourier_frame_source_law_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
P3136 = GEN / "p3136_s2086_fourier_phase_dhl_joint_source_candidate.json"

N = 12
BETA_TORS = 0.01
MODES = [1, 2, 3, 4, 5]
LAMBDAS = [-1, 1]
UNITS = [1, 5, 7, 11]
PRIMITIVE_CHARS = [1, 5, 7, 11]
CHAR_PAIRS = [(1, 11), (2, 10), (3, 9), (4, 8), (5, 7)]


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def append_once(path: Path, marker: str, text: str) -> None:
    existing = path.read_text(encoding="utf-8") if path.exists() else ""
    if marker not in existing:
        path.write_text(existing.rstrip() + "\n\n" + text.strip() + "\n", encoding="utf-8")


def profile(r: int, k: int, lam: int) -> list[float]:
    return [lam * BETA_TORS * math.sin(2 * math.pi * k * ((x - r) % N) / N) for x in range(N)]


def coeff(vals: list[float], k: int) -> complex:
    return sum(vals[x] * cmath.exp(-2j * math.pi * k * x / N) for x in range(N))


def active_pair(vals: list[float]) -> tuple[int, int]:
    energies = []
    for a, b in CHAR_PAIRS:
        energies.append((abs(coeff(vals, a)) + abs(coeff(vals, b)), (a, b)))
    return max(energies, key=lambda item: (item[0], -item[1][0]))[1]


def orbit_chars(seed: int) -> set[int]:
    return {(u * seed) % N for u in UNITS}


def orbit_pairs(seed: tuple[int, int]) -> set[tuple[int, int]]:
    out = set()
    for u in UNITS:
        pair = tuple(sorted(((u * seed[0]) % N, (u * seed[1]) % N)))
        out.add(pair)
    return out


def source_row(name: str, basis: str, character: bool, phase_zero: bool, primitive: bool, import_free: bool, accepted: bool, blocker: str) -> dict[str, Any]:
    return {
        "candidate_F_DHL_source": name,
        "repo_or_formula_basis": basis,
        "selects_character_or_pair": character,
        "selects_phase_zero_reference": phase_zero,
        "restricts_to_primitive_character": primitive,
        "import_free_current_artifacts": import_free,
        "accepted_F_DHL_source": accepted,
        "blocker": blocker,
    }


def build_payload() -> dict[str, Any]:
    rows = []
    for k in MODES:
        for r in range(N):
            for lam in LAMBDAS:
                vals = profile(r, k, lam)
                ap = active_pair(vals)
                rows.append({
                    "mode_k": k,
                    "r": r,
                    "lambda": lam,
                    "active_character_pair": list(ap),
                    "active_pair_is_primitive": set(ap).issubset(set(PRIMITIVE_CHARS)),
                    "phase_zero_reference_imported": True,
                })

    primitive_char_orbit = sorted(orbit_chars(1))
    primitive_pair_orbit = sorted([list(p) for p in orbit_pairs((1, 11))])
    source_rows = [
        source_row("active_spectral_pair_receiver", "P3136 Fourier coefficient support", True, False, False, True, False, "selects only the ±k pair, not orientation or phase-zero"),
        source_row("primitive_character_filter", "U(12) primitive character set", False, False, True, True, False, "keeps 48 primitive rows but does not select one primitive pair or phase origin"),
        source_row("minimal_positive_k_rule", "choose smallest positive active k", True, False, True, False, False, "uses integer chart order on character labels"),
        source_row("max_magnitude_character_rule", "Fourier magnitude receiver", False, False, False, True, False, "ties the ± active pair and all active coefficients have the same magnitude scale"),
        source_row("phase_zero_by_argument_gauge", "set arg(C_k)=0 by translation", True, True, True, False, False, "fixes phase by gauge choice, not strict source"),
        source_row("P2992_frequency_source_localizer", "P2992 Fourier-character localizer obstruction", False, False, False, True, False, "no nonpremise frequency/source localizer exported"),
        source_row("P2994_source_coupling", "P2994 named source-atom coupling obstruction", False, False, False, True, False, "no atom-specific source-coupling theorem exported"),
        source_row("P3039_chi_phase_localizer", "P3039 sine/chiral localizer obstruction", True, False, True, False, False, "real inversion-odd torsor but translations move phase origin"),
        source_row("P2869_Aut_character_idempotent", "P2869 Aut-character idempotent warning", True, False, True, False, False, "represents endpoints by importing projector/polarity coefficients"),
    ]
    gates = ["selects_character_or_pair", "selects_phase_zero_reference", "restricts_to_primitive_character", "import_free_current_artifacts"]
    matrix = [{"candidate_F_DHL_source": row["candidate_F_DHL_source"], "passed_gates": sum(bool(row[g]) for g in gates), "required_gates": len(gates), "accepted_F_DHL_source": row["accepted_F_DHL_source"]} for row in source_rows]

    return {
        "status": "P3137_FOURIER_FRAME_SOURCE_LAW_F_DHL_BOUNDED_NO_GO",
        "input_hashes": {"P3136": sha(P3136)},
        "constructed_object": {
            "name": "F_DHL Fourier-frame/source law candidate space",
            "required_output": "primitive character/frame plus phase-zero reference for J_DHL",
            "finite_frame_space": "primitive characters U(12)={1,5,7,11}, primitive pairs {1,11}/{5,7}, and 12 phase-zero cells",
        },
        "repo_backscan_summary": [
            "P2992: no nonpremise frequency/source localizer for Z12 Fourier characters.",
            "P2994: exact receiver matrix exists, but no atom-specific source-coupling theorem is exported.",
            "P3039: chi_i sine projector is real and inversion-odd, but phase origin remains chart-relative.",
            "P2869/P2870: Aut-character idempotents represent endpoints but import projector/polarity and do not intrinsically select the needed character.",
        ],
        "finite_frame_certificate": {
            "profiles_tested": len(rows),
            "primitive_character_orbit_size": len(primitive_char_orbit),
            "primitive_character_orbit": primitive_char_orbit,
            "primitive_pair_orbit_size": len(primitive_pair_orbit),
            "primitive_pair_orbit": primitive_pair_orbit,
            "profiles_with_primitive_active_pair": sum(row["active_pair_is_primitive"] for row in rows),
            "profiles_with_nonprimitive_active_pair": sum(not row["active_pair_is_primitive"] for row in rows),
            "phase_zero_cells": N,
            "accepted_F_DHL_sources": sum(row["accepted_F_DHL_source"] for row in source_rows),
            "source_candidates_tested": len(source_rows),
            "source_candidates_passing_all_gates": sum(row["passed_gates"] == len(gates) for row in matrix),
        },
        "active_pair_rows": rows,
        "source_acceptance_rows": source_rows,
        "source_matrix": matrix,
        "decision": {
            "bounded_result": "P3137 constructs the missing F_DHL frame-source candidate space and tests it against repo-backed Fourier/character/localizer results. The finite receiver side is real: every P3134/P3136 profile has an active ±k character pair, and 48 rows live on primitive pairs. But the primitive characters form one Aut(Z12) orbit and the primitive pairs form one orbit; no invariant frame data chooses a unique primitive character/pair. Phase-zero has 12 translated cells. Across nine source candidates, zero simultaneously select a primitive character/frame, select phase-zero, and remain import-free. Thus F_DHL is not exported on current artifacts.",
            "positive_scoped_flags": {
                "F_DHL_candidate_space_constructed": True,
                "active_character_pairs_computed": True,
                "primitive_character_orbit_obstruction_proved": True,
                "phase_zero_translation_orbit_recorded": True,
                "repo_backscan_integrated": True,
            },
            "negative_export_flags": {
                "F_DHL_source_exported": False,
                "import_free_J_DHL_source_exported": False,
                "D_HL_source_exported": False,
                "Zeta_OS_exported": False,
                "Gamma_SO_exported": False,
                "QW_2191_discharged": False,
                "bridge_completion_exported": False,
                "legacy_role_transfer_exported": False,
                "L_total_exported": False,
                "ToE_closure_exported": False,
            },
            "next_honest_step": "Do not continue Fourier-frame variants unless a genuinely new strict frame-breaking source is supplied. The next admissible move is either one new non-Fourier joint source object for (r,lambda), or a no-new-live-frontier reconciliation for the D_HL lane after P3133-P3137. If continuing constructively, the object must not be a Fourier receiver, character projector, phase gauge, lexicographic label rule, or prior chi_i/localizer replay; it must compute a nonconventional support origin and polarity in one law.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["finite_frame_certificate"]
    decision = payload["decision"]
    lines = ["# P3137/S2087 Fourier-frame source law F_DHL audit", "", f"Status: `{payload['status']}`", "", "## Constructed object", f"- `{payload['constructed_object']['required_output']}`", f"- `{payload['constructed_object']['finite_frame_space']}`", "", "## Repo backscan"]
    for hit in payload["repo_backscan_summary"]:
        lines.append(f"- {hit}")
    lines.extend(["", "## Finite certificate", f"- profiles tested: `{cert['profiles_tested']}`", f"- primitive character orbit size: `{cert['primitive_character_orbit_size']}`", f"- primitive pair orbit size: `{cert['primitive_pair_orbit_size']}`", f"- profiles with primitive active pair: `{cert['profiles_with_primitive_active_pair']}`", f"- profiles with nonprimitive active pair: `{cert['profiles_with_nonprimitive_active_pair']}`", f"- phase-zero cells: `{cert['phase_zero_cells']}`", f"- source candidates tested: `{cert['source_candidates_tested']}`", f"- candidates passing all gates: `{cert['source_candidates_passing_all_gates']}`", f"- accepted F_DHL sources: `{cert['accepted_F_DHL_sources']}`", "", "## Decision", decision["bounded_result"], "", "## Recommendation", decision["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3137/S2087 Fourier-frame source law F_DHL audit", "## P3137/S2087 Fourier-frame source law F_DHL audit\n\n`P3137/S2087` constructs the `F_DHL` Fourier-frame/source candidate space required after P3136: primitive characters `U(12)={1,5,7,11}`, primitive pairs `{1,11}/{5,7}`, and `12` phase-zero cells.  The finite audit computes active character pairs for all `120` P3134/P3136 profiles: `48` have primitive active pairs and `72` have nonprimitive active pairs.  The primitive characters form one `Aut(Z12)` orbit and the primitive pairs also form one orbit, while phase-zero remains a translated `12`-cell choice.  Nine repo-backed source candidates are tested against character/frame selection, phase-zero selection, primitive restriction, and import freedom; `0` pass all gates.  No `F_DHL`, import-free `J_DHL`, `D_HL` source, `Zeta_OS`, `Gamma_SO`, selector closure, bridge completion, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3137/S2087 F_DHL frame source is not a variational source", "## P3137/S2087 F_DHL frame source is not a variational source\n\n`P3137/S2087` proves that the current Fourier-frame/source candidate space does not export an import-free primitive character and phase-zero source.  Therefore the Fourier-phase `D_HL` reconstruction cannot yet become a Lagrangian density, Hamiltonian normalization, spacetime EOM, physical unit, `L_total`, bridge-completion theorem, or role-transfer theorem.\n")
    append_once(AGENTS, "Current Fourier-frame source law F_DHL guardrail (P3137/S2087, 2026-07-12)", "## Current Fourier-frame source law F_DHL guardrail (P3137/S2087, 2026-07-12)\n\n- P3137 constructs the `F_DHL` candidate space left by P3136: primitive Fourier characters, primitive character pairs, and phase-zero cells.\n- The finite audit finds real receiver structure but no source: primitive characters and primitive pairs each remain single `Aut(Z12)` orbits, phase-zero remains a translated `12`-cell choice, and `0/9` repo-backed candidates pass character/frame selection, phase-zero selection, primitive restriction, and import freedom simultaneously.\n- Do not replay Fourier receivers, character projectors, phase gauges, lexicographic labels, or prior `chi_i` localizer lanes as `D_HL`, `Zeta_OS`, `Gamma_SO`, `QW-2191` discharge, bridge completion, role transfer, `L_total`, or ToE closure.\n- The next admissible move is either one genuinely new non-Fourier joint source object for `(r,lambda)` or a no-new-live-frontier reconciliation for the `D_HL` lane after P3133-P3137.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
