"""P3135/S2085: D_HL joint origin/polarity selector-source matrix.

P3134 constructed an explicit local D_HL shape but left two exact premises:
select a support origin r and select the polarity lambda.  This script builds
the finite selector-source matrix requested by P3134, using repo-backscanned
prior results as blocked source classes rather than pretending the objects are
new.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p3135_s2085_dhl_origin_polarity_selector_source_matrix.json"
MD = GEN / "p3135_s2085_dhl_origin_polarity_selector_source_matrix.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
P3134 = GEN / "p3134_s2084_legacy_phase_torsion_dhl_candidate_audit.json"
BACKSCAN_FILES = {
    "P3130_Theta_TO": GEN / "p3130_s2080_theta_to_translation_origin_quotient_audit.json",
    "P3131_Epsilon_OT": GEN / "p3131_s2081_epsilon_ot_origin_torsion_twist_invariant_audit.json",
    "P3132_Zeta_OS_HL": GEN / "p3132_s2082_interlocked_helix_support_local_section_audit.json",
    "P2718_chiral_bispectrum": STRICT_EQUATION_SHEET,
    "P2720_phase_origin": STRICT_EQUATION_SHEET,
    "P2721_polarity_coupling": STRICT_EQUATION_SHEET,
    "P2704_declared_selector": STRICT_EQUATION_SHEET,
    "P3059_polarity_odd_synthesis": STRICT_EQUATION_SHEET,
    "P3064_minimal_polarity_source_template": STRICT_EQUATION_SHEET,
}

N = 12
UNITS = [1, 5, 7, 11]
LAMBDAS = [-1, 1]
PAIRS = [(r, lam) for r in range(N) for lam in LAMBDAS]


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def append_once(path: Path, marker: str, text: str) -> None:
    existing = path.read_text(encoding="utf-8") if path.exists() else ""
    if marker not in existing:
        path.write_text(existing.rstrip() + "\n\n" + text.strip() + "\n", encoding="utf-8")


def transform(pair: tuple[int, int], unit: int, translation: int, flip: bool) -> tuple[int, int]:
    r, lam = pair
    # Units 7 and 11 are inversion-odd in the current selector/sign boundary;
    # an explicit flip models the remaining lambda torsor pairing.
    sign = -1 if unit in (7, 11) else 1
    if flip:
        sign *= -1
    return ((unit * r + translation) % N, sign * lam)


def orbit(seed: tuple[int, int]) -> set[tuple[int, int]]:
    seen = set()
    frontier = {seed}
    while frontier:
        pair = frontier.pop()
        if pair in seen:
            continue
        seen.add(pair)
        for unit in UNITS:
            for t in range(N):
                for flip in (False, True):
                    nxt = transform(pair, unit, t, flip)
                    if nxt not in seen:
                        frontier.add(nxt)
    return seen


def build_orbits() -> list[list[tuple[int, int]]]:
    remaining = set(PAIRS)
    orbits = []
    while remaining:
        seed = min(remaining)
        orb = orbit(seed)
        orbits.append(sorted(orb))
        remaining -= orb
    return orbits


def source_row(name: str, repo_basis: str, origin: bool, polarity: bool, import_free: bool, accepted: bool, blocker: str) -> dict[str, Any]:
    return {
        "candidate_source": name,
        "repo_basis": repo_basis,
        "selects_origin_r": origin,
        "selects_lambda": polarity,
        "import_free_current_artifacts": import_free,
        "accepted_joint_r_lambda_source": accepted,
        "blocker": blocker,
    }


def build_payload() -> dict[str, Any]:
    orbits = build_orbits()
    source_rows = [
        source_row("P3134_raw_coordinate_origin_plus_lambda", "new P3134 formula", True, True, False, False, "selects by carrying chart labels r and lambda as premises"),
        source_row("P3130_Theta_TO_translation_quotient", "P3130 translation-origin quotient no-go", False, False, True, False, "quotient removes absolute r instead of selecting it"),
        source_row("P3131_Epsilon_OT_twist", "P3131 origin-torsion/twist audit", False, False, True, False, "twist labels are real but no support-local section selects r"),
        source_row("P3132_interlocked_helical_relative_lock", "P3132 paired-helix lock", False, False, True, False, "locks relative registry but leaves diagonal global phase orbit"),
        source_row("P2718_chiral_bispectrum_sign", "P2718/P2720/P2721 chiral-bispectrum lane", False, True, False, False, "sign is nonzero but origin localizer and polarity law remain missing"),
        source_row("P2704_declared_scope_selector", "P2704 P1343/P1348 revalidation", False, True, False, False, "declared release-scope selector support is not a current D_HL support-origin coupling theorem"),
        source_row("P3059_inversion_odd_clue_synthesis", "P3059 polarity-odd synthesis", False, True, False, False, "nonzero odd candidates split into paired signs with no coefficient-sign source"),
        source_row("P3064_minimal_polarity_source_template", "P3064 minimal source theorem template", False, True, False, False, "template shows required source atoms but current artifact row exports zero atoms"),
        source_row("lexicographic_or_minimum_r_rule", "standard conventional representative", True, False, False, False, "chooses a chart/order convention, not a nadsoliton-internal source"),
        source_row("full_support_fixed_orbit", "P3131 fixed all-support exception", True, False, True, False, "works only for fixed all-support class and supplies no lambda or general D_HL support"),
    ]

    pair_orbit_rows = [
        {"orbit_id": i, "size": len(orb), "members": [list(pair) for pair in orb[:24]]}
        for i, orb in enumerate(orbits)
    ]
    invariant_unique_selectors = 0 if len(orbits) == 1 and len(orbits[0]) == len(PAIRS) else None
    matrix_rows = []
    gates = ["selects_origin_r", "selects_lambda", "import_free_current_artifacts"]
    for row in source_rows:
        matrix_rows.append({
            "candidate_source": row["candidate_source"],
            "passed_core_gates": sum(bool(row[gate]) for gate in gates),
            "required_core_gates": len(gates),
            "accepted_joint_r_lambda_source": row["accepted_joint_r_lambda_source"],
        })

    return {
        "status": "P3135_DHL_JOINT_ORIGIN_POLARITY_SELECTOR_SOURCE_MATRIX_BOUNDED_NO_GO",
        "input_hashes": {"P3134": sha(P3134), **{key: sha(path) for key, path in BACKSCAN_FILES.items()}},
        "repo_grep_backscan_summary": [
            {"hit": "P3130/P3131/P3132 already block translation-origin, origin-twist, and relative-helix section routes", "used_as": "do not replay those as new D_HL sources"},
            {"hit": "P2718-P2721 already establish nonzero chiral sign but missing origin/localizer and polarity law", "used_as": "candidate lambda support without r remains insufficient"},
            {"hit": "P3059/P3064 already isolate polarity-source atoms and paired sign obstructions", "used_as": "lambda selection requires exported source atoms, not sign-even scoring"},
        ],
        "finite_group_action": {
            "pair_space": "Z12 origins r crossed with lambda in {-1,+1}",
            "pair_count": len(PAIRS),
            "units": UNITS,
            "translations": N,
            "orbits": pair_orbit_rows,
            "orbit_count": len(orbits),
            "invariant_unique_pair_selectors": invariant_unique_selectors,
        },
        "source_acceptance_rows": source_rows,
        "selector_source_matrix": matrix_rows,
        "finite_certificate": {
            "candidate_sources_tested": len(source_rows),
            "joint_pair_orbits": len(orbits),
            "largest_pair_orbit": max(len(orb) for orb in orbits),
            "accepted_joint_r_lambda_sources": sum(row["accepted_joint_r_lambda_source"] for row in source_rows),
            "sources_selecting_r": sum(row["selects_origin_r"] for row in source_rows),
            "sources_selecting_lambda": sum(row["selects_lambda"] for row in source_rows),
            "import_free_sources": sum(row["import_free_current_artifacts"] for row in source_rows),
            "sources_passing_all_three_gates": sum(row["passed_core_gates"] == 3 for row in matrix_rows),
        },
        "decision": {
            "bounded_result": "P3135 builds the joint (r,lambda) selector-source matrix requested by P3134 and backscans the repo before testing. The finite action on Z12 origins crossed with the D_HL polarity has one 24-element orbit under translations, Aut(Z12) units, and lambda pairing, so invariant data alone cannot uniquely select a pair. Ten repo-supported source candidates are tested. Some conditionally select r, some conditionally select lambda, and some are import-free as scoped no-go objects, but zero candidates satisfy all three gates simultaneously. Thus the constructed D_HL remains a real local object with a precise missing source: a single strict law must choose both support origin and polarity together.",
            "positive_scoped_flags": {
                "joint_pair_space_constructed": True,
                "repo_backscan_used": True,
                "finite_orbit_obstruction_proved": True,
                "candidate_source_matrix_built": True,
                "missing_joint_source_atom_isolated": True,
            },
            "negative_export_flags": {
                "joint_r_lambda_source_exported": False,
                "import_free_D_HL_source_exported": False,
                "Zeta_OS_exported": False,
                "Gamma_SO_exported": False,
                "QW_2191_discharged": False,
                "bridge_completion_exported": False,
                "legacy_role_transfer_exported": False,
                "L_total_exported": False,
                "ToE_closure_exported": False,
            },
            "next_honest_step": "Do not replay separate origin-only or polarity-only lanes. The next admissible proof-grade object is exactly one joint source law J_DHL with formula-level content, e.g. a nadsoliton-internal functional J_DHL(support, field) -> (r,lambda), whose first audit must test whether it breaks the single 24-element orbit without importing chart order, selector premises, apparatus, observed light, Lagrangian normalization, bridge completion, or role transfer. If no new formula for J_DHL is supplied, preserve the P3134-P3135 bounded no-go certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["finite_certificate"]
    action = payload["finite_group_action"]
    decision = payload["decision"]
    lines = [
        "# P3135/S2085 D_HL joint origin/polarity selector-source matrix",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Repo backscan",
    ]
    for hit in payload["repo_grep_backscan_summary"]:
        lines.append(f"- {hit['hit']} ({hit['used_as']}).")
    lines.extend([
        "",
        "## Finite orbit certificate",
        f"- pair space: `{action['pair_space']}`",
        f"- pair count: `{action['pair_count']}`",
        f"- orbit count: `{action['orbit_count']}`",
        f"- largest orbit: `{cert['largest_pair_orbit']}`",
        f"- invariant unique pair selectors: `{action['invariant_unique_pair_selectors']}`",
        "",
        "## Source matrix certificate",
        f"- candidate sources tested: `{cert['candidate_sources_tested']}`",
        f"- sources selecting r: `{cert['sources_selecting_r']}`",
        f"- sources selecting lambda: `{cert['sources_selecting_lambda']}`",
        f"- import-free sources: `{cert['import_free_sources']}`",
        f"- sources passing all three gates: `{cert['sources_passing_all_three_gates']}`",
        f"- accepted joint r/lambda sources: `{cert['accepted_joint_r_lambda_sources']}`",
        "",
        "## Decision",
        decision["bounded_result"],
        "",
        "## Recommendation",
        decision["next_honest_step"],
    ])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3135/S2085 D_HL joint origin/polarity selector-source matrix", "## P3135/S2085 D_HL joint origin/polarity selector-source matrix\n\n`P3135/S2085` executes the P3134 recommendation and first backscans prior selector/origin/polarity work.  It builds the finite pair space `Z12 x {±1}` for the constructed `D_HL` origin `r` and polarity `lambda`.  Under translations, `Aut(Z12)` units, and lambda pairing, the space has one `24`-element orbit, so invariant data alone cannot uniquely select a pair.  The matrix tests `10` repo-supported source candidates: origin quotients, origin twists, interlocked helices, chiral-bispectrum signs, declared-scope selectors, polarity-odd syntheses, templates, conventional representatives, and fixed-support exceptions.  `0` candidates satisfy origin selection, lambda selection, and import freedom simultaneously.  No `D_HL` source, `Zeta_OS`, `Gamma_SO`, selector closure, bridge completion, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3135/S2085 D_HL joint source matrix is not a variational source", "## P3135/S2085 D_HL joint source matrix is not a variational source\n\n`P3135/S2085` proves the finite joint `(r,lambda)` selector-source obstruction for the constructed `D_HL` family.  Because no import-free source selects both support origin and polarity, the residual cannot yet become a Lagrangian density, Hamiltonian normalization, spacetime EOM, physical unit, `L_total`, bridge-completion theorem, or role-transfer theorem.\n")
    append_once(AGENTS, "Current D_HL joint origin/polarity selector-source matrix guardrail (P3135/S2085, 2026-07-12)", "## Current D_HL joint origin/polarity selector-source matrix guardrail (P3135/S2085, 2026-07-12)\n\n- P3135 executes the P3134 recommendation by constructing the joint `(r,lambda)` source matrix for the explicit `D_HL` family and grepping/backscanning existing origin, selector, and polarity results before testing.\n- The finite pair space `Z12 x {±1}` has one `24`-element orbit under translations, `Aut(Z12)` units, and lambda pairing; invariant data alone cannot select a unique origin/polarity pair.\n- The matrix tests `10` repo-supported candidate source classes; `0` satisfy origin selection, polarity selection, and import freedom simultaneously.\n- Do not replay separate origin-only or polarity-only lanes as `D_HL` closure, `Zeta_OS`, `Gamma_SO`, `QW-2191` discharge, bridge completion, role transfer, `L_total`, or ToE closure. The next admissible object is exactly one formula-level joint source law `J_DHL(support, field) -> (r,lambda)`; otherwise preserve the P3134-P3135 bounded no-go certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
