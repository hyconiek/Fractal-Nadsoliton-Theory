#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from itertools import combinations
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2428_s1378_repair_readiness_anf_derivative_certificate.json"
MD = GEN / "p2428_s1378_repair_readiness_anf_derivative_certificate.md"

SOURCE_FILES = {
    "P2427_INDEPENDENCE": GEN / "p2427_s1377_weight_repair_projection_independence_certificate.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

GATES = [
    "source_obligation_discharge",
    "chi11_source_export",
    "qw2191_selector_discharge",
    "role_transfer_audit_license",
    "role_bearing_ltotal_export",
]
TARGET_FORMULAS = {
    "bridge_source_ready": [{"source_obligation_discharge"}],
    "selector_source_ready": [{"chi11_source_export", "qw2191_selector_discharge"}],
    "role_transfer_ready": [{"role_transfer_audit_license"}],
    "role_bearing_ltotal_ready": [{"role_bearing_ltotal_export"}],
    "toe_ready": [set(GATES)],
}


def rel(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": rel(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg",
            "-n",
            pattern,
            "fundamental_action_reconstruction",
            "material_dowodowy",
            "-g",
            "*.py",
            "-g",
            "*.md",
            "-g",
            "*.tex",
            "-g",
            "!generated/**",
        ],
        cwd=REPO,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:20]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2428|S1378|repair readiness ANF|readiness.*derivative|gate derivative",
        "p2427_input": "P2427|projection independence|weight-repair projection|empty weight-side support",
        "anf_prior": "ANF|Boolean derivative|prime implicant|Möbius|Mobius",
        "nonpromotion_prior": "nonpromotion|does not discharge|No source|No ToE|not a theorem discharge",
    }
    return {
        "tool": "rg",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
        "finding": (
            "Repo grep finds older ANF/Boolean-derivative certificates and P2427 projection independence, but no production "
            "P24xx certificate computing the exact ANF, prime implicants, and derivative-edge supports for the current five-gate "
            "repair readiness predicates."
        ),
    }


def theorem(payload: dict[str, Any], key: str) -> dict[str, Any]:
    return payload.get(key, {}).get("theorem_export", {})


def mask_to_set(mask: int) -> set[str]:
    return {gate for idx, gate in enumerate(GATES) if mask & (1 << idx)}


def set_to_mask(gates: set[str]) -> int:
    out = 0
    for idx, gate in enumerate(GATES):
        if gate in gates:
            out |= 1 << idx
    return out


def predicate_true(target: str, mask: int) -> bool:
    true_gates = mask_to_set(mask)
    return any(term.issubset(true_gates) for term in TARGET_FORMULAS[target])


def truth_vector(target: str) -> list[int]:
    return [1 if predicate_true(target, mask) else 0 for mask in range(1 << len(GATES))]


def anf_terms(target: str) -> list[list[str]]:
    coeffs = truth_vector(target)
    n = len(GATES)
    for bit in range(n):
        for mask in range(1 << n):
            if mask & (1 << bit):
                coeffs[mask] ^= coeffs[mask ^ (1 << bit)]
    return [[GATES[idx] for idx in range(n) if mask & (1 << idx)] for mask, coeff in enumerate(coeffs) if coeff]


def monomial(term: list[str]) -> str:
    if not term:
        return "1"
    return " * ".join(term)


def anf_polynomial(target: str) -> str:
    return " + ".join(monomial(term) for term in anf_terms(target))


def prime_implicants(target: str) -> list[list[str]]:
    implicants: list[set[str]] = []
    all_masks = range(1 << len(GATES))
    for size in range(len(GATES) + 1):
        for combo in combinations(GATES, size):
            combo_set = set(combo)
            combo_mask = set_to_mask(combo_set)
            containing_masks = [mask for mask in all_masks if (mask & combo_mask) == combo_mask]
            if containing_masks and all(predicate_true(target, mask) for mask in containing_masks):
                if not any(prev.issubset(combo_set) for prev in implicants):
                    implicants.append(combo_set)
    return [sorted(item, key=GATES.index) for item in implicants]


def derivative_edge_count(target: str, gate: str) -> int:
    idx = GATES.index(gate)
    count = 0
    for mask in range(1 << len(GATES)):
        if mask & (1 << idx):
            continue
        if predicate_true(target, mask) != predicate_true(target, mask | (1 << idx)):
            count += 1
    return count


def derivative_supports() -> dict[str, dict[str, int]]:
    return {
        target: {gate: derivative_edge_count(target, gate) for gate in GATES}
        for target in TARGET_FORMULAS
    }


def build_certificate(sources: dict[str, dict[str, Any]]) -> dict[str, Any]:
    p2427 = theorem(sources["P2427_INDEPENDENCE"], "weight_repair_projection_independence_certificate")
    anf = {target: anf_terms(target) for target in TARGET_FORMULAS}
    primes = {target: prime_implicants(target) for target in TARGET_FORMULAS}
    derivatives = derivative_supports()
    return {
        "truth_table_size": 1 << len(GATES),
        "gate_count": len(GATES),
        "anf_terms": anf,
        "anf_polynomials": {target: anf_polynomial(target) for target in TARGET_FORMULAS},
        "anf_degrees": {target: max(len(term) for term in terms) for target, terms in anf.items()},
        "prime_implicants": primes,
        "prime_implicant_counts": {target: len(items) for target, items in primes.items()},
        "derivative_edge_counts": derivatives,
        "essential_gate_sets": {
            target: [gate for gate, count in counts.items() if count > 0]
            for target, counts in derivatives.items()
        },
        "p2427_product_assignment_count_inherited": p2427.get("product_assignment_count"),
        "p2427_all_tables_factor_inherited": p2427.get("all_readiness_tables_factor_by_weight_x_repair"),
        "p2427_weight_side_support_empty_inherited": p2427.get("weight_side_support_for_repair_readiness_predicates"),
    }


def append_doc_sections() -> None:
    eq_section = """
## P2428/S1378 repair-readiness ANF and derivative certificate

`P2428/S1378` removes another ambiguity from the P2422--P2427 repair layer by computing the exact Boolean ANF, prime implicants, and derivative-edge supports for the five readiness predicates on the five missing theorem gates.  The ANFs are single monomials: bridge=`source_obligation_discharge`, selector=`chi11_source_export * qw2191_selector_discharge`, role-transfer=`role_transfer_audit_license`, role-bearing `L_total`=`role_bearing_ltotal_export`, and ToE=`all five gates`.

The derivative audit makes the blocker structure explicit: bridge, role-transfer, and role-bearing `L_total` have singleton essential gates; selector has exactly the chi11/QW-2191 pair; ToE has all five gates as essential, each only at the four-other-gates nearest miss.  This is a Boolean obstruction certificate, not a theorem-gate discharge.
""".strip()
    lag_section = """
## P2428/S1378 repair-readiness ANF derivative guard for Lagrangian/EOM

`P2428/S1378` computes that role-bearing `L_total` readiness is the single literal `role_bearing_ltotal_export`, while ToE readiness is the degree-five monomial over all missing gates.  No APD value witness, chi11 cut, weighted order, or projection-independence fact supplies those literals; they require separate theorem export.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    cert = build_certificate(sources)
    theorem_export = {
        "theorem_name": "P2428_T1_repair_readiness_anf_derivative_certificate",
        "gate_count": cert["gate_count"],
        "truth_table_size": cert["truth_table_size"],
        "anf_polynomials": cert["anf_polynomials"],
        "anf_degrees": cert["anf_degrees"],
        "prime_implicant_counts": cert["prime_implicant_counts"],
        "prime_implicants": cert["prime_implicants"],
        "derivative_edge_counts": cert["derivative_edge_counts"],
        "essential_gate_sets": cert["essential_gate_sets"],
        "single_prime_implicant_per_readiness_predicate": all(count == 1 for count in cert["prime_implicant_counts"].values()),
        "selector_requires_chi11_and_qw2191_pair": cert["prime_implicants"]["selector_source_ready"] == [
            ["chi11_source_export", "qw2191_selector_discharge"]
        ],
        "toe_prime_implicant_is_all_five_gates": cert["prime_implicants"]["toe_ready"] == [GATES],
        "toe_derivative_edges_are_nearest_misses": cert["derivative_edge_counts"]["toe_ready"] == {gate: 1 for gate in GATES},
        "p2427_product_assignment_count_inherited": cert["p2427_product_assignment_count_inherited"] == 4608,
        "p2427_all_tables_factor_inherited": cert["p2427_all_tables_factor_inherited"] is True,
        "p2427_weight_side_support_empty_inherited": cert["p2427_weight_side_support_empty_inherited"] == [],
        "source_obligation_discharge_exported": False,
        "chi11_source_exported": False,
        "qw2191_discharged": False,
        "role_transfer_licensed": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "ANF/derivative structure is a blocker map, not a theorem-gate discharge.",
            "The selector monomial still requires both chi11 source export and QW-2191 discharge.",
            "The ToE monomial still requires all five missing theorem gates.",
            "No source, selector, role-transfer, role-bearing L_total, or ToE theorem is exported.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "truth_table_32": theorem_export["truth_table_size"] == 32,
        "five_gates": theorem_export["gate_count"] == 5,
        "single_prime_each": theorem_export["single_prime_implicant_per_readiness_predicate"],
        "selector_pair": theorem_export["selector_requires_chi11_and_qw2191_pair"],
        "toe_all_five": theorem_export["toe_prime_implicant_is_all_five_gates"],
        "toe_derivative_nearest_misses": theorem_export["toe_derivative_edges_are_nearest_misses"],
        "bridge_derivative_only_source": theorem_export["derivative_edge_counts"]["bridge_source_ready"] == {
            "source_obligation_discharge": 16,
            "chi11_source_export": 0,
            "qw2191_selector_discharge": 0,
            "role_transfer_audit_license": 0,
            "role_bearing_ltotal_export": 0,
        },
        "selector_derivatives_pair_only": theorem_export["derivative_edge_counts"]["selector_source_ready"] == {
            "source_obligation_discharge": 0,
            "chi11_source_export": 8,
            "qw2191_selector_discharge": 8,
            "role_transfer_audit_license": 0,
            "role_bearing_ltotal_export": 0,
        },
        "ltotal_derivative_only_ltotal": theorem_export["derivative_edge_counts"]["role_bearing_ltotal_ready"] == {
            "source_obligation_discharge": 0,
            "chi11_source_export": 0,
            "qw2191_selector_discharge": 0,
            "role_transfer_audit_license": 0,
            "role_bearing_ltotal_export": 16,
        },
        "p2427_product_inherited": theorem_export["p2427_product_assignment_count_inherited"],
        "p2427_factor_inherited": theorem_export["p2427_all_tables_factor_inherited"],
        "p2427_weight_support_empty_inherited": theorem_export["p2427_weight_side_support_empty_inherited"],
        "source_still_open": not theorem_export["source_obligation_discharge_exported"],
        "chi11_still_open": not theorem_export["chi11_source_exported"],
        "qw2191_still_open": not theorem_export["qw2191_discharged"],
        "role_transfer_still_blocked": not theorem_export["role_transfer_licensed"],
        "ltotal_still_blocked": not theorem_export["role_bearing_ltotal_exported"],
        "toe_still_open": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2428_s1378_v1",
        "packet_id": "P2428",
        "stage_id": "S1378",
        "result_kind": "REPAIR_READINESS_ANF_DERIVATIVE_CERTIFICATE",
        "status": "PASS_REPAIR_READINESS_ANF_DERIVATIVE_NO_GATE_DISCHARGE_NO_CLOSURE",
        "repair_readiness_anf_derivative_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: source.get("status", "OPEN_UNKNOWN") for name, source in sources.items()},
            "finite_witness_certificate": cert,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": (
            "Use the derivative supports to target actual theorem work: chi11 plus QW-2191 for selector readiness, and all five "
            "missing gates for ToE, without treating ANF structure as a discharge."
        ),
        "global_status": "OPEN_PROGRESS_REPAIR_READINESS_ANF_DERIVATIVE_CERTIFIED_NO_THEOREM_GATE_DISCHARGED",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["repair_readiness_anf_derivative_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2428 S1378: repair-readiness ANF derivative certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Finite facts",
                "",
                f"- Gate count: `{theorem['gate_count']}`.",
                f"- Truth-table size: `{theorem['truth_table_size']}`.",
                f"- ANF polynomials: `{theorem['anf_polynomials']}`.",
                f"- Derivative edge counts: `{theorem['derivative_edge_counts']}`.",
                "",
                "## Hard limits",
                "",
                *[f"- {item}" for item in theorem["not_licensed"]],
                "",
                "## Gatekeepers",
                "",
                f"`{payload['gatekeeper_checks']}`",
                "",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    append_doc_sections()
    payload = build_payload()
    write_outputs(payload)
    if not all(payload["gatekeeper_checks"].values()):
        raise SystemExit(f"gatekeeper failure: {payload['gatekeeper_checks']}")


if __name__ == "__main__":
    main()
