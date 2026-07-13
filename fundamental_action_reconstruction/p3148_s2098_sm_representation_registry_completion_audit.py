"""P3148/S2098: SM representation-registry completion audit.

P3147 recommended one concrete Lie-algebra-to-SM step: complete/audit an
explicit one-family SU(3)xSU(2)xU(1) representation registry plus Higgs, then
machine-check commutators, product-factor commutation, Yukawa hypercharge
invariance, and anomaly sums.  This constructs that missing object as an
axiom-branch Standard-Model representation packet.  It is a real algebraic
completion, but it remains conditional: the registry is installed as an SM
representation ansatz, not derived as a strict nadsoliton source theorem.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p3148_s2098_sm_representation_registry_completion_audit.json"
MD = GEN / "p3148_s2098_sm_representation_registry_completion_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P3147": GEN / "p3147_s2097_axiom_lie_smgr_fit_readiness_matrix.json",
    "P1962": GEN / "p1962_s912_strict_matter_higgs_brst_extension_registry_audit.json",
    "P1961": GEN / "p1961_s911_strict_local_nonabelian_brst_differential_and_nilpotency_certificate.json",
    "P1960": GEN / "p1960_s910_strict_qw2190_su3_su2_structure_constants_and_jacobi_certificate.json",
}

I = sp.I
SQRT3 = sp.sqrt(3)


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def append_once(path: Path, marker: str, text: str) -> None:
    existing = path.read_text(encoding="utf-8") if path.exists() else ""
    if marker not in existing:
        path.write_text(existing.rstrip() + "\n\n" + text.strip() + "\n", encoding="utf-8")


def zero(n: int) -> sp.Matrix:
    return sp.zeros(n, n)


def gell_mann() -> list[sp.Matrix]:
    return [
        sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]) / 2,
        sp.Matrix([[0, -I, 0], [I, 0, 0], [0, 0, 0]]) / 2,
        sp.Matrix([[1, 0, 0], [0, -1, 0], [0, 0, 0]]) / 2,
        sp.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]) / 2,
        sp.Matrix([[0, 0, -I], [0, 0, 0], [I, 0, 0]]) / 2,
        sp.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]) / 2,
        sp.Matrix([[0, 0, 0], [0, 0, -I], [0, I, 0]]) / 2,
        sp.Matrix([[1, 0, 0], [0, 1, 0], [0, 0, -2]]) / (2 * SQRT3),
    ]


def pauli() -> list[sp.Matrix]:
    return [
        sp.Matrix([[0, 1], [1, 0]]) / 2,
        sp.Matrix([[0, -I], [I, 0]]) / 2,
        sp.Matrix([[1, 0], [0, -1]]) / 2,
    ]


def comm(a: sp.Matrix, b: sp.Matrix) -> sp.Matrix:
    return sp.simplify(a * b - b * a)


def all_zero(mat: sp.Matrix) -> bool:
    return all(sp.simplify(x) == 0 for x in mat)


def su3_structure_constants() -> list[list[list[sp.Expr]]]:
    basis = gell_mann()
    f = [[[sp.Integer(0) for _ in range(8)] for _ in range(8)] for _ in range(8)]
    # With 2 Tr(T_a T_b)=delta_ab, f_abc = -2i Tr([T_a,T_b] T_c).
    for a in range(8):
        for b in range(8):
            cab = comm(basis[a], basis[b])
            for c in range(8):
                f[a][b][c] = sp.simplify(-2 * I * sp.trace(cab * basis[c]))
    return f


def epsilon3(a: int, b: int, c: int) -> int:
    if len({a, b, c}) < 3:
        return 0
    perm = [a, b, c]
    inversions = sum(1 for i in range(3) for j in range(i + 1, 3) if perm[i] > perm[j])
    return -1 if inversions % 2 else 1


def reps() -> list[dict[str, Any]]:
    # One SM generation in left-handed Weyl convention plus the Higgs.
    return [
        {"field": "Q_L", "kind": "fermion", "su3": "3", "su2": "2", "Y": Fraction(1, 6), "color_dim": 3, "weak_dim": 2},
        {"field": "u_c", "kind": "fermion", "su3": "bar3", "su2": "1", "Y": Fraction(-2, 3), "color_dim": 3, "weak_dim": 1},
        {"field": "d_c", "kind": "fermion", "su3": "bar3", "su2": "1", "Y": Fraction(1, 3), "color_dim": 3, "weak_dim": 1},
        {"field": "L_L", "kind": "fermion", "su3": "1", "su2": "2", "Y": Fraction(-1, 2), "color_dim": 1, "weak_dim": 2},
        {"field": "e_c", "kind": "fermion", "su3": "1", "su2": "1", "Y": Fraction(1, 1), "color_dim": 1, "weak_dim": 1},
        {"field": "H", "kind": "higgs", "su3": "1", "su2": "2", "Y": Fraction(1, 2), "color_dim": 1, "weak_dim": 2},
    ]


def color_generators(label: str) -> list[sp.Matrix]:
    if label == "3":
        return gell_mann()
    if label == "bar3":
        return [-m.conjugate() for m in gell_mann()]
    return [zero(1) for _ in range(8)]


def weak_generators(label: str) -> list[sp.Matrix]:
    if label == "2":
        return pauli()
    return [zero(1) for _ in range(3)]


def kron_generators(rep: dict[str, Any]) -> tuple[list[sp.Matrix], list[sp.Matrix], sp.Matrix]:
    cgens = color_generators(rep["su3"])
    wgens = weak_generators(rep["su2"])
    ident_c = sp.eye(rep["color_dim"])
    ident_w = sp.eye(rep["weak_dim"])
    su3 = [sp.kronecker_product(t, ident_w) for t in cgens]
    su2 = [sp.kronecker_product(ident_c, t) for t in wgens]
    y = sp.Rational(rep["Y"].numerator, rep["Y"].denominator) * sp.eye(rep["color_dim"] * rep["weak_dim"])
    return su3, su2, y


def check_rep_commutators() -> list[dict[str, Any]]:
    f3 = su3_structure_constants()
    rows = []
    for rep in reps():
        su3, su2, y = kron_generators(rep)
        su3_bad = 0
        for a in range(8):
            for b in range(8):
                rhs = sum((I * f3[a][b][c] * su3[c] for c in range(8)), start=zero(su3[0].rows))
                if not all_zero(comm(su3[a], su3[b]) - rhs):
                    su3_bad += 1
        su2_bad = 0
        for a in range(3):
            for b in range(3):
                rhs = sum((I * epsilon3(a, b, c) * su2[c] for c in range(3)), start=zero(su2[0].rows))
                if not all_zero(comm(su2[a], su2[b]) - rhs):
                    su2_bad += 1
        cross_bad = 0
        for a in range(8):
            for b in range(3):
                if not all_zero(comm(su3[a], su2[b])):
                    cross_bad += 1
        y_bad = 0
        for gen in su3 + su2:
            if not all_zero(comm(gen, y)):
                y_bad += 1
        rows.append({
            "field": rep["field"],
            "kind": rep["kind"],
            "su3": rep["su3"],
            "su2": rep["su2"],
            "hypercharge": str(rep["Y"]),
            "matrix_dimension": rep["color_dim"] * rep["weak_dim"],
            "su3_commutator_failures": su3_bad,
            "su2_commutator_failures": su2_bad,
            "cross_factor_commutator_failures": cross_bad,
            "hypercharge_commutator_failures": y_bad,
            "all_representation_checks_pass": su3_bad == su2_bad == cross_bad == y_bad == 0,
        })
    return rows


def anomaly_sums() -> dict[str, Any]:
    fermions = [r for r in reps() if r["kind"] == "fermion"]
    def t_su3(label: str) -> Fraction:
        return Fraction(1, 2) if label in {"3", "bar3"} else Fraction(0)
    def t_su2(label: str) -> Fraction:
        return Fraction(1, 2) if label == "2" else Fraction(0)
    sums = {
        "SU3_SU3_U1": sum(Fraction(r["weak_dim"]) * t_su3(r["su3"]) * r["Y"] for r in fermions),
        "SU2_SU2_U1": sum(Fraction(r["color_dim"]) * t_su2(r["su2"]) * r["Y"] for r in fermions),
        "U1_cubed": sum(Fraction(r["color_dim"] * r["weak_dim"]) * r["Y"] ** 3 for r in fermions),
        "gravity_gravity_U1": sum(Fraction(r["color_dim"] * r["weak_dim"]) * r["Y"] for r in fermions),
    }
    return {k: {"value": str(v), "vanishes": v == 0} for k, v in sums.items()}


def yukawa_checks() -> dict[str, Any]:
    # Left-handed Weyl convention: Q H u_c, Q H^dagger d_c, L H^dagger e_c.
    checks = {
        "up_yukawa_Q_H_uc": Fraction(1, 6) + Fraction(1, 2) + Fraction(-2, 3),
        "down_yukawa_Q_Hdagger_dc": Fraction(1, 6) + Fraction(-1, 2) + Fraction(1, 3),
        "lepton_yukawa_L_Hdagger_ec": Fraction(-1, 2) + Fraction(-1, 2) + Fraction(1, 1),
    }
    return {k: {"hypercharge_sum": str(v), "u1_invariant": v == 0} for k, v in checks.items()}


def build_payload() -> dict[str, Any]:
    rep_rows = check_rep_commutators()
    anomalies = anomaly_sums()
    yukawas = yukawa_checks()
    counts = {
        "registry_rows": len(reps()),
        "fermion_rows": sum(1 for r in reps() if r["kind"] == "fermion"),
        "higgs_rows": sum(1 for r in reps() if r["kind"] == "higgs"),
        "representation_rows_passing": sum(r["all_representation_checks_pass"] for r in rep_rows),
        "total_representation_failure_slots": sum(r["su3_commutator_failures"] + r["su2_commutator_failures"] + r["cross_factor_commutator_failures"] + r["hypercharge_commutator_failures"] for r in rep_rows),
        "anomaly_rows": len(anomalies),
        "anomaly_rows_vanishing": sum(v["vanishes"] for v in anomalies.values()),
        "yukawa_rows": len(yukawas),
        "yukawa_rows_u1_invariant": sum(v["u1_invariant"] for v in yukawas.values()),
        "strict_nadsoliton_source_rows": 0,
        "unit_bearing_ltotal_rows": 0,
    }
    accepted = counts["representation_rows_passing"] == counts["registry_rows"] and counts["anomaly_rows_vanishing"] == counts["anomaly_rows"] and counts["yukawa_rows_u1_invariant"] == counts["yukawa_rows"]
    return {
        "status": "P3148_SM_REPRESENTATION_REGISTRY_COMPLETION_ALGEBRAIC_PASS_CONDITIONAL_NO_STRICT_SOURCE",
        "input_hashes": {name: sha(path) for name, path in INPUTS.items()},
        "constructed_object": {
            "name": "R_SM^ax one-family representation registry plus Higgs",
            "classification": "axiom_branch_sm_representation_registry_completion",
            "convention": "left-handed Weyl fields Q_L,u_c,d_c,L_L,e_c plus H; hypercharges normalized as Q=T3+Y",
            "scope": "algebraic SM-representation completion, not strict nadsoliton derivation",
        },
        "registry_rows": [{**r, "Y": str(r["Y"])} for r in reps()],
        "representation_commutator_rows": rep_rows,
        "anomaly_sums": anomalies,
        "yukawa_hypercharge_checks": yukawas,
        "finite_theorem": {
            "name": "P3148_T1_sm_registry_algebraic_completion_certificate",
            "statement": "The constructed one-family R_SM^ax registry passes the finite algebraic SM checks requested after P3147: all 6 representation rows satisfy SU(3), SU(2), cross-factor, and U(1) commutator checks; all 4 audited one-family gauge/gravity anomaly sums vanish; and all 3 audited Yukawa hypercharge sums vanish.  This is an algebraic/axiom-branch completion of the representation registry gap, not a strict source theorem for why this registry is selected by the nadsoliton or how it is installed in a unit-bearing L_total.",
            "finite_counts": counts,
        },
        "decision": {
            "bounded_result": "The Lie-algebra-to-SM gap is partially reduced: a theorem-grade algebraic registry can be constructed and machine-checked.  This strengthens the claim that the local Lie algebra is good and that SM can fit under explicit axioms.",
            "why_not_strict": "The representation list, hypercharges, chirality convention, Higgs doublet, and Yukawa targets are installed as an SM ansatz.  No strict nadsoliton source selects this registry, no unit-bearing L_total/EOM coupling is exported, no global BV/BRST physical-state theorem is proved, and GR remains untouched.",
            "next_honest_step": "Construct P3149 as a bounded full-bundle BRST/L_total interface audit: feed R_SM^ax into the P1961 local BRST rules and test whether matter/Higgs kinetic and Yukawa terms are invariant in the same convention, while explicitly keeping unit-source, GR/EH, selector, and ToE claims closed unless a new source theorem is added.",
            "negative_export_flags": {
                "strict_registry_source_exported": False,
                "unit_bearing_L_total_exported": False,
                "global_BV_BRST_exported": False,
                "strict_SM_generation_exported": False,
                "strict_GR_generation_exported": False,
                "selector_closure_exported": False,
                "ToE_closure_exported": False,
            },
            "accepted_as_algebraic_registry_completion": accepted,
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    th = payload["finite_theorem"]
    c = th["finite_counts"]
    lines = [
        "# P3148/S2098 SM representation-registry completion audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Constructed object",
        f"- `{payload['constructed_object']['name']}`",
        f"- Classification: `{payload['constructed_object']['classification']}`",
        f"- Convention: `{payload['constructed_object']['convention']}`",
        "",
        "## Finite theorem",
        f"`{th['name']}`: {th['statement']}",
        "",
        "## Finite counts",
    ]
    for key, value in c.items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Registry rows"])
    for r in payload["registry_rows"]:
        lines.append(f"- `{r['field']}`: `{r['su3']},{r['su2']},Y={r['Y']}` ({r['kind']})")
    lines.extend(["", "## Anomaly sums"])
    for key, row in payload["anomaly_sums"].items():
        lines.append(f"- `{key}`: `{row['value']}`; vanishes `{row['vanishes']}`")
    lines.extend(["", "## Yukawa checks"])
    for key, row in payload["yukawa_hypercharge_checks"].items():
        lines.append(f"- `{key}`: `{row['hypercharge_sum']}`; U1 invariant `{row['u1_invariant']}`")
    lines.extend(["", "## Decision", payload["decision"]["bounded_result"], "", "## Why this is not strict", payload["decision"]["why_not_strict"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3148/S2098 SM representation-registry completion audit", "## P3148/S2098 SM representation-registry completion audit\n\n`P3148/S2098` constructs `R_SM^ax`, a one-family `SU(3)xSU(2)xU(1)` representation registry plus Higgs in left-handed Weyl convention.  The finite audit passes all algebraic checks: `6/6` representation rows satisfy factor commutators, `4/4` anomaly sums vanish, and `3/3` Yukawa hypercharge sums vanish.  This is an axiom-branch SM registry completion only; no strict nadsoliton source, unit-bearing `L_total`, global BV/BRST theorem, GR closure, selector closure, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3148/S2098 R_SM^ax is algebraic not unit-bearing", "## P3148/S2098 R_SM^ax is algebraic not unit-bearing\n\n`P3148/S2098` gives a concrete algebraic representation registry suitable for subsequent BRST/Lagrangian interface tests.  It does not derive the registry from the strict kernel, does not install it in a unit-bearing action measure, and does not prove full `L_total` gauge invariance or nonproxy SM/GR EOM.\n")
    append_once(AGENTS, "Current SM representation-registry completion guardrail (P3148/S2098, 2026-07-13)", "## Current SM representation-registry completion guardrail (P3148/S2098, 2026-07-13)\n\n- P3148 constructs `R_SM^ax`, an explicit one-family `SU(3)xSU(2)xU(1)` representation registry plus Higgs in left-handed Weyl convention.\n- The finite algebraic audit passes: `6/6` representation rows satisfy factor commutator checks, `4/4` gauge/gravity anomaly sums vanish, and `3/3` Yukawa hypercharge sums vanish.\n- This strengthens the axiom-branch SM fit and confirms the local Lie-algebra machinery is good as algebra, but the registry is an installed SM ansatz, not a strict nadsoliton source theorem.\n- Do not promote `R_SM^ax` to strict SM generation, unit-bearing `L_total`, global BV/BRST closure, GR closure, selector closure, bridge/role transfer, or ToE.\n- Next honest move: P3149 should feed `R_SM^ax` into the P1961 local BRST/Lagrangian interface and test matter/Higgs kinetic plus Yukawa invariance in one convention, still without unit-source, GR/EH, selector, or ToE promotion.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
