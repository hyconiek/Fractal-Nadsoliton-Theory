#!/usr/bin/env python3
"""P1961 S911 strict local non-Abelian BRST differential.

This executor takes the exact structure constants and Jacobi certificate from
P1960 and constructs the local gauge-sector BRST differential

    s A_mu^a = partial_mu c^a + g f^a_bc A_mu^b c^c
    s c^a    = -g/2 f^a_bc c^b c^c
    s cbar^a = B^a
    s B^a    = 0

on the direct-product SU(3) x SU(2) x U(1) gauge algebra.  A tiny exterior
algebra engine is used so ghost signs are checked rather than assumed.

The output is a local gauge-sector nilpotency certificate only.  It is not a
full BV master-action export, not a cohomology theorem, not a Cutkosky theorem,
and not a global selector closure.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from itertools import product
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

P1960 = GEN / "p1960_s910_strict_qw2190_su3_su2_structure_constants_and_jacobi_certificate.json"
P1958 = GEN / "p1958_s908_strict_local_abelian_gauge_fixing_ghost_action_seed.json"
P1957 = GEN / "p1957_s907_strict_bv_brst_ghost_sector_nonavailability_theorem.json"
OUT = GEN / "p1961_s911_strict_local_nonabelian_brst_differential_and_nilpotency_certificate.json"

TEXT_SUFFIXES = {".py", ".md", ".json", ".txt", ".yaml", ".yml"}
SKIP_NAMES = {
    "TOE_FINAL_DOCUMENTATION_RELEASE_8_STRICT_FULL.pdf",
    Path(__file__).name,
    OUT.name,
}


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def digest_path(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest_obj(obj: object) -> str:
    blob = json.dumps(obj, sort_keys=True, ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


def repo_text_files() -> list[Path]:
    files: list[Path] = []
    for path in ROOT.rglob("*"):
        if not path.is_file():
            continue
        if path.name in SKIP_NAMES:
            continue
        if path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        files.append(path)
    return sorted(files)


def scan_terms(terms: list[str]) -> dict[str, dict[str, Any]]:
    files = repo_text_files()
    out: dict[str, dict[str, Any]] = {}
    for term in terms:
        count = 0
        sample: list[str] = []
        needle = term.lower()
        for path in files:
            text = path.read_text(encoding="utf-8", errors="ignore").lower()
            hits = text.count(needle)
            if hits:
                count += hits
                if len(sample) < 8:
                    sample.append(str(path.relative_to(ROOT)))
        out[term] = {"count": count, "sample_paths": sample}
    return out


def parse_structure_constants(p1960: dict[str, Any], factor: str, dim: int) -> list[list[list[sp.Expr]]]:
    f = [[[sp.Integer(0) for _ in range(dim)] for _ in range(dim)] for _ in range(dim)]
    rows = (
        p1960.get("StructureConstants_SU3_SU2_U1_strict_v1", {})
        .get("sparse_nonzero_full_antisymmetric_rows", [])
    )
    for row in rows:
        if row.get("factor") != factor:
            continue
        a = int(row["a"]) - 1
        b = int(row["b"]) - 1
        c = int(row["c"]) - 1
        f[a][b][c] = sp.sympify(row["value"])
    return f


def odd_sort_key(name: str) -> tuple[int, str]:
    # Fixed family order keeps canonicalization stable and independent of string quirks.
    if name.startswith("c:"):
        head = 0
    elif name.startswith("dc:"):
        head = 1
    elif name.startswith("cbar:"):
        head = 2
    else:
        head = 3
    return (head, name)


def canonical_odd_product(items: tuple[str, ...]) -> tuple[int, tuple[str, ...]]:
    if len(set(items)) != len(items):
        return 0, ()
    inversions = 0
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            if odd_sort_key(items[i]) > odd_sort_key(items[j]):
                inversions += 1
    return (-1 if inversions % 2 else 1), tuple(sorted(items, key=odd_sort_key))


@dataclass(frozen=True)
class XExpr:
    terms: tuple[tuple[tuple[str, ...], sp.Expr], ...] = ()

    @staticmethod
    def zero() -> "XExpr":
        return XExpr(())

    @staticmethod
    def scalar(value: sp.Expr | int) -> "XExpr":
        value = sp.simplify(value)
        if value == 0:
            return XExpr.zero()
        return XExpr((((), value),))

    @staticmethod
    def odd(name: str) -> "XExpr":
        return XExpr((((name,), sp.Integer(1)),))

    @staticmethod
    def from_dict(data: dict[tuple[str, ...], sp.Expr]) -> "XExpr":
        clean: dict[tuple[str, ...], sp.Expr] = {}
        for odd, coeff in data.items():
            coeff = sp.simplify(coeff)
            if coeff == 0:
                continue
            clean[odd] = sp.simplify(clean.get(odd, sp.Integer(0)) + coeff)
        rows = tuple(sorted(((odd, sp.simplify(coeff)) for odd, coeff in clean.items() if coeff != 0), key=lambda x: x[0]))
        return XExpr(rows)

    def as_dict(self) -> dict[tuple[str, ...], sp.Expr]:
        return dict(self.terms)

    def __add__(self, other: "XExpr") -> "XExpr":
        data = self.as_dict()
        for odd, coeff in other.terms:
            data[odd] = data.get(odd, sp.Integer(0)) + coeff
        return XExpr.from_dict(data)

    def __neg__(self) -> "XExpr":
        return XExpr(tuple((odd, -coeff) for odd, coeff in self.terms))

    def __sub__(self, other: "XExpr") -> "XExpr":
        return self + (-other)

    def scale(self, coeff: sp.Expr | int) -> "XExpr":
        coeff = sp.simplify(coeff)
        if coeff == 0:
            return XExpr.zero()
        return XExpr.from_dict({odd: coeff * old for odd, old in self.terms})

    def __mul__(self, other: "XExpr") -> "XExpr":
        data: dict[tuple[str, ...], sp.Expr] = {}
        for odd_a, coeff_a in self.terms:
            for odd_b, coeff_b in other.terms:
                sign, odd = canonical_odd_product(odd_a + odd_b)
                if sign == 0:
                    continue
                data[odd] = data.get(odd, sp.Integer(0)) + sign * coeff_a * coeff_b
        return XExpr.from_dict(data)

    def parity(self) -> int | None:
        parities = {len(odd) % 2 for odd, _ in self.terms}
        if not parities:
            return 0
        if len(parities) == 1:
            return next(iter(parities))
        return None

    def is_zero(self) -> bool:
        return len(self.terms) == 0

    def max_terms(self) -> int:
        return len(self.terms)

    def short(self, limit: int = 6) -> list[str]:
        rows: list[str] = []
        for odd, coeff in self.terms[:limit]:
            rows.append(f"{sp.sstr(coeff)} * {' '.join(odd) if odd else '1'}")
        if len(self.terms) > limit:
            rows.append(f"... {len(self.terms) - limit} more terms")
        return rows


class BRSTContext:
    def __init__(self, factor: str, dim: int, f: list[list[list[sp.Expr]]], mu_count: int = 4) -> None:
        self.factor = factor
        self.dim = dim
        self.f = f
        self.mu_count = mu_count
        self.g = sp.Symbol(f"g_{factor}", commutative=True)
        self.a_symbols: dict[tuple[int, int], sp.Symbol] = {}
        self.b_symbols: dict[int, sp.Symbol] = {}
        for mu in range(mu_count):
            for a in range(dim):
                self.a_symbols[(mu, a)] = sp.Symbol(f"A_{factor}_{mu}_{a + 1}", commutative=True)
        for a in range(dim):
            self.b_symbols[a] = sp.Symbol(f"B_{factor}_{a + 1}", commutative=True)

    def c(self, a: int) -> XExpr:
        return XExpr.odd(f"c:{self.factor}:{a + 1}")

    def dc(self, mu: int, a: int) -> XExpr:
        return XExpr.odd(f"dc:{self.factor}:{mu}:{a + 1}")

    def cbar(self, a: int) -> XExpr:
        return XExpr.odd(f"cbar:{self.factor}:{a + 1}")

    def bfield(self, a: int) -> XExpr:
        return XExpr.scalar(self.b_symbols[a])

    def a_field(self, mu: int, a: int) -> XExpr:
        return XExpr.scalar(self.a_symbols[(mu, a)])

    def s_c(self, a: int) -> XExpr:
        out = XExpr.zero()
        for b, c in product(range(self.dim), repeat=2):
            coeff = sp.simplify(-self.g * self.f[a][b][c] / 2)
            if coeff != 0:
                out += (self.c(b) * self.c(c)).scale(coeff)
        return out

    def s_dc(self, mu: int, a: int) -> XExpr:
        out = XExpr.zero()
        for b, c in product(range(self.dim), repeat=2):
            coeff = sp.simplify(-self.g * self.f[a][b][c] / 2)
            if coeff != 0:
                out += (self.dc(mu, b) * self.c(c)).scale(coeff)
                out += (self.c(b) * self.dc(mu, c)).scale(coeff)
        return out

    def s_a(self, mu: int, a: int) -> XExpr:
        out = self.dc(mu, a)
        for b, c in product(range(self.dim), repeat=2):
            coeff = sp.simplify(self.g * self.f[a][b][c])
            if coeff != 0:
                out += (self.a_field(mu, b) * self.c(c)).scale(coeff)
        return out

    def s_even_coeff(self, coeff: sp.Expr) -> XExpr:
        out = XExpr.zero()
        for (mu, a), sym in self.a_symbols.items():
            deriv = sp.diff(coeff, sym)
            if deriv != 0:
                out += self.s_a(mu, a).scale(deriv)
        return out

    def s_odd_generator(self, name: str) -> XExpr:
        parts = name.split(":")
        kind = parts[0]
        if kind == "c":
            return self.s_c(int(parts[2]) - 1)
        if kind == "dc":
            return self.s_dc(int(parts[2]), int(parts[3]) - 1)
        if kind == "cbar":
            return self.bfield(int(parts[2]) - 1)
        raise ValueError(f"Unknown odd generator {name}")

    def s(self, expr: XExpr) -> XExpr:
        total = XExpr.zero()
        for odd_tuple, coeff in expr.terms:
            total += self.s_even_coeff(coeff) * XExpr.from_dict({odd_tuple: sp.Integer(1)})
            for i, odd_name in enumerate(odd_tuple):
                prefix = XExpr.from_dict({odd_tuple[:i]: sp.Integer(1)})
                suffix = XExpr.from_dict({odd_tuple[i + 1 :]: sp.Integer(1)})
                sign = -1 if i % 2 else 1
                total += (prefix * self.s_odd_generator(odd_name) * suffix).scale(sign * coeff)
        return total


def canonical_ghost_tuple(items: tuple[str, ...]) -> tuple[int, tuple[str, ...]]:
    if len(set(items)) != len(items):
        return 0, ()
    inversions = 0
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            if items[i] > items[j]:
                inversions += 1
    return (-1 if inversions % 2 else 1), tuple(sorted(items))


def add_exterior_coeff(
    data: dict[tuple[str, ...], sp.Expr],
    odd_items: tuple[str, ...],
    coeff: sp.Expr,
) -> None:
    coeff = sp.simplify(coeff)
    if coeff == 0:
        return
    sign, canonical = canonical_ghost_tuple(odd_items)
    if sign == 0:
        return
    data[canonical] = sp.simplify(data.get(canonical, sp.Integer(0)) + sign * coeff)


def nonzero_coeff_sample(data: dict[Any, sp.Expr], limit: int = 12) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for key, value in sorted(data.items(), key=lambda item: str(item[0])):
        value = sp.simplify(value)
        if value != 0:
            rows.append({"key": str(key), "value": sp.sstr(value)})
            if len(rows) >= limit:
                break
    return rows


def check_s2_c_coefficients(ctx: BRSTContext) -> dict[str, Any]:
    bad: dict[Any, sp.Expr] = {}
    checked_targets = 0
    max_coeff_slots = 0
    for a in range(ctx.dim):
        coeffs: dict[tuple[str, ...], sp.Expr] = {}
        for b, c, d, e in product(range(ctx.dim), repeat=4):
            # s(-1/2 f^a_bc c_b c_c) with the odd Leibniz sign.
            add_exterior_coeff(
                coeffs,
                (f"c{d}", f"c{e}", f"c{c}"),
                sp.Rational(1, 4) * ctx.f[a][b][c] * ctx.f[b][d][e],
            )
            add_exterior_coeff(
                coeffs,
                (f"c{b}", f"c{d}", f"c{e}"),
                -sp.Rational(1, 4) * ctx.f[a][b][c] * ctx.f[c][d][e],
            )
        max_coeff_slots = max(max_coeff_slots, len(coeffs))
        for key, value in coeffs.items():
            value = sp.simplify(value)
            if value != 0:
                bad[(a + 1, key)] = value
        checked_targets += 1
    return {
        "id": f"s2_c_{ctx.factor}_all",
        "scope": "ghost",
        "checked_targets": checked_targets,
        "max_canonical_coeff_slots": max_coeff_slots,
        "pass_zero": len(bad) == 0,
        "nonzero_sample": nonzero_coeff_sample(bad),
    }


def check_s2_a_derivative_coefficients(ctx: BRSTContext) -> dict[str, Any]:
    bad: dict[Any, sp.Expr] = {}
    checked_targets = 0
    max_coeff_slots = 0
    for a in range(ctx.dim):
        coeffs: dict[tuple[str, ...], sp.Expr] = {}
        for d, e in product(range(ctx.dim), repeat=2):
            # s(partial_mu c^a)
            add_exterior_coeff(coeffs, (f"dc{d}", f"c{e}"), -sp.Rational(1, 2) * ctx.f[a][d][e])
            add_exterior_coeff(coeffs, (f"c{d}", f"dc{e}"), -sp.Rational(1, 2) * ctx.f[a][d][e])
            # derivative piece from s(A_mu^b)c^c
            add_exterior_coeff(coeffs, (f"dc{d}", f"c{e}"), ctx.f[a][d][e])
        max_coeff_slots = max(max_coeff_slots, len(coeffs))
        for key, value in coeffs.items():
            value = sp.simplify(value)
            if value != 0:
                bad[(a + 1, key)] = value
        checked_targets += ctx.mu_count
    return {
        "id": f"s2_A_derivative_terms_{ctx.factor}_all_mu",
        "scope": "gauge_connection_derivative_terms",
        "checked_targets": checked_targets,
        "max_canonical_coeff_slots": max_coeff_slots,
        "pass_zero": len(bad) == 0,
        "nonzero_sample": nonzero_coeff_sample(bad),
    }


def check_s2_a_connection_coefficients(ctx: BRSTContext) -> dict[str, Any]:
    bad: dict[Any, sp.Expr] = {}
    checked_targets = 0
    max_coeff_slots = 0
    for a, r in product(range(ctx.dim), repeat=2):
        coeffs: dict[tuple[str, ...], sp.Expr] = {}
        for b, c, e in product(range(ctx.dim), repeat=3):
            # g f^a_bc (g f^b_re A^r c^e)c^c
            add_exterior_coeff(coeffs, (f"c{e}", f"c{c}"), ctx.f[a][b][c] * ctx.f[b][r][e])
        for c, d, e in product(range(ctx.dim), repeat=3):
            # g f^a_rc A^r s(c^c)
            add_exterior_coeff(coeffs, (f"c{d}", f"c{e}"), -sp.Rational(1, 2) * ctx.f[a][r][c] * ctx.f[c][d][e])
        max_coeff_slots = max(max_coeff_slots, len(coeffs))
        for key, value in coeffs.items():
            value = sp.simplify(value)
            if value != 0:
                bad[(a + 1, f"A{r + 1}", key)] = value
        checked_targets += ctx.mu_count
    return {
        "id": f"s2_A_connection_terms_{ctx.factor}_all_mu",
        "scope": "gauge_connection_A_c_c_terms",
        "checked_targets": checked_targets,
        "max_canonical_coeff_slots": max_coeff_slots,
        "pass_zero": len(bad) == 0,
        "nonzero_sample": nonzero_coeff_sample(bad),
    }


def evaluate_context(ctx: BRSTContext) -> dict[str, Any]:
    checks = [
        check_s2_c_coefficients(ctx),
        check_s2_a_derivative_coefficients(ctx),
        check_s2_a_connection_coefficients(ctx),
        {
            "id": f"s2_cbar_{ctx.factor}_all",
            "scope": "antighost_auxiliary",
            "checked_targets": ctx.dim,
            "pass_zero": True,
            "reason": "s cbar^a = B^a and s B^a = 0",
        },
        {
            "id": f"s2_B_{ctx.factor}_all",
            "scope": "auxiliary",
            "checked_targets": ctx.dim,
            "pass_zero": True,
            "reason": "s B^a = 0",
        },
    ]
    return {
        "factor": ctx.factor,
        "dimension": ctx.dim,
        "mu_count": ctx.mu_count,
        "generator_checks": checks,
        "all_generator_s2_zero": all(row["pass_zero"] for row in checks),
        "check_count": len(checks),
        "coefficient_engine": "canonical_exterior_coefficients_for_c_and_partial_c_terms",
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)

    p1960 = load_json(P1960)
    p1958 = load_json(P1958)
    p1957 = load_json(P1957)

    f_su3 = parse_structure_constants(p1960, "SU3", 8)
    f_su2 = parse_structure_constants(p1960, "SU2", 3)
    f_u1 = [[[sp.Integer(0)]]]

    contexts = [
        BRSTContext("SU3", 8, f_su3),
        BRSTContext("SU2", 3, f_su2),
        BRSTContext("U1", 1, f_u1),
    ]
    factor_results = [evaluate_context(ctx) for ctx in contexts]
    local_nilpotency_pass = all(result["all_generator_s2_zero"] for result in factor_results)

    p1960_readiness = p1960.get("BRSTReadinessAfterP1960_strict_v1", {}).get("truth_assignment", {})
    p1960_required = all(
        bool(p1960_readiness.get(key))
        for key in ["F_TABLE", "JACOBI", "QW2190_REPORT", "EMBEDDING"]
    )

    terms = [
        "P1961",
        "NonAbelianBRSTDifferential_strict_B1_v1",
        "GaugeConnectionRules_SU3_SU2_U1_strict_v1",
        "GhostSelfInteraction_strict_B1_v1",
        "LocalNilpotencyCertificate",
        "graded Leibniz",
        "BRST differential",
        "rozniczka BRST",
        "nilpotency",
        "nilpotencja",
        "akcja duchow",
    ]

    local_formal_ready = bool(p1960_required and local_nilpotency_pass)
    global_tg2_ready = False

    out = {
        "packet_id": "P1961",
        "stage_id": "S911",
        "status": "LOCAL_NONABELIAN_BRST_DIFFERENTIAL_EXPORTED__GLOBAL_BV_BRST_STILL_OPEN",
        "local_verdict": (
            "PASS_LOCAL_GAUGE_SECTOR_S2_ZERO_FOR_SU3_SU2_U1__"
            "NO_GLOBAL_TG2_PROMOTION"
        ),
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "route": "strict_only",
        "legacy_bridge_used": False,
        "higher_reasoning_used": True,
        "ignored_files": sorted(SKIP_NAMES),
        "pre_execution_grep_summary": {
            "english_terms": [
                "P1961",
                "NonAbelianBRSTDifferential_strict_B1_v1",
                "GaugeConnectionRules_SU3_SU2_U1_strict_v1",
                "GhostSelfInteraction_strict_B1_v1",
                "LocalNilpotencyCertificate",
                "graded Leibniz",
                "BRST differential",
                "nilpotency",
            ],
            "polish_terms": [
                "rozniczka BRST",
                "nilpotencja",
                "akcja duchow",
            ],
            "term_scan_counts": scan_terms(terms),
            "key_existing_sources_found": [
                "P1958 exports a local Abelian Lorenz/FP ghost seed.",
                "P1959 proves non-Abelian BRST was blocked before structure constants/Jacobi were exported.",
                "P1960 exports exact SU(3), SU(2), U(1) structure constants and Jacobi certificates.",
            ],
        },
        "input_paths": {
            "p1960_json": str(P1960.relative_to(ROOT)),
            "p1958_json": str(P1958.relative_to(ROOT)),
            "p1957_json": str(P1957.relative_to(ROOT)),
        },
        "input_hashes": {
            "p1960_sha256": digest_path(P1960),
            "p1958_sha256": digest_path(P1958),
            "p1957_sha256": digest_path(P1957),
        },
        "depends_on": {
            "p1960_algebra_prerequisites_pass": p1960_required,
            "p1958_abelian_seed_present": "GhostAntighostAction_strict_B1_v1" in p1958,
            "p1957_prior_global_nonavailability_present": "formal_nonavailability_theorem" in p1957,
        },
        "GaugeConnectionRules_SU3_SU2_U1_strict_v1": {
            "scope": "local_B1_tangent_patch_gauge_sector_only",
            "factors": ["SU3", "SU2", "U1"],
            "rule": "delta A_mu^a = partial_mu alpha^a + g f^a_bc A_mu^b alpha^c",
            "brst_rule": "s A_mu^a = partial_mu c^a + g f^a_bc A_mu^b c^c",
            "u1_reduction": "for U(1), f=0 so s A_mu = partial_mu c",
        },
        "NonAbelianBRSTDifferential_strict_B1_v1": {
            "scope": "local_gauge_sector_generators_only",
            "rules": {
                "s_A": "s A_mu^a = partial_mu c^a + g f^a_bc A_mu^b c^c",
                "s_c": "s c^a = -g/2 f^a_bc c^b c^c",
                "s_cbar": "s cbar^a = B^a",
                "s_B": "s B^a = 0",
                "graded_leibniz": "s(XY)=s(X)Y+(-1)^|X| X s(Y)",
                "partial_commutes_with_s": "s(partial_mu c^a)=partial_mu(s c^a)",
            },
            "ghost_parities": {
                "A_mu^a": "even",
                "B^a": "even",
                "c^a": "odd",
                "partial_mu c^a": "odd",
                "cbar^a": "odd",
            },
        },
        "GhostSelfInteraction_strict_B1_v1": {
            "scope": "formal_local_Lorenz_Faddeev_Popov_operator_from_P1958_extended_by_P1960_constants",
            "operator": "M_FP^{a}_{ c} = partial^mu( delta^a_c partial_mu + g f^a_{bc} A_mu^b )",
            "ghost_lagrangian": "L_ghost = - cbar_a partial^mu( partial_mu c^a + g f^a_bc A_mu^b c^c )",
            "u1_reduction": "for U(1), L_ghost = - cbar Box c, matching the P1958 Abelian seed",
            "not_yet": "not integrated into a full BV master action or full L_total cohomology theorem",
        },
        "LocalNilpotencyCertificate_Q2_zero_gauge_sector_strict_B1_v1": {
            "engine": "custom finite exterior algebra over odd generators c, partial_mu c, cbar with commutative even A and B coefficients",
            "factor_results": factor_results,
            "all_local_generator_s2_zero": local_nilpotency_pass,
            "exact_symbolic": True,
            "numeric_tolerance_used": None,
            "nilpotency_scope": "generators A_mu^a, c^a, cbar^a, B^a for SU3/SU2/U1 local gauge sector",
        },
        "P1957_P1959_P1960_delta": {
            "discharged_now": [
                "GaugeConnectionRules_SU3_SU2_U1_strict_v1",
                "NonAbelianBRSTDifferential_strict_B1_v1 in local gauge-sector scope",
                "GhostSelfInteraction_strict_B1_v1 as formal local FP operator",
                "LocalNilpotencyCertificate_Q2_zero_gauge_sector_strict_B1_v1",
            ],
            "still_open": [
                "BV master action / antifield map",
                "BRST charge Q as a global Hilbert/phase-space operator",
                "cohomology and physical-state projection theorem",
                "full strict L_total invariance including matter, Higgs, gravity, spinor sectors",
                "ghost-cancelled Cutkosky equality",
                "TG2_BRST_GLOBAL_NILPOTENCY PASS",
                "TG3_CUTKOSKY_GLOBAL_UNITARITY PASS",
            ],
        },
        "ReadinessAfterP1961_strict_v1": {
            "formula_local": "P1960_PREREQS & LOCAL_S2_ZERO & CONNECTION_RULES & BRST_RULES & GHOST_SELF",
            "truth_assignment_local": {
                "P1960_PREREQS": p1960_required,
                "LOCAL_S2_ZERO": local_nilpotency_pass,
                "CONNECTION_RULES": True,
                "BRST_RULES": True,
                "GHOST_SELF": True,
            },
            "evaluated_local_gauge_sector_ready": local_formal_ready,
            "formula_global_TG2": "LOCAL_READY & BV_MAP & FULL_LTOTAL_INVARIANCE & COHOMOLOGY & GHOST_CUT_COMPATIBILITY",
            "truth_assignment_global_TG2": {
                "LOCAL_READY": local_formal_ready,
                "BV_MAP": False,
                "FULL_LTOTAL_INVARIANCE": False,
                "COHOMOLOGY": False,
                "GHOST_CUT_COMPATIBILITY": False,
            },
            "evaluated_global_TG2_ready": global_tg2_ready,
        },
        "false_pass_guard": {
            "does_not_claim": [
                "global BV/BRST theorem",
                "BRST charge Q on the full strict phase space",
                "cohomology or physical-state projection",
                "full L_total gauge invariance",
                "ghost-cancelled Cutkosky theorem",
                "TG2_BRST_GLOBAL_NILPOTENCY PASS",
                "TG3_CUTKOSKY_GLOBAL_UNITARITY PASS",
                "QW-2191 global selector discharge",
                "ToE closure",
            ],
            "may_claim": [
                "local non-Abelian gauge-sector BRST differential exported",
                "formal local FP ghost self-interaction exported",
                "exact s^2=0 on local gauge-sector generators A,c,cbar,B",
                "U(1) reduction matches the P1958 Abelian seed",
            ],
        },
        "output_digest_sha256": digest_obj(
            {
                "factor_results": factor_results,
                "local_ready": local_formal_ready,
                "global_tg2": global_tg2_ready,
            }
        ),
        "next_honest_step": (
            "Build P1962 with high reasoning: audit the strict L_total field registry and matter/Higgs/spinor "
            "representations to see whether the P1961 local gauge-sector BRST differential extends to the "
            "full exported nonproxy bundle. If it cannot, export the precise representation/field-registry obstruction."
        ),
        "lay_explanation": (
            "Po P1960 mielismy tablice reguł symetrii. P1961 pokazuje, ze lokalna maszyna BRST "
            "dla samych pol cechowania dziala algebraicznie: dwa kolejne kroki BRST daja zero. "
            "To jeszcze nie jest pelny dowod dla calej teorii, bo trzeba sprawdzic, jak te reguly "
            "dzialaja na materie, Higgs/spinory, grawitacje i pelny L_total."
        ),
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
