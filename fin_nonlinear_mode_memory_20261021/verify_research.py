#!/usr/bin/env python3
"""Exact normalization/scope controls, not a mechanized nonlinear PDE proof."""
import argparse
import ast
from fractions import Fraction as F
import hashlib
from math import factorial
from pathlib import Path
import json
import re
import subprocess
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
INPUTS = [
    "AGENTS.md", "research_tasks/00_EXECUTION_RULES.md",
    "research_tasks/RP004_SHARED_CONTRACT.json",
    "fin_coherent_mode_radiation_review_20261019/RETARDED_FLUX_PROOF_ROUTE.md",
    "fin_coherent_mode_radiation_review_20261019/PROOF_COMPLETIONS_AND_SCOPE.md",
    "fin_coherent_mode_radiation_campaign_20261017/MODE_AND_FLUX_CONTRACT.json",
    "fin_coherent_mode_radiation_campaign_20261017/ON_SHELL_CERTIFICATE_OR_FAILURE.md",
    "fin_rp004_flux_review_20261019/BOOTSTRAP_COMPLETION.md",
    "fin_rp004_flux_review_20261019/MOVING_CONE_THEOREM.md",
    "fin_rp004_flux_review_20261019/OUTWARD_FLUX_THEOREM.md",
    "fin_rp004_full_input_self_audit_20261020/CONTRACT.md",
    "fin_rp004_full_input_self_audit_20261020/CAUSAL_AND_SPECTRAL_AUDIT.md",
    "fin_nonlinear_stability_transfer_review_20261017/PROOF_COMPLETIONS_AND_SCOPE.md",
]
PROTECTED = ["fin_coherent_mode_radiation_review_20261019", "fin_rp004_flux_review_20261019",
             "fin_rp004_full_input_self_audit_20261020", "fin_nonlinear_stability_transfer_review_20261017"]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def cn(x):
    if isinstance(x, tuple):
        assert len(x) == 2
        return F(x[0]), F(x[1])
    assert not isinstance(x, (float, complex))
    return F(x), F(0)


ZERO = cn(0)


def ca(a, b):
    a, b = cn(a), cn(b)
    return a[0]+b[0], a[1]+b[1]


def cm(a, b):
    a, b = cn(a), cn(b)
    return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]


def conj(a):
    a = cn(a)
    return a[0], -a[1]


def norm(a):
    a = cn(a)
    return a[0]*a[0]+a[1]*a[1]


def cd(a, b):
    a, b = cn(a), cn(b)
    denominator = norm(b)
    assert denominator > 0
    out = cm(a, conj(b))
    return out[0]/denominator, out[1]/denominator


def polynomial(values):
    out = [cn(x) for x in values]
    while len(out) > 1 and out[-1] == ZERO:
        out.pop()
    return out or [ZERO]


def pa(a, b):
    a, b = polynomial(a), polynomial(b)
    return polynomial([ca(a[i] if i < len(a) else ZERO, b[i] if i < len(b) else ZERO)
                       for i in range(max(len(a), len(b)))])


def ps(a, scale):
    return polynomial([cm(x, scale) for x in polynomial(a)])


def pm(a, b):
    a, b = polynomial(a), polynomial(b)
    out = [ZERO]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] = ca(out[i+j], cm(x, y))
    return polynomial(out)


def pd(a):
    a = polynomial(a)
    return polynomial([cm(i, a[i]) for i in range(1, len(a))])


def dm(a, m):
    return pa(pm([1, 0, -1], pd(a)), ps(pm([0, m], a), -1))


def dx(a, m, wave):
    return pa(dm(a, m), ps(a, (F(0), F(wave))))


def h_operator(a, m, wave=F(0)):
    return pa(ps(dx(dx(a, m, wave), m, wave), -1), pm([-2, 0, 6], a))


def integrate_even_sech(m, a):
    assert m >= 2 and m % 2 == 0
    weight = polynomial([1])
    for _ in range(m//2-1):
        weight = pm(weight, [1, 0, -1])
    integrand = pm(weight, a)
    out = ZERO
    for i, coefficient in enumerate(integrand):
        if i % 2 == 0:
            out = ca(out, cm(coefficient, F(2, i+1)))
    assert out[1] == 0
    return out[0]


def gamma_half_coefficient(n):
    out = F(1)
    for j in range(n):
        out *= F(2*j+1, 2)
    return out


def manifest(folder, coverage=False):
    seen = set()
    for line in (folder/"MANIFEST.sha256").read_text().splitlines():
        h, name = line.split("  ", 1)
        assert name not in seen and not Path(name).is_absolute() and ".." not in Path(name).parts
        assert sha(folder/name) == h, (folder, name)
        seen.add(name)
    if coverage:
        actual = {p.relative_to(folder).as_posix() for p in folder.rglob("*") if p.is_file() and "__pycache__" not in p.parts}
        assert actual == seen | {"MANIFEST.sha256"}
    return {"entries": len(seen), "sha256": sha(folder/"MANIFEST.sha256"),
            "coverage": "complete_noncache" if coverage else "listed_members_only"}


def links():
    for p in HERE.glob("*.md"):
        assert all(line == line.rstrip() for line in p.read_text().splitlines()), p
        text = re.sub(r"```.*?```", "", p.read_text(), flags=re.S)
        text = re.sub(r"`[^`\n]*`", "", text)
        for link in re.findall(r"\[[^\]\n]+\]\(([^)\n]+)\)", text):
            if not link.startswith(("https:", "http:")):
                assert (p.parent/link.split("#")[0]).is_file(), (p, link)


def main():
    ap = argparse.ArgumentParser()
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--record", action="store_true")
    mode.add_argument("--check", action="store_true")
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    if args.check:
        manifest(HERE, True)
    inputs = {n: sha(ROOT/n) for n in INPUTS}
    protected = {n: manifest(ROOT/n) for n in PROTECTED}
    checks = []
    def ck(name, condition, **data):
        assert bool(condition), name
        checks.append({"name": name, "passed": True, **data})

    source = json.loads((ROOT/INPUTS[2]).read_text())
    ck("actual predeclared source equations and mass interval bound",
       source["equations"] == ["phi_tt-phi_xx+2phi(phi^2-1)+6phi psi^2=0",
                               "psi_tt-psi_xx+(a+6phi^2)psi+g4 psi^3+g6 psi^5=0"]
       and list(map(F, source["a_interval"])) == [F(-3, 4), F(-1, 4)])
    ck("prepared bound-mode versus shape frequency separation", 4*(F(-3, 4)+2)-3 == 2
       and F(-3, 4)+2 == F(5, 4) and F(-1, 4)+2 == F(7, 4))

    ck("true translation profile eigen identity", h_operator([1], 2) == polynomial([0]))
    ck("true odd shape eigen identity", h_operator([0, 1], 1) == polynomial([0, 3]))
    for wave in (F(0), F(1), F(-1), F(3, 2)):
        jost = [(-1-wave*wave, F(0)), (F(0), -3*wave), cn(3)]
        ck(f"exact Jost differential residual k{wave}", h_operator(jost, 0, wave) == ps(jost, 4+wave*wave),
           scope="pointwise differential control; completeness is the inherited Darboux theorem")

    norm_f = F(3, 4)*integrate_even_sech(4, [1])
    norm_shape = F(3, 2)*integrate_even_sech(2, [0, 0, 1])
    c4 = F(9, 16)*integrate_even_sech(8, [1])
    c6 = F(27, 64)*integrate_even_sech(12, [1])
    norm_f_feedback = F(9, 16)*integrate_even_sech(8, [0, 0, 1])
    ck("bound and shape modes correctly normalized", norm_f == norm_shape == 1)
    ck("signed self terms keep exact fourth and sixth projections", c4 == F(18, 35) and c6 == F(24, 77))
    ck("feedback norm and translation orthogonality", norm_f_feedback == F(2, 35)
       and integrate_even_sech(6, [0, 1]) == 0)
    beta_pi = gamma_half_coefficient(1)*gamma_half_coefficient(2)/factorial(3)
    shape_mass_pi2 = F(9, 16)*F(3, 2)*beta_pi**2
    ck("odd shape overlap coefficient", beta_pi == F(1, 16) and shape_mass_pi2 == F(27, 8192))
    ck("continuum spectral mass is strictly positive", norm_f_feedback-shape_mass_pi2*F(22, 7)**2 > 0)
    inverse_profile = ps([0, 1], F(1, 8))
    expected_same_factor = ps(pm([1, 0, -1], [0, 1]), F(3, 4))
    ck("static inverse is exact on the odd sector", h_operator(inverse_profile, 2) == expected_same_factor)
    chi0 = F(3, 32)*integrate_even_sech(6, [0, 0, 1])
    ck("static susceptibility and positive continuum remainder", chi0 == F(1, 70)
       and chi0-shape_mass_pi2*F(22, 7)**2/3 > 0)

    hf = h_operator([0, 1], 4)
    hf_inner = F(9, 16)*integrate_even_sech(8, pm([0, 1], hf))
    gradient_f = dm([0, 1], 4)
    energy_integrand = pa(pm(gradient_f, gradient_f), pm([-2, 0, 6], [0, 0, 1]))
    hf_energy = F(9, 16)*integrate_even_sech(8, energy_integrand)
    ck("independent integration-by-parts feedback moment", hf_inner == hf_energy == F(24, 77))

    # General Fourier coefficient identity, compared as coefficients of k^2.
    coefficient_k2 = F(1, 4)-F(5, 2)/20
    coefficient_const = F(10, 4)-F(5, 2)*F(16, 20)
    ck("Fourier transform reduction polynomial identity", coefficient_k2 == F(1, 8)
       and coefficient_const == F(1, 2))
    for wave in (F(1), F(3, 2)):
        denominator = cm((F(1), wave), (F(2), wave))
        m_feedback = cd((F(0), 3*wave*(wave*wave+4)), cm(32, denominator))
        m_source = cd((F(0), 9*wave*(wave*wave+4)), cm(32, denominator))
        ck(f"source S=3F and Gamma factor nine k{wave}", m_source == cm(3, m_feedback)
           and norm(m_source) == 9*norm(m_feedback))
    threshold_slope = F(3, 32)*4*integrate_even_sech(4, [1])/2
    a2_pi = threshold_slope**2/4
    b2_pi = a2_pi*4*2
    bpp_pi = 2*b2_pi
    ck("threshold and stationary-phase constants retain both roots", threshold_slope == F(1, 4)
       and a2_pi == F(1, 64) and bpp_pi == F(1, 4) and bpp_pi/4 == F(1, 16))
    ck("exact quadratic phase change coefficients", pa(pm([2, 0, 1], [2, 0, 1]), [-4]) == pm([0, 0, 1], [4, 0, 1]),
       scope="k^2=4u^2+u^4 and nu=2+u^2; full asymptotic remainder is analytical")

    # Reciprocal forces and causal turn-on; no field integration is simulated.
    kk, eta, psi, a0, g4, g6 = F(2, 3), F(1, 10), F(2, 7), F(-1, 2), F(-2, 3), F(5, 7)
    def grad_phi(phi, u):
        return 2*phi*(phi*phi-1)+6*phi*u*u
    def grad_psi(phi, u):
        return (a0+6*phi*phi)*u+g4*u**3+g6*u**5
    ck("exact full phi nonlinear remainder", grad_phi(kk+eta, psi)-grad_phi(kk, 0)-(6*kk*kk-2)*eta
       == 6*kk*eta**2+2*eta**3+6*kk*psi**2+6*eta*psi**2)
    ck("exact full psi nonlinear remainder with signed self terms", grad_psi(kk+eta, psi)-(a0+6*kk*kk)*psi
       == 12*kk*eta*psi+6*eta*eta*psi+g4*psi**3+g6*psi**5)
    ck("retarded feedback sign and DC source normalization", (-12)*(-6) == 72
       and F(6, 2)*F(3, 4) == F(9, 4)
       and 72*norm_f_feedback/2 == 36*norm_f_feedback)
    epsilon = F(1, 16)
    initial_acceleration = -(a0+2)*epsilon-g4*c4*epsilon**3-g6*c6*epsilon**5
    ck("quintic coefficient remains in the initial residual", initial_acceleration+(a0+2)*epsilon+g4*c4*epsilon**3
       == -g6*c6*epsilon**5 != 0)

    # Algebra-only two-oscillator control, not a continuum discretization or a new FIN law.
    z, zt, w, wt, overlap, omega2, lam = F(1, 7), F(-2, 5), F(3, 11), F(4, 13), F(3, 5), F(3, 2), F(3)
    ztt = -omega2*z-g4*c4*z**3-12*overlap*z*w
    wtt = -lam*w-6*overlap*z*z
    energy_derivative = zt*(ztt+omega2*z+g4*c4*z**3)+wt*(wtt+lam*w)
    energy_derivative += 12*overlap*z*zt*w+6*overlap*z*z*wt
    ck("reciprocal mode-field energy control", energy_derivative == 0,
       scope="rational oscillator fixture only, not a certified continuum or Gaussian approximation")
    wrong_wtt = -lam*w-3*overlap*z*z
    wrong_energy = zt*(ztt+omega2*z+g4*c4*z**3)+wt*(wrong_wtt+lam*w)
    wrong_energy += 12*overlap*z*zt*w+6*overlap*z*z*wt
    ck("halved feedback-force repair rejected", wrong_energy != 0)

    for sigma, tau in ((F(1, 2), F(2)), (F(2, 3), F(-3, 2))):
        for lam in (F(3), F(5)):
            s = (sigma, tau)
            den = ca(cm(s, s), lam)
            response = cd(s, den)
            expected = sigma*(sigma*sigma+tau*tau+lam)/norm(den)
            ck(f"positive-real response sigma{sigma} tau{tau} lambda{lam}", response[0] == expected > 0)
    retarded = cd(1, (F(-1), F(-1, 7)))
    advanced = cd(1, (F(-1), F(1, 7)))
    ck("retarded and advanced susceptibility signs separated", retarded[1] > 0 > advanced[1])
    orders = [(5, 3), (5, 3), (5, 2), (5, 2), (5, 0)]
    ck("projected residual order ledger", all(pe >= 5 and pt <= 3 for pe, pt in orders),
       scope="Sobolev products and absorption inequalities are proved analytically, not by this ledger")
    ck("growing-window versus damping-scale distinction", 2-3*F(1, 2) > 0 and 2-3*2 < 0)

    named = [n for n in INPUTS if n.endswith(".md") and n != "AGENTS.md"]
    pattern = r"Volterra|memory kernel|self.force|Darboux|shape|Gamma|epsilon.{0,20}tau"
    search = subprocess.run(["rg", "-n", "-i", "--max-columns", "300", "--", pattern, *named],
                            cwd=ROOT, capture_output=True, text=True, timeout=10)
    assert search.returncode in (0, 1), search.stderr
    assert inputs == {n: sha(ROOT/n) for n in INPUTS}
    assert protected == {n: manifest(ROOT/n) for n in PROTECTED}
    for p in HERE.glob("*.py"):
        ast.parse(p.read_text())
    for p in HERE.glob("*.json"):
        json.loads(p.read_text())
    result = {"date": "2026-10-21", "status": "SOURCE_BOUND_NONLINEAR_VOLTERRA_REDUCTION_WITH_EXPLICIT_MEMORY",
              "bounded_checks": len(checks), "checks": checks, "input_hashes": inputs,
              "protected_manifests": protected,
              "exact_moments": {"f4": str(c4), "f6": str(c6), "F_norm2": str(norm_f_feedback),
                                "shape_mass_pi2_coefficient": str(shape_mass_pi2),
                                "static_susceptibility": str(chi0), "F_HF": str(hf_inner)},
              "feedback_coefficient": 72, "Gamma_relation": "Gamma_old=9 Im chi_F",
              "continuum_tail": "sin(2t+3pi/4)/(16sqrt(pi)t^(3/2))+O(t^(-5/2))",
              "prediction_bound": "C epsilon^5(1+t)^4 on epsilon^2(1+T)^3<=kappa; analytical inherited/new proof",
              "finite_control_scope": "exact bounded identities, not PDE/continuum/quadrature or independent proof",
              "duplicate_search": {"files": named, "pattern": pattern, "matching_lines": sorted(search.stdout.splitlines()),
                                   "scope": "bounded named-source search, not exhaustive novelty"},
              "source_model": source["source_family"], "primary_HB_changed": False,
              "new_noise_or_incoming_field": False, "epsilon_inverse_square_decay_theorem": False,
              "F01": "OPEN_SOURCE", "independent_review": False, "AGENTS_or_originals_modified": False}
    intext = "".join(f"{h}  {n}\n" for n, h in sorted(inputs.items()))
    if args.record:
        (HERE/"RESULTS.json").write_text(json.dumps(result, indent=2)+"\n")
        (HERE/"INPUTS.sha256").write_text(intext)
        files = sorted(p for p in HERE.rglob("*") if p.is_file() and p.name != "MANIFEST.sha256" and "__pycache__" not in p.parts)
        (HERE/"MANIFEST.sha256").write_text("".join(f"{sha(p)}  {p.relative_to(HERE).as_posix()}\n" for p in files))
    else:
        saved = json.loads((HERE/"RESULTS.json").read_text())
        assert result == saved, [key for key in result if result[key] != saved.get(key)]
        assert intext == (HERE/"INPUTS.sha256").read_text()
    links()
    if args.output:
        assert Path("/tmp") in args.output.resolve().parents
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    literature = json.loads((HERE/"LITERATURE_REGISTER.json").read_text())["sources"][0]
    cache = Path(literature["optional_local_cache"])
    if cache.exists():
        assert cache.stat().st_size == literature["pdf_bytes"] and sha(cache) == literature["pdf_sha256"]
    print(json.dumps({"status": result["status"], "bounded_checks": len(checks), "inputs": len(inputs),
                      "protected_manifest_entries": sum(x["entries"] for x in protected.values()),
                      "static_susceptibility": str(chi0), "feedback_coefficient": 72,
                      "new_noise_or_incoming_field": False, "F01": "OPEN_SOURCE",
                      "optional_pdf_cache": "PRESENT_AND_HASH_MATCHED" if cache.exists() else "ABSENT_OPTIONAL_CACHE",
                      "AGENTS_or_originals_modified": False}, indent=2))


if __name__ == "__main__":
    main()
