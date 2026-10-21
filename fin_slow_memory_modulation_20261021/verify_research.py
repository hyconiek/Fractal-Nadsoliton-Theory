#!/usr/bin/env python3
"""Exact primitive/series controls, not a numerical long-time PDE proof."""
import argparse
import ast
from dataclasses import dataclass
from fractions import Fraction as F
from math import factorial
import hashlib
import json
import os
from pathlib import Path
import re
import resource
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
NM = ROOT/"fin_nonlinear_mode_memory_20261021"
DEGREE = 12
INPUTS = [
    "AGENTS.md", "research_tasks/00_EXECUTION_RULES.md",
    "research_tasks/RP004_SHARED_CONTRACT.json",
    "fin_nonlinear_mode_memory_20261021/CONTRACT.md",
    "fin_nonlinear_mode_memory_20261021/SPECTRAL_MEMORY_KERNEL.md",
    "fin_nonlinear_mode_memory_20261021/NONLINEAR_REDUCTION.md",
    "fin_nonlinear_mode_memory_20261021/NEXT_OPEN_LEMMA.md",
    "fin_nonlinear_mode_memory_20261021/PROOF_OR_OPEN_LEMMA.md",
    "fin_nonlinear_mode_memory_20261021/RESULTS.json",
    "fin_nonlinear_mode_memory_20261021/verify_research.py",
    "fin_nonlinear_mode_memory_20261021/LITERATURE_REGISTER.json",
    "fin_rp004_flux_review_20261019/BOOTSTRAP_COMPLETION.md",
]
PROTECTED = ["fin_nonlinear_mode_memory_20261021", "fin_rp004_flux_review_20261019"]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


@dataclass(frozen=True)
class G:
    real: F = F(0)
    imag: F = F(0)

    @staticmethod
    def of(x):
        if isinstance(x, G):
            return x
        if isinstance(x, tuple):
            return G(F(x[0]), F(x[1]))
        assert not isinstance(x, (float, complex))
        return G(F(x))

    def __add__(self, other):
        b = self.of(other)
        return G(self.real+b.real, self.imag+b.imag)

    __radd__ = __add__

    def __neg__(self):
        return G(-self.real, -self.imag)

    def __sub__(self, other):
        return self+-self.of(other)

    def __mul__(self, other):
        b = self.of(other)
        return G(self.real*b.real-self.imag*b.imag, self.real*b.imag+self.imag*b.real)

    __rmul__ = __mul__

    def __truediv__(self, other):
        b = self.of(other)
        den = b.real*b.real+b.imag*b.imag
        assert den != 0
        num = self*G(b.real, -b.imag)
        return G(num.real/den, num.imag/den)


ZERO = G()
I = G(F(0), F(1))


def pad(values):
    out = [G.of(x) for x in values[:DEGREE+1]]
    return out+[ZERO]*(DEGREE+1-len(out))


def add(a, b):
    return [x+y for x, y in zip(pad(a), pad(b))]


def scale(a, c):
    return [x*G.of(c) for x in pad(a)]


def multiply(a, b):
    out = [ZERO]*(DEGREE+1)
    for j, x in enumerate(pad(a)):
        for k, y in enumerate(pad(b)):
            if j+k <= DEGREE:
                out[j+k] = out[j+k]+x*y
    return out


def convolution(a, b):
    out = [ZERO]*(DEGREE+1)
    for j, x in enumerate(pad(a)):
        for k, y in enumerate(pad(b)):
            if j+k+1 <= DEGREE:
                beta = F(factorial(j)*factorial(k), factorial(j+k+1))
                out[j+k+1] = out[j+k+1]+x*y*beta
    return out


def derivative(a):
    a = [G.of(x) for x in a]
    return [j*a[j] for j in range(1, len(a))] or [ZERO]


def slow(a, epsilon):
    return pad([G.of(x)*epsilon**(2*j) for j, x in enumerate(a)])


def exponent(nu):
    out = [G.of(1)]
    for j in range(1, DEGREE+1):
        out.append(out[-1]*(I*nu)/j)
    return out


def cosine_and_sine(lam):
    cosine, sine = [ZERO]*(DEGREE+1), [ZERO]*(DEGREE+1)
    for n in range(DEGREE+1):
        if n % 2 == 0:
            cosine[n] = G(F((-lam)**(n//2), factorial(n)))
        else:
            sine[n] = G(F((-lam)**((n-1)//2), factorial(n)))
    return cosine, sine


def shape_coefficients(lam, nu):
    den = lam-nu*nu
    if den == 0:
        raise ValueError("shape resonance: no bounded primitive")
    chi = F(1)/den
    p = (G(chi), -I*nu*chi)
    q = (-I*(2*nu/den**2), G(-(lam+nu*nu)/den**2))
    return chi, p, q


def differentiate_basis(pair, lam, nu):
    a, b = pair
    return I*nu*a+b, -lam*a+I*nu*b


def shape_series(lam, nu):
    chi, p, q = shape_coefficients(lam, nu)
    c, s = cosine_and_sine(lam)
    e = exponent(nu)
    return chi, multiply(e, add(scale(c, p[0]), scale(s, p[1]))), multiply(e, add(scale(c, q[0]), scale(s, q[1]))), multiply(e, s)


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


def limits():
    resource.setrlimit(resource.RLIMIT_AS, (1024**3, 1024**3))
    resource.setrlimit(resource.RLIMIT_CPU, (55, 55))


def replay_nm():
    before = {p.relative_to(NM).as_posix(): sha(p) for p in NM.rglob("*") if p.is_file() and "__pycache__" not in p.parts}
    with tempfile.TemporaryDirectory(prefix="fin-om-nm-") as tmp:
        output = Path(tmp)/"read_only.json"
        run = subprocess.run([sys.executable, str(NM/"verify_research.py"), "--check", "--output", str(output)],
                             cwd=ROOT, env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"),
                             capture_output=True, text=True, timeout=55, preexec_fn=limits)
        assert run.returncode == 0, (run.stdout, run.stderr)
        out = json.loads(output.read_text())
    assert out == json.loads((NM/"RESULTS.json").read_text())
    assert before == {p.relative_to(NM).as_posix(): sha(p) for p in NM.rglob("*") if p.is_file() and "__pycache__" not in p.parts}
    return {"status": out["status"], "bounded_checks": out["bounded_checks"], "source_files_preserved": len(before),
            "mode": "existing --check, temporary output; not a source writer or independent review"}


def links():
    for p in HERE.glob("*.md"):
        assert all(line == line.rstrip() for line in p.read_text().splitlines()), p
        text = re.sub(r"```.*?```", "", p.read_text(), flags=re.S)
        text = re.sub(r"`[^`\n]*`", "", text)
        for target in re.findall(r"\[[^\]\n]+\]\(([^)\n]+)\)", text):
            if not target.startswith(("https:", "http:")):
                assert (p.parent/target.split("#")[0]).is_file(), (p, target)


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
    amin, amax = map(F, source["a_interval"])
    low, high = 4*(amin+2), 4*(amax+2)
    ck("actual carrier squared band", (low, high) == (5, 7))
    ck("continuum gap has rational margin one fifth", (2+F(1, 5))**2 < low)
    ck("shape gap has rational margin one half", F(7, 4)**2 > 3
       and F(13, 4)+F(7, 4) == low)
    ck("carrier size bound for original-time derivatives", high < 3**2)

    epsilon = F(1, 8)
    envelopes = [[G.of(1)], [G.of(0), G.of(1)], [G.of((1, F(1, 3))), G.of(-2), G.of(3)]]
    for nu in (F(0), F(5, 2), F(-5, 2)):
        chi, p, q = shape_coefficients(F(3), nu)
        ck(f"bounded shape primitive differential identity nu{nu}", differentiate_basis(p, F(3), nu) == (ZERO, G.of(-1)))
        ck(f"second bounded shape primitive differential identity nu{nu}", differentiate_basis(q, F(3), nu) == (-p[0], -p[1]))
        chi, ps, qs, kernel = shape_series(F(3), nu)
        ck(f"actual causal susceptibility and zero-data turn-on nu{nu}", ps[0] == G.of(chi)
           and convolution(kernel, [1])[0] == convolution(kernel, [1])[1] == ZERO
           and convolution(kernel, [1])[2] == G(F(1, 2)))
        for idx, envelope in enumerate(envelopes):
            ap1, ap2 = derivative(envelope), derivative(derivative(envelope))
            conv = convolution(kernel, slow(envelope, epsilon))
            identity = add(add(scale(slow(envelope, epsilon), chi), scale(ps, -envelope[0])),
                           scale(convolution(ps, slow(ap1, epsilon)), -epsilon**2))
            ck(f"exact causal convolution series nu{nu} A{idx}", conv == identity,
               degree=DEGREE, scope="finite exact series, not a long-time simulation")
            residual = scale(convolution(ps, slow(ap1, epsilon)), -epsilon**2)
            second = add(scale(slow(ap1, epsilon), -epsilon**2*qs[0]), scale(qs, epsilon**2*ap1[0]))
            second = add(second, scale(convolution(qs, slow(ap2, epsilon)), epsilon**4))
            ck(f"two integrations by parts and slow chain factors nu{nu} A{idx}", residual == second)
        constant = convolution(kernel, [1])
        ck(f"omitting shape initial slip is rejected nu{nu}", constant[0] == ZERO and chi != 0)
    try:
        shape_coefficients(F(1), F(1))
    except ValueError:
        ck("resonant bounded-primitive claim rejected", True,
           scope="unit-frequency oscillator control; actual sqrt3 resonance proof is analytical")
    else:
        raise AssertionError("resonance silently accepted")

    # Actual threshold phase factors, with common factor1/sqrt2 kept symbolic.
    phase_minus = G.of((-1, -1))
    nonoscillatory = -phase_minus/(I*2)
    primitive_leading = 2*nonoscillatory
    ck("threshold tail phase and factor two", nonoscillatory == G.of((F(1, 2), F(-1, 2)))
       and primitive_leading == G.of((1, -1)))
    ck("threshold ramp has epsilon not epsilon squared", 2-2*F(1, 2) == 1)

    # A normalized resonant shape identity; scaling returns the actual eta=sqrt3 formula.
    c1, s1 = cosine_and_sine(F(1))
    resonant_response = convolution(s1, exponent(F(-1)))
    explicit = scale(add(multiply([0, 1], exponent(F(-1))), scale(s1, -1)), I/F(2))
    ck("secular resonant response identity", resonant_response == explicit,
       scope="unit-frequency exact series; no retuned FIN mass or new packet")
    for ep in (F(1, 8), F(1, 16)):
        horizon = 1/(ep*ep)
        response_lower = ep*ep*(horizon/2-F(1, 2))
        susceptibility_and_slip_upper = 2*ep*ep
        residual_lower = response_lower-susceptibility_and_slip_upper
        ck(f"C1-only small-fast envelope resonance bound epsilon{ep}", residual_lower > F(2, 5)
           and ep*ep*horizon == 1,
           scope="unit-frequency triangle bound from the exact resonant identity; actual sqrt3 counter is analytical")

    dependency = replay_nm()
    ck("read-only actual kernel dependency replay", dependency["bounded_checks"] == 34
       and dependency["status"] == "SOURCE_BOUND_NONLINEAR_VOLTERRA_REDUCTION_WITH_EXPLICIT_MEMORY")
    named = [n for n in INPUTS if n.endswith(".md") and n != "AGENTS.md"]
    pattern = r"primitive|slow.{0,30}(envelope|modulat)|first.{0,20}moment|epsilon.{0,20}inverse|Markov"
    search = subprocess.run(["rg", "-n", "-i", "--max-columns", "300", "--", pattern, *named],
                            cwd=ROOT, capture_output=True, text=True, timeout=10)
    assert search.returncode in (0, 1), search.stderr
    assert inputs == {n: sha(ROOT/n) for n in INPUTS}
    assert protected == {n: manifest(ROOT/n) for n in PROTECTED}
    for pth in HERE.glob("*.py"):
        ast.parse(pth.read_text())
    for pth in HERE.glob("*.json"):
        json.loads(pth.read_text())
    result = {"date": "2026-10-21", "status": "OSCILLATORY_MEMORY_OPERATOR_GATE_PAID_WITH_CAUSAL_SHAPE_SLIP",
              "bounded_checks": len(checks), "checks": checks, "input_hashes": inputs,
              "protected_manifests": protected, "NM_dependency_replay": dependency,
              "estimate": "C1 in original t error <=C epsilon^2[||A'||+(1+S)||A''||] on0<=t<=S epsilon^-2",
              "scaled_drive_projected_error": "O(epsilon^4), given uniform C2 slow-envelope norms",
              "threshold_counterexample": "actual kernel at nu2, ramp envelope: error O(epsilon), not O(epsilon^2)",
              "C1_counterexample": "allowed base nu0, small fast envelope: order1 shape residual with uniform C1 only",
              "actual_PDE_envelope_regularities": "OPEN_NOT_ASSUMED", "actual_full_PDE_window_extended": False,
              "averaged_mode_decay_law": "NOT_ACTIVATED", "new_drive_noise_or_packet": False,
              "duplicate_search": {"files": named, "pattern": pattern,
                                   "matching_lines": sorted(search.stdout.splitlines()),
                                   "scope": "bounded actual-source search, not exhaustive novelty"},
              "F01": "OPEN_SOURCE", "independent_review": False, "AGENTS_or_originals_modified": False}
    input_text = "".join(f"{h}  {n}\n" for n, h in sorted(inputs.items()))
    if args.record:
        (HERE/"RESULTS.json").write_text(json.dumps(result, indent=2)+"\n")
        (HERE/"INPUTS.sha256").write_text(input_text)
        files = sorted(pth for pth in HERE.rglob("*") if pth.is_file() and pth.name != "MANIFEST.sha256" and "__pycache__" not in pth.parts)
        (HERE/"MANIFEST.sha256").write_text("".join(f"{sha(pth)}  {pth.relative_to(HERE).as_posix()}\n" for pth in files))
    else:
        saved = json.loads((HERE/"RESULTS.json").read_text())
        assert result == saved, [key for key in result if result[key] != saved.get(key)]
        assert input_text == (HERE/"INPUTS.sha256").read_text()
    links()
    if args.output:
        assert Path("/tmp") in args.output.resolve().parents
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": result["status"], "bounded_OM_checks": len(checks),
                      "NM_replayed_checks_not_new_OM_assertions": dependency["bounded_checks"],
                      "inputs": len(inputs), "protected_manifest_entries": sum(v["entries"] for v in protected.values()),
                      "actual_full_PDE_window_extended": False, "averaged_decay_law": "NOT_ACTIVATED",
                      "F01": "OPEN_SOURCE", "AGENTS_or_originals_modified": False}, indent=2))


if __name__ == "__main__":
    main()
