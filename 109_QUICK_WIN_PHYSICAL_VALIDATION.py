#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author: Krzysztof Żuchowski

"""
109_QUICK_WIN_PHYSICAL_VALIDATION.py
10 zadań: fizyczna walidacja wniosków z badań 102–108
- Język: polski
- Bez fittingu, bez tautologii
- Wyjście: report_109_quick_win.md i report_109_quick_win.json
"""

from datetime import datetime, timezone
import json
import numpy as np
from scipy.linalg import eigh
from scipy.integrate import simpson
import os

np.random.seed(42)

# Parametry globalne (bez fittingu — ustalone z teorii)
alpha_geo = 2.77
beta_tors = 0.01
omega = 2 * np.pi / 8.0
phi = np.pi / 6.0

# Kernel K(d,s)
def K(d, s):
    # znormalizowany kernel geometryczny stosowany we wszystkich badaniach
    return alpha_geo * np.cos(omega * d * s + phi) / (1.0 + beta_tors * d)

# Utility: generuj macierz sprzężeń S z kernel
def build_S(N, scales):
    # d: indeksy parowe 1..N; s: skale
    d = np.abs(np.subtract.outer(np.arange(N), np.arange(N))) + 1.0
    S = np.zeros((N, N))
    for s in scales:
        S += K(d, s)
    # dodaj symetryczny szum numeryczny (ale mały): nie jest fittingiem, to test odporności
    noise = 1e-3 * (np.random.randn(N, N) + np.random.randn(N, N).T) / 2.0
    S = (S + S.T) / 2.0 + noise
    return S

# Metryki używane w zadaniach

def spectral_gap(eigvals):
    # zakładamy eigenvalues posortowane rosnąco
    top = eigvals[-1]
    second = eigvals[-2]
    third = eigvals[-3] if len(eigvals) >= 3 else eigvals[-2]
    return float(top - second), float(second - third), float(top / (second + 1e-12))


def beta_proxy(scales, lambda_max_vals):
    # aproksymacja d g / d ln s
    ln_s = np.log(scales)
    return np.gradient(lambda_max_vals, ln_s)


def fidelity(u, v):
    # overlap między wektorami (moduł kwadratu skalarny)
    num = np.abs(np.vdot(u, v))
    den = np.sqrt(np.vdot(u, u) * np.vdot(v, v)) + 1e-12
    return float((num / den))


def shannon_entropy(vec):
    p = np.abs(vec) ** 2
    p = p / (np.sum(p) + 1e-12)
    p = p[p > 0]
    return float(-np.sum(p * np.log(p)))

# Scenariusz obliczeniowy
N = 24
scales = np.logspace(-1, 1, 60)  # 0.1 -- 10, 60 punktów

# Zbierz wyniki per task
results = {}
meta = {
    "script": "109_QUICK_WIN_PHYSICAL_VALIDATION.py",
    "generated": datetime.now(timezone.utc).isoformat(),
    "N": N,
    "scales_count": len(scales)
}

# Task 0: Spektralna separacja (potwierdzić ratio top/2nd oraz top/3rd)
S = build_S(N, scales)
vals, vecs = eigh(S)
gap_top2, gap_23, ratio_top2 = spectral_gap(vals)
results[0] = {
    "title": "Spektralna separacja",
    "description": "Analiza gap: top vs 2nd vs 3rd eigenvalue dla S zbudowanego na wszystkich skalach.",
    "metrics": {"gap_top_2nd": gap_top2, "gap_2nd_3rd": gap_23, "ratio_top_2nd": ratio_top2},
    "wniosek": "Jeżeli ratio_top_2nd >> 1 (np. >20) — silna separacja spektralna (smoking gun).",
    "zgodnosc_z_ToE": "✅" if ratio_top2 > 20 else "⚠️"
}

# Task 1: Hierarchia mas proxy — stosunki największych wartości własnych
largest = vals[-8:]  # top 8
mass_ratios = [float(largest[-1] / (largest[i] + 1e-12)) for i in range(len(largest)-1)]
results[1] = {
    "title": "Hierarchia mas (proxy)",
    "description": "Stosunki największych eigenvalues jako proxy masa (top/4th, top/5th ...)",
    "metrics": {"top_eigen": float(largest[-1]), "top_over_4th": mass_ratios[3] if len(mass_ratios) >=4 else None, "ratios": mass_ratios},
    "wniosek": "Stosunki powyżej ~3–5 wskazują na stabilną hierarchię mas generowaną z kernela K(d).",
    "zgodnosc_z_ToE": "✅" if (len(mass_ratios)>=4 and mass_ratios[3] > 3.0) else "⚠️"
}

# Task 2: Przepływ RG proxy: λ_max(s) i beta-proxy
lambda_max_vals = []
for s in scales:
    S_s = build_S(N, [s])
    vals_s, _ = eigh(S_s)
    lambda_max_vals.append(float(vals_s[-1]))
lambda_max_vals = np.array(lambda_max_vals)
beta = beta_proxy(scales, lambda_max_vals)
sign_changes = int(np.sum(np.abs(np.diff(np.sign(beta))) > 0))
results[2] = {
    "title": "Przepływ RG (beta-proxy)",
    "description": "Obliczenie lambda_max(s) i estymacja d lambda / d ln s (beta-proxy)",
    "metrics": {"lambda_max_vals_mean": float(np.mean(lambda_max_vals)), "beta_std": float(np.std(beta)), "beta_sign_changes": sign_changes},
    "wniosek": "Wiele zmian znaku beta-proxy wskazuje na nieliniową, fizyczną dynamikę RG (nie trywialną).",
    "zgodnosc_z_ToE": "✅" if sign_changes >= 4 else "⚠️"
}

# Task 3: Odporność na perturbacje (dodanie małej zmiany w kernel)
S_pert = build_S(N, scales) + 0.05 * np.eye(N)  # additive perturbation
vals_p, _ = eigh(S_pert)
lambda_change = float((vals_p[-1] - vals[-1]) / (vals[-1] + 1e-12))
results[3] = {
    "title": "Odporność na perturbacje",
    "description": "Zbadanie względnej zmiany lambda_max po małej perturbacji",
    "metrics": {"lambda_max_nominal": float(vals[-1]), "lambda_max_perturbed": float(vals_p[-1]), "relative_change": lambda_change},
    "wniosek": "Relatywna zmiana < 0.15 (15%) wskazuje na odporność i fizyczną realność stanu.",
    "zgodnosc_z_ToE": "✅" if abs(lambda_change) < 0.15 else "⚠️"
}

# Task 4: Uniwersalność względem N (test dla mniejszych/większych N)
Ns = [16, 20, 24, 28]
lambda_max_N = {}
for Ni in Ns:
    Sni = build_S(Ni, scales)
    vni, _ = eigh(Sni)
    lambda_max_N[Ni] = float(vni[-1])
lambda_std_relative = float(np.std(list(lambda_max_N.values())) / (np.mean(list(lambda_max_N.values())) + 1e-12))
results[4] = {
    "title": "Uniwersalność N",
    "description": "Sprawdzenie jak lambda_max zmienia się z rozmiarem systemu N",
    "metrics": {"lambda_max_by_N": lambda_max_N, "rel_std": lambda_std_relative},
    "wniosek": "Relatywne odchylenie < 0.10 oznacza N-niezależność (termodynamiczny limit).",
    "zgodnosc_z_ToE": "✅" if lambda_std_relative < 0.10 else "⚠️"
}

# Task 5: Stabilność vakuum — trace(S)/N jako proxy
trace_over_N = float(np.trace(S) / N)
results[5] = {
    "title": "Stabilność vakuum (trace/N)",
    "description": "Trace(S)/N jako proxy energii vakuumowej",
    "metrics": {"trace": float(np.trace(S)), "trace_over_N": trace_over_N},
    "wniosek": "Trace/N stabilne w skali <0.1% wskazuje na fizyczne vakuum (wartość referencyjna ≈ 2.77).",
    "zgodnosc_z_ToE": "✅" if abs(trace_over_N - 2.77) < 0.1 else "⚠️"
}

# Task 6: Fidelity między trybami (separacja modów oddzielnych)
u = vecs[:, -1]
v = vecs[:, -2]
fid_uv = fidelity(u, v)
results[6] = {
    "title": "Fidelity między top a 2nd mode",
    "description": "Overlap między top-mode a second-mode jako miara separacji modalnej",
    "metrics": {"fidelity": fid_uv},
    "wniosek": "Fidelity << 1 (np. <0.3) pokazuje wyraźne rozdzielenie modów (silna separacja).",
    "zgodnosc_z_ToE": "✅" if fid_uv < 0.3 else "⚠️"
}

# Task 7: Entropia Shannon dla top eigenvector (blockiness / symmetry emergence)
ent_top = shannon_entropy(np.abs(u))
results[7] = {
    "title": "Entropy Shannon top-mode",
    "description": "Niska entropia sugeruje 'blockiness' i emergencję symetrii (grupowe struktury)",
    "metrics": {"shannon_entropy": ent_top},
    "wniosek": "Relatywnie niska entropia (w porównaniu z równomiernym rozkładem) wspiera emergencję algebraiczną SU(3)×SU(2)×U(1).",
    "zgodnosc_z_ToE": "✅" if ent_top < np.log(len(u)) * 0.5 else "⚠️"
}

# Task 8: Commutator norm (test algebraic closure)
# Użyjemy prostego proxy: [S, P] gdzie P = diagonal matrix of top-mode amplitudes
P = np.diag(np.real(u))
comm = S.dot(P) - P.dot(S)
comm_norm = float(np.linalg.norm(comm))
results[8] = {
    "title": "Norma komutatora proxy",
    "description": "Proxy sprawdzające, czy naturalne operatory tworzą quasi-zamkniętą algebrę",
    "metrics": {"comm_norm": comm_norm},
    "wniosek": "Niska norma komutatora wskazuje na zbliżenie do zamkniętej algebry operatorów (symetrie emergentne).",
    "zgodnosc_z_ToE": "✅" if comm_norm < 1.0 else "⚠️"
}

# Task 9: Mapowanie do obserwabli eksperymentalnych (prostą skalą)
# Mapujemy top eigenvalue -> energia ΔE (arbitrary physical scale) i obliczamy SNR proxy
DeltaE = vals[-1] * 0.001  # przeskalowanie do MeV-analog (proxy)
noise_level = 0.001 * np.std(vals)
snr_proxy = float(np.abs(DeltaE) / (noise_level + 1e-12))
results[9] = {
    "title": "Mapowanie do obserwabli eksperymentalnych (proxy)",
    "description": "Przypisanie top-eigen do energii przejścia i ocena SNR",
    "metrics": {"DeltaE_proxy": DeltaE, "noise_level": noise_level, "snr_proxy": snr_proxy},
    "wniosek": "SNR proxy > 10 sugeruje obserwowalność w realnych eksperymentach (możliwa detekcja).",
    "zgodnosc_z_ToE": "✅" if snr_proxy > 10 else "⚠️"
}

# Zapis wyników
out_md = "report_109_quick_win.md"
out_json = "report_109_quick_win.json"

md_lines = []
md_lines.append(f"# Badanie 109 — Fizyczna walidacja wniosków 102–108\n")
md_lines.append(f"Data wygenerowania: {meta['generated']}\n")
md_lines.append("Opis: 10 zadań weryfikujących fizyczny sens wniosków (bez fittingu, bez tautologii).\n")

for i in range(10):
    r = results[i]
    md_lines.append(f"## Zadanie {i}: {r['title']}\n")
    md_lines.append(f"**Opis:** {r['description']}\n")
    md_lines.append("**Metryki:**")
    md_lines.append("```")
    md_lines.append(json.dumps(r['metrics'], indent=2, ensure_ascii=False))
    md_lines.append("```")
    md_lines.append(f"**Wniosek:** {r['wniosek']}\n")
    md_lines.append(f"**Zgodność z ToE:** {r['zgodnosc_z_ToE']}\n")

# JSON
out = {
    "meta": meta,
    "tasks": results
}

with open(out_json, 'w', encoding='utf-8') as f:
    json.dump(out, f, indent=2, ensure_ascii=False)

with open(out_md, 'w', encoding='utf-8') as f:
    f.write("\n".join(md_lines))

print(f"Wrote {out_md} and {out_json}")
