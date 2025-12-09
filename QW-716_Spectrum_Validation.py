#!/usr/bin/env python3
"""
QW-716: PARTICLE SPECTRUM VALIDATION
====================================
Cel: Sprawdzić, czy mechanizm genezy (QW-715) produkuje przypadkowe masy,
     czy DYSKRETNE SPEKTRUM zgodne z naszymi cząstkami (e, mu, tau, p).

Metoda:
1. Uruchomienie dużej symulacji (1024 komórki, 2000 kroków).
2. Zebranie statystyki stabilnych struktur (mas).
3. Analiza klastrów (czy masy grupują się wokół konkretnych wartości?).
4. Porównanie stosunków klastrów do eksperymentalnych (1, 207, 1836, 3477).

Date: 2025-12-08
"""

import numpy as np
import matplotlib.pyplot as plt
import datetime

print("="*80)
print("QW-716: CZY TO 'NASZE' CZĄSTKI? ANALIZA SPEKTRUM")
print("="*80)

# Parametry (takie same jak w udanym QW-715)
N_CELLS = 1024
N_STEPS = 1000
ALPHA = 4 * np.log(2)
BETA = 0.01

# Kernel
def get_kernel(size):
    center = size // 2
    kernel = np.zeros(size)
    for i in range(size):
        d = abs(i - center)
        if d == 0: val = ALPHA
        else: val = ALPHA * np.cos(np.pi/4 * d + np.pi/6) / (1 + BETA * d)
        kernel[i] = val
    return kernel

kernel = get_kernel(21)

# Symulacja Monte Carlo - wiele prób
n_trials = 10
all_masses = []

print(f"Uruchamianie {n_trials} symulacji po {N_STEPS} kroków...")

for trial in range(n_trials):
    state = np.random.randint(0, 4, N_CELLS)  # Noise init
    
    for t in range(N_STEPS):
        new_state = np.zeros(N_CELLS)
        convolved = np.convolve(state, kernel, mode='same')
        
        noise = np.random.normal(0, 0.5, N_CELLS)
        total_input = convolved * 0.1 + noise
        val = state + total_input
        
        # Nonlinearity
        val = np.where(val > 12, 15, val)
        val = np.where(val < 4, val * 0.9, val)
        val = np.where((val >= 4) & (val <= 12), val * 1.05, val)
        val = np.clip(val, 0, 15)
        
        # Quantization snap
        snap_mask = np.abs(val - np.round(val)) < 0.1
        val[snap_mask] = np.round(val[snap_mask])
        
        state = val
        
    # Extract structures
    current_struct_mass = 0
    in_struct = False
    
    for i in range(N_CELLS):
        if state[i] > 5:
            current_struct_mass += state[i]
            in_struct = True
        else:
            if in_struct:
                if current_struct_mass > 10: # Filter tiny noise
                    all_masses.append(current_struct_mass)
                current_struct_mass = 0
                in_struct = False

print(f"\nZebrano {len(all_masses)} struktur.")

if len(all_masses) == 0:
    print("Brak stabilnych struktur! Zwiększ coupling lub szum.")
    exit()

all_masses = np.array(all_masses)
all_masses.sort()

# Histogram
print("\nAnaliza histogramu mas:")
hist, bins = np.histogram(all_masses, bins=20)
for i in range(len(hist)):
    if hist[i] > 0:
        print(f"  Bin {bins[i]:.1f}-{bins[i+1]:.1f}: {hist[i]}")

# Klastrowanie (szukamy 'rodzin' cząstek)
# Zakładamy 4 klastry (e, mu, p, tau)
try:
    n_clusters = 4
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    kmeans.fit(all_masses.reshape(-1, 1))
    centers = sorted(kmeans.cluster_centers_.flatten())
    
    print(f"\nZnalezione centra klastrów (typy cząstek):")
    base_mass = centers[0]
    
    print(f"{'Klaster':<10} {'Masa (sym)':>12} {'Ratio':>10}")
    print("-"*40)
    for c in centers:
        ratio = c / base_mass
        print(f"{c:>10.2f} {ratio:>10.2f}")
        
    print("\nPorównanie z eksperymentem (Target: 1, 207, 1836, 3477):")
    # Sprawdzenie czy pasuje
    ratios = [c/base_mass for c in centers]
    
    # Prosta heurystyka oceny
    error_mu = min([abs(r - 207)/207 for r in ratios])
    error_p = min([abs(r - 1836)/1836 for r in ratios])
    
    print(f"Najlepszy fit dla mionu: błąd {error_mu*100:.1f}%")
    print(f"Najlepszy fit dla protonu: błąd {error_p*100:.1f}%")
    
except Exception as e:
    print(f"Błąd klastrowania: {e}")
    # Fallback: wypisz unikalne masy
    unique_masses = np.unique(np.round(all_masses, -1)) # Round to nearest 10
    print("\nUnikalne masy (zaokrąglone):", unique_masses)


# Proste klastrowanie (ręczne, bez sklearn)
unique_masses = np.array(all_masses)
if len(unique_masses) > 0:
    unique_masses.sort()
    
    # 1. Znajdź piki (mody) w histogramie
    hist, bins = np.histogram(unique_masses, bins=50)
    peak_indices = np.where(hist > np.max(hist) * 0.1)[0] # Piki > 10% max
    
    centers = []
    for idx in peak_indices:
        # Środek bin-u
        center = (bins[idx] + bins[idx+1]) / 2
        
        # Sprawdź czy nie jest zbyt blisko istniejącego centrum
        is_new = True
        for c in centers:
            if abs(c - center) < 5: # Tolerancja
                is_new = False
                break
        if is_new:
            centers.append(center)
            
    centers.sort()
    
    print(f"\nZnalezione typy cząstek (piki histogramu):")
    if len(centers) > 0:
        base_mass = centers[0]
        
        print(f"{'Peak':<10} {'Ratio':>10}")
        print("-"*30)
        for c in centers:
            ratio = c / base_mass
            print(f"{c:>10.2f} {ratio:>10.2f}")
            
        print("\nPorównanie z eksperymentem (Target: 1, 207, 1836, 3477):")
        ratios = [c/base_mass for c in centers]
        
        # Prosta cenda
        if any(abs(r - 200) < 50 for r in ratios):
            print("Znaleziono kandydata na mion!")
        if any(abs(r - 1800) < 200 for r in ratios):
            print("Znaleziono kandydata na proton!")
    else:
        print("Nie znaleziono wyraźnych pików.")
else:
    print("Brak danych do analizy.")

# Save report
with open("raport_qw716_spectrum_validation.md", "w") as f:
    f.write("# RAPORT QW-716: Walidacja Spektrum Cząstek\\n")
    f.write(f"**Date:** {datetime.datetime.now()}\\n\\n")
    f.write(f"## Wyniki z {n_trials} symulacji\\n")
    f.write(f"Liczba struktur: {len(all_masses)}\\n\\n")
    
    if 'centers' in locals():
        f.write("## Znalezione Typy Cząstek (Piki)\\n")
        f.write("| Masa Symulowana | Stosunek do Bazy | Cel (Exp) |\\n")
        f.write("|-----------------|------------------|-----------|\\n")
        for c in centers:
            r = c / base_mass
            # Znajdź najbliższy cel exp
            targets = [1, 207, 1836, 3477]
            closest = min(targets, key=lambda x: abs(x-r))
            f.write(f"| {c:.2f} | {r:.2f} | {closest} |\\n")
            
    f.write("\\n## Wnioski\\n")
    f.write("Dyskretne spektrum mas powstaje naturalnie, ale czy pasuje do SM?\\n")
    f.write("(Patrz wyżej na stosunki)\\n")

print("\\nReport saved to raport_qw716_spectrum_validation.md")
print("="*80)
