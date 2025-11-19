#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author: Krzysztof Żuchowski

"""
run_optimization_v1.py

Finalny skrypt do przeprowadzenia wielkoskalowej, globalnej optymalizacji
parametrów uniwersalnego jądra sprzężeń hydrodynamicznych w celu odtworzenia
hierarchii mas leptonów z Modelu Standardowego.

- Używa `scipy.differential_evolution` do robustnego, globalnego poszukiwania.
- Loguje postępy do pliku CSV w czasie rzeczywistym.
- Zapisuje finalne wyniki, parametry i wykresy.
- Zaprojektowany do długich, nieprzerwanych sesji na platformach takich jak Kaggle.
"""

import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution
from scipy.linalg import eigh
import matplotlib.pyplot as plt
import os
import time

print("="*80)
print("START: GLOBALNA OPTYMALIZACJA HIERARCHII MAS NADSOLITONU")
print("="*80)

# --- 1. KONFIGURACJA SYMULACJI I OPTYMALIZACJI ---

# Parametry siatki i pola (zgodne z udaną analizą)
N_OCTAVES = 12
NX, NY = 64, 64
LX, LY = 20.0, 20.0
M0_SQUARED = 1.0  # Bazowa (goła) masa kwadratowa

# Parametry optymalizatora
MAX_ITER = 500  # Zwiększona liczba iteracji dla pełnej zbieżności
POP_SIZE = 20    # Rozmiar populacji (15-20 jest standardem)
TOLERANCE = 0.01
SEED = 42

# Ścieżki do plików wyjściowych
OUTPUT_DIR = "optimization_results"
os.makedirs(OUTPUT_DIR, exist_ok=True)
LOG_FILE = os.path.join(OUTPUT_DIR, "optimization_log.csv")
BEST_PARAMS_FILE = os.path.join(OUTPUT_DIR, "best_parameters.txt")
CONVERGENCE_PLOT = os.path.join(OUTPUT_DIR, "convergence_plot.png")
FINAL_SPECTRUM_PLOT = os.path.join(OUTPUT_DIR, "final_mass_spectrum.png")


# --- 2. DEFINICJE FIZYCZNE (JĄDRO, HAMILTONIAN) ---

def initialize_vortex_fields(n_octaves, nx, ny, Lx, Ly):
    """Inicjalizuje 12 pól oktawowych z różnymi liczbami 'krętności'."""
    psi_octaves = []
    winding_numbers = []
    x = np.linspace(-Lx/2, Lx/2, nx)
    y = np.linspace(-Ly/2, Ly/2, ny)
    X, Y = np.meshgrid(x, y, indexing='ij')
    r_grid = np.sqrt(X**2 + Y**2)
    phi_grid = np.arctan2(Y, X)

    for o in range(n_octaves):
        n_wind = o % 3
        winding_numbers.append(n_wind)
        core_rad = 2.0 + 0.2 * (o % 4)
        amp = 0.3 + 0.05 * (o % 3)

        f_r = amp * np.tanh(r_grid / core_rad)
        phase = n_wind * phi_grid
        psi = f_r * np.exp(1j * phase)
        psi_octaves.append(psi)

    print(f"✓ Zainicjalizowano {n_octaves} pól oktawowych.")
    print(f"  Winding numbers: {winding_numbers}")
    return psi_octaves, winding_numbers

def K_universal(i, j, psi_i, psi_j, n_i, n_j, params):
    """Uniwersalne jądro sprzężeń."""
    A, omega, phi_tors, alpha_geo, alpha_res, beta_topo = params

    # 1. Geometryczne
    K_geo = 2**(-alpha_geo * abs(i - j))

    # 2. Rezonansowe
    corr = np.abs(np.corrcoef(np.abs(psi_i.flatten()), np.abs(psi_j.flatten()))[0, 1])
    K_res = 1.0 + alpha_res * (corr if np.isfinite(corr) else 0.0)

    # 3. Skrętne (Torsional)
    r_char = 1.0  # Uproszczenie: reprezentatywna skala
    K_tors = A * np.cos(omega * r_char + phi_tors)

    # 4. Topologiczne
    K_topo = np.exp(-beta_topo * abs(n_i - n_j))

    return K_geo * K_res * (1.0 + K_tors) * K_topo

def construct_hamiltonian(psi_octaves, winding_numbers, params):
    """Konstruuje macierz oddziaływań (Hamiltonian)."""
    n_oct = len(psi_octaves)
    H = np.zeros((n_oct, n_oct))

    for i in range(n_oct):
        for j in range(n_oct):
            if i == j:
                H[i, j] = M0_SQUARED
            else:
                H[i, j] = K_universal(i, j, psi_octaves[i], psi_octaves[j],
                                      winding_numbers[i], winding_numbers[j], params)
    return H

# Inicjalizacja pól (robimy to raz na początku)
PSI_OCTAVES, WINDING_NUMBERS = initialize_vortex_fields(N_OCTAVES, NX, NY, LX, LY)


# --- 3. FUNKCJA KOSZTU DLA OPTYMALIZACJI ---

def cost_function(params_vec):
    """Funkcja kosztu: logarytmiczna różnica od stosunków mas leptonów."""
    try:
        H = construct_hamiltonian(PSI_OCTAVES, WINDING_NUMBERS, params_vec)
        eigenvalues = eigh(H, eigvals_only=True)

        positive_evals = eigenvalues[eigenvalues > 0]

        if len(positive_evals) < 3:
            return 1e6  # Duża kara, jeśli nie ma 3 stabilnych stanów

        masses = np.sqrt(positive_evals)
        m1, m2, m3 = masses[0], masses[1], masses[2] # 3 najlżejsze stabilne

        ratio_21 = m2 / m1
        ratio_31 = m3 / m1

        target_21 = 206.77  # m_mu / m_e
        target_31 = 3477.15 # m_tau / m_e

        cost = (np.log(ratio_21) - np.log(target_21))**2 + \
               (np.log(ratio_31) - np.log(target_31))**2

        return cost if np.isfinite(cost) else 1e6
    except Exception:
        return 1e7 # Kara za błędy numeryczne

# --- 4. PĘTLA OPTYMALIZACYJNA Z LOGOWANIEM ---

# Granice dla parametrów [A, omega, phi_tors, alpha_geo, alpha_res, beta_topo]
bounds = [
    (0.1, 2.0),      # A
    (0.1, 2.0),      # omega
    (0.0, 2 * np.pi),# phi_tors
    (0.01, 0.5),     # alpha_geo
    (0.1, 5.0),      # alpha_res
    (0.1, 3.0)       # beta_topo
]

# Przygotowanie pliku logów
with open(LOG_FILE, "w", newline="") as f:
    f.write("iteration,cost,A,omega,phi_tors,alpha_geo,alpha_res,beta_topo\n")

# Callback do logowania postępów
iteration_count = 0
def callback(xk, convergence):
    global iteration_count
    cost = cost_function(xk)
    with open(LOG_FILE, "a", newline="") as f:
        f.write(f"{iteration_count},{cost},{','.join(map(str, xk))}\n")

    if iteration_count % 10 == 0:
        print(f"  Iteracja {iteration_count:4d}: Najlepszy koszt = {cost:.6f}")
    iteration_count += 1

print("\nRozpoczynanie globalnej optymalizacji...")
print(f"Maksymalna liczba iteracji: {MAX_ITER}")
print(f"Logi będą zapisywane w: {LOG_FILE}")
print("="*60)

start_time = time.time()

result = differential_evolution(
    cost_function,
    bounds,
    strategy='best1bin',
    maxiter=MAX_ITER,
    popsize=POP_SIZE,
    tol=TOLERANCE,
    seed=SEED,
    callback=callback,
    disp=False,
    workers=-1  # Użyj wszystkich dostępnych rdzeni
)

end_time = time.time()

print("\n" + "="*60)
print("OPTYMALIZACJA ZAKOŃCZONA")
print("="*60)
print(f"Czas trwania: {end_time - start_time:.1f} sekund")
print(f"Liczba iteracji: {result.nit}")
print(f"Sukces: {result.success}")
print(f"Wiadomość: {result.message}")
print(f"\nFinalny, najlepszy koszt: {result.fun:.6f}")

optimal_params = result.x
param_names = ['A', 'omega', 'phi_tors', 'alpha_geo', 'alpha_res', 'beta_topo']
with open(BEST_PARAMS_FILE, "w") as f:
    f.write("OPTIMAL PARAMETERS:\n")
    for name, val in zip(param_names, optimal_params):
        print(f"  {name:12s}: {val:8.4f}")
        f.write(f"  {name}: {val}\n")
print(f"✓ Optymalne parametry zapisane do: {BEST_PARAMS_FILE}")


# --- 5. FINALNA ANALIZA I WIZUALIZACJA ---

print("\n" + "="*80)
print("FINALNA ANALIZA DLA OPTYMALNYCH PARAMETRÓW")
print("="*80)

# Wykres zbieżności
log_df = pd.read_csv(LOG_FILE)
plt.figure(figsize=(10, 6))
plt.plot(log_df['iteration'], log_df['cost'])
plt.yscale('log')
plt.title('Zbieżność Funkcji Kosztu w Trakcie Optymalizacji')
plt.xlabel('Iteracja')
plt.ylabel('Koszt (logarytmiczny)')
plt.grid(True, which='both', linestyle='--')
plt.savefig(CONVERGENCE_PLOT)
plt.close()
print(f"✓ Wykres zbieżności zapisany do: {CONVERGENCE_PLOT}")

# Finalne widmo mas
H_final = construct_hamiltonian(PSI_OCTAVES, WINDING_NUMBERS, optimal_params)
eigenvalues_final = eigh(H_final, eigvals_only=True)
positive_evals_final = eigenvalues_final[eigenvalues_final > 0]
masses_final = np.sqrt(positive_evals_final)

m1, m2, m3 = masses_final[0], masses_final[1], masses_final[2]
ratio_21_final = m2 / m1
ratio_31_final = m3 / m1

print("\nFinalne, zoptymalizowane stosunki mas:")
print(f"  m_μ / m_e: Model = {ratio_21_final:.2f}, SM = 206.77 (Błąd: {abs(ratio_21_final/206.77-1)*100:.1f}%)")
print(f"  m_τ / m_e: Model = {ratio_31_final:.2f}, SM = 3477.15 (Błąd: {abs(ratio_31_final/3477.15-1)*100:.1f}%)")

# Wykres widma mas
plt.figure(figsize=(10, 6))
plt.bar(range(len(masses_final)), masses_final, color='skyblue', edgecolor='black')
plt.yscale('log')
plt.title('Finalne, Zoptymalizowane Widmo Mas')
plt.xlabel('Indeks Stanu Własnego')
plt.ylabel('Masa (jednostki modelu, skala log)')
plt.grid(True, which='both', linestyle='--')
for i in range(3):
    plt.text(i, masses_final[i]*1.2, f'Gen {i+1}\n{masses_final[i]:.2f}', ha='center', color='red')
plt.savefig(FINAL_SPECTRUM_PLOT)
plt.close()
print(f"✓ Wykres finalnego widma mas zapisany do: {FINAL_SPECTRUM_PLOT}")

print("\n" + "="*80)
print("ANALIZA ZAKOŃCZONA")
print("="*80)
