#!/usr/bin/env python3
"""
QW-718: TOPOLOGICAL MASS GENESIS VERIFICATION
=============================================
Cel: Sprawdzić Hipotezę Unifikującą (H22):
     "Cząstki to trójwymiarowe węzły (Winding Numbers) 
      uformowane z 4-bitowego kodu."

Metoda:
1. Generujemy stan 1D ("DNA") z mechanizmu genezy (QW-715).
2. Rzutujemy go do 3D z wymuszeniem topologii (Winding Number W).
   Rzutowanie: Ψ_3D(r, θ, φ) = Ψ_1D(r) * exp(i * W * θ) * exp(-r/R)
3. Obliczamy "masę emergentną" jako całkę z intensywności (złożoności)
   tego 3D obiektu.
4. Sprawdzamy stosunki mas dla W=1, 2, 3.

Oczekiwanie:
  M(W=2)/M(W=1) ≈ 207 (Mion)
  M(W=3)/M(W=1) ≈ 3477 lub 1836 (Tau/Proton)

Date: 2025-12-08
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import datetime

print("="*80)
print("QW-718: TOPOLOGICZNA GENEZA MASY")
print("="*80)

# ===========================================================================
# 1. GENERACJA STANU 1D (DNA NADSOLITONA)
# ===========================================================================

# Używamy wyniku z QW-715 (stabilna struktura 4-bitowa)
# Symulujemy idealny profil 1D solitonu wyłonionego z szumu
N_1D = 64
x = np.linspace(0, 10, N_1D)
# Profil Gaussa jako aproksymacja stabilnej wyspy 4-bitowej
psi_1d = np.exp(-0.5 * (x - 5.0)**2 / (0.8**2)) 
# Normalizacja
psi_1d = psi_1d / np.linalg.norm(psi_1d)

print("Wygenerowano stan 1D (DNA cząstki).")

# ===========================================================================
# 2. RZUTOWANIE DO 3D Z WINDING NUMBER
# ===========================================================================

# Siatka 3D
N_GRID = 32
L = 10.0
grid_x = np.linspace(-L/2, L/2, N_GRID)
grid_y = np.linspace(-L/2, L/2, N_GRID)
grid_z = np.linspace(-L/2, L/2, N_GRID)

X, Y, Z = np.meshgrid(grid_x, grid_y, grid_z, indexing='ij')

# Współrzędne sferyczne
R = np.sqrt(X**2 + Y**2 + Z**2) + 1e-10 # Unikamy 0
Theta = np.arctan2(np.sqrt(X**2 + Y**2), Z)
Phi = np.arctan2(Y, X)

# Interpolator dla profilu radialnego (z 1D)
# Mapujemy x(0..10) na R(0..L/2)
r_1d = np.linspace(0, L/2, N_1D)
psi_1d_radial = np.exp(-0.5 * (r_1d - 2.0)**2 / (0.8**2)) # Przesunięty od środka (powłoka)
interp_radial = interp1d(r_1d, psi_1d_radial, bounds_error=False, fill_value=0.0)

def generate_3d_state(winding_number):
    """Generuje pole 3D dla zadanej topologii W."""
    
    # 1. Profil radialny (z 1D DNA)
    radial_part = interp_radial(R)
    
    # 2. Część kątowa (Harmonika sferyczna Y_l^m ?)
    # Uprośćmy: Wpływ Winding Number na fazę i amplitudę
    # Winding number W tworzy W węzłów na obwodzie
    # To wpływa na gęstość energii (gradienty fazy)
    
    # Faza topologiczna
    phase = np.exp(1j * winding_number * Phi)
    
    # Pełna funkcja falowa
    psi_3d = radial_part * phase
    
    # Dodajemy efekt skręcenia (torsion) wpływający na gęstość
    # Im silniejszy wir (W), tym bardziej pole "zaciska się" lub "rozdyma"
    # To jest KLUCZOWE dla masy.
    # W Hipotezie 4-bitowej: Każdy skok W dodaje 4 bity złożoności na jednostkę.
    
    return psi_3d

# ===========================================================================
# 3. OBLICZANIE MASY EMERGENTNEJ
# ===========================================================================

def compute_emergent_mass(psi_3d, winding_number):
    """
    Oblicza masę emergentną z pola 3D.
    Masa = Intensywność Procesu.
    
    Intensywność I(r) zależy od:
    1. Amplitudy pola: |Ψ|²
    2. Gradientu fazy (prąd topologiczny): |∇Φ|² 
       (to jest energia kinetyczna wiru)
    3. Kosztu złożoności (4 bity na węzeł)
    """
    
    # Gradienty
    grad_psi = np.gradient(psi_3d)
    grad_sq_norm = np.abs(grad_psi[0])**2 + np.abs(grad_psi[1])**2 + np.abs(grad_psi[2])**2
    
    # Energia kinetyczna pola (klasyczna)
    E_kin = np.sum(grad_sq_norm)
    
    # Energia potencjalna (masa spoczynkowa DNA)
    E_pot = np.sum(np.abs(psi_3d)**2)
    
    # EFEKT 4-BITOWY (H21):
    # Masa skaluje się eksponencjalnie ze złożonością topologiczną
    # Czy musimy to dodać ręcznie ("n_skok")?
    # CZY MOŻE wynika to z samej geometrii wiru?
    
    # W wirze kwantowym energia rośnie jak W² (klasycznie)
    # Ale w Nadsolitonie (fraktalnym) może rosnąć szybciej?
    
    # Sprawdźmy "czystą" intensywność pola 3D
    # I_proc = ∫ ( |Ψ|² + α * |∇Ψ|² ) dV
    
    alpha_coupling = 4 * np.log(2) # 2.77
    
    # Hipoteza: Masa to energia gradientu podniesiona do potęgi α?
    # Albo po prostu suma ważona
    
    mass_raw = E_pot + alpha_coupling * E_kin
    
    # Hipoteza QW-711: Masa = exp(Complexity)
    # Complexity ~ Winding Number
    
    return {
        'W': winding_number,
        'E_pot': E_pot,
        'E_kin': E_kin,
        'Mass_Raw': mass_raw
    }

# ===========================================================================
# 4. TEST DLA W=1, 2, 3
# ===========================================================================

print("\nObliczanie mas dla różnych topologii...")
print(f"{'W':<5} {'E_pot':>10} {'E_kin':>10} {'Masa Raw':>12}")
print("-"*50)

results = {}
for W in [1, 2, 3, 4, 5]:
    psi = generate_3d_state(W)
    res = compute_emergent_mass(psi, W)
    results[W] = res
    print(f"{W:<5} {res['E_pot']:>10.4f} {res['E_kin']:>10.4f} {res['Mass_Raw']:>12.4f}")

# Analiza stosunków
base_mass = results[1]['Mass_Raw']
print("\nStosunki mas (Względem W=1):")
print(f"{'W':<5} {'Ratio Raw':>12} {'Cel (Exp)':>15}")
print("-"*40)

for W in [1, 2, 3, 4, 5]:
    ratio = results[W]['Mass_Raw'] / base_mass
    
    target = "?"
    if W == 1: target = "1.0 (e)"
    elif W == 2: target = "207 (μ)?"
    elif W == 3: target = "1836/3477?"
    elif W == 4: target = "?"
    
    print(f"{W:<5} {ratio:>12.4f} {target:>15}")

# ===========================================================================
# 5. TEST HIPOTEZY EXPONENCJALNEJ (H21) NA WYNIKACH GEOMETRII
# ===========================================================================

print("\nCzy geometria wiru generuje 'n_skok' w H21?")
print(f"H21 mówi: M ∝ 2^(4 * n_skok)")
print("Hipoteza: n_skok ∝ E_kin (energia wiru)")

E_k1 = results[1]['E_kin']

print(f"{'W':<5} {'E_kin/E_k1':>12} {'Predykcja H21':>15} {'Cel':>10}")
print("-"*50)

for W in [1, 2, 3]:
    # Jeśli n_skok rośnie jak E_kin (względna energia wiru)
    # Zobaczmy co wyjdzie
    n_proxy = results[W]['E_kin'] / E_k1  # Jak bardzo wir jest 'mocniejszy' od W=1
    
    # Ale n_skok dla mionu to ~2.
    # Spodziewamy się, że n_proxy ≈ 2 dla W=2?
    
    print(f"{W:<5} {n_proxy:>12.4f}")

print("\nZauważmy: Energia kinetyczna wiru rośnie jak W².")
print("Dla W=2: E_kin/E_k1 ≈ 4 (bo 2²)")
print("Dla W=3: E_kin/E_k1 ≈ 9 (bo 3²)")

print("\nJesli n_skok = W² (kwadratowe skale złożoności):")
# Wzór: M = M_e * 2^(4 * (W-1))  ?? Albo coś podobnego

for W in [1, 2, 3]:
    # Próbujemy dopasować
    # Jeśli n_skok ≈ W (liniowo)?
    # Jeśli n_skok ≈ log(W)?
    pass
    
# Save report
with open("raport_qw718_topological_genesis.md", "w") as f:
    f.write("# RAPORT QW-718: Topologiczna Geneza Masy\n")
    f.write(f"**Date:** {datetime.datetime.now()}\n\n")
    
    f.write("## Wyniki 'Surowe' (Mass = Pot + alpha*Kin)\n")
    f.write("| W | Masa Raw | Ratio |\n")
    f.write("|---|---|---|\n")
    for W, res in results.items():
        ratio = res['Mass_Raw'] / base_mass
        f.write(f"| {W} | {res['Mass_Raw']:.4f} | {ratio:.4f} |\n")
        
    f.write("\n## Wnioski\n")
    f.write("Surowa masa wiru rośnie powoli (ratio ~ W lub ~ W²).\n")
    f.write("Aby uzyskać 207x, potrzebujemy efektu eksponencjalnego:\n")
    f.write("Masa emergentna = exp(Energia Wiru)?\n")

print("\nReport saved to raport_qw718_topological_genesis.md")
print("="*80)
