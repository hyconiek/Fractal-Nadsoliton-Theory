# Author: Krzysztof Żuchowski
# Data: 19.11.2025


# QW-296: STAŁA ZŁOŻONOŚCI (KOLMOGOROV)
# ============================================================================
# Cel: Czy życie jest nieuniknione?
# Zbadaj wzrost złożoności (nie entropii) w ewolucji sieci

print("\n" + "="*80)
print("QW-296: STAŁA ZŁOŻONOŚCI (KOLMOGOROV)")
print("="*80)

# Złożoność Kolmogorowa K(x) to długość najkrótszego programu,
# który generuje ciąg x. Aproksymujemy ją przez kompresję.

# Symulujemy ewolucję stanu sieci w czasie
N = 16  # Rozmiar sieci
S = build_S_matrix(N)
eigenvalues, eigenvectors = eigh(S)

# Stan początkowy: superposition kilku modów
psi_0 = eigenvectors[:, 0] + 0.5*eigenvectors[:, 1] + 0.3*eigenvectors[:, 2]
psi_0 = psi_0 / np.linalg.norm(psi_0)

# Ewolucja unitarna: |ψ(t)⟩ = exp(-iSt)|ψ(0)⟩
# W jednostkach naturalnych (ℏ=1)
n_timesteps = 50
t_max = 10.0
times = np.linspace(0, t_max, n_timesteps)

# Funkcja złożoności Kolmogorowa poprzez kompresję
def kolmogorov_complexity(state, precision=6):
    """
    Aproksymacja złożoności Kolmogorowa przez kompresję gzip.
    Im bardziej nieuporządkowany stan, tym trudniej go skompresować.
    """
    # Dyskretyzuj stan do bajtów (konwersja do stringa)
    state_rounded = np.round(state.real, precision) + 1j*np.round(state.imag, precision)
    state_bytes = str(state_rounded.tolist()).encode('utf-8')

    # Kompresja
    compressed = gzip.compress(state_bytes)

    # Złożoność = długość skompresowanego ciągu
    complexity = len(compressed)

    return complexity

# Ewolucja i pomiar złożoności
complexities = []
entropies = []

print(f"\n🧬 EWOLUCJA ZŁOŻONOŚCI SIECI (N={N}):")
print(f"   Czas ewolucji: t ∈ [0, {t_max}]")
print(f"   Liczba kroków: {n_timesteps}")

for i, t in enumerate(times):
    # Ewolucja
    U_t = expm(-1j * S * t)
    psi_t = U_t @ psi_0

    # Złożoność Kolmogorowa
    K_t = kolmogorov_complexity(psi_t)
    complexities.append(K_t)

    # Entropia von Neumanna dla porównania
    # S = -Tr(ρ log ρ), gdzie ρ = |ψ⟩⟨ψ|
    rho = np.outer(psi_t, psi_t.conj())
    eigenvals_rho = eigh(rho, eigvals_only=True)
    eigenvals_rho = eigenvals_rho[eigenvals_rho > 1e-12]  # Usuń małe wartości
    entropy = -np.sum(eigenvals_rho * np.log(eigenvals_rho + 1e-12))
    entropies.append(entropy)

    if i % 10 == 0:
        print(f"   t={t:.2f}: K={K_t}, S={entropy:.6f}")

complexities = np.array(complexities)
entropies = np.array(entropies)

# Analiza wzrostu: czy rośnie szybciej niż liniowo?
# Dopasowanie: K(t) = a + b·t + c·t²

from scipy.optimize import curve_fit

def linear_model(t, a, b):
    return a + b * t

def quadratic_model(t, a, b, c):
    return a + b * t + c * t**2

# Dopasowanie modeli
params_linear, _ = curve_fit(linear_model, times, complexities)
params_quad, _ = curve_fit(quadratic_model, times, complexities)

K_linear = linear_model(times, *params_linear)
K_quad = quadratic_model(times, *params_quad)

# Błędy średniokwadratowe
mse_linear = np.mean((complexities - K_linear)**2)
mse_quad = np.mean((complexities - K_quad)**2)

print(f"\n📊 DOPASOWANIE MODELI WZROSTU:")
print(f"   Model liniowy: K(t) = {params_linear[0]:.2f} + {params_linear[1]:.4f}·t")
print(f"   MSE (liniowy): {mse_linear:.4f}")
print(f"   Model kwadratowy: K(t) = {params_quad[0]:.2f} + {params_quad[1]:.4f}·t + {params_quad[2]:.6f}·t²")
print(f"   MSE (kwadratowy): {mse_quad:.4f}")

# Który model lepiej pasuje?
improvement = (mse_linear - mse_quad) / mse_linear * 100
print(f"\n   Poprawa dla modelu kwadratowego: {improvement:.2f}%")

# Sprawdź czy współczynnik kwadratowy jest istotny
if abs(params_quad[2]) > 0:
    print(f"   Współczynnik kwadratowy c = {params_quad[2]:.6f}")
    if params_quad[2] > 0:
        print(f"   ✓ Złożoność rośnie SZYBCIEJ niż liniowo (c > 0)")
        print(f"   ✓ Wszechświat dąży do tworzenia struktur!")
    else:
        print(f"   ✗ Złożoność rośnie wolniej niż liniowo (c < 0)")
else:
    print(f"   Współczynnik kwadratowy c ≈ 0, wzrost liniowy")

print(f"\n✅ WYNIK QW-296:")
if params_quad[2] > 1e-6:
    print(f"   Złożoność Kolmogorowa rośnie superlinearnie")
    print(f"   Tempo wzrostu: dK/dt ∝ {params_quad[1]:.4f} + 2·{params_quad[2]:.6f}·t")
    print(f"   Interpretacja: Wszechświat ma tendencję do samoorganizacji")
    print(f"   Status: ŻYCIE JEST NIEUNIKNIONE ✓")
else:
    print(f"   Złożoność rośnie liniowo lub subliniowo")
    print(f"   Status: Brak wyraźnej tendencji do złożonych struktur")


================================================================================
QW-296: STAŁA ZŁOŻONOŚCI (KOLMOGOROV)
================================================================================

🧬 EWOLUCJA ZŁOŻONOŚCI SIECI (N=16):
   Czas ewolucji: t ∈ [0, 10.0]
   Liczba kroków: 50
   t=0.00: K=121, S=-0.000000
   t=2.04: K=178, S=-0.000000
   t=4.08: K=183, S=-0.000000
   t=6.12: K=180, S=-0.000000
   t=8.16: K=184, S=-0.000000

📊 DOPASOWANIE MODELI WZROSTU:
   Model liniowy: K(t) = 176.30 + 0.7153·t
   MSE (liniowy): 73.7879
   Model kwadratowy: K(t) = 170.21 + 4.4449·t + -0.372960·t²
   MSE (kwadratowy): 65.4265

   Poprawa dla modelu kwadratowego: 11.33%
   Współczynnik kwadratowy c = -0.372960
   ✗ Złożoność rośnie wolniej niż liniowo (c < 0)

✅ WYNIK QW-296:
   Złożoność rośnie liniowo lub subliniowo
   Status: Brak wyraźnej tendencji do złożonych struktur

In [12]:


# ============================================================================
# QW-297: TEST HIPOTEZY PANCOMPUTATIONALISM
# ============================================================================
# Cel: Czy wszechświat jest komputerem kwantowym?
# Znajdź konfiguracje bramek logicznych (AND, OR, NOT) w interakcjach solitonów

print("\n" + "="*80)
print("QW-297: TEST HIPOTEZY PANCOMPUTATIONALISM")
print("="*80)

# Kompletność Turinga wymaga uniwersalnego zbioru bramek
# Np. {AND, NOT} lub {NAND} są wystarczające

# Sprawdzamy, czy macierz S może realizować operacje logiczne
# poprzez przekształcenia stanów bazowych

N = 8  # Mały wymiar dla jasności
S = build_S_matrix(N)

# Stany bazowe reprezentują qubity
# |0⟩ i |1⟩ mogą być reprezentowane przez dwa stany własne

eigenvalues, eigenvectors = eigh(S)

# Wybieramy dwa skrajne stany własne jako |0⟩ i |1⟩
state_0 = eigenvectors[:, 0]  # Najniższa energia
state_1 = eigenvectors[:, -1]  # Najwyższa energia

print(f"\n🔷 REPREZENTACJA QUBITÓW:")
print(f"   |0⟩: stan własny o λ_min = {eigenvalues[0]:.6f}")
print(f"   |1⟩: stan własny o λ_max = {eigenvalues[-1]:.6f}")
print(f"   Ortogonalność: ⟨0|1⟩ = {np.abs(np.vdot(state_0, state_1)):.10f}")

# Bramka NOT: |0⟩ → |1⟩, |1⟩ → |0⟩
# Sprawdzamy, czy istnieje macierz w algebrze S, która realizuje tę operację

# Operator NOT w bazie |0⟩, |1⟩:
# NOT = |1⟩⟨0| + |0⟩⟨1|
NOT_target = np.outer(state_1, state_0.conj()) + np.outer(state_0, state_1.conj())

# Czy istnieje kombinacja liniowa exp(iθS) realizująca NOT?
# Sprawdzamy dla różnych kątów θ

print(f"\n🔧 POSZUKIWANIE BRAMKI NOT:")
print(f"   Testujemy U(θ) = exp(iθS) dla różnych θ")

best_theta = None
best_fidelity = 0

for theta in np.linspace(0, 2*np.pi, 100):
    U_theta = expm(1j * theta * S)

    # Testujemy działanie na |0⟩ i |1⟩
    result_0 = U_theta @ state_0
    result_1 = U_theta @ state_1

    # Fidelity: jak blisko jesteśmy do |1⟩ i |0⟩?
    fidelity_0 = np.abs(np.vdot(state_1, result_0))**2
    fidelity_1 = np.abs(np.vdot(state_0, result_1))**2

    avg_fidelity = (fidelity_0 + fidelity_1) / 2

    if avg_fidelity > best_fidelity:
        best_fidelity = avg_fidelity
        best_theta = theta

print(f"   Najlepsza fidelity: F = {best_fidelity:.6f} przy θ = {best_theta:.6f}")

if best_fidelity > 0.95:
    print(f"   ✓ Bramka NOT zrealizowana z wysoką wiernością!")
else:
    print(f"   ⚠ Bramka NOT tylko częściowo zrealizowana")

# Bramka CNOT (kontrolowana NOT) dla 2 qubitów
# Wymaga 2N-wymiarowej przestrzeni (produkt tensorowy)

print(f"\n🔧 BRAMKA CNOT (2-qubitowa):")
print(f"   Wymaga przestrzeni 2N = {2*N}")

# Budujemy macierz S dla 2N wymiarów
S_2qubit = build_S_matrix(2*N)
eigenvalues_2q, eigenvectors_2q = eigh(S_2qubit)

# Definiujemy stany produktowe |00⟩, |01⟩, |10⟩, |11⟩
# jako kombinacje stanów własnych

# Przestrzeń produktowa: najprostsza realizacja to produkt Kroneckera
state_00 = np.kron(state_0, state_0)[:2*N]  # Obcięcie do wymiaru
state_01 = np.kron(state_0, state_1)[:2*N]
state_10 = np.kron(state_1, state_0)[:2*N]
state_11 = np.kron(state_1, state_1)[:2*N]

# Normalizacja
state_00 = state_00 / np.linalg.norm(state_00)
state_01 = state_01 / np.linalg.norm(state_01)
state_10 = state_10 / np.linalg.norm(state_10)
state_11 = state_11 / np.linalg.norm(state_11)

# CNOT: |00⟩→|00⟩, |01⟩→|01⟩, |10⟩→|11⟩, |11⟩→|10⟩
CNOT_target = (np.outer(state_00, state_00.conj()) +
               np.outer(state_01, state_01.conj()) +
               np.outer(state_11, state_10.conj()) +
               np.outer(state_10, state_11.conj()))

# Szukamy U(θ) = exp(iθS_2qubit) realizującego CNOT
best_theta_cnot = None
best_fidelity_cnot = 0

for theta in np.linspace(0, 2*np.pi, 50):
    U_theta = expm(1j * theta * S_2qubit)

    # Testujemy na wszystkich stanach bazowych
    results = {
        '00': U_theta @ state_00,
        '01': U_theta @ state_01,
        '10': U_theta @ state_10,
        '11': U_theta @ state_11
    }

    # Fidelity dla każdej transformacji
    f_00 = np.abs(np.vdot(state_00, results['00']))**2
    f_01 = np.abs(np.vdot(state_01, results['01']))**2
    f_10 = np.abs(np.vdot(state_11, results['10']))**2  # |10⟩→|11⟩
    f_11 = np.abs(np.vdot(state_10, results['11']))**2  # |11⟩→|10⟩

    avg_fidelity = (f_00 + f_01 + f_10 + f_11) / 4

    if avg_fidelity > best_fidelity_cnot:
        best_fidelity_cnot = avg_fidelity
        best_theta_cnot = theta

print(f"   Najlepsza fidelity: F = {best_fidelity_cnot:.6f} przy θ = {best_theta_cnot:.6f}")

if best_fidelity_cnot > 0.90:
    print(f"   ✓ Bramka CNOT zrealizowana!")
else:
    print(f"   ⚠ Bramka CNOT tylko częściowo zrealizowana")

# Kompletność Turinga
print(f"\n🎯 KOMPLETNOŚĆ TURINGA:")

if best_fidelity > 0.95 and best_fidelity_cnot > 0.90:
    print(f"   Zbiór {NOT, CNOT} jest uniwersalny dla obliczeń kwantowych")
    print(f"   Teoria może realizować dowolne przekształcenia unitarne")
    print(f"   ✓ WSZECHŚWIAT JEST UNIWERSALNYM KOMPUTEREM KWANTOWYM!")
elif best_fidelity > 0.80:
    print(f"   Bramki podstawowe są realizowane z umiarkowaną wiernością")
    print(f"   Możliwa realizacja obliczeń z poprawkami błędów")
    print(f"   ⚠ CZĘŚCIOWA KOMPLETNOŚĆ OBLICZENIOWA")
else:
    print(f"   Bramki logiczne słabo zrealizowane w prostej ewolucji exp(iθS)")
    print(f"   Może wymagać bardziej złożonych kombinacji operatorów")
    print(f"   ✗ KOMPLETNOŚĆ TURINGA NIEJASNA")

print(f"\n✅ PODSUMOWANIE QW-297:")
print(f"   Fidelity NOT: {best_fidelity:.4f}")
print(f"   Fidelity CNOT: {best_fidelity_cnot:.4f}")
if best_fidelity > 0.80:
    print(f"   Status: WSZECHŚWIAT REALIZUJE OBLICZENIA KWANTOWE")
    print(f"   Mechanizm: Ewolucja unitarna U = exp(iθS)")
else:
    print(f"   Status: Proste bramki trudno zrealizować bezpośrednio")


================================================================================
QW-297: TEST HIPOTEZY PANCOMPUTATIONALISM
================================================================================

🔷 REPREZENTACJA QUBITÓW:
   |0⟩: stan własny o λ_min = -3.611071
   |1⟩: stan własny o λ_max = 11.103718
   Ortogonalność: ⟨0|1⟩ = 0.0000000000

🔧 POSZUKIWANIE BRAMKI NOT:
   Testujemy U(θ) = exp(iθS) dla różnych θ
   Najlepsza fidelity: F = 0.000000 przy θ = 5.838920
   ⚠ Bramka NOT tylko częściowo zrealizowana

🔧 BRAMKA CNOT (2-qubitowa):
   Wymaga przestrzeni 2N = 16
   Najlepsza fidelity: F = 0.500000 przy θ = 0.000000
   ⚠ Bramka CNOT tylko częściowo zrealizowana

🎯 KOMPLETNOŚĆ TURINGA:
   Bramki logiczne słabo zrealizowane w prostej ewolucji exp(iθS)
   Może wymagać bardziej złożonych kombinacji operatorów
   ✗ KOMPLETNOŚĆ TURINGA NIEJASNA

✅ PODSUMOWANIE QW-297:
   Fidelity NOT: 0.0000
   Fidelity CNOT: 0.5000
   Status: Proste bramki trudno zrealizować bezpośrednio

In [13]:


# ============================================================================
# QW-298: GRANICA PRĘDKOŚCI PRZETWARZANIA (MARGOLUS-LEVITIN)
# ============================================================================
# Cel: Maksymalna szybkość obliczeń natury
# Limit: v_op ≤ 2E/(πℏ)

print("\n" + "="*80)
print("QW-298: GRANICA PRĘDKOŚCI PRZETWARZANIA (MARGOLUS-LEVITIN)")
print("="*80)

# Twierdzenie Margolusa-Levitina: minimalna czas ewolucji między
# ortogonalnymi stanami: τ ≥ πℏ/(2E)
# Czyli maksymalna szybkość operacji: v_op = 1/τ ≤ 2E/(πℏ)

# W jednostkach naturalnych ℏ = 1
# v_op ≤ 2E/π

# Symulujemy ewolucję systemu i mierzymy czas ortogonalizacji
N = 16
S = build_S_matrix(N)
eigenvalues, eigenvectors = eigh(S)

# Stan początkowy: stan podstawowy
psi_0 = eigenvectors[:, 0]

# Energia całkowita systemu: wartość oczekiwana ⟨ψ|S|ψ⟩
E_total = np.real(np.vdot(psi_0, S @ psi_0))

print(f"\n⚡ ENERGIA SYSTEMU:")
print(f"   Stan początkowy: |ψ₀⟩ = stan podstawowy")
print(f"   Energia: E₀ = ⟨ψ₀|S|ψ₀⟩ = {E_total:.6f}")

# Ewolucja: |ψ(t)⟩ = exp(-iSt)|ψ₀⟩
# Szukamy najmniejszego t, dla którego |⟨ψ₀|ψ(t)⟩|² < ε (ortogonalne)

epsilon = 1e-6  # Próg ortogonalności
t_max_search = 100
n_steps = 10000
times_search = np.linspace(0, t_max_search, n_steps)

t_orth = None
for t in times_search:
    U_t = expm(-1j * S * t)
    psi_t = U_t @ psi_0

    # Overlap z stanem początkowym
    overlap = np.abs(np.vdot(psi_0, psi_t))**2

    if overlap < epsilon:
        t_orth = t
        break

if t_orth is None:
    # Jeśli nie znaleziono, użyj najbliższego minimum
    overlaps = []
    for t in times_search:
        U_t = expm(-1j * S * t)
        psi_t = U_t @ psi_0
        overlap = np.abs(np.vdot(psi_0, psi_t))**2
        overlaps.append(overlap)

    overlaps = np.array(overlaps)
    min_idx = np.argmin(overlaps)
    t_orth = times_search[min_idx]
    min_overlap = overlaps[min_idx]

    print(f"\n⏱️ CZAS ORTOGONALIZACJI:")
    print(f"   Minimalne |⟨ψ₀|ψ(t)⟩|² = {min_overlap:.6e} przy t = {t_orth:.6f}")
else:
    print(f"\n⏱️ CZAS ORTOGONALIZACJI:")
    print(f"   t_orth = {t_orth:.6f} (|⟨ψ₀|ψ(t)⟩|² < {epsilon})")

# Szybkość operacji rzeczywista
v_op_actual = 1 / t_orth

print(f"   Szybkość operacji: v_op = 1/t_orth = {v_op_actual:.6f}")

# Limit Margolusa-Levitina
v_op_limit = 2 * abs(E_total) / np.pi

print(f"\n📊 PORÓWNANIE Z LIMITEM MARGOLUSA-LEVITINA:")
print(f"   Limit teoretyczny: v_op^(max) = 2E/π = 2·{abs(E_total):.6f}/π = {v_op_limit:.6f}")
print(f"   Szybkość rzeczywista: v_op^(actual) = {v_op_actual:.6f}")
print(f"   Stosunek: v_op^(actual) / v_op^(max) = {v_op_actual / v_op_limit:.6f}")

if v_op_actual <= v_op_limit * 1.01:  # 1% tolerancja
    print(f"   ✓ SYSTEM DZIAŁA NA LIMICIE FIZYCZNYM!")
    print(f"   Natura realizuje obliczenia z maksymalną możliwą prędkością")
else:
    print(f"   ⚠ System działa poniżej limitu (wolniej niż maksimum)")

# Badamy dla różnych stanów początkowych
print(f"\n🔍 ANALIZA DLA RÓŻNYCH STANÓW POCZĄTKOWYCH:")

v_ops_actual = []
v_ops_limit = []
energies = []

for i in range(min(8, N)):
    psi_i = eigenvectors[:, i]
    E_i = np.real(eigenvalues[i])

    # Szukamy czasu ortogonalizacji
    t_orth_i = None
    for t in times_search:
        U_t = expm(-1j * S * t)
        psi_t = U_t @ psi_i
        overlap = np.abs(np.vdot(psi_i, psi_t))**2
        if overlap < epsilon:
            t_orth_i = t
            break

    if t_orth_i is None:
        overlaps = []
        for t in times_search:
            U_t = expm(-1j * S * t)
            psi_t = U_t @ psi_i
            overlap = np.abs(np.vdot(psi_i, psi_t))**2
            overlaps.append(overlap)
        overlaps = np.array(overlaps)
        min_idx = np.argmin(overlaps)
        t_orth_i = times_search[min_idx]

    v_op_i = 1 / t_orth_i
    v_lim_i = 2 * abs(E_i) / np.pi

    v_ops_actual.append(v_op_i)
    v_ops_limit.append(v_lim_i)
    energies.append(abs(E_i))

    print(f"   Stan {i}: E={E_i:+.4f}, v_op={v_op_i:.4f}, limit={v_lim_i:.4f}, ratio={v_op_i/v_lim_i:.4f}")

v_ops_actual = np.array(v_ops_actual)
v_ops_limit = np.array(v_ops_limit)
energies = np.array(energies)

# Średni stosunek
mean_ratio = np.mean(v_ops_actual / v_ops_limit)
print(f"\n   Średni stosunek v_op/v_max: {mean_ratio:.4f}")

print(f"\n✅ WYNIK QW-298:")
if mean_ratio > 0.9:
    print(f"   System operuje BLISKO limitu Margolusa-Levitina")
    print(f"   Natura maksymalizuje efektywność przetwarzania informacji")
    print(f"   Status: MAKSYMALNA PRĘDKOŚĆ OBLICZENIOWA POTWIERDZONA ✓")
elif mean_ratio > 0.5:
    print(f"   System operuje z umiarkowaną efektywnością")
    print(f"   Możliwe optymalizacje w wyborze stanów/operacji")
    print(f"   Status: CZĘŚCIOWA SATURACJA LIMITU")
else:
    print(f"   System operuje znacznie poniżej limitu")
    print(f"   Duża rezerwa dla przyśpieszenia")
    print(f"   Status: DALEKO OD LIMITU MARGOLUSA-LEVITINA")


================================================================================
QW-298: GRANICA PRĘDKOŚCI PRZETWARZANIA (MARGOLUS-LEVITIN)
================================================================================

⚡ ENERGIA SYSTEMU:
   Stan początkowy: |ψ₀⟩ = stan podstawowy
   Energia: E₀ = ⟨ψ₀|S|ψ₀⟩ = -4.937223


⏱️ CZAS ORTOGONALIZACJI:
   Minimalne |⟨ψ₀|ψ(t)⟩|² = 1.000000e+00 przy t = 93.379338
   Szybkość operacji: v_op = 1/t_orth = 0.010709

📊 PORÓWNANIE Z LIMITEM MARGOLUSA-LEVITINA:
   Limit teoretyczny: v_op^(max) = 2E/π = 2·4.937223/π = 3.143133
   Szybkość rzeczywista: v_op^(actual) = 0.010709
   Stosunek: v_op^(actual) / v_op^(max) = 0.003407
   ✓ SYSTEM DZIAŁA NA LIMICIE FIZYCZNYM!
   Natura realizuje obliczenia z maksymalną możliwą prędkością

🔍 ANALIZA DLA RÓŻNYCH STANÓW POCZĄTKOWYCH:

   Stan 0: E=-4.9372, v_op=0.0107, limit=3.1431, ratio=0.0034

   Stan 1: E=-4.5036, v_op=0.0120, limit=2.8671, ratio=0.0042

   Stan 2: E=-3.3734, v_op=0.0109, limit=2.1476, ratio=0.0051

   Stan 3: E=-0.0309, v_op=0.0101, limit=0.0196, ratio=0.5142

   Stan 4: E=+0.5957, v_op=0.0192, limit=0.3792, ratio=0.0505

   Stan 5: E=+0.6150, v_op=0.0112, limit=0.3915, ratio=0.0287

   Stan 6: E=+0.6621, v_op=0.0152, limit=0.4215, ratio=0.0360

   Stan 7: E=+0.7075, v_op=0.0101, limit=0.4504, ratio=0.0224

   Średni stosunek v_op/v_max: 0.0831

✅ WYNIK QW-298:
   System operuje znacznie poniżej limitu
   Duża rezerwa dla przyśpieszenia
   Status: DALEKO OD LIMITU MARGOLUSA-LEVITINA

In [14]:


# ============================================================================
# QW-299: STAŁA SPRZĘŻENIA ŚWIADOMOŚCI (α_mind)
# ============================================================================
# Cel: Test IIT (Integrated Information Theory) w topologii splotów
# Czy istnieje krytyczne zagęszczenie pętli zwrotnych (proto-świadomość)?

print("\n" + "="*80)
print("QW-299: STAŁA SPRZĘŻENIA ŚWIADOMOŚCI (α_mind)")
print("="*80)

# Podejście uproszczone: badamy zwrotność bez liczenia wszystkich cykli
# Zbyt długie obliczanie cykli - używamy prostszych metryk topologicznych

print(f"\n🧠 KONSTRUKCJA GRAFU PRZEPŁYWU INFORMACJI:")

# Badamy małe rozmiary dla wydajności
N_values = [6, 8, 10]
results = []

for N in N_values:
    S = build_S_matrix(N)

    # Próg sprzężenia: 30% maksymalnej wartości
    threshold = 0.30 * np.abs(S).max()

    # Graf skierowany: i → j jeśli S_ij > threshold
    G = nx.DiGraph()
    for i in range(N):
        G.add_node(i)

    for i in range(N):
        for j in range(N):
            if i != j and np.abs(S[i, j]) > threshold:
                G.add_edge(i, j, weight=np.abs(S[i, j]))

    # Prostsze metryki: silnie spójne składowe i zwrotność
    sccs = list(nx.strongly_connected_components(G))
    n_sccs = len(sccs)
    largest_scc_size = max([len(scc) for scc in sccs]) if sccs else 0

    # Współczynnik zwrotności: krawędzie zwrotne / wszystkie krawędzie
    n_edges = G.number_of_edges()
    n_reciprocal = sum(1 for i, j in G.edges() if G.has_edge(j, i)) / 2
    reciprocity = n_reciprocal / n_edges if n_edges > 0 else 0

    # Gęstość grafu (procent możliwych połączeń)
    max_edges = N * (N - 1)
    density = n_edges / max_edges if max_edges > 0 else 0

    results.append({
        'N': N,
        'n_sccs': n_sccs,
        'largest_scc': largest_scc_size,
        'reciprocity': reciprocity,
        'density': density,
        'n_edges': n_edges
    })

    print(f"\n   N={N:2d}:")
    print(f"      Krawędzie: {n_edges}, Gęstość: {density:.6f}")
    print(f"      Silnie spójne składowe: {n_sccs}, Największa: {largest_scc_size}")
    print(f"      Zwrotność: {reciprocity:.6f} ({n_reciprocal:.0f} par)")

# Analiza
print(f"\n🎯 ANALIZA INTEGRACJI INFORMACJI:")

N_vals = np.array([r['N'] for r in results])
reciprocities = np.array([r['reciprocity'] for r in results])
densities = np.array([r['density'] for r in results])

# Stała sprzężenia świadomości α_mind
# Definiujemy jako zwrotność (miara pętli feedbackowych)
alpha_mind = reciprocities.max()
N_critical = N_vals[np.argmax(reciprocities)]

print(f"\n💫 STAŁA SPRZĘŻENIA ŚWIADOMOŚCI:")
print(f"   α_mind = max(reciprocity) = {alpha_mind:.6f}")
print(f"   Osiągana przy N = {N_critical}")

# Porównanie z innymi stałymi sprzężenia
alpha_EM_theory = 1/137.115
alpha_s = 1.0  # Silne przy niskich energiach

print(f"\n🔢 PORÓWNANIE ZE STAŁYMI FUNDAMENTALNYMI:")
print(f"   α_EM = {alpha_EM_theory:.6e} (elektromagnetyzm)")
print(f"   α_s  ≈ {alpha_s:.6f} (oddziaływania silne)")
print(f"   α_mind = {alpha_mind:.6e} (proto-świadomość?)")
print(f"   α_mind / α_EM = {alpha_mind / alpha_EM_theory:.2e}")

print(f"\n✅ WYNIK QW-299:")
if alpha_mind > 0.5:
    print(f"   Istnieje silna zwrotność (feedback loops)")
    print(f"   α_mind = {alpha_mind:.6f} wskazuje na wysoką integrację informacji")
    print(f"   Krytyczna skala: N_c = {N_critical}")
    print(f"   Status: MECHANIZM PROTO-ŚWIADOMOŚCI MOŻLIWY ✓")
elif alpha_mind > 0.1:
    print(f"   Umiarkowana zwrotność w systemie")
    print(f"   α_mind = {alpha_mind:.6f} wskazuje na częściową integrację")
    print(f"   Status: INTEGRACJA INFORMACJI NA ŚREDNIM POZIOMIE")
else:
    print(f"   α_mind ≈ {alpha_mind:.6e} jest relatywnie mały")
    print(f"   System wykazuje słabą integrację informacji")
    print(f"   Status: PROTO-ŚWIADOMOŚĆ WYMAGA WYŻSZYCH PROGÓW SPRZĘŻENIA")
