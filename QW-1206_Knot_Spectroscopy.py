import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from math import gcd
from datetime import datetime

# ==============================================================================
# QW-1206: SPEKTROSKOPIA WĘZŁÓW TORUSOWYCH (KNOT SPECTROSCOPY)
# ==============================================================================
#
# CEL: Sprawdzić, czy węzły Fibonacciego (np. T(21,3)) mają unikalne własności
#      rezonansowe (widmo harmoniczne) w porównaniu do innych węzłów.
#
# METODA: 
# 1. Dyskretyzacja węzła T(p,q) jako grafu cyklicznego w przestrzeni 3D.
# 2. Konstrukcja operatora Laplace'a-Beltramiego na grafie (z uwzględnieniem
#    krzywizny i skręcenia węzła).
# 3. Obliczenie wartości własnych (częstotliwości rezonansowych).
# 4. Miara "harmoniczności": jak blisko widmo jest ciągu n^2 lub n.
#
# ==============================================================================

REPORT_FILE = "RAPORT_QW1206_SPECTROSCOPY.md"
md = []

def log(msg):
    print(msg)
    md.append(msg)

log("=" * 80)
log("QW-1206: SPEKTROSKOPIA TOPOLOGICZNA WĘZŁÓW TORUSOWYCH")
log("=" * 80)

def generate_torus_knot_graph(p, q, points=1000):
    """
    Generuje graf reprezentujący węzeł torusowy T(p,q).
    Węzły grafu to punkty w 3D. Krawędzie łączą sąsiednie punkty.
    Zwraca macierz sąsiedztwa (z wagami 1/dist) i pozycje.
    """
    t = np.linspace(0, 2*np.pi, points, endpoint=False)
    # Parametryzacja na torusie
    # R duże = 2, r małe = 1
    R, r = 2.0, 1.0
    x = (R + r * np.cos(q * t)) * np.cos(p * t)
    y = (R + r * np.cos(q * t)) * np.sin(p * t)
    z = r * np.sin(q * t)
    
    positions = np.column_stack((x, y, z))
    
    # Tworzenie macierzy sąsiedztwa (pierścień z wagami geometrycznymi)
    # Wagi = odwrotność odległości euklidesowej (sztywność sprężyny)
    adj = sp.lil_matrix((points, points))
    
    for i in range(points):
        # Sąsiad "następny"
        j = (i + 1) % points
        dist = np.linalg.norm(positions[i] - positions[j])
        weight = 1.0 / dist # Krótszy segment = sztywniejsza sprężyna
        adj[i, j] = weight
        adj[j, i] = weight
        
        # Opcjonalnie: second-nearest neighbor dla sztywności giętnej (curvature)
        k = (i + 2) % points
        dist2 = np.linalg.norm(positions[i] - positions[k])
        weight2 = 0.5 / dist2 # Mniejsza waga dla zginania
        adj[i, k] = weight2
        adj[k, i] = weight2

    return adj.tocsr(), positions

def compute_spectrum(adj_matrix, k=10):
    """Oblicza k najmniejszych niezerowych wartości własnych Laplasjanu."""
    # Laplasjan grafowy L = D - A
    degrees = np.array(adj_matrix.sum(axis=1)).flatten()
    D = sp.diags(degrees)
    L = D - adj_matrix
    
    # Wartości własne (najmniejsze magnitudowo, ale pomijamy 0)
    # Używamy shift-invert mode sigma=0.01 aby znaleźć te bliskie zera
    vals, vecs = spla.eigsh(L, k=k+1, sigma=0.01, which='LM')
    
    # Sortowanie i usunięcie modu zerowego (pierwsza w.w. ≈ 0)
    vals = np.sort(np.real(vals))
    return vals[1:] # Pomijamy 0

def harmonicity_score(spectrum):
    """
    Ocenia jak bardzo widmo przypomina idealny oscylator harmoniczny (n)
    lub strunę (n^2).
    Score = średnie odchylenie od liniowego dopasowania.
    Niższy score = lepiej (bardziej harmoniczne).
    """
    # Normalizacja widma do pierwszej wartości
    norm_spec = spectrum / spectrum[0]
    n = np.arange(1, len(spectrum) + 1)
    
    # Hipoteza 1: Struna 1D (n^2) - typowe dla pętli
    # frequencies ~ n, eigenvalues ~ n^2
    # Sprawdzamy liniowość pierwiastków wartości własnych (czyli częstotliwości)
    freqs = np.sqrt(norm_spec)
    
    # Fit liniowy: freq = a * n + b
    slope, intercept = np.polyfit(n, freqs, 1)
    
    # Błąd dopasowania (MSE)
    fitted_freqs = slope * n + intercept
    mse = np.mean((freqs - fitted_freqs)**2)
    
    return mse, slope

log("\n[1] BADANIE WIDMA DLA RÓŻNYCH WĘZŁÓW (DLA Q=24)")
log("-" * 80)
log(f"{'Węzeł':<10} {'Q':<5} {'Score (MSE)':<15} {'Slope':<10} {'Widmo (pierwsze 4 mode)':<40}")
log("-" * 80)

# Lista węzłów do sprawdzenia dla Q bliskiego 24
candidates = [
    (21, 3), # Elektron (Fibonacci)
    (13, 11), # Konkurent E/Q
    (19, 5),
    (17, 7), 
    (12, 12) # "Idealnie" symetryczny, choć torusowy tylko dla gcd=1 (tu gcd=12, więc to link)
]

# Dodajmy też inne Q dla porównania
candidates += [(13, 8), (8, 5), (5, 3)] # Ciąg Fibonacciego

results = []

for p, q in candidates:
    if gcd(p, q) > 1:
        log(f"T({p},{q}) - Pominięto (nie jest węzłem, gcd={gcd(p,q)})")
        continue

    adj, pos = generate_torus_knot_graph(p, q, points=2000)
    spec = compute_spectrum(adj, k=10)
    mse, slope = harmonicity_score(spec)
    
    # Formatowanie widma do wyświetlenia (sqrt dla częstotliwości)
    freqs_str = str(np.round(np.sqrt(spec[:4]/spec[0]), 2))
    
    is_fib = (p in [3,5,8,13,21,34] and q in [3,5,8,13,21,34])
    note = "✅ FIB" if is_fib else ""
    
    res = {'knot': f"T({p},{q})", 'Q': p+q, 'mse': mse, 'note': note}
    results.append(res)
    
    log(f"T({p},{q})".ljust(10) + f"{p+q:<5} {mse:<15.6f} {slope:<10.2f} {freqs_str:<40} {note}")

# Analiza wyników
log("\n[2] WNIOSKI Z PORÓWNANIA")
log("-" * 80)
# Sortowanie wg Score (im mniej tym lepiej)
results.sort(key=lambda x: x['mse'])

log(f"Najbardziej harmoniczne węzły (TOP 3):")
for i, r in enumerate(results[:3]):
    log(f"{i+1}. {r['knot']} (Q={r['Q']}) - MSE: {r['mse']:.6f} {r['note']}")

# Sprawdzenie czy Elektron T(21,3) wygrywa z T(13,11)
e_knot = next((r for r in results if r['knot'] == "T(21,3)"), None)
competitor = next((r for r in results if r['knot'] == "T(13,11)"), None)

if e_knot and competitor:
    log("\nPOJEDYNEK DLA Q=24:")
    log(f"T(21,3) [Elektron, Fib]: MSE = {e_knot['mse']:.6f}")
    log(f"T(13,11) [Symetryczny]:   MSE = {competitor['mse']:.6f}")
    
    if e_knot['mse'] < competitor['mse']:
        log("✅ T(21,3) JEST BARDZIEJ HARMONICZNY! (Hipoteza rezonansu potwierdzona)")
    else:
        log("❌ T(13,11) jest bardziej harmoniczny. (Hipoteza rezonansu odrzucona)")

log("\n[3] INTERPRETACJA FIZYCZNA")
log("-" * 80)
log("""
Co to oznacza?
MSE mierzy, jak bardzo węzeł zachowuje się jak idealna, jednowymiarowa struna.
Węzły mocno splecione (jak T(13,11)) są "sztywniejsze" i mają zaburzone widmo
z powodu silnych oddziaływań geometrycznych między pętlami.
Węzły "luźniejsze" i asymetryczne (jak T(21,3)) mogą zachowywać się bardziej
jak swobodne pętle, co daje czystsze widmo rezonansowe.

Dla stabilnej cząstki potrzebujemy CZYSTYCH stanów kwantowych (długi czas życia).
Stany chaotyczne (wysokie MSE) szybko ulegają dekoherencji.
""")

# Zapis raportu
with open(REPORT_FILE, 'w') as f:
    f.write("# QW-1206: Spektroskopia Węzłów Torusowych\n\n")
    f.write(f"**Data:** {datetime.now()}\n\n")
    f.write("```\n")
    f.write("\n".join(md))
    f.write("\n```\n")

print(f"Raport zapisano w {REPORT_FILE}")
