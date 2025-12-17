#!/usr/bin/env python3
"""
QW-1513: FALE GRAWITACYJNE JAKO FALE TORSYJNE
==============================================
Na podstawie analizy QW-1512, QW-1202, QW-1214.

KLUCZOWY WNIOSEK:
Fale GW to nie fonony (drgania pozycji) ale TORSJE (drgania kąta/fazy).
Analogia: w krysztale mamy fonony (dźwięk) i magnony (fale spinowe).
Fale GW są odpowiednikiem "magnonów topologicznych".

MECHANIZM:
- Defekt topologiczny ma fazę / kąt nawiniecia
- Oscylacja defektu = oscylacja fazy
- Fala fazy propaguje się przez sieć
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.fft import fft, fftfreq
import datetime

print("="*80)
print("QW-1513: FALE GRAWITACYJNE JAKO FALE TORSYJNE")
print("="*80)

# Parametry z najnowszych badań (QW-1202, QW-1214)
ALPHA_GEO = 4 * np.log(2)  # 2.7726
OMEGA_FUND = np.pi / 4     # Fundamentalna częstość z QW-1214
BETA_TORS = 0.01           # Torsja z hierarchii gauge
N_OCTAVES = 12

def K(d):
    """Kernel K(d) z parametrami frozen"""
    if d < 0.1:
        d = 0.1
    return ALPHA_GEO * np.cos(OMEGA_FUND * d + np.pi/6) / (1 + BETA_TORS * d)

def run_torsion_wave_simulation():
    """Symulacja fal torsyjnych (nie fononowych)"""
    
    print("\n[1] KONFIGURACJA MODELU TORSYJNEGO")
    print("-" * 60)
    
    # Sieć 1D z węzłami reprezentującymi FAZĘ (nie pozycję)
    N = 200  # Liczba węzłów
    L = 50.0  # Długość "łańcucha"
    dx = L / N
    
    # Każdy węzeł ma fazę θ[i] (zamiast pozycji x[i])
    # Równanie ruchu dla fazy: d²θ/dt² = c² × (θ[i+1] - 2θ[i] + θ[i-1]) / dx²
    
    # Prędkość fal torsyjnych:
    # c_tors wynika z K(d) - sztywność sieci
    K_0 = K(0)
    c_tors = np.sqrt(abs(K_0))
    
    print(f"Prędkość fal torsyjnych: c_tors = {c_tors:.4f}")
    print(f"Liczba węzłów: N = {N}")
    print(f"Rozdzielczość: dx = {dx:.4f}")
    
    # Źródło: dwa orbitujące defekty (symulacja mergera)
    source_pos = N // 4  # Pozycja źródła
    source_freq = 0.05   # Częstotliwość orbitalna (Hz)
    source_amplitude = 0.5  # Amplituda oscylacji fazy
    
    # Fale GW mają 2x częstotliwość orbitalną!
    gw_freq = 2 * source_freq
    
    print(f"\nŹródło (merger):")
    print(f"  Pozycja: węzeł {source_pos}")
    print(f"  Częstotliwość orbitalna: f_orbit = {source_freq}")
    print(f"  Częstotliwość fal GW: f_gw = 2 × f_orbit = {gw_freq}")
    
    # Detektory na różnych odległościach
    detector_positions = [N//2, 3*N//4, N-10]
    print(f"  Detektory na węzłach: {detector_positions}")
    
    print("\n[2] RÓWNANIE RUCH DLA FAL TORSYJNYCH")
    print("-" * 60)
    print("d²θ/dt² = c² × ∇²θ - γ × dθ/dt + Source(t)")
    print(f"c = {c_tors:.4f}, γ = {BETA_TORS}")
    
    # Symulacja
    def torsion_dynamics(t, state):
        """Równanie ruchu dla pola fazy θ."""
        theta = state[:N]
        dtheta_dt = state[N:]
        
        # Laplacian (dyskretny)
        d2theta = np.zeros(N)
        for i in range(1, N-1):
            d2theta[i] = (c_tors**2 / dx**2) * (theta[i+1] - 2*theta[i] + theta[i-1])
        
        # Warunki brzegowe (otwarte / absorpcyjne)
        d2theta[0] = 0
        d2theta[N-1] = 0
        
        # Źródło: oscylująca faza w miejscu źródła (chirp)
        # Dla mergera: amplituda rośnie, częstotliwość rośnie (chirp)
        # Uproszczenie: stała częstotliwość, stała amplituda
        d2theta[source_pos] += source_amplitude * np.sin(2 * np.pi * gw_freq * t)
        
        # Tłumienie torsyjne (z beta_tors)
        d2theta -= BETA_TORS * dtheta_dt
        
        return np.concatenate([dtheta_dt, d2theta])
    
    # Parametry czasowe
    t_max = 200.0
    dt = 0.5
    t_eval = np.arange(0, t_max, dt)
    
    print(f"\nCzas symulacji: {t_max} jednostek")
    
    # Warunki początkowe: wszystko w spoczynku (θ = 0, dθ/dt = 0)
    initial_state = np.zeros(2 * N)
    
    print("\n[3] INTEGRACJA RÓWNAŃ RUCHU")
    print("-" * 60)
    
    sol = solve_ivp(torsion_dynamics, [0, t_max], initial_state,
                    method='RK45', t_eval=t_eval, rtol=1e-6)
    
    theta_history = sol.y[:N, :]  # Historia pola fazy
    times = sol.t
    
    print(f"Integracja zakończona: {len(times)} punktów czasowych")
    
    print("\n[4] DETEKCJA FAL GRAWITACYJNYCH")
    print("-" * 60)
    
    print(f"Oczekiwana częstotliwość: f_gw = {gw_freq}")
    print("")
    
    wave_detected = False
    results = []
    
    for det_pos in detector_positions:
        # Sygnał na detektorze
        signal = theta_history[det_pos, :]
        signal_ac = signal - np.mean(signal)
        
        # FFT
        freqs = fftfreq(len(times), dt)
        spectrum = np.abs(fft(signal_ac))**2
        
        # Peak w dodatnich częstotliwościach
        pos_mask = freqs > 0
        if np.sum(pos_mask) > 0:
            peak_idx = np.argmax(spectrum[pos_mask])
            peak_freq = freqs[pos_mask][peak_idx]
            peak_power = spectrum[pos_mask][peak_idx]
            amplitude = np.std(signal_ac)
        else:
            peak_freq = 0
            peak_power = 0
            amplitude = 0
        
        # Sprawdź zgodność
        freq_match = abs(peak_freq - gw_freq) < 0.02
        
        distance = (det_pos - source_pos) * dx
        status = "✅ WYKRYTO" if freq_match else "❌ NIE WYKRYTO"
        
        print(f"Detektor r={distance:.1f}: f_peak={peak_freq:.4f}, A={amplitude:.6f} → {status}")
        
        if freq_match:
            wave_detected = True
        
        results.append({
            "position": det_pos,
            "distance": distance,
            "peak_freq": peak_freq,
            "amplitude": amplitude,
            "detected": freq_match
        })
    
    print("\n[5] WERDYKT")
    print("=" * 60)
    
    if wave_detected:
        print("✅ FALE GRAWITACYJNE (TORSYJNE) WYKRYTE!")
        print(f"   Mechanizm: Fale fazy θ propagują się przez sieć")
        print(f"   Prędkość: c_tors = {c_tors:.4f}")
        verdict = "SUKCES"
    else:
        print("❌ FALE NIE WYKRYTE")
        print("   Możliwe przyczyny: za silne tłumienie, za słabe źródło")
        verdict = "PORAŻKA"
    
    # Raport
    report = f"""# QW-1513: Fale Grawitacyjne jako Fale Torsyjne

**Data:** {datetime.datetime.now()}
**Status:** {verdict}

## Metodologia
Na podstawie analizy QW-1512, QW-1202, QW-1214:
- Fale GW to fale **fazy/torsji**, nie fale pozycji (fonony)
- Równanie: d²θ/dt² = c² ∇²θ - γ dθ/dt

## Parametry
- Prędkość fal: c_tors = {c_tors:.4f}
- Tłumienie: γ = {BETA_TORS}
- Częstotliwość źródła: f_gw = {gw_freq}

## Wyniki Detekcji
| Odległość | f_peak | Amplituda | Wykryto? |
|-----------|--------|-----------|----------|
"""
    for r in results:
        status = "✅" if r["detected"] else "❌"
        report += f"| {r['distance']:.1f} | {r['peak_freq']:.4f} | {r['amplitude']:.6f} | {status} |\n"
    
    report += f"""
## Wnioski
{"Fale torsyjne propagują się przez sieć Nadsolitonu z prędkością c_tors." if wave_detected else "Model wymaga dalszej analizy."}
"""
    
    with open("QW-1513_Torsion_Waves.md", "w") as f:
        f.write(report)
    
    print(f"\n[SAVED] QW-1513_Torsion_Waves.md")

if __name__ == "__main__":
    run_torsion_wave_simulation()
