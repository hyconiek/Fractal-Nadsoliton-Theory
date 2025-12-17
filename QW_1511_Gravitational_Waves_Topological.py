#!/usr/bin/env python3
"""
QW-1511: FALE GRAWITACYJNE W MODELU DEFEKTÓW TOPOLOGICZNYCH
============================================================
Bazując na QW-722 (statyczna grawitacja z n=-2.26), sprawdzamy
czy oscylujące defekty mogą generować fale grawitacyjne.

Metodologia z QW-722:
- Masy = defekty topologiczne (winding numbers)
- Siła F ∝ -K_eff × m1 × m2 / r²
- K_eff = K(d) zmodyfikowane przez defekty

Nowy element:
- Oscylujący defekt (symulacja mergera dwóch mas)
- Mierzymy propagację perturbacji do dalszych obserwatorów
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from scipy.fft import fft, fftfreq
import datetime

print("="*80)
print("QW-1511: FALE GRAWITACYJNE W MODELU DEFEKTÓW TOPOLOGICZNYCH")
print("="*80)

# Parametry z QW-722
ALPHA = 4 * np.log(2)  # 2.7726
OMEGA = np.pi / 4
PHI = np.pi / 6
BETA = 0.1  # Z QW-722 (nie 0.01!)
N_OCTAVES = 12

def K(d):
    """K(d) kernel z QW-722"""
    if d < 0.1:
        d = 0.1
    return ALPHA * np.cos(OMEGA * d + PHI) / (1 + BETA * d)

class TopologicalDefect:
    """Defekt topologiczny z QW-722"""
    def __init__(self, position, octave, winding_number=1):
        self.position = np.array(position, dtype=float)
        self.octave = octave
        self.winding = winding_number

def calculate_force(defect1, defect2):
    """Oblicza siłę między defektami (z QW-722)"""
    r = np.linalg.norm(defect1.position - defect2.position)
    d_octave = abs(defect1.octave - defect2.octave)
    
    # Mieszanie skali (z QW-722)
    d_eff = d_octave + 0.1 * r
    K_eff = K(d_eff)
    K_eff *= (1.0 + 0.5 * defect1.winding)
    K_eff *= (1.0 + 0.5 * defect2.winding)
    
    m1, m2 = defect1.winding, defect2.winding
    F = -K_eff * m1 * m2 / (r**2 + 0.1)
    return F

def run_gravitational_wave_test():
    """Główna symulacja fal grawitacyjnych"""
    
    print("\n[1] KONFIGURACJA SYMULACJI")
    print("-" * 60)
    
    # Oscylujące źródło (symulacja mergera)
    # Dwa defekty orbitujące wokół siebie
    source_orbit_radius = 2.0
    source_freq = 0.1  # Hz - częstotliwość orbitalna
    source_octave = 3
    source_winding = 2  # Ciężkie masy (jak czarne dziury)
    
    print(f"Źródło: 2 defekty orbitujące")
    print(f"  Promień orbity: {source_orbit_radius}")
    print(f"  Częstotliwość: f = {source_freq} Hz")
    print(f"  Winding number: {source_winding}")
    
    # Obserwatorzy (detektory) na różnych odległościach
    observer_distances = [5.0, 10.0, 15.0, 20.0, 30.0]
    observer_octave = 3  # Ta sama oktawa co źródło
    
    print(f"\nObserwatorzy na odległościach: {observer_distances}")
    
    # Symulacja czasowa
    t_max = 100.0
    dt = 0.5
    times = np.arange(0, t_max, dt)
    n_steps = len(times)
    
    print(f"Czas symulacji: {t_max} jednostek, dt = {dt}")
    
    print("\n[2] SYMULACJA DYNAMIKI")
    print("-" * 60)
    
    # Historia sił na obserwatorach
    force_history = {d: [] for d in observer_distances}
    
    for t in times:
        # Pozycje dwóch defektów źródłowych (orbita)
        angle = 2 * np.pi * source_freq * t
        pos1 = np.array([source_orbit_radius * np.cos(angle), 
                        source_orbit_radius * np.sin(angle), 0])
        pos2 = -pos1  # Po przeciwnej stronie orbity
        
        defect1 = TopologicalDefect(pos1, source_octave, source_winding)
        defect2 = TopologicalDefect(pos2, source_octave, source_winding)
        
        # Oblicz siłę na każdym obserwatorze
        for obs_distance in observer_distances:
            obs_pos = np.array([obs_distance, 0, 0])
            observer = TopologicalDefect(obs_pos, observer_octave, 1)
            
            # Siła od obu defektów źródłowych
            F1 = calculate_force(defect1, observer)
            F2 = calculate_force(defect2, observer)
            
            # Całkowita siła (składowa x - w kierunku radialnym)
            F_total = F1 + F2
            
            force_history[obs_distance].append(F_total)
    
    # Konwertuj na tablice
    for d in observer_distances:
        force_history[d] = np.array(force_history[d])
    
    print("Symulacja zakończona.")
    
    print("\n[3] ANALIZA FAL GRAWITACYJNYCH")
    print("-" * 60)
    
    # Dla każdego obserwatora, sprawdź czy widzi sygnał źródła
    wave_detected = False
    detected_freqs = []
    detected_amplitudes = []
    
    expected_freq = 2 * source_freq  # Fale GW mają 2x częstotliwość orbitalną
    
    print(f"Oczekiwana częstotliwość fal: f_gw = 2 × f_orbit = {expected_freq} Hz")
    print("")
    
    for obs_distance in observer_distances:
        signal = force_history[obs_distance]
        
        # Usuń DC
        signal_ac = signal - np.mean(signal)
        
        # FFT
        freqs = fftfreq(len(times), dt)
        spectrum = np.abs(fft(signal_ac))**2
        
        # Znajdź peak (tylko dodatnie częstotliwości)
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
        
        # Sprawdź zgodność z oczekiwaną częstotliwością
        freq_match = abs(peak_freq - expected_freq) < 0.05
        
        status = "✅ WYKRYTO" if freq_match else "❌ NIE WYKRYTO"
        print(f"r = {obs_distance:5.1f}: f_peak = {peak_freq:.4f}, A = {amplitude:.6f} → {status}")
        
        if freq_match:
            wave_detected = True
            detected_freqs.append(peak_freq)
            detected_amplitudes.append(amplitude)
    
    print("\n[4] ANALIZA AMPLITUDY vs ODLEGŁOŚĆ")
    print("-" * 60)
    
    if wave_detected and len(detected_amplitudes) > 2:
        # Fale grawitacyjne powinny maleć jak 1/r
        r_values = np.array([observer_distances[i] for i, f in enumerate(
            [force_history[d] for d in observer_distances]) 
            if abs(fftfreq(len(times), dt)[pos_mask][np.argmax(
                np.abs(fft(f - np.mean(f)))**2[pos_mask])] - expected_freq) < 0.05])
        
        amplitudes = np.array(detected_amplitudes)
        
        if len(r_values) >= 2 and len(amplitudes) >= 2:
            # Fit A = A0 / r^n
            def power_law(r, A0, n):
                return A0 * np.power(r, n)
            
            try:
                popt, _ = curve_fit(power_law, r_values[:len(amplitudes)], amplitudes, 
                                   p0=[1.0, -1.0], maxfev=5000)
                A0_fit, n_fit = popt
                
                print(f"Dopasowanie: A = {A0_fit:.4f} × r^({n_fit:.2f})")
                print(f"Oczekiwane dla fal GW: n = -1.0")
                print(f"Błąd: {abs(n_fit + 1.0):.2f}")
                
                if abs(n_fit + 1.0) < 0.3:
                    amplitude_result = "✅ ZGODNE z 1/r (fale GW)"
                else:
                    amplitude_result = f"🟡 n = {n_fit:.2f} ≠ -1.0"
            except:
                n_fit = 0
                amplitude_result = "❌ Nie udało się dopasować"
        else:
            n_fit = 0
            amplitude_result = "❌ Za mało punktów"
    else:
        n_fit = 0
        amplitude_result = "❌ Brak wykrytych fal"
    
    print(f"\n{amplitude_result}")
    
    print("\n[5] WERDYKT")
    print("=" * 60)
    
    if wave_detected:
        verdict = "✅ FALE GRAWITACYJNE WYKRYTE"
        conclusion = f"Oscylujące defekty topologiczne generują fale propagujące z f = {np.mean(detected_freqs):.3f} Hz"
    else:
        verdict = "❌ FALE GRAWITACYJNE NIE WYKRYTE"
        conclusion = "Model wymaga dalszych modyfikacji"
    
    print(f"\n{verdict}")
    print(f"{conclusion}")
    
    # Zapisz raport
    report_file = "QW-1511_Gravitational_Waves_Topological.md"
    with open(report_file, "w") as f:
        f.write("# QW-1511: Fale Grawitacyjne w Modelu Defektów Topologicznych\n\n")
        f.write(f"**Data:** {datetime.datetime.now()}\n\n")
        f.write("## 1. Metoda\n")
        f.write("Bazując na QW-722 (statyczna grawitacja n=-2.26), symulujemy\n")
        f.write("oscylujące defekty topologiczne (merger dwóch mas).\n\n")
        f.write("## 2. Wyniki Detekcji\n")
        f.write(f"| Odległość | f_peak | Amplituda | Status |\n")
        f.write("|-----------|--------|-----------|--------|\n")
        for obs_distance in observer_distances:
            signal = force_history[obs_distance]
            signal_ac = signal - np.mean(signal)
            freqs = fftfreq(len(times), dt)
            spectrum = np.abs(fft(signal_ac))**2
            pos_mask = freqs > 0
            if np.sum(pos_mask) > 0:
                peak_idx = np.argmax(spectrum[pos_mask])
                peak_freq = freqs[pos_mask][peak_idx]
                amplitude = np.std(signal_ac)
            else:
                peak_freq = 0
                amplitude = 0
            status = "✅" if abs(peak_freq - expected_freq) < 0.05 else "❌"
            f.write(f"| {obs_distance:.1f} | {peak_freq:.4f} | {amplitude:.6f} | {status} |\n")
        f.write(f"\n## 3. Werdykt\n")
        f.write(f"### {verdict}\n")
        f.write(f"{conclusion}\n")
    
    print(f"\n[SAVED] {report_file}")

if __name__ == "__main__":
    run_gravitational_wave_test()
