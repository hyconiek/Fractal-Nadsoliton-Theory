#!/usr/bin/env python3
import json
import numpy as np
from scipy.optimize import minimize
from pathlib import Path

# Definicje jądra
alpha_geo = 4.0 * np.log(2.0)  # ~2.7725887

def K_legacy_ont(d):
    omega = np.pi / 4.0
    phi = np.pi / 6.0
    beta_tors = 0.01
    return alpha_geo * np.cos(omega * d + phi) / (1.0 + beta_tors * d)

def K_strict_gate(d):
    omega = 0.18575
    phi = 0.16250
    beta = 1.0
    eta = 1.8
    return np.cos(omega * d + phi) / (1.0 + beta * (d ** eta))

def main():
    print("Rozpoczynanie badania mostu matematycznego...")
    
    # Punkty odległości d na pierścieniu Z_12
    d_vals = np.arange(1, 12, dtype=float)
    
    y_leg = K_legacy_ont(d_vals)
    y_str = K_strict_gate(d_vals)
    
    # 1. Sprawdźmy prostą hipotezę skróconego wymiaru / deformacji współrzędnej:
    # d_eff = d^gamma * scale
    # K_strict(d) = cos(omega_leg * d_eff + phi_leg) / (1 + beta_tors * d_eff) * (jakieś sprzężenie / absorpcja)
    
    # Zoptymalizujmy parametry transformacji współrzędnej d -> d_eff = scale * d^gamma
    # oraz ewentualnego czynnika fazowego lub tłumienia.
    def loss_func(params):
        scale, gamma, amp_scale, phase_shift = params
        d_eff = scale * (d_vals ** gamma)
        # Obliczamy modelowany K z parametrami legacy ale na zdeformowanej współrzędnej
        omega_leg = np.pi / 4.0
        phi_leg = np.pi / 6.0
        beta_tors = 0.01
        
        y_model = amp_scale * np.cos(omega_leg * d_eff + phi_leg + phase_shift) / (1.0 + beta_tors * d_eff)
        return np.sum((y_model - y_str) ** 2)
        
    initial_guess = [0.2, 1.0, 1.0/alpha_geo, 0.0]
    res = minimize(loss_func, initial_guess, method='L-BFGS-B')
    
    scale_opt, gamma_opt, amp_scale_opt, phase_shift_opt = res.x
    
    # Sprawdźmy dopasowanie
    d_eff_opt = scale_opt * (d_vals ** gamma_opt)
    omega_leg = np.pi / 4.0
    phi_leg = np.pi / 6.0
    beta_tors = 0.01
    y_fitted = amp_scale_opt * np.cos(omega_leg * d_eff_opt + phi_leg + phase_shift_opt) / (1.0 + beta_tors * d_eff_opt)
    residuals = np.sum((y_fitted - y_str) ** 2)
    
    # 2. Sprawdźmy hipotezę dodatkowego dephasingu (rozfazowania) jako nowej właściwości:
    # K_strict(d) = 1/alpha_geo * K_legacy(d) * dephasing_factor(d)
    # Wyznaczmy dephasing_factor(d) punktowo:
    dephasing_factors = []
    for d in d_vals:
        val_leg = K_legacy_ont(d)
        val_str = K_strict_gate(d)
        # Zapiszmy dephasing jako stosunek
        ratio = val_str / (val_leg / alpha_geo)
        dephasing_factors.append(ratio)
        
    # Spróbujmy dopasować dephasing_factor do fizycznego modelu nieliniowej dyspersji:
    # ratio(d) = cos(w_str * d + phi_str) / cos(w_leg * d + phi_leg) * (1 + beta_tors * d) / (1 + beta * d^eta)
    
    report = {
        "title": "Badanie Mostu Konforemnego i Deformacji Współrzędnych Jądra",
        "description": "Próba zrekonstruowania przejścia Legacy -> Strict za pomocą dodatkowego sprzężenia konforemnego / deformacji wymiarowej d na pierścieniu Z12.",
        "d_values": d_vals.tolist(),
        "legacy_values": y_leg.tolist(),
        "strict_values": y_str.tolist(),
        "hypothesis_1_coordinate_deformation": {
            "formula_d_eff": "d_eff = scale * d^gamma",
            "optimized_parameters": {
                "scale": scale_opt,
                "gamma": gamma_opt,
                "amplitude_absorption_factor": amp_scale_opt,
                "expected_amplitude_absorption": 1.0 / alpha_geo,
                "phase_shift": phase_shift_opt
            },
            "residual_sum_of_squares": residuals,
            "interpretation": "Pokazuje czy deformacja wymiarowa d (np. wymiar fraktalny sieci) wyjaśnia przejście częstotliwości i fazy."
        },
        "hypothesis_2_pointwise_dephasing_factors": {
            "d_values": d_vals.tolist(),
            "ratios_strict_to_normalized_legacy": dephasing_factors,
            "interpretation": "Stosunek wartości jądra ścisłego do znormalizowanego jądra legacy. Wskazuje na istnienie nieliniowego sprzężenia tłumiącego i przesunięcia fazowego, które nie były uwzględnione w uproszczonym modelu legacy."
        }
    }
    
    out_dir = Path("scratch")
    out_dir.mkdir(exist_ok=True)
    out_path = out_dir / "bridge_analysis_report.json"
    out_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Raport zapisany w {out_path}")

if __name__ == "__main__":
    main()
