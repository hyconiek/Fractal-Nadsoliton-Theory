#!/usr/bin/env python3
import json
import numpy as np
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
    print("Rozpoczynanie zoptymalizowanej rygorystycznej weryfikacji mostu fraktalnego...")
    
    d_vals = np.arange(1, 12, dtype=float)  # kształt (11,)
    y_str = K_strict_gate(d_vals)          # kształt (11,)
    
    omega_leg = np.pi / 4.0
    phi_leg = np.pi / 6.0
    beta_tors = 0.01
    
    # Generowanie siatek o mniejszym kroku ale zwektoryzowane
    # Zmniejszmy siatkę i zróbmy pełną wektoryzację
    lambdas = np.linspace(0.01, 3.0, 100) # (100,)
    alphas = np.linspace(0.05, 3.0, 100)   # (100,)
    phases = np.linspace(-np.pi, np.pi, 100) # (100,)
    
    # Rozszerzamy do kształtów dla broadcasting:
    # d_vals: (1, 1, 1, 11)
    # lambdas: (100, 1, 1, 1)
    # alphas: (1, 100, 1, 1)
    # phases: (1, 1, 100, 1)
    
    L = lambdas[:, np.newaxis, np.newaxis, np.newaxis]
    A = alphas[np.newaxis, :, np.newaxis, np.newaxis]
    P = phases[np.newaxis, np.newaxis, :, np.newaxis]
    D = d_vals[np.newaxis, np.newaxis, np.newaxis, :]
    
    # d_eff shape: (100, 100, 1, 11)
    d_eff = L * (D ** A)
    
    # y_model shape: (100, 100, 100, 11)
    y_model = np.cos(omega_leg * d_eff + phi_leg + P) / (1.0 + beta_tors * d_eff)
    
    # Błąd: suma kwadratów różnic po ostatniej osi (d)
    # err shape: (100, 100, 100)
    err = np.sum((y_model - y_str[np.newaxis, np.newaxis, np.newaxis, :]) ** 2, axis=3)
    
    # Znajdujemy indeks minimalnego błędu
    min_idx = np.unravel_index(np.argmin(err), err.shape)
    
    l_opt = lambdas[min_idx[0]]
    a_opt = alphas[min_idx[1]]
    p_opt = phases[min_idx[2]]
    ss_res = err[min_idx]
    
    d_eff_opt = l_opt * (d_vals ** a_opt)
    y_fitted = np.cos(omega_leg * d_eff_opt + phi_leg + p_opt) / (1.0 + beta_tors * d_eff_opt)
    
    # Obliczmy współczynnik R^2 (determinacji)
    ss_tot = np.sum((y_str - np.mean(y_str)) ** 2)
    r2 = 1.0 - (ss_res / ss_tot)
    
    report = {
        "title": "Rygorystyczny Raport Weryfikacyjny Mostu Fraktalnego (Zoptymalizowany)",
        "status": "ZWERYFIKOWANO" if r2 > 0.95 else ("KORELACJA_DOBRA" if r2 > 0.8 else "SLABA_KORELACJA"),
        "r2_score": r2,
        "best_fit_parameters": {
            "scaling_lambda": l_opt,
            "fractal_dimension_scaling_alpha": a_opt,
            "phase_shift": p_opt,
            "residual_sum_of_squares": ss_res
        },
        "points": []
    }
    
    for i, d in enumerate(d_vals):
        report["points"].append({
            "d": d,
            "d_eff": d_eff_opt[i],
            "strict_actual": y_str[i],
            "model_fitted": y_fitted[i],
            "error": y_str[i] - y_fitted[i]
        })
        
    out_path = Path("scratch/bridge_verification_report.json")
    out_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Zoptymalizowana weryfikacja zakończona. R^2 score: {r2:.4f}")
    print(f"Raport zapisany w {out_path}")

if __name__ == "__main__":
    main()
