"""
Phase 51: SNR Scan — At what signal strength does H=0.23 background
produce Cross-H ≈ 0.31 in damped channels?

If there exists an SNR where Cross-H passes through 0.31,
then the measured 0.31 COULD be a projection of true 0.23.
If Cross-H jumps directly from ~0.50 to ~0.23 without passing through 0.31,
then 0.31 is NOT a projection of 0.23.
"""

import numpy as np
import json
import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

def generate_fgn(n, H, seed=None):
    if seed is not None:
        np.random.seed(seed)
    N = n
    f = np.fft.fftfreq(N)
    f[0] = 1e-10
    power = np.abs(f) ** (-(2*H - 1))
    power[0] = 0
    phases = np.random.uniform(0, 2*np.pi, N)
    spectrum = np.sqrt(power) * np.exp(1j * phases)
    if N % 2 == 0:
        spectrum[N//2] = np.abs(spectrum[N//2])
    for i in range(1, N//2):
        spectrum[N-i] = np.conj(spectrum[i])
    spectrum[0] = 0
    x = np.real(np.fft.ifft(spectrum))
    x = (x - np.mean(x)) / (np.std(x) + 1e-30)
    return x

def cross_mfdfa_q0(x, y, scales):
    z = np.cumsum((x - np.mean(x)) * (y - np.mean(y)))
    F, valid_scales = [], []
    N = len(z)
    for s in scales:
        n = N // s
        if n == 0: continue
        rms = []
        for i in range(n):
            seg = z[i*s:(i+1)*s]
            p = np.polyfit(np.arange(s), seg, 1)
            var = np.mean((seg - np.polyval(p, np.arange(s)))**2)
            if var > 0: rms.append(var)
        rms = np.array(rms)
        if len(rms) > 0:
            F.append(np.exp(0.5 * np.mean(np.log(rms + 1e-300))))
            valid_scales.append(s)
    if len(F) < 3: return np.nan
    return float(np.polyfit(np.log(valid_scales), np.log(F), 1)[0])

def main():
    N = 2**19  # 512K samples (faster, still statistically valid)
    N_TRIALS = 10
    H_BG = 0.23
    H_LOCAL = 0.04

    # SNR scan: ratio of local_noise_amplitude / background_amplitude
    # SNR=0 means pure background (no local noise)
    # SNR=inf means pure local noise (no background)
    snr_values = [0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 50.0]

    scales = np.logspace(3, np.log10(N//4), 15).astype(int)
    scales = np.unique(scales)

    results = {}

    for snr in snr_values:
        log.info(f"--- SNR = {snr} (local/background amplitude ratio) ---")
        cross_vals = []

        for trial in range(N_TRIALS):
            bg = generate_fgn(N, H_BG, seed=trial*1000)
            n1 = generate_fgn(N, H_LOCAL, seed=trial*1000 + 1)
            n2 = generate_fgn(N, H_LOCAL, seed=trial*1000 + 2)

            x = bg + snr * n1
            y = bg + snr * n2

            hc = cross_mfdfa_q0(x, y, scales)
            cross_vals.append(hc)

        mean_h = float(np.mean(cross_vals))
        std_h = float(np.std(cross_vals))
        log.info(f"  Cross-H = {mean_h:.4f} ± {std_h:.4f}")

        results[str(snr)] = {
            "Cross_H_mean": mean_h,
            "Cross_H_std": std_h
        }

    with open("QW_1660_v51_SNR_Scan.json", "w") as f:
        json.dump(results, f, indent=2)

    print("\n--- SNR SCAN RESULTS ---")
    print(f"{'SNR':>6} | {'Cross-H':>12} | {'Std':>8}")
    print("-" * 35)
    for snr in snr_values:
        r = results[str(snr)]
        marker = " <-- 0.31?" if abs(r["Cross_H_mean"] - 0.31) < 0.02 else ""
        print(f"{snr:6.1f} | {r['Cross_H_mean']:12.4f} | {r['Cross_H_std']:8.4f}{marker}")

if __name__ == "__main__":
    main()
