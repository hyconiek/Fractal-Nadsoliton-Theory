import numpy as np
import matplotlib.pyplot as plt

def main():
    print("QW-586: Multi-Layer Fractal Gravity (MOND Fix)")
    print("----------------------------------------------")
    print("Hypothesis: Gravity filters through fractal layers")
    print("           Different layers dominate at different scales")
    
    # Parameters from odkrycie_fraktalne_warstwy.md
    BETA_TORS = 0.01  # Layer damping
    N_LAYERS = 20     # Fractal depth (Planck to Macro)
    
    # From QW-480: G(N) = G_0 * beta^N
    # At our scale (N=20): G_macro = G_0 * 0.01^20 = 10^-40 G_0
    
    def force_multi_layer(r, M):
        """
        Force from multiple fractal layers.
        Each layer N contributes with:
        - Strength: G_N = G_0 * beta^N
        - Range: Effective at r ~ r_N = r_0 * scale^N
        
        Hypothesis: At small r, shallow layers (N~0) dominate -> 1/r^2
                   At large r, deep layers (N~20) dominate -> modified law
        """
        G_0 = 1.0
        r_0 = 1.0  # Base scale
        scale_factor = 2.0  # Each layer is 2x larger
        
        F_total = 0.0
        
        for N in range(N_LAYERS + 1):
            # Layer strength (exponential damping)
            G_N = G_0 * (BETA_TORS ** N)
            
            # Layer effective range
            r_N = r_0 * (scale_factor ** N)
            
            # Layer contribution (Gaussian or similar envelope)
            # Layer N is most effective near r ~ r_N
            # We use a resonance-like contribution
            
            # Model 1: Each layer has 1/r^2 force but different coupling
            # Weight = exp(-(r - r_N)^2 / (2 * r_N^2))
            # This makes layer N "resonate" at radius r_N
            
            weight = np.exp(-((r - r_N)**2) / (2 * r_N**2))
            
            # Force contribution
            F_N = weight * G_N * M / (r**2 + 0.1)  # Regularize at r->0
            
            F_total += F_N
            
        return F_total
    
    def orbital_velocity_fractal(r, M):
        """v^2 / r = F -> v = sqrt(r * F)"""
        F = force_multi_layer(r, M)
        return np.sqrt(r * F)
    
    # Test Multiple Galaxy Masses
    masses = [10, 50, 100, 200, 500]
    R_MAX = 50.0
    radii = np.linspace(1, R_MAX, 50)
    
    plt.figure(figsize=(12, 6))
    
    # Subplot 1: Rotation Curves
    plt.subplot(1, 2, 1)
    
    v_flats = []
    
    for M in masses:
        velocities = [orbital_velocity_fractal(r, M) for r in radii]
        plt.plot(radii, velocities, label=f'M={M}')
        v_flats.append(velocities[-1])
        
    plt.xlabel('Radius r')
    plt.ylabel('Rotation Velocity v(r)')
    plt.title('Multi-Layer Fractal Rotation Curves')
    plt.legend()
    plt.grid(True)
    
    # Subplot 2: Tully-Fisher
    plt.subplot(1, 2, 2)
    
    log_M = np.log10(masses)
    log_v = np.log10(v_flats)
    
    # Fit
    slope, intercept = np.polyfit(log_v, log_M, 1)
    
    plt.scatter(log_v, log_M, c='r', s=50)
    plt.plot(log_v, slope*log_v + intercept, 'k--', label=f'Slope = {slope:.2f}')
    
    plt.xlabel('Log10(v_flat)')
    plt.ylabel('Log10(Mass)')
    plt.title(f'Tully-Fisher: M ~ v^{slope:.2f}')
    plt.legend()
    plt.grid(True)
    
    plt.tight_layout()
    plt.savefig('QW-586_Result.png')
    
    print(f"\nTully-Fisher Slope: {slope:.2f}")
    print(f"Expected (Newtonian): 2.0")
    print(f"Expected (MOND): 4.0")
    
    if 3.5 <= slope <= 4.5:
        print("✅ SUCCESS: MOND-like Tully-Fisher reproduced!")
    elif 1.8 <= slope <= 2.2:
        print("⚠️ NEWTONIAN: Standard gravity (no dark matter)")
    else:
        print(f"❓ UNEXPECTED: Slope {slope:.2f}")
    
    # Additional Analysis: Check if curves flatten
    # Compute dv/dr at large r
    v_last = [orbital_velocity_fractal(radii[-1], M) for M in masses]
    v_mid = [orbital_velocity_fractal(radii[-10], M) for M in masses]
    
    avg_slope_outer = np.mean([(v_last[i] - v_mid[i])/(radii[-1] - radii[-10]) for i in range(len(masses))])
    
    print(f"\nOuter slope dv/dr: {avg_slope_outer:.4f}")
    if abs(avg_slope_outer) < 0.01:
        print("✅ Curves are FLAT (Dark Matter signature)")
    else:
        print("❌ Curves are NOT flat")

if __name__ == "__main__":
    main()
