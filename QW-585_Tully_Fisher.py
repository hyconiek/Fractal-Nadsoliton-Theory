
import numpy as np
import matplotlib.pyplot as plt

def main():
    print("QW-585: Tully-Fisher Relation Test (Galaxy Rotation)")
    print("----------------------------------------------------")
    print("Hypothesis: Entropic Gravity produces Flat Rotation Curves (Dark Matter mimicry)")
    
    # Simulation Parameters
    R_MAX = 50.0
    DR = 1.0
    
    # We will simulate multiple Galaxy Masses
    # In our model, Mass = Central Connectivity Hub Strength (M_hub)
    # Effect: Force F(r) ~ M_hub / r^alpha
    
    # Standard Gravity: alpha = 2.0 -> v ~ 1/sqrt(r) (Keplerian)
    # Entropic/Fractal Gravity?
    # From QW-583 we saw density P(r) ~ 1/r^0.56.
    # Force is gradient of effective potential ~ ln(P). 
    # If P(r) ~ r^-n, then U(r) ~ -ln(r^-n) ~ n ln(r).
    # Force F = -dU/dr ~ -n/r.
    # CONSTANT FORCE? Or 1/r force?
    # If F ~ 1/r, then v^2/r ~ 1/r => v^2 ~ const => v ~ const.
    # THIS WOULD EXPLAIN FLAT ROTATION CURVES!
    
    # Let's test this derivation numerically.
    # If our density was indeed P(r) ~ r^-n, does it imply F ~ 1/r?
    
    # We assume 'n' (exponent of density) depends on Mass? Or is it universal?
    # QW-583 gave n=0.56. Let's assume F_entropic(r) = C_coupling * Mass / r (in 2D/Fractal limit).
    
    def orbital_velocity(r, M, model='newton'):
        if model == 'newton':
            # F = G M / r^2 -> v = sqrt(G M / r)
            G = 1.0
            if r == 0: return 0
            return np.sqrt(G * M / r)
        elif model == 'entropic_fin':
            # F = epsilon * M / r (Entropic 1/r decay due to fractal dimension D~2)
            # v^2 / r = F = eps M / r
            # v^2 = eps M
            # v = sqrt(eps M) -> CONSTANT !
            epsilon = 0.1
            if r < 2.0: # Core region (Newtonian dominate?)
                return np.sqrt(1.0 * M / r) 
            return np.sqrt(epsilon * M) # Flat regime
            
    # Actually, let's SIMULATE particle motion to measure v, not just assume formula.
    # Use the Force Law inferred from QW-583: F ~ 1/r (roughly, from n=0.56 density slope implies 1/r force in entropic logic).
    
    masses = [10, 50, 100, 200, 500]
    results_v_flat = []
    
    plt.figure(figsize=(10, 6))
    
    for M in masses:
        radii = np.linspace(1, R_MAX, 50)
        velocities = []
        
        # We model the force explicitly as a mix of Newtonian (short range) and Entropic (long range)
        # F_total = G*M/r^2 + lambda*M/r (MOND-like / Entropic)
        
        G = 1.0
        Lambda = 0.05 # Entropic coupling strength
        
        for r in radii:
            # Centripetal balance: v^2 / r = F_total
            F = (G * M) / (r**2) + (Lambda * M) / r
            v = np.sqrt(r * F)
            velocities.append(v)
            
        plt.plot(radii, velocities, label=f'Mass={M}')
        
        # Estimate "Flat" velocity (v_flat) at R_MAX
        results_v_flat.append(velocities[-1])
        
    plt.xlabel('Radius r')
    plt.ylabel('Rotation Velocity v')
    plt.title('QW-585: Galaxy Rotation Curves (Entropic + Newtonian)')
    plt.legend()
    plt.grid(True)
    plt.savefig('QW-585_Rotation_Curves.png')
    
    # Tully-Fisher Analysis
    # Relation: M ~ v^alpha
    # We plot log(M) vs log(v_flat)
    
    plt.figure(figsize=(6, 6))
    log_M = np.log10(masses)
    log_v = np.log10(results_v_flat)
    
    # Fit line
    slope, intercept = np.polyfit(log_v, log_M, 1) # M against v
    
    plt.scatter(log_v, log_M, c='r', s=50)
    plt.plot(log_v, slope*log_v + intercept, 'k--', label=f'Slope alpha={slope:.2f}')
    
    plt.xlabel('Log10(Velocity)')
    plt.ylabel('Log10(Mass)')
    plt.title(f'Tully-Fisher Relation: M ~ v^{slope:.2f}')
    plt.legend()
    plt.grid(True)
    
    plt.savefig('QW-585_Tully_Fisher.png')
    print(f"Tully-Fisher Slope: {slope:.2f}")
    
    if 3.5 <= slope <= 4.5:
        print("SUCCESS: Matches Observed Tully-Fisher Relation (alpha ~ 4)")
    else:
        print(f"FAIL: Slope {slope:.2f} differs from M ~ v^4")

if __name__ == "__main__":
    main()
