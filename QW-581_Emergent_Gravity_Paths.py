
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def main():
    print("QW-581: Emergent Gravity (Action Weighted)")
    print("------------------------------------------")
    
    # Constants
    N_PATHS = 50000        # High count for statistics
    PATH_LENGTH = 60
    R_BOX = 10.0
    MASS_CENTER = np.array([5.0, 5.0])
    
    # Physics Parameters
    # Refractive Index n(r) = 1 + alpha / r
    # Action S = Sum (n(r) * ds)
    # Weight w = exp(-beta * S)
    ALPHA = 0.5
    BETA = 2.0  # Inverse "Planck constant" - controls quantumness
                # High Beta -> Classical limit (Selection of Geodesic)
                # Low Beta -> High fluctuations
                
    print(f"Simulating {N_PATHS} paths. Alpha={ALPHA}, Beta={BETA}")
    
    # Grid for density
    GRID_RES = 20 # 20x20 grid
    density_map = np.zeros((GRID_RES, GRID_RES))
    
    # Launch configuration: Beam from left
    start_x = 1.0
    end_x = 9.0
    y_starts = np.linspace(2, 8, 30) # Emitters
    
    all_paths = []
    all_weights = []
    
    BATCH_SIZE = N_PATHS // len(y_starts)
    
    for start_y in y_starts:
        start_pos = np.array([start_x, start_y])
        end_pos = np.array([end_x, start_y]) # Target (not forced, just bias guidance)
        
        for i in range(BATCH_SIZE):
            path = [start_pos]
            current_pos = start_pos
            action = 0.0
            
            valid = True
            
            for step in range(PATH_LENGTH):
                # 1. Biased Walk (Guiding Function to ensure paths reach other side)
                # We bias slightly towards +x to simulate propagation
                angle = np.random.normal(0, 1.0) # Angle relative to motion?
                # No, absolute steps.
                # Bias towards target direction roughly
                dir_vec = end_pos - current_pos
                dist_rem = np.linalg.norm(dir_vec)
                if dist_rem > 0.1:
                    dir_vec /= dist_rem
                else:
                    dir_vec = np.array([1.0, 0.0])
                    
                # Mix random direction with target direction
                # v = v_random + 0.5 * v_target
                rng_angle = np.random.uniform(0, 2*np.pi)
                v_rnd = np.array([np.cos(rng_angle), np.sin(rng_angle)])
                
                # Step vector
                step_vec = (v_rnd * 0.5 + dir_vec * 0.5) 
                step_vec /= np.linalg.norm(step_vec) # Normalize
                
                step_len = np.random.exponential(0.5) # Stepsize
                
                dx = step_vec * step_len
                proposed_pos = current_pos + dx
                
                # Check bounds
                if not (0 <= proposed_pos[0] <= R_BOX and 0 <= proposed_pos[1] <= R_BOX):
                    valid = False
                    break
                
                # Compute Action Increment
                r_mid = np.linalg.norm(current_pos + dx/2 - MASS_CENTER)
                # Avoid singularity
                if r_mid < 0.2:
                    n_index = 5.0 # Cap
                else:
                    n_index = 1.0 + ALPHA / r_mid
                    
                ds = step_len
                dS = n_index * ds
                action += dS
                
                path.append(proposed_pos)
                current_pos = proposed_pos
                
                # Stop if reached x > 9.5
                if current_pos[0] > 9.5:
                    break
            
            if valid:
                all_paths.append(np.array(path))
                all_weights.append(np.exp(-BETA * action))
                
    # Normalize weights
    all_weights = np.array(all_weights)
    max_log_w = np.max(-BETA * np.array([p.shape[0] for p in all_paths])) # Wait, Action is roughly proportional to length
    # Action varies.
    # To avoid underflow, shift weights?
    # Actually, let's just use what we have. If beta is high, weights might drop fast.
    # Check max weight
    max_w = np.max(all_weights)
    print(f"Max weight: {max_w}, Min weight: {np.min(all_weights)}")
    
    # Map accumulation
    print("Accumulating weighted density...")
    for p, w in zip(all_paths, all_weights):
        for pos in p:
            ix = int(pos[0] / (R_BOX/GRID_RES))
            iy = int(pos[1] / (R_BOX/GRID_RES))
            if 0 <= ix < GRID_RES and 0 <= iy < GRID_RES:
                density_map[ix, iy] += w

    # Analyze Density vs Radius from Mass
    print("Profiling density...")
    radii = []
    densities = []
    
    for ix in range(GRID_RES):
        for iy in range(GRID_RES):
            x = (ix + 0.5) * (R_BOX/GRID_RES)
            y = (iy + 0.5) * (R_BOX/GRID_RES)
            r = np.sqrt((x - 5.0)**2 + (y - 5.0)**2)
            if r > 0.5 and r < 4.5:
                radii.append(r)
                densities.append(density_map[ix, iy])
                
    radii = np.array(radii)
    densities = np.array(densities)
    
    # Fit
    def power_law(r, a, b):
        return a * r**(-b)
        
    try:
        # Ignore zeros
        mask = densities > 0
        popt, _ = curve_fit(power_law, radii[mask], densities[mask], p0=[1.0, 1.0], maxfev=1000)
        exponent = popt[1]
        print(f"Exponent n = {exponent:.4f}")
    except:
        exponent = 0.0
        print("Fit failed")
        
    # Visualization
    plt.figure(figsize=(12, 6))
    
    plt.subplot(1, 2, 1)
    plt.imshow(density_map.T, origin='lower', extent=[0, 10, 0, 10], cmap='inferno')
    plt.colorbar(label='Weighted Probability')
    plt.scatter([5], [5], c='cyan', marker='x', label='Mass')
    plt.title(f"QW-581: Quantum Probability Density\nBeta={BETA}, Alpha={ALPHA}")
    
    plt.subplot(1, 2, 2)
    plt.scatter(radii, densities, alpha=0.5, c='k', s=5, label='Grid Points')
    if exponent != 0:
        r_line = np.linspace(0.5, 4.5, 100)
        plt.plot(r_line, power_law(r_line, *popt), 'r-', label=f'Fit: $r^{{-{exponent:.2f}}}$')
    plt.xlabel('Radius from Mass')
    plt.ylabel('Density')
    plt.yscale('log')
    plt.legend()
    plt.grid(True, which="both", ls="-")
    plt.title("Radial Density Falloff")
    
    plt.savefig("QW-581_Result_Weighted.png")
    print("Saved QW-581_Result_Weighted.png")

if __name__ == "__main__":
    main()
