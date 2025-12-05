
import numpy as np
import matplotlib.pyplot as plt

def main():
    print("QW-584: Grand Unified Simulation (FIN Theory)")
    print("---------------------------------------------")
    print("Combining: Topological Gravity + Recursive Time Dilation")
    
    # Configuration
    GRID_SIZE = 60
    STEPS = 2000
    N_PARTICLES = 50
    
    # 1. Spacetime Topology (with Mass)
    # We use a weighted graph approach implicitly on a grid.
    # Connectivity density rho(r) ~ 1/r
    
    CENTER = np.array([GRID_SIZE/2, GRID_SIZE/2])
    MASS_RADIUS = 5.0
    
    # The "Singularity" Portal
    PORTAL_LOC = (int(CENTER[0]), int(CENTER[1]))
    
    # Inner Universe (for Time Dilation)
    # Simple delay model: T_inner = T_outer * Dilaton_Factor
    DILATION_FACTOR = 100 
    
    print("Simulating Geodesics in Fractal Spacetime...")
    
    # Store trajectories
    trajectories = []
    time_delays = []
    
    for p_idx in range(N_PARTICLES):
        # Start from left, aim for right (Scattering experiment)
        start_y = np.random.normal(GRID_SIZE/2, GRID_SIZE/6)
        pos = np.array([2.0, start_y])
        
        path = [pos.copy()]
        proper_time = 0
        coordinate_time = 0
        
        # Velocity (drift to right)
        vel = np.array([1.0, 0.0])
        
        active = True
        
        while active and coordinate_time < STEPS:
            coordinate_time += 1
            
            # 1. Topological Gravity (P(r) ~ connectivity)
            # Probability to move towards center is enhanced by "available links"
            r_vec = CENTER - pos
            r_dist = np.linalg.norm(r_vec)
            
            # Effective Force / Drift towards higher connectivity
            # Gradient of Entropy S ~ ln(Connectivity).
            # If Connectivity ~ 1/r, then Force ~ 1/r^2.
            
            if r_dist > 0.5:
                # Stochastic bias towards center
                bias_strength = 20.0 / (r_dist**2 + 1.0) # 1/r^2 attraction
                # Limit bias
                bias_strength = min(bias_strength, 0.5)
                
                # Update velocity (acceleration)
                vel += bias_strength * (r_vec / r_dist) * 0.1
                
                # Drag/Friction (Viscosity of Ether) - limits speed
                vel *= 0.95
                
                # Noise (Quantum fluctuations)
                vel += np.random.normal(0, 0.1, 2)
            
            # Normalize speed to roughly 1 unless in potential well?
            # Actually, near mass speed should increase (falling) or decrease (refraction)?
            # In QW-580 we used refraction (slower).
            # In Entropic Gravity, particles accelerate.
            # Let's simple allow acc.
            
            # 2. Portal Check (Recursive Nesting)
            # If particle hits the center pixel
            curr_grid_pos = (int(pos[0]), int(pos[1]))
            if curr_grid_pos == PORTAL_LOC:
                # Enter Inner Universe
                # Add delay to Proper Time but NOT Coordinate X/Y change
                # Wait, Proper Time increases, but Coordinate Time also passes?
                # Observer at infinity sees particle "stuck" there.
                # So Coordinate Time ticks increase, Pos stays same.
                
                delay = int(np.random.exponential(DILATION_FACTOR))
                coordinate_time += delay
                # Proper time (for particle) might be just "a short trip inside"
                proper_time += 10 # Subjective
                
                # Log
                # print(f"Particle {p_idx} trapped in portal for {delay} steps.")
                # Teleport out randomly nearby?
                pos += np.random.normal(0, 1.0, 2)
            else:
                proper_time += 1
                
                # Update Position
                pos += vel
            
            path.append(pos.copy())
            
            # Exit conditions
            if pos[0] > GRID_SIZE - 2 or pos[0] < 0 or pos[1] > GRID_SIZE or pos[1] < 0:
                active = False
                
        trajectories.append(np.array(path))
        time_delays.append(proper_time / coordinate_time)
        
    print("Simulation Complete.")
    
    # Visualization
    plt.figure(figsize=(10, 10))
    
    # Plot mass density/gravity well
    X, Y = np.meshgrid(np.linspace(0, GRID_SIZE, 100), np.linspace(0, GRID_SIZE, 100))
    R = np.sqrt((X-CENTER[0])**2 + (Y-CENTER[1])**2)
    Potential = -1.0 / (R + 1.0)
    plt.contourf(X, Y, Potential, levels=20, cmap='Greys', alpha=0.3)
    
    # Plot trajectories
    for tr in trajectories:
        plt.plot(tr[:, 0], tr[:, 1], 'b-', alpha=0.3, linewidth=1)
        
    # Mark Portal
    plt.plot(CENTER[0], CENTER[1], 'ro', markersize=10, label='Singularity (Recursive Hub)')
    
    plt.xlim(0, GRID_SIZE)
    plt.ylim(0, GRID_SIZE)
    plt.title("QW-584: Grand Unified Simulation\nEntropic Gravity + Recursive Time Dilation")
    plt.xlabel("Space X")
    plt.ylabel("Space Y")
    plt.legend()
    
    output_png = "QW-584_Result.png"
    plt.savefig(output_png)
    print(f"Saved {output_png}")
    
    # Stats
    avg_dilation = 1.0 / np.mean(time_delays)
    print(f"Average Time Dilation Factor: {avg_dilation:.2f}")

if __name__ == "__main__":
    main()
