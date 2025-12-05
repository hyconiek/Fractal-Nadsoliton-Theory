
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def main():
    print("QW-583: Statistical Entropic Gravity")
    print("------------------------------------")
    print("Hypothesis: Mass = Region of High Connectivity (Information Density)")
    print("Prediction: Walker Probability P(r) ~ Degree k(r) ~ 1/r (Entropic Attraction)")
    
    # Constants
    GRID_SIZE = 50
    N_NODES = GRID_SIZE * GRID_SIZE
    N_WALKERS = 2000
    STEPS = 5000
    
    # 1. Initialize 2D Grid Graph (Vacuum)
    adj = {}
    
    print(f"Initializing {GRID_SIZE}x{GRID_SIZE} Lattice...")
    
    def get_idx(x, y):
        return x * GRID_SIZE + y
        
    def get_pos(idx):
        return (idx // GRID_SIZE, idx % GRID_SIZE)
        
    for x in range(GRID_SIZE):
        for y in range(GRID_SIZE):
            u = get_idx(x, y)
            adj[u] = []
            # 4-neighbors 
            for dx, dy in [(-1,0), (1,0), (0,-1), (0,1)]:
                nx, ny = x + dx, y + dy
                if 0 <= nx < GRID_SIZE and 0 <= ny < GRID_SIZE:
                    v = get_idx(nx, ny)
                    adj[u].append(v)
                    
    # 2. Add Mass (Deformation of Topology)
    # Model: Probability of extra link decreases with distance from center
    CENTER_X = GRID_SIZE // 2
    CENTER_Y = GRID_SIZE // 2
    
    # Warping Factor
    # We add M extra links.
    # Source node u is random (or weighted by centrality?)
    # Target node v is random (or weighted?)
    # Let's effectively increase DEGREE k(r) ~ 1/r + const.
    
    print("Warping Topology (Adding Links ~ 1/r)...")
    
    MAX_DEGREE_BOOST = 20 # Max extra links at center
    
    added_links = 0
    
    for x in range(GRID_SIZE):
        for y in range(GRID_SIZE):
            r = np.sqrt((x - CENTER_X)**2 + (y - CENTER_Y)**2) + 0.5
            # Target degree boost
            boost = int(MAX_DEGREE_BOOST / r)
            
            u = get_idx(x, y)
            
            # Add 'boost' random links to local neighborhood or global?
            # Metric metric is local. 
            # If we add links to far away nodes, that's wormholes.
            # If we add links to neighbors, we just make multi-edges (weights).
            # Let's simulating "Finer mesh" by adding links to 2nd nearest neighbors or random local.
            # Or just assume the graph is weighted.
            
            # Let's add shortcuts to random nodes within radius R_LOCAL
            R_LOCAL = 5
            
            for _ in range(boost):
                # Pick random offset
                dx = np.random.randint(-R_LOCAL, R_LOCAL+1)
                dy = np.random.randint(-R_LOCAL, R_LOCAL+1)
                nx, ny = x + dx, y + dy
                if 0 <= nx < GRID_SIZE and 0 <= ny < GRID_SIZE:
                    v = get_idx(nx, ny)
                    if v != u:
                        adj[u].append(v)
                        adj[v].append(u) # Undirected
                        added_links += 1

    print(f"Added {added_links} topological defects (links).")

    # 3. Simulation: Random Walkers to Equilibrium
    walkers = np.random.randint(0, N_NODES, N_WALKERS)
    visits = np.zeros(N_NODES)
    
    print(f"Simulating {N_WALKERS} walkers for {STEPS} steps...")
    
    for t in range(STEPS):
        new_walkers = []
        for w in walkers:
            if len(adj[w]) > 0:
                next_node = np.random.choice(adj[w])
                new_walkers.append(next_node)
                # Count visits only in last 1000 steps (Equilibrium)
                if t > STEPS - 1000:
                    visits[next_node] += 1
            else:
                new_walkers.append(w)
        walkers = np.array(new_walkers)
        
        if t % 1000 == 0:
            print(f"Step {t}...")
            
    # 4. Analysis
    print("Analyzing Distribution...")
    
    radii = []
    densities = [] # Normalized visits / Area (or just raw visits per node?)
    # Visits per node is proportional to P(r).
    
    for x in range(GRID_SIZE):
        for y in range(GRID_SIZE):
            r = np.sqrt((x - CENTER_X)**2 + (y - CENTER_Y)**2)
            u = get_idx(x, y)
            if r > 1.0 and r < 20.0:
                radii.append(r)
                densities.append(visits[u])
                
    radii = np.array(radii)
    densities = np.array(densities)
    
    # Binning
    bins = np.linspace(1.0, 20.0, 30)
    digitized = np.digitize(radii, bins)
    
    bin_r = []
    bin_rho = []
    
    for i in range(1, len(bins)):
        mask = digitized == i
        if np.sum(mask) > 0:
            # Density per NODE is what we measured (visits[u]).
            # We want to know if P(u) ~ 1/r_u.
            bin_r.append(np.mean(radii[mask]))
            bin_rho.append(np.mean(densities[mask]))
            
    bin_r = np.array(bin_r)
    bin_rho = np.array(bin_rho)
    
    # Fit
    def power_law(r, a, b):
        return a * r**(-b) + 50 # Background noise
        
    try:
        popt, pcov = curve_fit(power_law, bin_r, bin_rho, p0=[1000, 1.0], maxfev=2000)
        exponent = popt[1]
        print(f"Fit Result: Density ~ r ^ -{exponent:.4f}")
    except Exception as e:
        print(f"Fit failed: {e}")
        exponent = 0.0
        popt = [0, 0]

    # Visualization
    plt.figure(figsize=(12, 6))
    
    plt.subplot(1, 2, 1)
    visits_grid = visits.reshape((GRID_SIZE, GRID_SIZE))
    plt.imshow(visits_grid, origin='lower', cmap='plasma')
    plt.colorbar(label='Walker Density')
    plt.title(f"QW-583: Statistical Gravity (Connectivity)\nExponent n = {exponent:.2f}")

    plt.subplot(1, 2, 2)
    plt.scatter(bin_r, bin_rho, c='k', label='Simulated Data')
    if popt[0] != 0:
        plt.plot(bin_r, power_law(bin_r, *popt), 'r--', label=f'Fit: $r^{{-{exponent:.2f}}}$')
    plt.xlabel('Radius r')
    plt.ylabel('Probability Density P(r)')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.grid(True, which="both", ls="-")
    plt.title("Radial Density Distribution (Log-Log)")
    
    output_file = "QW-583_Result.png"
    plt.savefig(output_file)
    print(f"Saved {output_file}")


if __name__ == "__main__":
    main()
