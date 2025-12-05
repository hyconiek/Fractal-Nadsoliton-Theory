
import numpy as np
import matplotlib.pyplot as plt


def main():
    print("QW-582: Recursive Nesting & Time Dilation")
    print("-----------------------------------------")
    
# 1. Define Outer Universe (Macro Scale)
    # A simple path or cycle
    N_OUTER = 20
    # Outer graph: Node i connects to (i-1)%N and (i+1)%N
    
    # Define the "Portal" node that contains the Micro Universe
    PORTAL_NODE = 10
    
    # 2. Define Inner Universe (Micro Scale)
    N_INNER = 100
    # Inner graph: Random connectivity
    # inner_adj = {i: [list of neighbors]}
    inner_adj = {}
    for i in range(N_INNER):
        # Watts-Strogatz like: connect to neighbors + random
        # Simple Circle + Random shortcuts
        neighbors = [(i-1)%N_INNER, (i+1)%N_INNER]
        # Random shortcut
        if np.random.rand() < 0.1:
            neighbors.append(np.random.randint(0, N_INNER))
        inner_adj[i] = neighbors
    
    print(f"Outer Graph: {N_OUTER} nodes. Portal at Node {PORTAL_NODE}.")
    print(f"Inner Graph: {N_INNER} nodes (hidden inside Portal).")
    
    # 3. Simulation: Walker
    # Walker moves on Outer graph. If it hits Portal, it enters Inner graph.
    # It must traverse Inner graph from Entry to Exit to return to Outer.
    
    path_log = [] # List of strings describing state
    steps_outer = 0
    steps_total = 0
    
    current_node = 0
    in_portal = False
    
    inner_node = 0
    inner_target = N_INNER // 2 # Target to escape
    
    MAX_STEPS = 200
    
    print("\n--- Starting Simulation ---")
    
    history_x = []
    history_y = [] # For plotting "Effective Time" vs "Proper Time"
    
    while steps_outer < MAX_STEPS:
        steps_total += 1
        
        if not in_portal:
            # OUTER DYNAMICS
            steps_outer += 1
            path_log.append(f"Outer Step {steps_outer}: Node {current_node}")
            
            # Simple random walk on outer graph
            # neighbors = list(G_outer.neighbors(current_node))
            neighbors = [(current_node - 1) % N_OUTER, (current_node + 1) % N_OUTER]
            
            # Drift towards portal for demonstration
            # Find neighbor closest to portal
            best_n = neighbors[0]
            min_dist = abs(best_n - PORTAL_NODE)
            for n in neighbors:
                d = abs(n - PORTAL_NODE)
                if d < min_dist:
                    min_dist = d
                    best_n = n
            
            # 80% chance to move towards portal, 20% random
            if np.random.rand() < 0.8:
                current_node = best_n
            else:
                current_node = np.random.choice(neighbors)
            
            # Check Portal Entry
            if current_node == PORTAL_NODE:
                print(f"  -> Observer entered Portal at Outer Step {steps_outer}!")
                in_portal = True
                inner_node = 0 # Start at 0 of inner graph
        else:
            # INNER DYNAMICS
            # Outer time is FROZEN (or ticks very slowly relative to inner)
            # We count Total Steps (Proper Time of Observer) vs Outer Steps (Coordinate Time)
            
            # Random walk on inner graph
            neighbors = inner_adj[inner_node]
            
            # Drift towards exit (inner_target)
            # BFS distance check would be better but expensive loops.
            # Let's just do random walk + slight bias if we knew topology.
            # Watts-Strogatz is small world, random walk finds target fast enough?
            
            # Simple bias: if any neighbor is closer to target? No topology info
            inner_node = np.random.choice(neighbors)
            
            if inner_node == inner_target:
                print(f"  <- Observer escaped Portal! Spent {steps_total - steps_outer} ticks inside.")
                in_portal = False
                # Eject to a neighbor of Portal in outer graph
                current_node = (PORTAL_NODE + 1) % N_OUTER
        
        history_x.append(steps_outer) # Coordinate Time (Observer at Infinity)
        history_y.append(steps_total) # Proper Time (Observer)

    print("--- Simulation Complete ---")
    
    # Analysis
    # The slope of Proper Time vs Coordinate Time indicates Time Dilation
    # Slope = 1 (Normal space)
    # Slope >> 1 (Inside Portal/Black Hole)
    
    history_x = np.array(history_x)
    history_y = np.array(history_y)
    
    # Visualization
    plt.figure(figsize=(10, 6))
    
    plt.plot(history_x, history_y, 'b-', linewidth=2, label='Worldline')
    plt.plot([0, MAX_STEPS], [0, MAX_STEPS], 'k--', alpha=0.5, label='Flat Space (1:1)')
    
    plt.xlabel('Outer Coordinate Time (Steps)')
    plt.ylabel('Observer Proper Time (Ticks)')
    plt.title('QW-582: Recursive Nesting Time Dilation')
    plt.legend()
    plt.grid(True)
    
    # Annotate region
    # Find jump
    diff = history_y - history_x
    if np.max(diff) > 5:
        plt.annotate('Trapped in Inner Universe', xy=(MAX_STEPS/2, np.max(history_y)/2), 
                     xytext=(MAX_STEPS/2 + 20, np.max(history_y)/2 - 50),
                     arrowprops=dict(facecolor='red', shrink=0.05))

    plt.savefig("QW-582_Result.png")
    print("Saved QW-582_Result.png")

if __name__ == "__main__":
    main()
