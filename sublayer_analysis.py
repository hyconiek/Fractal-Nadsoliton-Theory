import numpy as np

# 1. Data Setup
# Masses in MeV
particles = {
    "Top": {"mass": 172760, "type": "quark", "gen": 3},
    "Bottom": {"mass": 4180, "type": "quark", "gen": 3},
    "Tau": {"mass": 1776.86, "type": "lepton", "gen": 3},
    "Charm": {"mass": 1270, "type": "quark", "gen": 2},
    "Muon": {"mass": 105.66, "type": "lepton", "gen": 2},
    "Strange": {"mass": 93, "type": "quark", "gen": 2},
    "Down": {"mass": 4.7, "type": "quark", "gen": 1},
    "Up": {"mass": 2.2, "type": "quark", "gen": 1},
    "Electron": {"mass": 0.511, "type": "lepton", "gen": 1}
}

M_top = 172760
GAMMA = 1.52

def calculate_Q(mass):
    ratio = mass / M_top
    if ratio <= 0: return 0
    log_val = np.log(ratio) / np.log(4)
    Q = - (4 / GAMMA) * log_val
    return Q

# 2. Analyze Q Sublayers (Modulo 4)
print(f"{'Particle':<10} {'Mass':<10} {'Calc Q':<8} {'Round Q':<8} {'k (Mod 4)':<10} {'Fib Match'}")
print("-" * 70)

fibs = [0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89]

results = []

for name, data in particles.items():
    Q_calc = calculate_Q(data['mass'])
    Q_round = round(Q_calc)
    
    k = Q_round % 4
    is_fib = Q_round in fibs
    
    print(f"{name:<10} {data['mass']:<10.2f} {Q_calc:<8.2f} {Q_round:<8} {k:<10} {is_fib}")
    
    results.append({
        "name": name,
        "type": data['type'],
        "Q": Q_round,
        "k": k
    })

print("\n=== SUBSLAYER ANALYSIS ===")
# Group by k manually
for k in range(4):
    print(f"\nSublayer k={k} (Fraction 0.{k*25}):")
    # Filter list
    group = [x for x in results if x['k'] == k]
    if group:
        print(f"{'Name':<10} {'Type':<10} {'Q':<5}")
        print("-" * 30)
        for item in group:
            print(f"{item['name']:<10} {item['type']:<10} {item['Q']:<5}")
    else:
        print("No particles found.")

print("\n=== FIBONACCI MODULO 4 ANALYSIS ===")
print("Fibonacci Integers and their Sublayer positions:")
for f in fibs:
    print(f"F={f} -> k={f%4}")
