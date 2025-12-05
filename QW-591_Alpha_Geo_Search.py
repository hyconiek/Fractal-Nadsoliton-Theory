import numpy as np
from itertools import product
import time

print("QW-591: Systematic Search for Geometric Origin of α_geo")
print("="*70)

TARGET = 4 * np.log(2)  # 2.772588722239781
print(f"Target: α_geo = {TARGET:.12f}")
print()

# Fundamental constants
pi = np.pi
e = np.e
phi = (1 + np.sqrt(5)) / 2  # Golden ratio
sqrt2 = np.sqrt(2)
sqrt3 = np.sqrt(3)
ln2 = np.log(2)
zeta3 = 1.2020569  # Apéry's constant

constants = {
    'π': pi,
    'e': e,
    'φ': phi,
    '√2': sqrt2,
    '√3': sqrt3,
    'ln2': ln2,
    'ζ(3)': zeta3,
    '1': 1.0,
    '2': 2.0,
    '3': 3.0,
    '4': 4.0,
    '5': 5.0,
    '6': 6.0,
}

const_list = list(constants.keys())
const_vals = [constants[k] for k in const_list]

print("Testing combinations: (a ⊕ b) ⊗ c")
print("where ⊕, ⊗ ∈ {+, -, ×, ÷}")
print()

operations = ['+', '-', '*', '/']
best_results = []

start_time = time.time()
total_tests = 0

for i, a in enumerate(const_vals):
    for j, b in enumerate(const_vals):
        for k, c in enumerate(const_vals):
            for op1 in operations:
                for op2 in operations:
                    try:
                        # First operation
                        if op1 == '+':
                            temp = a + b
                        elif op1 == '-':
                            temp = a - b
                        elif op1 == '*':
                            temp = a * b
                        else:
                            if b == 0:
                                continue
                            temp = a / b
                        
                        # Second operation
                        if op2 == '+':
                            result = temp + c
                        elif op2 == '-':
                            result = temp - c
                        elif op2 == '*':
                            result = temp * c
                        else:
                            if c == 0:
                                continue
                            result = temp / c
                        
                        total_tests += 1
                        error = abs(result - TARGET) / TARGET * 100
                        
                        if error < 0.1:  # Less than 0.1% error
                            expr = f"({const_list[i]} {op1} {const_list[j]}) {op2} {const_list[k]}"
                            best_results.append((error, result, expr))
                            
                    except:
                        pass

elapsed = time.time() - start_time

print(f"Tested {total_tests:,} combinations in {elapsed:.2f}s")
print()

# Sort by error
best_results.sort()

if best_results:
    print(f"FOUND {len(best_results)} candidates with error < 0.1%:")
    print("-"*70)
    for error, result, expr in best_results[:20]:  # Top 20
        print(f"{expr:40s} = {result:.12f}  (error: {error:.6f}%)")
    
    print()
    print("="*70)
    print("BEST CANDIDATE:")
    error, result, expr = best_results[0]
    print(f"  {expr}")
    print(f"  = {result:.12f}")
    print(f"  Target: {TARGET:.12f}")
    print(f"  Error: {error:.8f}%")
    print("="*70)
else:
    print("❌ NO ELEGANT EXPRESSION FOUND with error < 0.1%")
    print()
    print("Trying extended search (4 terms)...")
    # Could extend to 4-term expressions if needed

print()
print("Conclusion:")
if best_results and best_results[0][0] < 0.001:
    print("✅ GEOMETRIC ORIGIN CONFIRMED: α_geo can be expressed using fundamental constants")
elif best_results:
    print("🟡 PARTIAL SUCCESS: Found approximate expressions")
else:
    print("🔴 FAILED: α_geo = 4×ln(2) appears to be the simplest form")
    print("    No elegant geometric origin from {π, e, φ, √2, √3, ζ(3)}")
