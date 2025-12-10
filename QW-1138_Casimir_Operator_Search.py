
import numpy as np
from itertools import product, combinations

print('SEARCHING FOR CASIMIR OPERATOR D')
print('=' * 60)

# Target values (d) and their "Quarter-Bit" integers (Q = 4*d)
targets = {
    'Top': {'d': 0.0, 'Q': 0},
    'Bottom': {'d': 1.75, 'Q': 7},
    'Tau': {'d': 2.25, 'Q': 9},
    'Muon': {'d': 3.5, 'Q': 14},
    'Electron': {'d': 6.0, 'Q': 24}
}

target_Q = [0, 7, 9, 14, 24]
print(f'Target Integers (Q = 4*d): {target_Q}')
print('Differences:', [target_Q[i+1]-target_Q[i] for i in range(len(target_Q)-1)])

# HYPOTHESIS 1: Information Operator on 4-bit state
# State |psi> is defined by (Generation, Winding, Spin?) or pure bit sequence?
# Nadsoliton is a KNOT.
# Knots are classified by Crossing Number? Or Braiding?

print('\nTEST 1: 4-BIT HAMMING OPERATOR')
# Can we construct Q from popcounts of 5-bit or n-bit integers?
# We found binary patterns earlier:
# 0 (00000), 7 (00111), 9 (01001), 14 (01110), 24 (11000)

ints = [0, 7, 9, 14, 24]
for x in ints:
    print(f'  {x:2d}: {x:05b} (Hamming: {bin(x).count("1")})')

# Hamming weights: 
# 0 -> 0
# 7 -> 3
# 9 -> 2
# 14 -> 3
# 24 -> 2
# Pattern: 0, 3, 2, 3, 2. Not unique.

print('\nTEST 2: FRACTAL ADDRESSING OPERATOR')
# Fractal address: similar to Sierpinski triangle or Cantor set?
# Cantor set removes middle thirds. Base-4 removes?
# Base-4 representation of integers?

for x in ints:
    # Convert to base 4 digits
    s = ""
    n = x
    if n == 0: s = "0"
    while n > 0:
        s = str(n % 4) + s
        n //= 4
    print(f'  {x:2d} in Base-4: {s:>3}')

# 0 -> 0
# 7 -> 13
# 9 -> 21
# 14 -> 32
# 24 -> 120 (Wait: 1*16 + 2*4 + 0 = 24)
#
# Digits sum?
# 0 -> 0
# 7 -> 1+3 = 4
# 9 -> 2+1 = 3
# 14 -> 3+2 = 5
# 24 -> 1+2+0 = 3
# No obvious pattern.

print('\nTEST 3: CASIMIR OF SYMMETRY GROUP')
# Assume symmetry group G.
# C = Casimir eigenvalues.
# Maybe SU(2) x SU(2)? Or SO(4)? Or Fractal Group?

# Try polynomial generator
# Find P(x) such that P(0)=0, P(1)=7, P(2)=9, P(3)=14, P(4)=24?
# x could be "generation index" or "winding number"

# Index i: 0 (Top), 1 (Bottom), 2 (Tau), 3 (Muon), 4 (Electron)
# Wait, ordering by mass is Top -> Bottom -> Tau -> Muon -> Electron
# But ordering by d is Top (0) -> Bottom (1.75) -> Tau (2.25) -> Muon (3.5) -> Electron (6.0)
# This is monotonic.
# Let x = 0, 1, 2, 3, 4.
# Search for P(x) = ax^2 + bx + c
from scipy.optimize import curve_fit

x_data = np.array([0, 1, 2, 3, 4])
y_data = np.array([0, 7, 9, 14, 24])

def poly2(x, a, b, c): return a*x**2 + b*x + c
popt, _ = curve_fit(poly2, x_data, y_data)
a, b, c = popt
print(f'Fit Quadratic: {a:.2f}x^2 + {b:.2f}x + {c:.2f}')
residuals = y_data - poly2(x_data, *popt)
print(f'Residuals: {residuals}')

# 1.36 x^2 + 0.69 x + 1.2 ??? Bad fit (residuals ~ 1-2)

print('\nTEST 4: TOPOLOGICAL KNOT INVARIANT?')
# 3D Knot theory (Nadsoliton).
# Jones Polynomial? Alexander Polynomial?
# Maybe Q = Crossing Number?
# Known particles as knots (e.g. Trefoil = 3 crossings).
# 
# Hypothesis:
# Top = Unknot? (0 crossings) -> Q=0
# Bottom = ? (7 crossings?)
# Tau = ? (9 crossings?)
# Muon = ? (14 crossings?)
# Electron = ? (24 crossings?)
#
# Is there a knot sequence 0, 7, 9, 14, 24?
# Torus knots T(p,q)? Crossing = min(p(q-1), q(p-1))
# T(2,3) = 3 (Trefoil)
# T(2,5) = 5
# T(2,7) = 7 (Wow, Bottom?)
# T(3,4) = 8
# T(3,5) = 11?
#
# Check Torus Knots T(p,q) crossing numbers: c = p*q - p - q (?) No.
# For T(p,q), c = p*q - p - q (if p,q > 0) ? No, crossings are usually p(q-1).
# Formula: c = p(q-1) assuming p < q.
#
# Let's search for (p,q) that give 7, 9, 14, 24.

print('Searching Torus Knots T(p,q) with c = p(q-1) for c in [7, 9, 14, 24]')
found_knots = {c: [] for c in target_Q}
found_knots[0] = ["Unknot T(1,k)"]

for p in range(2, 10):
    for q in range(p+1, 20):
        c = p * (q - 1)
        if c in target_Q:
            found_knots[c].append(f'T({p},{q})')

for c in target_Q:
    print(f'  Q={c}: {found_knots[c]}')
    
# 7: T(7,2) but p<q -> T(?,?)
# If p=1, c = 1*(q-1) = q-1.
# So T(1,8) -> 7. But T(1,k) is unknot.
#
# Formula for crossing number of T(p,q) is actually p*q - p? No.
# Standard: c(T(p,q)) = min(p(q-1), q(p-1)).
#
# Let's check T(p,q) that match.

print('\nTEST 5: FRACTAL LEVEL OPERATOR')
# Q = sum of digits in base phi? Or base 2?
# Let's look at the "Fractal Information" aspect.
#
# Fractal Information => Depth + Complexity.
#
# What if d = N * (N-1) ?
# 0 -> 0
# 1 -> 0
# 2 -> 2
# 3 -> 6
# 4 -> 12
# 5 -> 20
#
# Sequence: 0, 7, 9, 14, 24.
# Diffs: 7, 2, 5, 10.
#
# Note: 7 = 2^3 - 1
#       2 = 2^1
#       5 = 2^2 + 1
#       10 = 2^3 + 2
#
# Is there a recursive relation?
# a(n) = a(n-1) + ...

print('\nTEST 6: EIGENVALUES OF FRACTAL LAPLACIAN?')
# On a gasket or carpet?
# Just checking if specific spectral numbers match 0, 1.75...

print('\nTEST 7: "FIBONACCI" BIT OPERATOR')
# User mentioned "Fractal Information Nadsoliton".
# Maybe related to the Golden Mean (phi) in base 2?
# Or just Fibonacci numbers: 0, 1, 1, 2, 3, 5, 8, 13, 21, 34...
#
# Targets: 0, 7, 9, 14, 24.
#
# 0 = F_0
# 7 = 5 + 2 = F_5 + F_3
# 9 = 8 + 1 = F_6 + F_2
# 14 = 13 + 1 = F_7 + F_2
# 24 = 21 + 3 = F_8 + F_4
#
# Looks like sums of Fibonacci!
# Q_n = F_{n+k} + F_{m}?
#
# Let's check indices.
# Top (0): 0
# Bottom (1): 7 = F_5 + F_3  (Indices 5,3)
# Tau (2):    9 = F_6 + F_2  (Indices 6,2)
# Muon (3):  14 = F_7 + F_2  (Indices 7,2)
# Electron (4): 24 = F_8 + F_4 (Indices 8,4)
#
# Regularity:
# Primary term: F_5, F_6, F_7, F_8... The index increases by 1!
# Secondary term: F_3, F_2, F_2, F_4... Irregular?
#
# Wait, F_2 = 1.
# 7 = 5 + 2 = F_5 + F_3 (2 is F_3)
# 9 = 8 + 1 = F_6 + F_2 (1 is F_2)
# 14 = 13 + 1 = F_7 + F_2 (1 is F_2)
# 24 = 21 + 3 = F_8 + F_4 (3 is F_4)
#
# Let's index particles by n = 0, 1, 2, 3, 4.
# Q(0) = 0?
# Q(1) (Bottom) = F_5 + F_3
# Q(2) (Tau) = F_6 + F_2
# Q(3) (Muon) = F_7 + F_2
# Q(4) (Electron) = F_8 + F_4
#
# Hypothesis: Q(n) = F_{n+4} + Correction?
# F_4=3, F_5=5, F_6=8, F_7=13, F_8=21.
#
# Observed: 0, 7, 9, 14, 24.
# Base F:   3, 5, 8, 13, 21. (Indices 4,5,6,7,8)
# Diff:    -3, +2, +1, +1, +3.
#
# This looks VERY CLOSE.
# The "diff" is small integers.
#
# Let's print this cleanly.

print('FIBONACCI HYPOTHESIS')
fib = [0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89]
# F_0=0, F_1=1, F_2=1, F_3=2, F_4=3, F_5=5...

for i, q in enumerate(target_Q):
    # Try to match q = F_k + F_j
    matches = []
    for k in range(12):
        for j in range(k+1):
            if fib[k] + fib[j] == q:
                matches.append(f'F_{k} + F_{j}')
            if fib[k] - fib[j] == q:
                matches.append(f'F_{k} - F_{j}')
    print(f'  Q={q} ({["Top", "Bottom", "Tau", "Muon", "Electron"][i]}): {matches}')

