import math
def h(K,L=30):
    q=math.log(K)
    for _ in range(L): q=2*q+math.log1p(math.exp(-q))-math.log(2)
    return q/(2**L)
for K in (4,8,12): print(K,h(K))
