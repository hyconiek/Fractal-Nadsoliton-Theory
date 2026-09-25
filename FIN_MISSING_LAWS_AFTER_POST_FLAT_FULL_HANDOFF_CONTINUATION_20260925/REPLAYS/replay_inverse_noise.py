import math
# Algebraic example of the sufficient bound.
def eps_max(sigma,wmin): return wmin*sigma*sigma/(4+wmin*sigma)
def deltaR(eps,sigma): return 2*eps/(sigma*(sigma-eps))
for sigma,w in [(0.5,0.2),(1.0,0.1),(2.0,0.05)]:
    em=eps_max(sigma,w); e=0.9*em
    print('sigma',sigma,'wmin',w,'eps_sufficient',em,'test_deltaR',deltaR(e,sigma),'target',w/2)
    assert deltaR(e,sigma)<w/2
