from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).parents[1]/"PHYS-002"))
import safe_seed as s
_,_,_,_,A=s.build_rank7()
for N in (1,2,3):
 st=s.compositions(N); pi=s.count_pi(st,0.2,A)
 print(N,{k:s.static_S(st,pi,N,k) for k in range(1,7)})
