import itertools,re
from pathlib import Path
BASE=Path(__file__).resolve().parent.parent/'continuation_03'
# orbit result map from all 29 D12 representatives already computed
rep=[]
for line in open(BASE/'d12_four_subset_orbits.tsv'):
    mask,size,pts=line.strip().split('\t')
    mask=int(mask);size=int(size);pts=tuple(map(int,pts.split(',')))
    txt=open(BASE/f'impurity_{mask}.txt').read()
    gap=float(re.search(r'gap ([0-9eE+.-]+)',txt).group(1))
    d=4 if gap>12.55 else 3
    rep.append((set(pts),d,size))

def transforms(S):
    out=[]
    for eps in (1,-1):
      for a in range(12): out.append(frozenset((eps*x+a)%12 for x in S))
    return out
label={}
for S,d,size in rep:
    orb=set(transforms(S))
    assert len(orb)==size
    for T in orb:
        if T in label: assert label[T]==d
        label[T]=d
assert len(label)==495

def modular_exception(S):
    counts=[sum(x%4==r for x in S) for r in range(4)]
    transversal=(counts==[1,1,1,1])
    coset3=(len({x%3 for x in S})==1)
    return transversal and not coset3

bad=[]; n4=0
for S in map(frozenset,itertools.combinations(range(12),4)):
    pred=4 if modular_exception(S) else 3
    got=label[S]
    n4 += (got==4)
    if pred!=got: bad.append((sorted(S),pred,got))
print('total',len(label),'d4',n4,'predicted_d4',sum(modular_exception(S) for S in label))
print('mismatches',len(bad),bad[:5])
print('count_formula',3**4-3)
assert not bad and n4==78==3**4-3
