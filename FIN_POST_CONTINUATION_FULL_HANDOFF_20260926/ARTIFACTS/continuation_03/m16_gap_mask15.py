import re, ast, functools, numpy as np, math
src=open('/mnt/data/FIN_NEXT_RESEARCH_20260925/continuation_02/fin12_m32_one_fast.cpp').read()
m=re.search(r'static const double G\[12\]\[12\]=\{(.*?)\};',src,re.S)
s='['+m.group(1).replace('{','[').replace('}',']')+']'
G=np.array(ast.literal_eval(s),float)
root=tuple(2 if i<4 else 1 for i in range(12))
TOL=1e-10
@functools.lru_cache(None)
def splits(c):
 h=sum(c)//2;out=[];a=[0]*12
 def rec(i,rem):
  if i==12:
   if rem==0:
    aa=tuple(a);bb=tuple(c[j]-aa[j] for j in range(12))
    if aa<=bb:out.append((aa,bb))
   return
  tail=sum(c[i+1:]);
  for x in range(max(0,rem-tail),min(c[i],rem)+1):a[i]=x;rec(i+1,rem-x)
  a[i]=0
 rec(0,h);return out
@functools.lru_cache(None)
def norm(c):
 x=np.array(c,float);return float(x@G@x)
def delta(a,b,c):return (2*norm(a)+2*norm(b)-norm(c))/sum(c)
def merge_levels(items):
 items=sorted(items,key=lambda x:x[0]);out=[]
 for e,g in items:
  if out and abs(e-out[-1][0])<TOL:out[-1]=(out[-1][0],out[-1][1]+g)
  else:out.append((e,g))
  if len(out)>=3 and e>out[1][0]+TOL:break
 return tuple(out[:2])
@functools.lru_cache(None)
def levels(c):
 if sum(c)==1:return ((0.0,1),)
 cand=[]
 for a,b in splits(c):
  A=levels(a);B=levels(b);d=delta(a,b,c)
  if a!=b:
   for ea,ga in A:
    for eb,gb in B:cand.append((d+2*(ea+eb),ga*gb))
  else:
   for i,(ea,ga) in enumerate(A):
    cand.append((d+4*ea,ga*(ga+1)//2))
    for j in range(i+1,len(A)):
     eb,gb=A[j];cand.append((d+2*(ea+eb),ga*gb))
 return merge_levels(cand)
L=levels(root)
print('root',root,'levels',L,'gap',L[1][0]-L[0][0],'states',levels.cache_info(),'splits',splits.cache_info())
