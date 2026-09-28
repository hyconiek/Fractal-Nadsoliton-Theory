import numpy as np, json, time

def analyze(N):
 z=np.load(f'/mnt/data/defect_research_inputs/base_N{N}.npz',allow_pickle=True);st=z['states'];lab=z['labels'];pi=z['pi']; idx={tuple(map(int,s)):i for i,s in enumerate(st)};j0=np.where(lab==0)[0];loc={int(g):a for a,g in enumerate(j0)};parent=np.arange(len(j0));sz=np.ones(len(j0),int)
 def find(a):
  while parent[a]!=a:parent[a]=parent[parent[a]];a=parent[a]
  return a
 def union(a,b):
  a=find(a);b=find(b)
  if a==b:return
  if sz[a]<sz[b]:a,b=b,a
  parent[b]=a;sz[a]+=sz[b]
 for aa,gidx in enumerate(j0):
  n=st[gidx]
  for i in np.where(n>0)[0]:
   for jj in range(12):
    if jj==i:continue
    n2=n.copy();n2[i]-=1;n2[jj]+=1;h=idx[tuple(map(int,n2))]
    if lab[h]==0:union(aa,loc[int(h)])
 roots=np.array([find(i) for i in range(len(j0))]);uniq=np.unique(roots);seedg=np.where((st[:,0]==N)&(st[:,1:].sum(1)==0))[0][0];seedroot=roots[loc[int(seedg)]]
 res={'N':N,'components':len(uniq),'component_sizes':sorted([int(np.sum(roots==r)) for r in uniq],reverse=True)[:10]}
 for kap in [0.,12.]:
  lw=np.log(pi[j0])+kap*st[j0,0]/N;lw-=lw.max();w=np.exp(lw);w/=w.sum();masses=[float(w[roots==r].sum()) for r in uniq];res[f'kappa_{kap:g}']={'seed_component_mass':float(w[roots==seedroot].sum()),'largest_component_mass':max(masses),'outside_seed_component_mass':1-float(w[roots==seedroot].sum()),'top_component_masses':sorted(masses,reverse=True)[:10]}
 return res
out=[]
for N in range(3,11):
 t=time.time();r=analyze(N);r['seconds']=time.time()-t;out.append(r);print(json.dumps(r),flush=True)
open('/mnt/data/defect_research_inputs/PREP_COMPONENTS_330.json','w').write(json.dumps({'rows':out},indent=2))
