#!/usr/bin/env python3
"""Self-contained exact Fraction checks; NOT q12 all-N certification."""
from fractions import Fraction as F
from itertools import combinations
from collections import deque
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parent

def all_sets(n):
    for mask in range(1<<n):
        yield tuple(i for i in range(n) if (mask>>i)&1)

def source(n, edge_list):
    s=[F(0)]*n
    for u,v,c in edge_list:
        assert u!=v and c>=0
        s[u]+=F(c);s[v]-=F(c)
    return s

def cut_bound(n, h_edges, s):
    assert sum(s)==0
    largest=F(0); bad=[]; maximizing=[]
    for S in all_sets(n):
        b=sum(s[i] for i in S)
        out=sum(c for u,v,c in h_edges if u in S and v not in S)
        if b>0 and not out:
            bad.append({'S':S,'source':str(b),'out_native':'0'})
            continue
        if out and b/out>largest:
            largest=b/out;maximizing=[S]
        elif out and b/out==largest:
            maximizing.append(S)
    return (None if bad else largest),bad,maximizing

def exact_maxflow(n, h_edges, s, K):
    # Edmonds-Karp on finite directed graph; source and sink are external.
    sup=n;snk=n+1; z=n+2
    cap=[[F(0) for _ in range(z)] for _ in range(z)]
    for u,v,c in h_edges:cap[u][v]+=K*c
    for x,t in enumerate(s):
        if t>0:cap[sup][x]+=t
        if t<0:cap[x][snk]+=-t
    target=sum(t for t in s if t>0);total=F(0)
    original=[row[:] for row in cap]
    while total<target:
        prev=[None]*z;prev[sup]=sup; q=deque([sup])
        while q and prev[snk] is None:
            u=q.popleft()
            for v in range(z):
                if prev[v] is None and cap[u][v]>0:
                    prev[v]=u;q.append(v)
        if prev[snk] is None:break
        path=[];v=snk
        while v!=sup:
            u=prev[v];path.append((u,v));v=u
        df=min(cap[u][v] for u,v in path)
        for u,v in path:cap[u][v]-=df;cap[v][u]+=df
        total+=df
    flows={f'{u}->{v}':str(original[u][v]-cap[u][v]) for u,v,c in h_edges}
    return total==target,flows,{'flow':str(total),'target':str(target)}

def q_sample(t):
    # Levels: 0 at level 0, 1 and 2 at level 1, 3 at level 2.
    # Phase: 0,2,3 A; 1 B. Allowed H 0->2, 2->3;
    # deleted native F 0->1, 1->3; both equal weight 1.
    h=[(0,2,F(t)),(2,3,F(t))]
    f=[(0,1,F(1)),(1,3,F(1))]
    return h,f,source(4,h+f)

def no_up_path(start,target,n,edges,levels):
    seen={start};stack=[start]
    while stack:
        u=stack.pop()
        for x,v,c in edges:
            assert levels[v]-levels[x]==1
            if x==u and v not in seen:
                if v==target:return False
                seen.add(v);stack.append(v)
    return True

def run():
    flags={}
    h,f,s=q_sample(F(1))
    levels=(0,1,1,2)
    flags['toy_native_source']=(s==[F(2),F(0),F(0),F(-2)])
    k,bad,sets=cut_bound(4,h,s)
    flags['toy_sharp_relative_cap']=(k==2 and not bad)
    feasible,flows,meta=exact_maxflow(4,h,s,k)
    flags['toy_primal_at_optimal_cap']=feasible and flows=={'0->2':'2','2->3':'2'}
    flags['toy_below_cap_fails']=not exact_maxflow(4,h,s,F(199,100))[0]
    flags['toy_no_up_detour_for_each_forbidden']=all(no_up_path(x,y,4,h,levels) for x,y,c in f)
    flags['toy_phase_balance']=s[0]+s[2]+s[3]==0 and s[1]==0
    flags['toy_level_cuts']=(sum(s[i] for i in [0])==2 and sum(s[i] for i in [0,1,2])==2)
    zero_edges=[];ss=source(2,[(0,1,F(1))]);kf,fail,sets0=cut_bound(2,zero_edges,ss)
    flags['forbidden_only_obstruction']=(kf is None and len(fail)==1 and fail[0]['S']==(0,))
    for N in (2,4,10,100):
        hN,fN,sN=q_sample(F(1,N))
        kN,badN,_=cut_bound(4,hN,sN)
        ok,flowN,_=exact_maxflow(4,hN,sN,kN)
        flags[f'varying_relative_cap_N{N}']=(kN==N+1 and ok and flowN['0->2']==str(1+F(1,N)))
        flags[f'bounded_absolute_exit_N{N}']=all(abs(z)<=F(2) for z in sN)
    # Exact demonstration: a legal final H-flow may need to DECREASE a retained native edge.
    h8=[(1,2,F(2)),(2,4,F(2)),(4,6,F(2)),(0,3,F(1)),(3,5,F(1)),(5,7,F(1))]
    f8=[(0,2,F(1)),(4,7,F(1))]
    s8=source(8,h8+f8)
    k8,bad8,_=cut_bound(8,h8,s8)
    feasible8,flow8,_=exact_maxflow(8,h8,s8,k8)
    f8source=source(8,f8)
    ka,bada,_=cut_bound(8,h8,f8source)
    flags['signed_repair_has_feasible_final_flow']= feasible8 and k8==2
    flags['additive_only_repair_is_impossible']=ka is None and any(b['S']==(4,6) for b in bada)
    flags['signed_repair_decreases_one_native_edge']= feasible8 and flow8['2->4']=='1'
    flags['empty_zero_supply']= cut_bound(2,[],[F(0),F(0)])[0]==0
    assert all(flags.values()),flags
    doc={'status':'PASS','check_count':len(flags),'flags':flags,'exact_fixtures':{
      'four_vertex':{'h':[[u,v,str(c)] for u,v,c in h], 'f':[[u,v,str(c)] for u,v,c in f], 's':[str(x) for x in s], 'K_star':str(k), 'maximizing_cuts':sets,'constructive_flow':flows},
      'negative':{'H':[],'native_deleted':[[0,1,'1']],'s':['1','-1'],'K_star':'INFINITY','failed_cuts':fail},
      'eight_vertex_signed_repair':{'h':[[u,v,str(c)] for u,v,c in h8], 'f':[[u,v,str(c)] for u,v,c in f8], 's':[str(x) for x in s8], 'K_star':str(k8), 'constructive_flow':flow8,'additive_only_K_star':'INFINITY','additive_failed_cut':[b for b in bada if b['S']==(4,6)]},
      'nonnecessary_uniform_K':{'K_star_N':'N+1','constructed_exit_upper_bound': '1+1/N <= 2 at each non-isolated vertex for weights W=1; artificial graded toy, not q12'}}}
    (ROOT/'verification.json').write_text(json.dumps(doc,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({'status':'PASS','flags':len(flags),'toy_optimal_K':str(k),'actual_FIN':'NOT_TESTED'}))
if __name__=='__main__':run()
