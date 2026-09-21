from pathlib import Path
import sys
sys.set_int_max_str_digits(0)
import json,time,hashlib
from fractions import Fraction as F
ROOT=Path('/mnt/data/fin_rank7_next_campaign')
H=ROOT/'inputs/FR223_20260916'
sys.path.insert(0,str(H/'src'))
import frontier_shifted_boxes as fsb
ledger=json.load(open(H/'results/FR223_ACTIVE_MASK_LEDGER.json'))
fail={'FR32','FR48','FR52','FR54','FR56','FR58','FR60'}
by={m['name']:m for m in ledger['shifted_masks']}

def box_from(m):
    return [(F(str(m['x'][0])),F(str(m['x'][1]))),(F(str(m['u'][0])),F(str(m['u'][1]))),(F(str(m['v'][0])),F(str(m['v'][1]))),(F(0),F(str(m['e'])))]

def jsI(I): return [str(I.lo),str(I.hi)]
def res_json(r):
    return {'status':r['status'],'reason':r['reason'],'rstar':jsI(r['rstar']),'P':jsI(r['P']),'P1':jsI(r['P1']),'c2':jsI(r['c2'])}
def endpoint(b):return [[str(a),str(c)] for a,c in b]

def certify_recursive(name,b,maxdepth=4):
    calls=0
    def rec(cell,depth,path):
        nonlocal calls
        calls+=1;t=time.time();r=fsb.raw_shifted_box(cell);dt=time.time()-t
        print(name,path,'depth',depth,r['status'],r['reason'],'secs',round(dt,3),flush=True)
        if r['status']=='INTERVAL_CERTIFIED':
            return {'kind':'leaf','path':path,'box':endpoint(cell),'checker_result':res_json(r),'runtime_s':dt}
        if depth>=maxdepth:
            return {'kind':'unresolved','path':path,'box':endpoint(cell),'checker_result':res_json(r),'runtime_s':dt}
        # choose largest normalized width among x,u,v; e is relaxation slab and kept intact
        widths=[float(c-a) for a,c in cell[:3]]
        ax=max(range(3),key=lambda i:widths[i])
        a,c=cell[ax];mid=(a+c)/2
        left=list(cell);right=list(cell);left[ax]=(a,mid);right[ax]=(mid,c)
        return {'kind':'split','path':path,'axis':ax,'split':str(mid),'box':endpoint(cell),
                'failed_parent':res_json(r),'children':[rec(left,depth+1,path+'0'),rec(right,depth+1,path+'1')]}
    tree=rec(b,0,'')
    return tree,calls

def leaves(tree):
    if tree['kind']!='split': return [tree]
    return sum((leaves(c) for c in tree['children']),[])

existing_path=ROOT/'results/fresh_repair_trees.json'
if existing_path.exists():
    out=json.load(open(existing_path))
else:
    out={'threshold':'sigma (historical strict interval provider)','checker':'inputs/FR223_20260916/src/frontier_shifted_boxes.py::raw_shifted_box','source_ledger':'inputs/FR223_20260916/results/FR223_ACTIVE_MASK_LEDGER.json','repairs':{}}
selected=[x for x in sys.argv[1:] if x in fail] or sorted(fail,key=lambda s:int(s[2:]))
for name in selected:
    tree,calls=certify_recursive(name,box_from(by[name]),4)
    ls=leaves(tree); cert=[x for x in ls if x['kind']=='leaf']; unr=[x for x in ls if x['kind']=='unresolved']
    out['repairs'][name]={'tree':tree,'calls':calls,'certified_leaf_count':len(cert),'unresolved_leaf_count':len(unr)}
    tmp=ROOT/'certificates'/f'{name}_fresh_repair.json';tmp.parent.mkdir(exist_ok=True)
    json.dump(out['repairs'][name],open(tmp,'w'),indent=2)
json.dump(out,open(ROOT/'results/fresh_repair_trees.json','w'),indent=2)
print('SUMMARY',[(k,v['certified_leaf_count'],v['unresolved_leaf_count'],v['calls']) for k,v in out['repairs'].items()],flush=True)
