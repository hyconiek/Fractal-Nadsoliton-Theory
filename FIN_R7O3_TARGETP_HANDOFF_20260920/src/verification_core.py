from __future__ import annotations
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,copy

def sha_bytes(b:bytes): return hashlib.sha256(b).hexdigest()
def sha_file(p): return sha_bytes(Path(p).read_bytes())
def canonical_hash(obj): return sha_bytes(json.dumps(obj,sort_keys=True,separators=(',',':')).encode())
def load_jsonl(p): return [json.loads(x) for x in Path(p).read_text().splitlines() if x.strip()]
def C(cell): return tuple((F(a),F(b)) for a,b in cell)

def cert_payload_hash(c): return canonical_hash({k:v for k,v in c.items() if k!='certificate_sha256'})
def tree_record_hash(r): return canonical_hash({k:v for k,v in r.items() if k!='tree_sha256'})

def _det3(A):
 return A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])-A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])+A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0])
def _gram(c):
 den=int(c['basis_den']); Bn=c['basis_num'];B=[[F(int(Bn[i][j]),den) for j in range(3)] for i in range(4)]
 return [[str(sum(B[r][a]*B[r][b] for r in range(4))) for b in range(3)] for a in range(3)]

def validate_cert(c):
 e=[]
 req=['leaf_id','original_index','original_parent_sha256','cell','threshold','basis_num','basis_den','center_c','rank_minor_rows','rank_minor_exact','gram_matrix_exact','pd_method','pd_bounds_exact','certificate_sha256']
 for k in req:
  if k not in c:e.append('missing:'+k)
 if e:return e
 try:
  if c['threshold']!='67/250':e.append('threshold')
  if len(c['basis_num'])!=4 or any(len(r)!=3 for r in c['basis_num']):e.append('basis_shape')
  if int(c['basis_den'])<=0:e.append('basis_den')
  if len(c['center_c'])!=3 or any(F(x).denominator<=0 for x in c['center_c']):e.append('center_c')
  rows=tuple(c['rank_minor_rows'])
  if len(rows)!=3 or len(set(rows))!=3 or any(r not in range(4) for r in rows):e.append('rank_rows')
  else:
   Bn=c['basis_num'];A=[[int(Bn[r][j]) for j in range(3)] for r in rows];detq=F(_det3(A),int(c['basis_den'])**3)
   if detq==0:e.append('rank_zero')
   if str(detq)!=c['rank_minor_exact']:e.append('rank_evidence')
  if _gram(c)!=c['gram_matrix_exact']:e.append('gram')
  pb=c['pd_bounds_exact']
  if c['pd_method']=='SYLVESTER':
   if not all(F(pb[k][0])>0 for k in ['d1','d2','d3']):e.append('pd_sign')
  elif c['pd_method']=='GERSHGORIN':
   if not all(F(x)>0 for x in pb['gersh_lower_bounds']):e.append('pd_sign')
  else:e.append('pd_method')
  if cert_payload_hash(c)!=c['certificate_sha256']:e.append('certificate_hash')
 except Exception as ex:e.append('exception:'+type(ex).__name__+':'+str(ex))
 return e

def _expected_children(cell,axis,split):
 c=list(C(cell));lo,hi=c[axis];m=F(split)
 if not lo<m<hi: raise ValueError('split_outside')
 L=list(c);R=list(c);L[axis]=(lo,m);R[axis]=(m,hi);return tuple(L),tuple(R)
def check_tree_node(node,parent_cell,leafmap,orig_index):
 errs=[];lids=[]
 try:
  cell=C(node['cell'])
  if cell!=C(parent_cell):errs.append('node_cell_mismatch:'+str(node.get('node_id')))
  kind=node.get('kind')
  if kind=='SAFE':
   lid=node.get('leaf_id')
   if lid not in leafmap:return errs+['missing_leaf:'+str(lid)],lids
   c=leafmap[lid];lids.append(lid)
   if c['original_index']!=orig_index:errs.append('leaf_parent_index')
   if C(c['cell'])!=cell:errs.append('leaf_cell_mismatch')
  elif kind=='SPLIT':
   ax=int(node['axis']);sp=node['split'];L,R=_expected_children(cell,ax,sp)
   el,ll=check_tree_node(node['left'],L,leafmap,orig_index);er,rr=check_tree_node(node['right'],R,leafmap,orig_index)
   errs+=el+er;lids+=ll+rr
   if set(ll)&set(rr):errs.append('duplicate_leaf_across_children')
  else:errs.append('bad_kind:'+str(kind))
 except Exception as ex:errs.append('exception:'+type(ex).__name__+':'+str(ex))
 return errs,lids

def check_tree_record(rec,parent,leafmap):
 errs=[]
 if rec['original_index']!=parent['original_index']:errs.append('index')
 if rec.get('original_path')!=parent['original_path']:errs.append('path')
 if tree_record_hash(rec)!=rec.get('tree_sha256'):errs.append('tree_hash')
 e,lids=check_tree_node(rec['tree'],C(parent['original_cell']),leafmap,parent['original_index']);errs+=e
 return errs,lids

def validate_registry(reg):
 e=[];ps=reg.get('parents',[]);ids=[p.get('original_index') for p in ps]
 if len(ps)!=5432:e.append('parent_count')
 if len(set(ids))!=len(ids):e.append('duplicate_parent_index')
 if set(ids)!=set(range(5432)):e.append('missing_or_extra_parent_index')
 if reg.get('certified_closed')!=5432 or any(reg.get(k,0)!=0 for k in ['partial','unprocessed','invalid','counterexample']):e.append('state_counts')
 return e

def formula_record_match(cert,out):
 return out.get('ok') and all(cert.get(k)==out.get(k) for k in ['gram_matrix_exact','internal_chart_bounds','moment_entry_enclosures','pd_bounds_exact'])
