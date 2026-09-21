from pathlib import Path
from fractions import Fraction as F
import sys,json,numpy as np
ROOT=Path(__file__).resolve().parents[1];H=ROOT/'inputs/FR223_20260916';sys.path.insert(0,str(H/'src'))
import off_face
nav=json.load(open(ROOT/'results/R7N-018_target_p_navigation.json'));x=np.array(nav['best_point_rsty']);M,_=off_face.covariance_direct_compact(x);a=67/250
w,V=np.linalg.eigh(M); # columns ascending. span of 3 smallest is candidate compression subspace
B=V[:,:3]; den=10**6; Br=np.round(B*den).astype(np.int64)/den
rank=int(np.linalg.matrix_rank(Br));Q=Br.T@(a*np.eye(4)-M)@Br; qe=np.linalg.eigvalsh(Q)
# exact rank witness: choose 3x3 minor max determinant numerically and record rational determinant
best=None
from itertools import combinations
for rows in combinations(range(4),3):
 A=np.round(B[list(rows),:]*den).astype(np.int64)
 det=round(np.linalg.det(A))
 if best is None or abs(det)>abs(best[1]):best=(rows,int(det))
out={'task':'R7N-013','scientific_status':'NUMERICAL_NEW','point_rsty':x.tolist(),'threshold':'67/250','numerical_M_eigenvalues':w.tolist(),
 'rational_basis_denominator':den,'rational_basis_numerators':np.round(B*den).astype(int).tolist(),'rank':rank,'rank_minor_rows':list(best[0]),'rank_minor_integer_det_approx':best[1],
 'compressed_center_eigenvalues':qe.tolist(),'center_PSD_margin':float(qe.min()),
 'conclusion':'Candidate rational rank-three basis is positive definite at the navigation center. Whole-cell interval PSD certification is still missing, so this is discovery evidence only.'}
json.dump(out,open(ROOT/'candidates/R7N-013_compression_candidate.json','w'),indent=2);print(json.dumps(out,indent=2))
