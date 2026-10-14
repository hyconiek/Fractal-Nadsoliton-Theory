"""Numerical full hidden memory integral M0,M1 and reduced comparison.
No formal tail certificate; explicit CG residual and spectral caveat.
"""
import numpy as np
from scipy.sparse.linalg import cg,LinearOperator
from scipy.linalg import solve
from model import *

def run():
    outputs={}
    for name,A in models().items():
        m=build(A);n=N;S=m['S'];U=m['U'];c=np.asarray(m['C'][:,0].toarray()).ravel()
        def q(v):return v-U@(U.T@v)
        def mv(v):return -q(S@q(v))
        H=LinearOperator((N,N),matvec=mv,dtype=float)
        it=[0]
        x,info=cg(H,c,rtol=1e-12,atol=0.,maxiter=3000,callback=lambda x:it.__setitem__(0,it[0]+1))
        res=float(np.linalg.norm(mv(x)-c)); M0=float(c@x);M1=float(x@x)
        K0=float(c@c);Kprime=float(c@(S@c))
        # comparison selected 27-dimensional channel/preparation-specific subspace
        B,_=channel_conditioned_basis(m,'conditional_00',6)
        Hhat=-B.T@(S@B);cc=B.T@c
        xhat=np.linalg.solve(Hhat,cc)
        outputs[name]={'unscaled_resolved_channel':'a_00 = 1 in orthonormal pi-resolved basis','K0_00':K0,'Kprime0_00':Kprime,'M0_integral_00_numeric':M0,'M1_first_integral_moment_00_numeric':M1,'CG_iterations':it[0],'CG_status':int(info),'CG_residual_l2':res,'CG_solution_orthogonal_to_U':float(np.linalg.norm(U.T@x)),'reduced_hidden_rank':int(B.shape[1]),'reduced_M0_integral_00_numeric':float(cc@xhat),'reduced_M1_integral_00_numeric':float(xhat@xhat),'reduced_H_min_eig':float(np.linalg.eigvalsh(Hhat)[0]),'numerical_only':True,'missing':'rigorous global H_min gap bound and numerical tail enclosure for integrated memory'}
    z=build(A7(),0.);outputs['FIN_A7_g0']={'K0_00':float(np.linalg.norm(z['C'][:,0].toarray())**2),'exact_expected':0.}
    write('RG-D03/MEMORY_REDUCTION.json',outputs)
    print(outputs)
if __name__=='__main__':run()
