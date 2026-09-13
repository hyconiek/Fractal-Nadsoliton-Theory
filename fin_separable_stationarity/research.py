"""New goal: separable stationary strict marginals and quantum discord.

Preparation is explicit but uses a supplied target-encoded program state.
No physical apparatus, intrinsic source, or fundamental closure is claimed.
"""
from fractions import Fraction as F
import importlib.util
import itertools
import json
from pathlib import Path

import numpy as np


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('quantum_lift',ROOT/'fin_quantum_lift/research.py')
q=importlib.util.module_from_spec(spec);spec.loader.exec_module(q)


def hadamard12():
    residues={k*k%11 for k in range(1,11)}
    H=np.ones((12,12),dtype=int)
    for i in range(11):
        for j in range(11):H[i+1,j+1]=-1 if i==j else (1 if (j-i)%11 in residues else -1)
    assert np.array_equal(H@H.T,12*np.eye(12,dtype=int))
    assert np.array_equal(H[:,1:].sum(axis=0),np.zeros(11,dtype=int))
    return H


def cut_certificate(row_permutation=None):
    H=hadamard12()
    if row_permutation is not None:
        if sorted(row_permutation)!=list(range(12)):raise ValueError('A permutation of all twelve labels is required')
        H=H[row_permutation]
    cuts=np.concatenate([H[:,1:].T,-H[:,1:].T],axis=0)
    incidence=(cuts+1)//2
    assert np.all(incidence.sum(axis=1)==6)
    assert np.all(incidence.sum(axis=0)==11)
    assert np.array_equal(incidence.T@incidence,6*np.eye(12,dtype=int)+5*np.ones((12,12),dtype=int))
    return cuts,dict(oriented_branches=22,unoriented_cuts=11,block_size=6,
        point_replication=11,pair_replication=5,coherence_retention='5/11',
        exact_hadamard_and_incidence_identities=True,
        optimality_scope='At least 11 rank-one cuts are needed for this isotropic balanced-cut second moment; no universal preparation-cost optimality.')


def stationary_separable(C,row_permutation=None):
    n=len(C)
    if (n!=12 or not np.allclose(C,C.conj().T) or not np.isfinite(C).all()
            or not np.allclose(np.diag(C),np.ones(n)/n)):
        raise ValueError('The declared Hermitian 12-label uniform-diagonal class is required')
    retention=5/11;program=np.eye(n)/n+(C-np.eye(n)/n)/retention
    if np.linalg.eigvalsh(program)[0]<-1e-12:raise ValueError('The amplified program is not positive')
    cuts,certificate=cut_certificate(row_permutation);R=np.zeros((n*n,n*n),complex);branches=[]
    for cut in cuts:
        a=(cut+1)/2;b=1-a
        left=2*(a[:,None]*program*a[None,:]);right=2*(b[:,None]*program*b[None,:])
        branches.append((left,right));R+=np.kron(left,right)/len(cuts)
    return R,program,branches,certificate


def positivity_certificate():
    spec=importlib.util.spec_from_file_location('strict_spectrum',ROOT/'fin_projected_learning/research.py')
    previous=importlib.util.module_from_spec(spec);spec.loader.exec_module(previous)
    source=previous.certify_strict_spectrum();eig=[tuple(map(F,p)) for p in source['eigenvalue_intervals']]
    gamma=F(1,20);amplified=gamma*F(11,5)
    minimum=min(F(1,12)+amplified*lo for lo,hi in eig)
    assert minimum>F(1,125)
    assert eig[0][1]<F(5,3)
    ell=(-eig[6][1],-eig[6][0])
    guaranteed_loading=(F(5,11)/(12*ell[1]),F(5,11)/(12*ell[0]))
    return dict(output_loading=str(gamma),program_loading=str(amplified),
        exact_program_eigenvalue_lower=str(minimum),simple_program_floor='1/125',
        off_coincidence_state_eigenvalue_lower=str(F(12,11)*minimum**2),
        guaranteed_separable_loading_interval=[str(x) for x in guaranteed_loading],
        safe_loading_upper=str(guaranteed_loading[0]),
        loading_threshold_scope='Sufficient for this construction, not a separability boundary.',
        exact_joint_coincidence_probability='0',independent_joint_coincidence_probability='1/12')


def discord_block(n=12):
    V,S,D,P=q.structures(n)
    u=np.ones(n)/np.sqrt(n);v=(-1.)**np.arange(n)/np.sqrt(n)
    block=np.einsum('i,iajb,j->ab',u,V.reshape(n,n,n,n),v)
    J=np.diag((-1.)**np.arange(n))
    expected=(np.outer(u,v)+np.outer(v,u))/2-J/n
    assert np.linalg.norm(block-expected)<1e-13
    return block,u,v


def frame_ambiguity(C):
    permutation=list(range(12));permutation[1],permutation[2]=permutation[2],permutation[1]
    cuts0,_=cut_certificate();cuts1,_=cut_certificate(permutation)
    inc0=(cuts0+1)//2;inc1=(cuts1+1)//2
    witness=None
    for i,j,k,l in itertools.permutations(range(12),4):
        count0=int(np.sum(inc0[:,i]*inc0[:,j]*(1-inc0[:,k])*(1-inc0[:,l])))
        count1=int(np.sum(inc1[:,i]*inc1[:,j]*(1-inc1[:,k])*(1-inc1[:,l])))
        if count0!=count1:
            witness=(i,j,k,l,count0,count1);break
    assert witness is not None
    moment0=np.einsum('ri,rj,rk,rl->ijkl',cuts0,cuts0,cuts0,cuts0,optimize=True)
    moment1=np.einsum('ri,rj,rk,rl->ijkl',cuts1,cuts1,cuts1,cuts1,optimize=True)
    equivalent=False
    for sign in [1,-1]:
        for shift in range(12):
            p=(sign*np.arange(12)+shift)%12
            equivalent|=np.array_equal(moment0[np.ix_(p,p,p,p)],moment1)
    assert not equivalent
    R0,program,_,_=stationary_separable(C);R1,_,_,_=stationary_separable(C,permutation)
    i,j,k,l,count0,count1=witness
    expected=4*(count1-count0)/22*program[i,j]*program[k,l]
    assert abs((R1-R0)[12*i+k,12*j+l]-expected)<1e-14
    return dict(second_frame_row_permutation=permutation,four_label_witness=[i,j,k,l],
        oriented_counts=[count0,count1],exact_entry_difference_coefficient=str(F(2*(count1-count0),11)),
        example_joint_entry_difference=float(expected.real),
        dihedral_fourth_moment_equivalent=bool(equivalent),
        different_fourth_moments_proved=True,
        scope='Same labelled kernel, cut count and first/two-point cut statistics; joint predictions still require the frame/higher moments.')


def run():
    n=12;W=q.strict();C=np.eye(n)/n+.05*W
    R,program,branches,cuts=stationary_separable(C);V,S,D,P=q.structures(n)
    marginal=q.marginal(R,n);partial=q.partial_transpose(R,n)
    assert np.linalg.norm(marginal-C)<1e-13
    assert np.linalg.norm(V@R-R@V)<1e-13
    assert np.linalg.norm(D@R)<1e-13 and np.linalg.norm(S@R-R@S)<1e-13
    assert np.linalg.norm(partial-R)<1e-13
    # Explicit local projective instrument: the success probability is 1/2.
    a=(hadamard12()[:,1]+1)/2;b=1-a
    K1=np.diag(np.kron(a,b));K2=np.diag(np.kron(b,a))
    source=np.kron(program,program)
    success=K1@source@K1+K2@source@K2
    assert abs(np.trace(success)-.5)<1e-13
    expected=(np.kron(branches[0][0],branches[0][1])+np.kron(branches[11][0],branches[11][1]))/2
    assert np.linalg.norm(2*success-expected)<1e-13
    block,u,v=discord_block()
    conditional_comm=R@np.kron(C,np.eye(n))-np.kron(C,np.eye(n))@R
    distance=sum(abs(np.linalg.eigvalsh(R-np.kron(C,C))))/2
    state_eig=np.linalg.eigvalsh(R)
    return dict(status='Explicit separable stationary construction and conditional LOCC preparation; canonical zero-discord obstruction proved separately.',
        exact_cut_certificate=cuts,positivity_certificate=positivity_certificate(),
        numerical_checks=dict(marginal_residual=float(np.linalg.norm(marginal-C)),
            stationary_residual=float(np.linalg.norm(V@R-R@V)),
            partial_transpose_equality_residual=float(np.linalg.norm(partial-R)),
            state_rank=int(np.linalg.matrix_rank(R,tol=1e-10)),
            smallest_positive_state_eigenvalue=float(min(x for x in state_eig if x>1e-10)),
            trace_distance_to_product=float(distance),
            mutual_information_nats=2*q.entropy(C)-q.entropy(R),
            local_marginal_commutator_frobenius=float(np.linalg.norm(conditional_comm)),
            preparation_success_probability=float(np.trace(success).real),
            invertible_cross_block_singular_values=np.linalg.svd(block,compute_uv=False).tolist()),
        cut_frame_ambiguity=frame_ambiguity(C),
        interpretation=['Explicit sum of 22 product density matrices proves separability, not just PPT.',
            'The amplified program already contains W; the construction is not source-independent emergence.',
            'Stationarity concerns the unconditioned pair after the preparation label is discarded.',
            'Separable does not mean classical-quantum or zero discord.'])


if __name__=='__main__':
    data=run();Path(__file__).with_name('results.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(data,indent=2))
