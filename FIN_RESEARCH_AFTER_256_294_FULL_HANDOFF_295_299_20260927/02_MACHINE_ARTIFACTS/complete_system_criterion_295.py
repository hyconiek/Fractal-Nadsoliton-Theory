import itertools, json, math
from collections import Counter


def H(prob):
    s=0.0
    for p in prob:
        if p>0:
            s -= p*math.log2(p)
    return s


def conditional_entropies(K, mu=None):
    n=len(K)
    if mu is None:
        mu=[1/n]*n
    joint=[[mu[i]*K[i][j] for j in range(n)] for i in range(n)]
    px=[sum(row) for row in joint]
    py=[sum(joint[i][j] for i in range(n)) for j in range(n)]
    hxy=H([joint[i][j] for i in range(n) for j in range(n)])
    hx=H(px); hy=H(py)
    return {"H_out_given_in": hxy-hx, "H_in_given_out": hxy-hy}


def reset_kernel(q):
    return [[1/q]*q for _ in range(q)]


def identity_kernel(n):
    return [[1.0 if i==j else 0.0 for j in range(n)] for i in range(n)]


def swap_kernel(q):
    n=q*q
    K=[[0.0]*n for _ in range(n)]
    for x in range(q):
        for y in range(q):
            i=x*q+y
            j=y*q+x
            K[i][j]=1.0
    return K


def one_site_swap_kernel(q):
    # Neighbor is unobserved and initially uniform.
    return reset_kernel(q)


def tape_step(state,q,L):
    x,tape,p=state
    tape=list(tape)
    y=tape[p]
    tape[p]=x
    return (y,tuple(tape),(p+1)%L)


def tape_inverse(state,q,L):
    x,tape,p=state
    p0=(p-1)%L
    tape=list(tape)
    old_x=tape[p0]
    tape[p0]=x
    return (old_x,tuple(tape),p0)


def verify_tape_bijection(q,L):
    states=[(x,t,p) for x in range(q) for t in itertools.product(range(q), repeat=L) for p in range(L)]
    images=[tape_step(s,q,L) for s in states]
    inv_ok=all(tape_inverse(tape_step(s,q,L),q,L)==s for s in states)
    return len(states), len(set(images)), inv_ok


def tape_visible_trajectory(q,L,tape,x0=0,horizon=None):
    if horizon is None:
        horizon=L+1
    state=(x0,tuple(tape),0)
    out=[]
    for _ in range(horizon):
        state=tape_step(state,q,L)
        out.append(state[0])
    return tuple(out)


def trajectory_distribution(q,L,h,x0=0):
    c=Counter()
    tapes=list(itertools.product(range(q), repeat=L))
    for tape in tapes:
        c[tape_visible_trajectory(q,L,tape,x0,h)] += 1
    total=len(tapes)
    return {k:v/total for k,v in c.items()}


def iid_distribution(q,h):
    p=1/(q**h)
    return {seq:p for seq in itertools.product(range(q), repeat=h)}


def tv(d1,d2):
    keys=set(d1)|set(d2)
    return 0.5*sum(abs(d1.get(k,0)-d2.get(k,0)) for k in keys)


def main():
    q=3; L=5
    reset=conditional_entropies(reset_kernel(q))
    swap_full=conditional_entropies(swap_kernel(q))
    swap_one=conditional_entropies(one_site_swap_kernel(q))
    nstates,nimages,inv_ok=verify_tape_bijection(q,L)
    horizon=[]
    for h in range(1,L+2):
        d=trajectory_distribution(q,L,h)
        iid=iid_distribution(q,h)
        horizon.append({
            "h":h,
            "support_tape":len(d),
            "support_iid":len(iid),
            "tv_to_iid_reset":tv(d,iid),
            "max_abs_prob_error":max(abs(d.get(k,0)-iid.get(k,0)) for k in set(d)|set(iid)),
        })
    alpha_rows=[]
    for alpha in [0,0.01,0.1,0.5,1.0]:
        alpha_rows.append({
            "alpha":alpha,
            "fresh_bits_rate_per_site_per_rho":alpha*math.log2(q),
            "two_sided_reversibility_defect_rate_per_site_per_rho":2*alpha*math.log2(q),
        })
    out={
        "q":q,
        "entropy_defects_bits":{
            "single_site_full_reset":reset,
            "full_two_site_swap":swap_full,
            "one_site_projection_of_swap_with_uniform_neighbor":swap_one,
        },
        "finite_reversible_tape":{
            "L":L,
            "state_count":nstates,
            "image_count":nimages,
            "bijective": nstates==nimages and inv_ok,
            "trajectory_tests":horizon,
        },
        "alpha_family_conditioned_on_event_type":alpha_rows,
    }
    print(json.dumps(out,indent=2,sort_keys=True))

if __name__=='__main__':
    main()
