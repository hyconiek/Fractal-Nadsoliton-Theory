# IF-003: an actual entropy-improving same-mean occupation direction

Date: 2026-10-20. The comparison is in the enlarged diagnostic constraint
class, NOT a same-Gibbs counterexample to HYP or a change to primary FIN.

## 1. Exact feasible direction, not a multiplier-only obstruction

At N2, q=12, let c_j=cos(2pi j/12), and define

    p_j(t)=p_j-t c_j,
    pi_(a,b)(t)=p_(b-a)(t)/12.

For 0<t<min p this is positive and sums to one. Fourier orthogonality
gives sum_j E_j c_j=0, so its stationary mean pair interaction is EXACTLY
the same as that of p, not just equal to first order. It preserves common
label translations, reflection and slot exchange. Legal one-slot support
and the one-self reference stay unchanged. The new maximizing process
P*(pi(t)) exists by the already paid finite scaling theorem.

Let h(t)=h*(pi(t)) be its full LABELLED entropy rate per intrinsic step.
The derivative in IF-001 gives

    h'(0)=sum_j c_j f_j=6 C_1(f)>0.                  (1)

Thus a real legal stationary history has greater entropy at exactly the
same mean interaction. This proves failure of the proposed one-mean
reduction without supposing that observed histories were Markov beforehand.
The occupation itself changes: this is not the old same-equilibrium
activity freedom and is not advertised as such.

## 2. A paid finite displacement, not only an infinitesimal direction

Freeze g=12/Lambda6 as in IF-002. Since |gE_j|<10/3,

    min p >=1/[12 exp(20/3)]>1/(12*29^2)>p0,
    p0=1/12000.

For 0<=t<=p0/2, min p(t)>p0/2, so every labelled adjacent occupation
ratio is at most 2/p0=24000. Apply the actual IS Hamming theorem for
N=2,q=12 to these perturbed weights, not only to Gibbs weights:

    B=1+22[24000+23999*11]=6335759,
    every P*(pi(t))(x,x)>=1/B.

For v=pi'(t), ||v||_1=sum_j|c_j|<=12 and
||v/pi(t)||_infinity<=2/p0. Equation IF-001(6) therefore gives on the
whole interval

    |h''(t)|<=24(B+1)/p0=1824698880000<2*10^12.      (2)

By IF-002(8), h'(0)>1/100. Set the exact rational

    t_*=1/10^16.

This lies inside [0,p0/2]. Taylor's integral remainder and (2) imply

    h(t_*)-h(0)>t_*/100-10^12 t_*^2
                  =99/10^20>1/(2*10^18)>0.          (3)

The deliberately tiny bound is a finite conditional existence certificate,
not a practical effect size or a measured entropy difference. The full
transition matrix of this positive history is defined exactly by scaling
pi(t_*); no approximate drift equality or numerical root is substituted.

## 3. What the constructed object does and does not provide

We have one explicit type-compatible competitor with the same supplied
mean interaction and legal reference, and a proved strict entropy advantage.
Therefore that scalar mean cannot independently recover the entire Gibbs
occupation premise of the authorized history principle. Any successful
reduction would need genuinely additional justified constraints, not a
retuned multiplier or a field reconstructed backwards from target pi.

This does not disprove Gibbs statistics, maximum entropy with the FULL
occupation constraint, arbitrary richer physical principles, or the
existence of an intrinsic source. It does not supply a physical conservation
law, units, clock, space, observer or QM/GR/ToE. No comparison law is adopted.
