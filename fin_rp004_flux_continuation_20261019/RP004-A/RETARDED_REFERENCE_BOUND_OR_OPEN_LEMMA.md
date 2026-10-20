# RP004-A — exact retarded global-energy growth, exterior balance, moving-ray localization OPEN

Date: 2026-10-19. Disposition: **PARTIAL_REFERENCE_GLOBAL_ENERGY_AND_EXTERIOR_BALANCE_ONLY; MOVING_RAY_OPEN_LEMMA**. Classical conditional result in the exact frozen F_RP sparse family; no new wave law.

## 1. Frozen reference and spectral formula

On the odd real subspace, L=-∂_x²+4-6 sech²(x) >=3; the even translation eigenfunction at 0 is excluded by parity, while the odd shape eigenfunction at eigenvalue 3 is retained. Let Ω=2 sqrt(a+2), a∈[-3/4,-1/4], so Ω²∈[5,7], q=sqrt(Ω²-4)∈[1,sqrt(3)]. Set S=(9/4)tanh(x) sech⁴(x). Take precisely the prescribed retarded solution:

    u_tt+Lu=-S (1+cos Ωt),  u(0)=u_t(0)=0.

By the spectral theorem, with the removable value at λ=Ω² taken continuously,

    u(t) = -[(1-cos(t sqrt(L)))/L
            +(cos(Ωt)-cos(t sqrt(L)))/(L-Ω²)] S .        (A.1)

Here (A.1) is an operator function of L; the individually displayed second quotient is continuous at λ=Ω² for each finite t. At the resonance its value is +t sin(Ωt)/(2Ω), and the overall minus sign gives the correctly driven resonant oscillator. S is odd and exponentially localized. The continuous spectrum is [4,∞); the odd discrete shape eigenvalue 3 is included, while the translation mode is absent for this specific preparation. Formula (A.1) keeps the DC/turn-on terms and BOTH q and -q scattering channels; it does not insert a free radiation packet.

## 2. Paid strict growth of the *global* retarded reference energy

Define the quadratic (full-line) reference energy

    E_u(t)=1/2 (||u_t(t)||_2²+<u(t),Lu(t)>).

The inherited normalized generalized Fourier states and source transform M(q), with dq/(2π), yield

    Γ=|M(q)|²/(2q)>125/73728.

The spectral measure μ_S of S in the frequency variable ν=sqrt(λ) has a continuum density smooth near Ω,

    ρ_S(Ω) = (Ω/(π q)) |M(q)|² = (2Ω/π) Γ.

For the cosine-driven component alone, the exact scalar-oscillator response at frequency ν is

    b_ν(t)= -(cos Ωt-cos νt)/(ν²-Ω²),
    b'_ν(t)= (Ω sin Ωt-ν sin νt)/(ν²-Ω²).

Accordingly its energy is the integral of [b'_ν(t)²+ν² b_ν(t)²]/2 against dμ_S(ν). Substitute ν=Ω+z/t in a fixed resonant neighborhood. The leading rescaled energy integrand is

    t² sin²(z/2)/(2z²) ,

with dν=dz/t. Since ∫_{-∞}^∞ sin²(z/2)/z² dz = π/2, the resonant contribution is

    (π/4) ρ_S(Ω) t = (Ω Γ/2) t.

For completeness: split |ν-Ω|<r and its complement with a fixed r smaller than the uniform spectral-gap separation to ν=2 and ν=sqrt(3). On the complement b and b' have uniformly bounded integrable squared amplitudes (smooth rapidly decaying scattering transform at infinity). In the resonant neighborhood, the bound C min(t²,|ν-Ω|^{-2}) and smoothness of ρ justify dominated convergence after z scaling, uniformly in a in the stated compact interval. Odd shape-mode and DC responses have uniformly bounded energy (L_odd>=3 and Ω²-3>=2); by Cauchy-Schwarz their cross energy with the harmonic part is O(sqrt(t)). Hence

    E_u(t)/t -> Ω Γ/2 >0                 (A.2)

uniformly in a in [-3/4,-1/4] (the stated compact smoothness bounds are uniform). This is finite-time reference energy *growth*, not a theorem on its spatial location and not an exact nonlinear flux theorem.

## 3. Paid exterior identity and almost-nonnegative signed flux

For any R>=1 set

    E_ext(t;R)=1/2 ∫_{|x|>R} [u_t²+u_x²+(4-6 sech²x)u²] dx.

The exterior potential is strictly positive for R>=1. Direct multiplication of the actual retarded PDE by u_t, and spatial integration over both exterior half-lines, gives the exact signed formula

    F_u(R,T) := ∫_0^T [J_u(t,R)-J_u(t,-R)] dt
              = E_ext(T;R) + ∫_0^T∫_{|x|>R} S(x)(1+cos Ωt) u_t(t,x) dx dt,
    J_u=-u_tu_x.                                    (A.3)

This sign follows from e_t+J_x=-S(1+cosΩt)u_t; the initial exterior energy is zero. The global linear wave energy estimate, using L_odd>=3, yields ||u_t(t)||_2<=C t; ||S 1_{|x|>R}||_2<=C e^{-4R}. Thus

    F_u(R,T) >= -C e^{-4R} T²,   R>=1, T>=0.      (A.4)

At the frozen ray R=T/6, for T>=6,

    F_u(T/6,T) >= -C T² exp(-2T/3).              (A.5)

A.5 pays the sign of the outside forcing contribution without assuming compact source support. It **does not** give F_u>=cT, because E_ext(T;T/6) has not been bounded below proportionally to T. The positive full-line E_u(T) in (A.2) is not a spatial localization theorem.

## 4. The exact missing analytic estimate (A-GATE)

A sufficient unpaid statement is a uniform-in-a lower bound

    liminf_{T->infty} T^{-1} E_ext(T;T/6) >= c_ext>0.  (A-GATE)

A stronger and explicit target is a limiting-amplitude/turn-on remainder in a *moving cone* (x/T in [1/6,1/3], t=T) for u and first derivatives. The outgoing stationary solution has amplitude |A_+|=|M(q)|/(2q) on each half-line, but one must prove that its retarded finite-T remainder (including DC, the odd shape mode, the turn-on wave packet, both continuum roots and nonthreshold dispersive terms) is o(1) in a suitable integrated energy norm **uniformly on the moving cone**. Were this paid, the band energy limit would be

    T^{-1} E_{T/6<|x|<T/3}(T) -> Ω² Γ/(12q)>0.

This expression is a **conditional calculation for the proposed cone lemma**, not a proved bound. q/Ω >=1/sqrt(5)>1/3 puts the proposed spatial band behind the harmonic group-velocity ray; finite propagation alone and Γ>0 do not prove uniform limiting amplitude there. Fixed-x resolvent asymptotics are insufficient because x grows with T.

Do not choose a different ray, retune a or omit DC/initial-slip to obtain positivity. One controlled spectral route and the alternative exterior identity have been exhausted here. **Return A-GATE as OPEN_PROOF.**
