# MP7-025 — simple-fold normal-form data from the accepted interval root

Scientific state: **PROVED_INTERVAL_ASSISTED LOCAL NORMAL-FORM COEFFICIENTS**.

This task does not use the decimal gain `3.51564471684` as an exact input.  It
imports the accepted 9-variable Krawczyk box `R7P-031`, whose variables are

```
(s3,s4,s5,s6,g,v3,v4,v5,v6),
```

with equations

```
grad_s Phi(s,g)=0,
H4(s,g) v=0,
||v||^2=1.
```

Thus the fold vector `v` is normalized by an equation inside the validated
root system, not normalized afterwards numerically.

The accepted fold gain interval is

```
3.5156447068395917... <= g_fold <= 3.5156447268395917...
```

and the normalized null-vector box is the last four coordinates of the same
Krawczyk box, centered near

```
v=(0.50736753967725,
   0.528685655313811,
   0.553760694068750,
   0.395498105244206).
```

A fixed 3-by-3 principal restriction of `H4` has three strictly positive LDL
pivots throughout the root box, while `H4 v=0`; the sine block has three
strictly positive pivots.  Hence the isolated full 7D stationary root has one
zero and six positive transverse directions.

Let `F=grad_s Phi`, `epsilon=g-g_fold`, and choose the sign of the normalized
fold coordinate `xi` so the imported `v` is used.  The Lyapunov--Schmidt
linear/cubic coefficients are

```
a = v^T partial_g F
  in [-0.2125716400840138741, -0.2125716260204590083],

b = D^3 Phi[v,v,v]
  in [0.1189825274504987752, 0.1189826671032236644].
```

Therefore `a<0` and `b>0` are certified.  After solving the positive transverse
variables locally, the reduced potential has the form

```
Psi(xi,epsilon)=Psi0(epsilon)+a epsilon xi+(b/6)xi^3+R(xi,epsilon),
```

where the existence of a smooth reduction follows from the certified positive
transverse Hessian.  This task certifies the simple-fold premises and the
leading coefficients; **MP7-026 still has to pay an explicit uniform bound on
R and its first two xi derivatives before the square-root and 3/2 scaling laws
are promoted with error bars.**

No claim is made that this is the first/global fold or a physical-time event.
