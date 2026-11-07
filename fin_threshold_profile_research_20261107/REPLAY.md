# Replay

From the repository root run with Python and NumPy.

```sh
python fin_threshold_profile_research_20261107/check_profile.py
```

The resonant comparison uses Fourier(G)(k)=exp(-k^2/2), momentum forcing
epsilon^3 cos(m t)G and m^2=16/3. The exact oscillator response is
(cos(m t)-cos(nu t))/(nu^2-m^2), evaluated with a sinc formula near k=0.
Frequency integration covers [0,8] and includes the reflected half-line.
Gaussian tails beyond 8 are truncated in this finite diagnostic.
Observation-kernel quadrature uses 64 points on [-1,1]. Critical times
are 64,256,1024. Time and frequency refinements are recorded in JSON.
Three symmetric translations test dependence on source position.

A second diagnostic integrates the squared difference of
cos(kL)exp(-k^2/2) from its threshold value with measure dk/k on (0,1/2).
It compares two quadrature steps and the derived logarithmic growth.
The analytical proof uses a fixed smooth compact frequency cutoff for
its low/high decomposition, and establishes uniform bounds separately.
The Gaussian runs check representative resonant responses.

The script records NumPy's version and predecessor hash checks.
MANIFEST.json binds package files other than itself.
