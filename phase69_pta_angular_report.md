# Raport Fazy 69: FIN Hellings-Downs Analog (PTA Angular Correlation Test)

Ten raport dokumentuje "killer test" teorii FIN (Fractal Information Nadsoliton) na danych z rejestru Pulsar Timing Arrays (PTA). 

Teoria sugerowała, że ułamek $H_{cross} \approx 0.23$ posiada charakter geometryczny, tzn. koreluje z układem przestrzennym pulsarów (zależność $H_{cross}(\theta)$). Klasyczny szum czerwony (common red noise) oraz artefakty wycieku z modelu czasowego (timing leakage) produkują płaską linię ($const$) niezależnie od geometrii.

## Kod Źródłowy: `FIN_Search_PTA_Angular_v3.py`

```python
"""
FIN Search in PTA Phase 69: FIN Hellings-Downs Analog Test
Testing if Fractal Memory (Cross-H) depends on Angular Separation.
Author: Antigravity & User | Date: Feb 2026
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress
from tqdm import tqdm
import tarfile, os
from pathlib import Path
from itertools import combinations
from astropy.coordinates import SkyCoord
import astropy.units as u
import warnings

warnings.filterwarnings('ignore')
sns.set_style("darkgrid")
plt.rcParams['figure.figsize']=(12,6)

DATA = Path("./nano15")
tar_path = DATA/"nano15.tar.gz"
extract_dir = DATA/"residuals"
par_dir = DATA/"parfiles"

# [EXTRACTION AND LOADING OF .RES DATA OMITTED HERE FOR BREVITY]

def load_pulsar_positions(parfile_dir):
    pos = {}
    for par in Path(parfile_dir).rglob("*.par"):
        name = par.stem
        base_name = name.split('_')[0]
        ra, dec, elong, elat = None, None, None, None
        with open(par) as f:
            for line in f:
                if line.startswith("RAJ"): ra = line.split()[1]
                if line.startswith("DECJ"): dec = line.split()[1]
                if line.startswith("ELONG"): elong = line.split()[1]
                if line.startswith("ELAT"): elat = line.split()[1]
        try:
            if ra and dec:
                pos[base_name] = SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame='icrs')
            elif elong and elat:
                pos[base_name] = SkyCoord(lon=float(elong)*u.deg, lat=float(elat)*u.deg, frame='barycentrictrueecliptic').transform_to('icrs')
        except: pass
    return pos

positions = load_pulsar_positions(par_dir)

def angular_sep(psr1, psr2):
    p1_base, p2_base = psr1.split('_')[0], psr2.split('_')[0]
    if p1_base in positions and p2_base in positions:
        return positions[p1_base].separation(positions[p2_base]).deg
    return None

def cross_dfa(x, y, min_scale=15):
    # [STANDARD CROSS-DFA IMPLEMENTATION OMITTED HERE FOR BREVITY]
    # [...]
    return slope

# 5. Execute Angular Cross-DFA on 500 Random Pairs
psr_list = list(residuals.keys())
pairs = list(combinations(psr_list, 2))
np.random.seed(42)
indices = np.random.choice(len(pairs), min(500, len(pairs)), replace=False)
pairs = [pairs[i] for i in indices]

theta_vals, Hcross_real, Hcross_phase = [], [], []

for p1, p2 in tqdm(pairs, desc="Cross DFA pairs"):
    sep = angular_sep(p1, p2)
    if sep is None: continue
    x, y = get_uniform_grid_and_interp(residuals[p1], residuals[p2], N_points=500)
    Hxy = cross_dfa(x, y)
    Hxy_p = cross_dfa(phase_randomize(x), phase_randomize(y))

    if not np.isnan(Hxy) and not np.isnan(Hxy_p):
        theta_vals.append(sep)
        Hcross_real.append(Hxy)
        Hcross_phase.append(Hxy_p)

# Plotting Binning and Statistical Dump Omitted
```

## Wyniki JSON:
Plik JSON wygenerowany przez kod zawiera test zgodności z płaskim modelem rozkładu ("Flat Line correlation"). Brak korelacji między separacją kątową pulsarów $\theta_{ij}$ a wyliczanym $H_{cross}$ oznacza losowy rozrzut punktów, co falsyfikuje fraktalne pola geometryczne na rzecz powszechnych szumów detektora.

```json
{
    "R2_real_vs_angle": 0.016891542846312124,
    "R2_phase_vs_angle": 0.10427312458891623,
    "real_mean": 0.7631922652981568,
    "phase_mean": 0.7606462155790492
}
```

### Konkluzja Ostateczna:

Z wyników jasno wynika:
1. Przeciętna siła $H_{cross}$ dla sygnałów (REAL = 0.763) nie wykazuje drastycznej preferencji wobec zrandomizowanych faz Surrogate (PHASE = 0.760). Rzeczywiste korelacje pamięciowe tkwią niemal całkowicie w modelu urojonego Power Spectral Density (PSD) pochodzącego z szumu "common red noise".
2. Badanie krzywej dla zróżnicowanego kąta dało wynik $R^2 \approx 0.017$. Kształt zależności $H(\theta)$ rysuje płaską linię, co całkowicie zaprzecza tezie geometrycznej. 

Zatem, struktura ta jest klasycznym **płaskim tłem (flat fluctuation)** i dowodzi bezsprzecznie statystycznej awarii FIN na dystansach obserwacyjnych Pulsar Timing Array.
