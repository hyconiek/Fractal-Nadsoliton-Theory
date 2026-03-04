# QW1660 Colab Runbook (Strict Methodology)

## Goal

Run full strict methodological rerun for QW-1660 on Colab with higher compute budget:

- v61 full null model (no hardcoded observation),
- v63 whitening with null calibration,
- v65 micro-timeshift with permutation significance,
- QW-2116 repair gate.

## 1) Colab setup

In a Colab cell:

```bash
!git clone https://github.com/hyconiek/Fractal-Nadsoliton-Theory.git
%cd Fractal-Nadsoliton-Theory
!python3 -m pip install --upgrade pip
!python3 -m pip install numpy scipy h5py gwpy
```

## 2) Strict full rerun (recommended baseline)

```bash
!QW1660_STRICT_N_SAMPLES=524288 \
  QW1660_STRICT_N_TRIALS=400 \
  QW1660_STRICT_N_PERM=400 \
  python3 run_qw1660_strict_colab.py
```

Optional stronger rerun:

```bash
!QW1660_STRICT_N_SAMPLES=1048576 \
  QW1660_STRICT_N_TRIALS=800 \
  QW1660_STRICT_N_PERM=800 \
  python3 run_qw1660_strict_colab.py
```

## 3) Output artifacts

- `QW_1660_v61_FullNullModel_strict.json`
- `QW_1660_v63_Whitening_strict.json`
- `QW_1660_v65_MicroTimeShift_strict.json`
- `report_qw2116_gw1660_method_repair_gate.json`
- `RAPORT_QW2116_GW1660_METHOD_REPAIR_GATE.md`
- Notebook export cell additionally creates: `/content/qw1660_export.zip` and attempts `files.download(...)`.

## 4) Strict interpretation

- Methodological pass requires `GW1660_METHOD_REPAIR_GATE_PASS`.
- Scientific anomaly claim remains unsupported unless v61/v65 strict significance and physical lag criteria both pass.
- If gate passes but v65 says `MICRO_TIMESHIFT_SIGNAL_NOT_SUPPORTED`, then the branch is methodologically repaired but does not support a robust GW cross-H anomaly claim.
