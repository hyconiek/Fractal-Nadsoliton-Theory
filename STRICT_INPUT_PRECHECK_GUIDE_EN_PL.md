# Strict Input Precheck Guide (EN + PL)

Date: 2026-03-04

## EN

This guide defines strict precheck requirements for the two currently non-closing channels:
- `h0` + `lambda_cosmological` (`QW-2099 -> QW-2102 -> QW-2090`)
- `G_newton` (`QW-2101 -> QW-2103 -> QW-2092`)

### A) H(z) channel strict-ready requirements

Input artifact: `h0_lambda_decoupling_input_qw2090.json` (built by `QW_2099_*`).

Hard thresholds (`QW-2102`):
- `n_nodes >= 5`
- `z_span >= 0.8`
- `e_span >= 1.0` where `E(z)=Omega_m(1+z)^3+Omega_r(1+z)^4`
- `cond([E,1]) < 8`

Metadata/provenance requirements:
- non-placeholder `source/citation/reference_url`
- integrity metadata (`source_sha256` or `source_version`)
- `provenance_anchor_free = true`

Run:

```bash
python3 QW_2099_HZ_EXTERNAL_DECOUPLING_AUTOCOLLECTOR.py \
  --nodes-csv external_hz_nodes_qw2099.csv \
  --citation "YOUR_DATASET_CITATION" \
  --reference-url "YOUR_DATASET_URL_OR_DOI" \
  --source-version "YOUR_DATASET_VERSION" \
  --require-strict-ready

python3 QW_2102_HZ_DECOUPLING_IDENTIFIABILITY_GATE.py --input h0_lambda_decoupling_input_qw2090.json
python3 QW_2090_H0_LAMBDA_DECOUPLING_GATE.py --input h0_lambda_decoupling_input_qw2090.json
```

Interpretation:
- if `QW-2099 --require-strict-ready` exits non-zero, input lever-arm is insufficient.
- if `QW-2102` is pending, `QW-2090` cannot be promoted to strict closure.

### B) G_newton channel strict-ready requirements

Input artifact: `gnewton_si_bridge_input_qw2092.json` (built by `QW_2101_*`).

Hard provenance requirements (`QW-2103`):
- `bridge_observable_origin = external_dimensionless_observable`
- `provenance_anchor_free = true`
- `seeded_from_registry = false`
- `g_si_input_optional` must be `null` (not SI-primary)

Run:

```bash
python3 QW_2101_GNEWTON_BRIDGE_EXTERNAL_AUTOCOLLECTOR.py \
  --source-file external_gnewton_bridge_qw2101.json \
  --citation "YOUR_DIMENSIONLESS_BRIDGE_CITATION" \
  --reference-url "YOUR_DIMENSIONLESS_BRIDGE_URL_OR_DOI" \
  --source-version "YOUR_DIMENSIONLESS_BRIDGE_VERSION" \
  --strict-dimensionless-only \
  --omit-g-si-optional \
  --require-strict-ready

python3 QW_2103_GNEWTON_DIMENSIONLESS_PROVENANCE_GATE.py --input gnewton_si_bridge_input_qw2092.json
python3 QW_2092_GNEWTON_SI_BRIDGE_GATE.py --input gnewton_si_bridge_input_qw2092.json
```

Interpretation:
- if strict mode fails in `QW-2101`, the bridge is still tautological or provenance-incomplete.
- if `QW-2103` is pending, `QW-2092` cannot be promoted to strict closure.

---

## PL

Ten przewodnik definiuje ścisłe wymagania wejściowe dla dwóch obecnie niedomkniętych kanałów:
- `h0` + `lambda_cosmological` (`QW-2099 -> QW-2102 -> QW-2090`)
- `G_newton` (`QW-2101 -> QW-2103 -> QW-2092`)

### A) Kanał H(z): wymagania strict-ready

Artefakt wejściowy: `h0_lambda_decoupling_input_qw2090.json` (budowany przez `QW_2099_*`).

Progi twarde (`QW-2102`):
- `n_nodes >= 5`
- `z_span >= 0.8`
- `e_span >= 1.0`, gdzie `E(z)=Omega_m(1+z)^3+Omega_r(1+z)^4`
- `cond([E,1]) < 8`

Wymagania metadanych/proweniencji:
- nie-placeholder `source/citation/reference_url`
- metadane integralności (`source_sha256` lub `source_version`)
- `provenance_anchor_free = true`

Uruchomienie:

```bash
python3 QW_2099_HZ_EXTERNAL_DECOUPLING_AUTOCOLLECTOR.py \
  --nodes-csv external_hz_nodes_qw2099.csv \
  --citation "YOUR_DATASET_CITATION" \
  --reference-url "YOUR_DATASET_URL_OR_DOI" \
  --source-version "YOUR_DATASET_VERSION" \
  --require-strict-ready

python3 QW_2102_HZ_DECOUPLING_IDENTIFIABILITY_GATE.py --input h0_lambda_decoupling_input_qw2090.json
python3 QW_2090_H0_LAMBDA_DECOUPLING_GATE.py --input h0_lambda_decoupling_input_qw2090.json
```

Interpretacja:
- jeśli `QW-2099 --require-strict-ready` zwraca błąd, dźwignia wejścia jest za słaba.
- jeśli `QW-2102` jest pending, `QW-2090` nie może być promowany do strict closure.

### B) Kanał G_newton: wymagania strict-ready

Artefakt wejściowy: `gnewton_si_bridge_input_qw2092.json` (budowany przez `QW_2101_*`).

Wymagania twarde proweniencji (`QW-2103`):
- `bridge_observable_origin = external_dimensionless_observable`
- `provenance_anchor_free = true`
- `seeded_from_registry = false`
- `g_si_input_optional` musi być `null` (nie SI-primary)

Uruchomienie:

```bash
python3 QW_2101_GNEWTON_BRIDGE_EXTERNAL_AUTOCOLLECTOR.py \
  --source-file external_gnewton_bridge_qw2101.json \
  --citation "YOUR_DIMENSIONLESS_BRIDGE_CITATION" \
  --reference-url "YOUR_DIMENSIONLESS_BRIDGE_URL_OR_DOI" \
  --source-version "YOUR_DIMENSIONLESS_BRIDGE_VERSION" \
  --strict-dimensionless-only \
  --omit-g-si-optional \
  --require-strict-ready

python3 QW_2103_GNEWTON_DIMENSIONLESS_PROVENANCE_GATE.py --input gnewton_si_bridge_input_qw2092.json
python3 QW_2092_GNEWTON_SI_BRIDGE_GATE.py --input gnewton_si_bridge_input_qw2092.json
```

Interpretacja:
- jeśli tryb strict w `QW-2101` pada, mostek jest nadal tautologiczny lub ma niepełną proweniencję.
- jeśli `QW-2103` jest pending, `QW-2092` nie może być promowany do strict closure.
