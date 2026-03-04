# External Data Sources (No Large Binaries in Git)

Date: 2026-03-04

## Purpose
This project keeps large raw archives out of git history.
Use this file as the canonical source list for downloading external public data.

## Policy
- Do not commit large binary datasets (`.tar.gz`, `.h5`, `.hdf5`) to this repository.
- Keep reproducibility via:
1. explicit source URLs,
2. deterministic collection scripts,
3. SHA256 manifests of generated analysis artifacts.

## Required Public Sources

### 1) NANOGrav 15yr timing archive
- Dataset: `NANOGrav15yr_PulsarTiming_v2.1.0.tar.gz`
- Public URL:
  - `https://zenodo.org/records/16051178/files/NANOGrav15yr_PulsarTiming_v2.1.0.tar.gz?download=1`
- Example download:

```bash
mkdir -p external_data
curl -L -o external_data/NANOGrav15yr_PulsarTiming_v2.1.0.tar.gz \
  "https://zenodo.org/records/16051178/files/NANOGrav15yr_PulsarTiming_v2.1.0.tar.gz?download=1"
sha256sum external_data/NANOGrav15yr_PulsarTiming_v2.1.0.tar.gz
```

### 2) GWOSC GWTC event catalog API
- API:
  - `https://www.gw-openscience.org/eventapi/json/GWTC/`
- Example cache download:

```bash
mkdir -p external_data
curl -L -o external_data/gwosc_gwtc_eventapi.json \
  "https://www.gw-openscience.org/eventapi/json/GWTC/"
sha256sum external_data/gwosc_gwtc_eventapi.json
```

## Optional Source (for raw strain reproductions)
- Provider: GWOSC (Gravitational Wave Open Science Center)
- Base portal:
  - `https://www.gw-openscience.org/`
- Note:
  - raw strain `.h5` caches should be created locally for each run and not committed.

## How to run the external autocollector without storing binaries in git

```bash
python3 QW_2014_TRUE_EXTERNAL_BETA_CHANNEL_AUTOCOLLECTOR_V2.py \
  --nanograv-archive external_data/NANOGrav15yr_PulsarTiming_v2.1.0.tar.gz \
  --gwosc-api "https://www.gw-openscience.org/eventapi/json/GWTC/"
```

## Independent bundle relation
- QW-2033 and QW-2050 freeze bundles track scripts/reports/manifests.
- External raw archives are referenced by source URL and acquisition procedure, not by committing payloads.

## Operational reminder
- Before sharing results:
1. verify manifests (`manifest_qw2033.json`, `manifest_qw2050.json`),
2. verify protocol freeze files,
3. provide this source document with runbooks.
