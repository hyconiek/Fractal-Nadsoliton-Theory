# Independent Check Guide (EN + PL)

## EN: How other people can verify this independently

## Goal
Run a truly independent confirmatory check with the locked protocol, without retuning.

## 1) Freeze reference
Use and publish:
- commit hash,
- lock file SHA256,
- exact Python version and OS.

Current lock:
- `external_confirmatory_v2/independent_multiteam_lock_qw2053/protocol_lock_qw2053.json`
- SHA256: `2b49385e38a9985b7cf031992da821c3da23a6c0390e4643d56b111adf07c880`

## 2) Independent team setup (must be separate)
- Separate machine/account (not your development environment).
- Fresh clone of repository.
- No local edits before run.

## 3) External data policy
- Do not use local private copies hidden from protocol.
- Download only from documented public sources:
  - `DATA_SOURCES_EXTERNAL_DOWNLOADS.md`
- Verify checksums before running gates.

## 4) Required execution order (locked)
Run exactly:
1. `python3 QW_2048_SPECTRAL_PHASE_LOCKED_POINTWISE_DERIVATION.py`
2. `python3 QW_2049_SPECTRAL_MICRO_STAGEC_INTERSECTION_GATE.py`
3. `python3 QW_2051_INDEPENDENT_REHEARSAL_GATE.py`
4. `python3 QW_2052_EXTERNAL_SOURCE_ONLY_GOVERNANCE_GATE.py`
5. `python3 QW_2054_PROTOCOL_LOCK_INTEGRITY_GATE.py`
6. `python3 QW_2067_STRICT_FIRST_PRINCIPLES_INTERNAL_CLOSURE_STRENGTHENED_GATE.py`

## 5) Pass criteria
Minimum independent pass package:
- `QW-2049`: pass with all hard flags true,
- `QW-2051`: independent rehearsal pass,
- `QW-2052`: external-source-only governance pass,
- `QW-2054`: lock integrity pass,
- `QW-2067`: strengthened strict internal closure pass.

## 6) What each team should publish
Each team should publish:
- commit hash,
- lock SHA256,
- all generated report JSON files,
- exact command log,
- environment manifest (OS, Python, package versions),
- statement: "no retune/no threshold edits/no post-hoc selection".

## 7) What would count as a real independent confirmation
At least 2-3 independent teams should reproduce the same locked verdicts from fresh environments.

---

## PL: Jak inni ludzie mają to niezależnie sprawdzić

## Cel
Wykonać prawdziwie niezależny check confirmatory na zablokowanym protokole, bez dostrajania.

## 1) Referencja freeze
Każdy zespół publikuje:
- hash commita,
- SHA256 lock-file,
- dokładną wersję Pythona i systemu.

Aktualny lock:
- `external_confirmatory_v2/independent_multiteam_lock_qw2053/protocol_lock_qw2053.json`
- SHA256: `2b49385e38a9985b7cf031992da821c3da23a6c0390e4643d56b111adf07c880`

## 2) Środowisko niezależnego zespołu
- Osobna maszyna/konto (nie Twoje środowisko deweloperskie).
- Świeży clone repo.
- Brak lokalnych zmian przed uruchomieniem.

## 3) Zasada danych zewnętrznych
- Nie używać ukrytych lokalnych kopii spoza protokołu.
- Pobierać tylko z udokumentowanych źródeł publicznych:
  - `DATA_SOURCES_EXTERNAL_DOWNLOADS.md`
- Sprawdzić checksumy przed uruchomieniem gate'ów.

## 4) Wymagana kolejność uruchomienia (lock)
Uruchomić dokładnie:
1. `python3 QW_2048_SPECTRAL_PHASE_LOCKED_POINTWISE_DERIVATION.py`
2. `python3 QW_2049_SPECTRAL_MICRO_STAGEC_INTERSECTION_GATE.py`
3. `python3 QW_2051_INDEPENDENT_REHEARSAL_GATE.py`
4. `python3 QW_2052_EXTERNAL_SOURCE_ONLY_GOVERNANCE_GATE.py`
5. `python3 QW_2054_PROTOCOL_LOCK_INTEGRITY_GATE.py`
6. `python3 QW_2067_STRICT_FIRST_PRINCIPLES_INTERNAL_CLOSURE_STRENGTHENED_GATE.py`

## 5) Kryteria zaliczenia
Minimalny pakiet niezależnego PASS:
- `QW-2049`: PASS + wszystkie hard flagi true,
- `QW-2051`: PASS próby niezależnej,
- `QW-2052`: PASS governance external-source-only,
- `QW-2054`: PASS integralności locka,
- `QW-2067`: PASS strengthened strict internal closure.

## 6) Co każdy zespół powinien opublikować
Każdy zespół publikuje:
- hash commita,
- SHA256 locka,
- wszystkie wygenerowane raporty JSON,
- log komend,
- manifest środowiska (OS, Python, wersje pakietów),
- oświadczenie: "no retune/no threshold edits/no post-hoc selection".

## 7) Co uznajemy za realne niezależne potwierdzenie
Co najmniej 2-3 niezależne zespoły muszą odtworzyć te same zablokowane werdykty na świeżych środowiskach.
