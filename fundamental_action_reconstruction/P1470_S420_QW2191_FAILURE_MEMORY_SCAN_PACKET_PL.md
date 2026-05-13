# P1470 — S4.20 QW-2191 Failure-Memory Scan (PL)

Status: `P1470_EXECUTED_QW2191_FAILURE_MEMORY_LOCAL_ONLY`
As of: `2026-05-13`

## Cel

Nie powtarzać ślepych uliczek: zeskanować repo pod kątem udokumentowanych dróg,
które **nie** rozwiązały `QW-2191`, i zapisać je jako jawne anty-wzorce.

## Zakres

- lane viscosity/proxy (`V1..V7`),
- lane grammar/certificate (`T6..T9`),
- lokalne checkpointy FAIL z serii `P1452+`.

## Zasada

Każda przyszła propozycja musi deklarować, że nie powiela żadnego anty-wzorca
z P1470 lub musi jawnie wyjaśnić, co jest nowe.
