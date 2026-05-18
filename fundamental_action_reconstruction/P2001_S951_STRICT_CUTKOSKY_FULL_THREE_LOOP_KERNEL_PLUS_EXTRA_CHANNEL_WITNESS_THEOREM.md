# P2001 / S951 — Strict Cutkosky full-three loop-kernel plus extra-channel witness

P2001 wykonuje następny uczciwy krok po P2000:

1. promuje `hh` z proxy do jawnej formy jądra pętlowego,
2. dodaje nową klasę stanu pośredniego `gx` jako test brakującego kanału.

- strict lane only,
- brak transferu roszczeń legacy,
- status pozostaje `OPEN_OBSTRUCTION_WITH_TRACE`.

To jest pierwszy etap, w którym wszystkie trzy główne kanały mają jawne formy
loop-kernel, a analiza `Delta_opt(s)` zawiera dodatkowy kanał kontrolny.
