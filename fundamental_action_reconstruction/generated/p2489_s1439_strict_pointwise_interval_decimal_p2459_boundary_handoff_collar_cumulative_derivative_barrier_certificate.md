# P2489/S1439 strict pointwise interval-Decimal P2459 boundary-handoff collar cumulative derivative barrier certificate

Status: `FINITE_COLLAR_CUMULATIVE_DERIVATIVE_BARRIER_CERTIFICATE_PROOF_COMPRESSION_NO_NEW_REPLAY_NO_ROOT_WINDOW_NO_GLOBAL_ANALYTIC_MONOTONICITY_NO_SELECTOR_SOURCE_NO_TOE`

## Finite cumulative derivative lower barrier

Audited chain: `P2487/S1437, P2488/S1438`.
Barrier rows reused from P2487/P2488: `134`.
All barrier preconditions met: `True`.
Minimum entry cumulative lower barrier: `5.06293652999729190722297749614506090018367516699953176293930422739711844797226919218097E-7`.
Minimum exit cumulative lower barrier: `5.43967212900842071743255837843696112852118027122229118583045690889003349902508414970795064E-7`.
Minimum positive derivative gain on one row: `2.94452516595360198366833367866235432932404026826144357149393971855114312428267327692163754E-10`.
Total certified derivative gain over collar: `2.63670284865778124097316884096847344624475419545110405715616914967156550658511793601823408E-7`.
Right endpoint cumulative lower barrier: `7.69963937865507314819614633711353434642842936245063582009547337706868395455738712819920408E-7`.
Transport gain matches P2488: `True`.

## Coverage semantics

New Decimal replay rows in P2489: `0`.
Reused barrier rows (not a P2459 coverage count): `134`.
Diagnostic row ratio, not coverage fraction: `134/99846`.
New P2459 unreplayed cells added by P2489: `0`.
P2489 is barrier proof compression, not new replay: `True`.

## Negative controls

P2489 does not export directed rounding, symbolic/global continuum root exclusion, root-window exclusion, global analytic monotonicity, selector/source/gauge closure, QW-2191 discharge, role-bearing L_total, legacy-role transfer, physical-value generation, or ToE closure.

## Lay summary

This packet integrates the certified positive derivative lower bounds from P2487 across the checked P2481 collar, using the P2488 monotonicity lemma as the finite proof context. The cumulative lower barrier remains positive at every row entry and exit, but this is still a collar-local finite certificate outside the excluded root window.

## Fingerprints

Barrier fingerprint: `46a06dfc5895a0adc6a803033a583df4805bcac803b79deb2e96bf8e19324b3c`.
Theorem fingerprint: `49e35751c74badd75c6831e72afe4d538352cae86c65f531cb4c4b7e504164d0`.
