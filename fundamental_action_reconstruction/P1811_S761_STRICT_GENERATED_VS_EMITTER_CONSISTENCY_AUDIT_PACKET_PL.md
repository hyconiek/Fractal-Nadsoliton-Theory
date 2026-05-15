# P1811 S761 Strict Generated-vs-Emitter Consistency Audit Packet (PL)

Status: `P1811_EXECUTED_STRICT_GENERATED_VS_EMITTER_CONSISTENCY_AUDIT_PACKET_NO_FALSE_PASS`

## Cel

Domknąć proceduralny bottleneck: zgodność między plikami `generated/*.json` a aktualną logiką emitterów.

Bez tego można mieć stare artefakty z innym statusem gate niż wynika z bieżących reguł kodu.

## Zakres

Audit obejmuje krytyczny fragment BW lane:

- `p1808` emitter vs `generated/p1808...json`,
- `p1810` merger inputs + resolved TG1.

## Reguła

Jeśli `generated` różni się od bieżącego emittera:
- status = `OPEN_OBSTRUCTION_WITH_TRACE`,
- wymagana regeneracja artefaktu,
- locki mają używać wyłącznie statusu po regeneracji.

## Co zostało dowiedzione

1. Potrzebny jest formalny check „artifact freshness”, nie tylko logic contracts.
2. W BW lane spójność emitter->generated jest teraz jawnie audytowana.

## Co pozostaje OPEN

1. Każda przyszła zmiana emitera wymaga rerun audytu.

## Ryzyka false-pass

1. Użycie przestarzałego `generated` przy nowych regułach locków.
2. Manualna edycja json bez ponownego uruchomienia emitera.

## Następny uczciwy krok

Wpiąć audit do stałego workflow: `run emitter -> run consistency audit -> dopiero state-vector update`.
