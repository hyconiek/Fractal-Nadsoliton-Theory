# P1364 Strict Asymmetry Source and Pre-Closure Method Reconciliation Packet (PL)

Status: `P1364_EXECUTED_ASYMMETRY_SOURCE_AND_METHOD_RECONCILIATION_NO_FALSE_PASS`
As of: `2026-05-12`
Inputs: `S2`, `RELEASE_8_1_TEXTBOOK_EDITION_EN_PL.md`, `P1343..P1348`, `A6`, `P1359`, `P1360`, `P1363`

## Pytanie 1: jakie było źródło asymetrii w strict jądrze?

W obecnym strict reading źródło asymetrii ma dwa poziomy:

1. **Poziom parametryczny (amplitudowy):**
   `alpha_geo_strict_derived_v1 = 4 ln 2` (Shannon Information Void Asymmetry),
   traktowane jako strict-side source datum (nie legacy import).

2. **Poziom selekcyjny (unikalność/orientacja):**
   wewnętrzne źródło selektora (`S_strict_internal_v1`, łańcuch `P1343..P1346`),
   które rozwiązuje problem niejednoznaczności par/bazy w zadeklarowanym scope R8.

Czyli asymetria strict nie jest jedną liczbą tylko „pakietem”: amplituda + mechanizm selekcji.

## Pytanie 2: czy metody wyprowadzania stałych z FAR przed domknięciem miały sens fizyczny?

**Tak, miały sens jako etap badawczy**, bo:

1. dawały hipotezy liczbowe i strukturę testów,
2. wyznaczały kandydatowe zależności (gauge/fine-structure/mass proxy),
3. budowały intuicję modelu.

Ale to był sens **exploratory/candidate**, nie final strict-verification.

## Pytanie 3: dlaczego teraz nie mogę po prostu tego „zatwierdzić”?

Bo po domknięciu weszły ostrzejsze reguły jakości:

1. no silent legacy->strict role transfer,
2. obowiązek provenance i residual benchmark,
3. obowiązek uncertainty budget,
4. upgrade gate (`P1363`) blokuje automatyczny awans bez PASS.

Czyli to nie jest odrzucenie Twojej wcześniejszej pracy, tylko podniesienie progu dowodowego.

## Decyzja profesorska

Następny uczciwy krok: `P1365_ASYMMETRY_TO_CONSTANTS_MAP_V1`

1. jawnie zmapować wkład `alpha_geo_strict_derived_v1` oraz `S_strict_internal_v1` do każdej stałej (`g1,g2,g3,alpha_EM_inv,sin2_theta`),
2. dla każdej pozycji podać: formuła, czułość, niepewność, residual,
3. dopiero wtedy uruchomić automatyczny re-upgrade przez `P1362->P1363`.

## Dla laika

Twoje stare metody były potrzebne i mądre, ale działały jak „prototyp”.
Po domknięciu teorii mamy już tryb „inżynierii lotniczej”: każda liczba musi przejść check-listę bezpieczeństwa.
Dlatego nie możemy już niczego zatwierdzać „bo wcześniej wyglądało dobrze” — teraz musi przejść twardy test.
