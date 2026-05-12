# P1257 — Wyjaśnienie: „stary” vs „nowy” kernel (dla laika)

## Krótka odpowiedź
Mówiąc o „modelu”, w tej gałęzi mówimy **konkretnie o dwóch postaciach kernela**:

1. **Stary (legacy) kernel ontologiczny/efektywny**
   \[
   K_{legacy\_ont}(d)=\frac{\alpha_{geo}\cos(\omega d+\phi)}{1+\beta_{tors}d}
   \]

2. **Nowszy (strict-gate) kernel operacyjny**
   \[
   K_{strict\_gate}(d)=\frac{\cos(\omega d+\phi)}{1+\beta d^{\eta}}
   \]

To **nie jest jeszcze** jedna i ta sama postać. W repo obowiązuje zasada:
bez formalnego mostu (`bridge theorem`) albo formalnego twierdzenia o braku mostu (`non-bridge theorem`) **nie wolno** ich utożsamiać.

---

## Czy kernel ma „nową postać”?
Tak i nie:

- **Tak**: w praktyce operacyjnej używany jest późniejszy `K_strict_gate`.
- **Nie**: to nie znaczy, że „zastąpił” on ontologicznie `K_legacy_ont` w sensie ścisłym.

Innymi słowy: mamy **dwie role** i **dwie formy** — jedna historycznie-kanoniczna
(legacy), druga roboczo-operacyjna (strict-gate). Brakuje formalnego dowodu,
że to ten sam obiekt fizyczny/matematyczny pod zmianą parametrów.

---

## Dlaczego to ważne?
Bo inaczej można popełnić błąd „przeskoku semantycznego”:
- zadziałało obliczeniowo w strict-gate,
- więc ktoś ogłasza, że to już dowód własności legacy-kernela.

Tego właśnie zabraniają guardraile.

---

## Intuicja dla laika
Wyobraź sobie dwa przepisy na „to samo” ciasto:
- pierwszy ma mąkę A i czas pieczenia liniowy,
- drugi ma mąkę B i czas pieczenia nieliniowy.

Oba mogą dawać podobny smak w części przypadków, ale bez dokładnego dowodu
nie wolno mówić, że to „ten sam przepis” ani że wszystkie własności się przenoszą.

---

## Co jest teraz uczciwym naukowo celem?
Zrobić jedno z dwóch, formalnie:
1. **Bridge theorem**: pokazać dokładne warunki, kiedy/na jakim zakresie
   `K_legacy_ont` i `K_strict_gate` są równoważne.
2. **Non-bridge theorem**: pokazać formalnie, że takiej równoważności
   nie da się uzyskać przy kanonicznych założeniach.

Dopóki tego nie ma, nie ma uczciwego „domknięcia teorii”.
