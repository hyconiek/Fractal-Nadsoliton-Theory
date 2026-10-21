# Trzy kanały oddawania energii i rzeczywista przeszkoda progowa

Data: 2026-10-21. Warunkowe dowody analityczne, audyt własny.

Wyprowadziłem pełny bilans trzech kwadratowych kanałów promieniowania
dwóch modów związanych. Współczynniki mają dodatnie, dokładne dolne
ograniczenia. To nowy obiekt potrzebny do wspólnej kontroli modów,
nie dokładniejsze cyfry starego Gamma.

Ważne ograniczenie: obliczony bilans dotyczy zamrożonej odpowiedzi
referencyjnej. Nie jest jeszcze prawem ewolucji działań rzeczywistego PDE.

Przez browser-harness przeczytałem pierwotną pracę o metodzie wirialnej.
Dotyczy jednego pola i perturbacji nieparzystych. Nasz dodatkowy kanał
parzysty ma konkretną przeszkodę: całka lokalnej energii liniowej może
rosnąć jak log T. Wyprowadziłem więc jawną, czasowo przyczynową składową
progową Theta i oszacowanie rzeczywistej reszty po jej lokalnym odjęciu.
Nie jest to dodatkowe pole fizyczne ani nowy postulat.

[Dowód kanałów](THREE_CHANNEL_RESPONSE.md),
[przeszkoda i naprawa progowa](VIRIAL_SCOPE_AND_THRESHOLD_GATE.md).
[Następny krok](NEXT_JOINT_VIRIAL_GATE.md) wymaga wspólnej koercywności
oraz małego błędu całkowanego na czasie epsilon^-2. Tego nie zastąpi
dodatniość samych współczynników. Źródła i AGENTS zachowano; F01 i fizyczne
bramki QM/GR/ToE pozostają otwarte.
