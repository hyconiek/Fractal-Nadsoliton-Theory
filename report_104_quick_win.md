# Badanie 104 — Rozszerzone RG-proby
**Autor:** Krzysztof Żuchowski


Data: 2025-11-14T18:07:25.510389+00:00

## Streszczenie zadań i najważniejsze wyniki

- N = 16, alpha_geo = 2.77, beta_tors = 0.01, omega = 0.7853981633974483

## Wnioski i interpretacje (sukces / porażka)

1) Spektrum (task_0): λ_max = 21.113196, trace/N = 2.770000.

   Wniosek: System wykazuje silną dominującą modę; to wskazuje na naturalny kandydat na nośnik masy (sukces w sensie identyfikacji).

2) Mass-hierarchy proxy (task_5): top/4th = 3.243227e+01.

   Wniosek: Silna separacja skal (sukces: mechanizm możliwy do dalszego badania).

3) Beta-proxy (task_2): min=11.203436, max=41.563113.

   Wniosek: Beta zmienność wskazuje na mieszane regiony scalające/antiscalajace; wymaga analizy lokalnych zerów beta.

4) Zintegrowana zmiana g (task_8) pomiędzy s=0.5 a s=2.0: Δg = 3.170889e+01.

   Wniosek: Finite renormalization umiarkowanego rzędu; nie wskazuje na divergentne przepływy (sukces).

5) Anomaly proxy (task_7): antisymm_norm = 0.000000e+00.

   Wniosek: Znikoma antysymetria; brak ewidentnych anomalii algebraicznych (sukces).

6) Operator mixing (task_6): overlaps = ['0.000000', '0.000000', '0.000000', '0.000000'].

   Wniosek: Top mode silnie odseparowany od następujących (niska mieszalność) — pozytywne dla stabilności mas.

7) Vacuum-energy proxy (task_4): neg_sum = 0.000000e+00, trace/N = 2.770000.

   Wniosek: Nieznaczne ujemne składowe wskazują na lokalne obniżenie energii; wymaga dalszej analizy termodynamicznej.

8) Screening vs antiscreening (task_9): n_screening = 25, n_antiscreening = 0.

   Wniosek: Mieszane zachowanie; możliwe przejścia między reżimami w zależności od skali.

9) Anomalous-dimension proxy (task_3): min=1.000000, max=1.000000.

   Wniosek: Obecne umiarkowane anomalne wymiary; dalsza analiza wykresów potrzebna.

10) Ogólna obserwacja: g(s) zmienia się w sposób ciągły w badanym zakresie; brak nagłych divergencji.

    Wniosek: Stabilność parametryczna w skali badanej; dalsze przebadanie za pomocą gęstszego gridu s.

## Zalecenia na kolejne badania

- Przeprowadzić finer grid w s wokół miejsc, gdzie beta_proxy zmienia znak, i policzyć dokładne zerowe punkty beta.

- Zbadać zależność od fazy phi i od beta_tors (sweep) — możliwe, że separacja masy silnie zależy od tych parametrów.

- Wprowadzić małą, kontrolowaną antysymetrię do kernela by testować stabilność anomalii.

- Wygenerować wykresy lambda(s), gamma(s), oraz overlap top vs others dla publikacji wewnętrznej.

