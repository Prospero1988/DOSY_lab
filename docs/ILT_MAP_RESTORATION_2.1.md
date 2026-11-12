# DOSY Lab 2.1 · Przywrócenie eksploracyjnych map ILT

## Co się zmieniło

W wersji 2.0.0 na mapie Global ILT można było zobaczyć cały zmierzony profil ROI przy maksimum wspólnego rozkładu, nawet jeśli maksimum nie opisywało zaniku danego sygnału. W AILT część sygnałów poza grupą była rysowana przy D z dopasowania monoeksponencjalnego. W commicie `e601223` oba przypadki trafiły pod kryteria rekonstrukcji. W projekcie strychniny cztery pozycje AILT zniknęły z mapy 2D, mimo że były obecne na widmie 1D. Dwa ROI miały zbieżny wynik grupowy, lecz nie przeszły testu zaniku dla D maksimum. Pozostałe dwa były poza zgodną grupą i nie miały stabilnego niezależnego rozdzielenia, które dopuszczał poprzedni widok.

Teraz w oknie **Mapa DOSY 2D** można przełączać trzy widoki bez ponownego obliczania ILT. **Eksploracyjna** rysuje profile ze zbieżnego rozwiązania, zachowując ostrzeżenia. **Filtrowana QC** zachowuje wcześniejsze kryteria przypisania. **Odniesienie monoeksponencjalne** pokazuje oryginalne profile przy mono D i ma własny tytuł w PNG. Wybrana polityka, liczba prób, zbieżnych i narysowanych ROI oraz liczby pominiętych i nierozdzielonych są widoczne w legendzie. Powody pominięcia, reszty i udziały wyłączonych składowych są w podpowiedzi ROI.

## Jak powstaje widok eksploracyjny

Dla każdej składowej bazą zaniku jest suma eksponent z jej pełnego przedziału siatki ILT, z wagami odzyskanymi dla danego ROI. Pozycję D na mapie wyznacza najsilniejsza komórka tego ROI w danym przedziale. W Global ILT może się ona różnić od maksimum zbiorczego. Dla rozkładu z kilkoma przedziałami dopasowanie NNLS wykorzystuje wszystkie ich bazy jednocześnie do odtworzenia widma ze znakiem w kolejnych warstwach. Każdy bin ma te same współrzędne ppm co w pomiarze. Osobno wyznaczony profil może mieć mały udział, jeżeli ma dla siebie dodatnią masę w rozwiązaniu ILT i w dopasowaniu widmowym. Kod nie wymaga, by każda składowa opisywała zanik całego ROI.

Pole każdej składowej ogranicza jej przewidywana całka na warstwie referencyjnej. Sumę intensywności wszystkich składowych w każdym binie ogranicza dodatnia część zmierzonego widma referencyjnego. Oba ograniczenia mogą tylko zmniejszyć intensywność. Wyliczona reszta widmowa obejmuje ten etap. Wspólne biny nakładających się ROI pozostają zamaskowane. Przy jednym przedziale kryterium zgodności zaniku z pojedynczym D pozostaje w diagnostyce QC, a eksploracyjna mapa może pokazać profil dopasowany przy wykorzystaniu całego przedziału. Jeżeli dekompozycja nie ma pełnego rzędu lub zabraknie dodatniego profilu, mierzony sygnał nadal można zobaczyć na nakładce 1D.

AILT wykonuje inwersję indywidualną dla ROI z użytecznymi danymi, które nie weszło do zgodnej grupy. Podpowiedź oznacza ten przypadek. W trybie **Indywidualne ILT** mapa wykorzystuje wyniki uzyskane przez przycisk **Analizuj ILT** w tabeli głównej. Wyniki z kolejnych sygnałów są dostępne do zmiany danych lub processingu.

## Porównanie na projekcie strychniny

We wszystkich wariantach użyłem tych samych 20 ROI, 16 gradientów, warstwy referencyjnej 2, fazy, korekcji bazy, progu 0% i wzmocnienia 1×. Dla Global ILT przyjąłem 100 punktów w siatce logarytmicznej, λ = 10⁻⁶ i μ = 10⁻⁴. Dla AILT użyłem tolerancji D 5%, balansu amplitud 10 i λ = 10⁻⁶. Indywidualne ILT korzystało z domyślnych parametrów okna analizy.

| Mapa | 2.0.0 | Przed zmianą w 2.1 | Eksploracyjna | Filtrowana QC |
| --- | ---: | ---: | ---: | ---: |
| Mono-fit | 20 | 20 | 20 jako odniesienie | — |
| Shared-D | 20 | 20 | — | — |
| Indywidualne ILT | — | — | 20 | 18 |
| Global ILT | 20 | 18 | 20 | 18 |
| Adaptive Local ILT | 20, w tym 2 profile z mono D | 16 | 20 | 16 |

W porównaniu z `e601223` wartości b, ppm, całki, wartości mono D, wagi i przewidywania wspólnego i indywidualnego ILT oraz wagi grupowego AILT były identyczne. Względem wydania 2.0.0 maksymalna różnica mono D wyniosła 1,70 × 10⁻¹⁷ m²/s. Maksymalne różnice bezwzględne w znormalizowanych wagach wyniosły 5,68 × 10⁻⁹ dla Global ILT, 2,92 × 10⁻⁹ dla grupowego AILT i 1,15 × 10⁻⁹ dla indywidualnego ILT. W tym projekcie różnice widoczności wynikają z kryteriów projekcji i dołożenia inwersji dla dwóch ROI poza grupami. [Wyniki benchmarku](../benchmarks/ilt_map_comparison.json)

**S2, strychnina i chloroform.** Mono-fit daje 2,427 × 10⁻⁹ m²/s dla wspólnej całki. Global ILT w 2.0.0 rysował całe ROI przy maksimum zbiorczym 2,128 × 10⁻⁹ m²/s. W widoku eksploracyjnym Global ILT pokazuje jedno maksimum przy 2,417 × 10⁻⁹ m²/s, a indywidualny ILT i niezależna inwersja AILT pokazują po jednym przy 2,501 × 10⁻⁹ m²/s. Eksploracyjna mapa przypisuje około 98% dodatniego pola S2. Filtrowana mapa pozostawia ten sygnał na nakładce 1D. Podane wartości opisują różne reprezentacje tego samego mieszanego zaniku. Dane i użyte ustawienia nie dają dwóch rozdzielonych składowych.

## Testy i ograniczenia

Przeszło 100 testów, w tym 63 testy obliczeń i 37 testów GUI.

- W syntetycznym przypadku z dwiema nakładającymi się składowymi i znaną bazą z rozwiązania ILT widok eksploracyjny zachował składnik o udziale około 2,91%, wraz z jego własnym polem i pozycją D. Filtr QC pozostawił ten przypadek bez przypisania.
- Test z ILT wyjaśniającym 75% pola potwierdził, że mapa nie uzupełnia brakujących 25%. Osobny test utrzymał zamaskowanie binów nakładających się ROI, także gdy tylko jeden ROI miał odtworzony profil.
- W teście rzeczywistego projektu w oknie Qt sprawdziłem wszystkie trzy polityki dla każdej metody ILT i dziewięć eksportów PNG. Przełączenie polityki nie zmieniło wersji obliczenia ani wyników monoeksponencjalnych.
- Maksima i pola zależą od próbkowania b, S/N, regularyzacji, fazy i odejmowania bazy. Widok eksploracyjny może zachować składową, której zanik przy pojedynczym D jest słabym przybliżeniem, albo maksimum na granicy siatki. Pozycja pasma odpowiada maksimum przedziału, a szerokość graficzna nie jest rozrzutem ani niepewnością D. Nakładających się chemicznie sygnałów nie można uznać za rozdzielone na podstawie liczby narysowanych pasm.

Porównanie można odtworzyć poleceniem poniżej. Skrypt uruchamia archiwalne źródła w osobnych katalogach tymczasowych. Bez flagi `--summary` zwraca również siatki, wagi, pola i pozycje pasm dla każdego ROI.

```powershell
& '.venv/Scripts/python.exe' benchmarks/compare_ilt_maps.py 'projects/Strychnina_CDCl3_bpp_ste_led_diffusion-1-1.proj' --summary
```

Zmiany w projekcji są w `bin/dosy_alignment.py` i `bin/adaptive_ilt.py`. Przełącznik, diagnostyka i eksport są w `bin/dosy_map.py` i `bin/dosy_alignment_panel.py`. Testy regresji są w `tests/test_qc_maps.py` i `tests/test_gui.py`.
