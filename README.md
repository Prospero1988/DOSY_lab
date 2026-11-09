# DOSY Lab

DOSY Lab to aplikacja desktopowa do analizy klasycznych eksperymentów dyfuzyjnych NMR z surowych plików JEOL JDF. Otwiera serię jednowymiarowych widm z różnymi gradientami, pozwala wspólnie je przetworzyć i wyznaczyć współczynnik dyfuzji **D** dla wskazanych sygnałów.

![Główne okno DOSY Lab z widmem, całkami, zanikami i resztami](docs/images/widmo-integracja.png)

*Dane syntetyczne. Widmo referencyjne u góry, całki i diagnostyka dopasowań poniżej.*

## Co można zrobić

- Wyświetlić dowolną warstwę lub nałożyć kilka widm. Aplikacja po otwarciu pokazuje warstwę o najniższym dostępnym **b**.
- Ustawić apodyzację, dopełnienie zerami, referencję ppm i korekcję filtra cyfrowego. Skorygować PH0, PH1 i pivot interaktywnie, wspólną fazą dla całej serii.
- Wybrać fragmenty tła do korekcji linii bazowej. Dodawać, przesuwać, dzielić i kalibrować całki na widmie.
- Sprawdzić pochodzenie wartości **b** lub wczytać własną listę. Obejrzeć zanik całki, model i reszty dla każdego sygnału. Ręcznie decydować o udziale poszczególnych warstw w dopasowaniu.
- Porównać wartości **D** w obrębie wskazanej grupy sygnałów oraz zestawić podsumowania osobnych eksperymentów.
- Zapisać projekt, wyeksportować tabele do CSV i XLSX, wykresy do PNG i SVG oraz dane zespolone do NPZ.

## Pierwsze uruchomienie na Windows

Wymagany jest Python 3.11 lub nowszy. Aplikacja działa lokalnie. Połączenie z internetem jest potrzebne przy pobraniu projektu i pierwszej instalacji bibliotek.

Pobierz repozytorium przez **Code → Download ZIP** i rozpakuj je albo użyj Gita.

~~~powershell
git clone https://github.com/Prospero1988/DOSY_lab.git
cd DOSY_lab
~~~

W otwartym folderze kliknij dwukrotnie **start.cmd**. Przy pierwszym uruchomieniu skrypt przygotuje środowisko `.venv` i zainstaluje zależności. W PowerShell ten sam skrypt uruchomisz przez `.\start.cmd`.

<details>
<summary>Instalacja ręczna z PowerShell</summary>

~~~powershell
py -3 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe app.py
~~~

Jeżeli instalator Pythona udostępnił polecenie `python` zamiast `py`, użyj `python -m venv .venv` w pierwszym kroku.

</details>

## Przykład bez własnego pomiaru

Po pierwszym uruchomieniu z katalogu aplikacji wygeneruj syntetyczny eksperyment z 16 warstwami.

~~~powershell
.\.venv\Scripts\python.exe -m examples.generate_demo
.\start.cmd DOSY_demo.jdf
~~~

Plik można też wybrać w aplikacji przez **Otwórz DOSY**. Dla szybkiej próby kliknij **Dodaj całkę ręcznie** i wpisz `3,60; 3,68`. Powtórz czynność dla `4,66; 4,74`. Granice są podane w ppm. Oczekiwane **D** wynoszą w przybliżeniu odpowiednio **3,8 × 10⁻¹⁰ m²/s** i **8,0 × 10⁻¹⁰ m²/s**. W tabeli głównej kolumna `D / 10⁻¹⁰` pokaże więc wartości w okolicy **3,8** i **8,0**. Generator odmawia nadpisania istniejącego JDF. Przy kolejnym generowaniu można podać inną nazwę jako argument.

## Podręcznik

[Instrukcja użytkownika](docs/INSTRUKCJA.md) prowadzi od otwarcia JDF do eksportu. Zawiera opis gestów, parametrów processingu, fazowania, bazy, całek, fitów, statystyki i interpretacji wyników.

- [Korekcja fazy i gesty myszy](docs/INSTRUKCJA.md#3-korekcja-fazy)
- [Pochodzenie b i ręczny import](docs/INSTRUKCJA.md#6-wartości-b)
- [Statystyka D i diagnostyka](docs/INSTRUKCJA.md#8-statystyka-d)
- [Projekty i pliki eksportu](docs/INSTRUKCJA.md#9-projekty-i-eksport)

## Format danych i zakres obsługi

- **Wejście** to surowy JEOL JDF z zespolonym FID w wymiarze bezpośrednim i kolejnymi warstwami w wymiarze gradientowym. Czytnik obsługuje układy `small_two_d` i `two_d`, liczby 32- i 64-bitowe oraz obie kolejności bajtów. Pojedynczy obsługiwany FID można oglądać jak widmo 1D.
- **Automatyczne b** pochodzi z rozpoznanej osi **b** w JDF albo z modelu rozpoznanych fabrycznych wariantów JEOL ¹H BPP-STE-LED i BPP-DSTE-LED. Model sprawdza zapisany `st_coef`, jeśli jest dostępny. Dla pozostałych sekwencji można podać niezależnie wyznaczone **b** dla każdej warstwy.
- **Zakres aplikacji** obejmuje klasyczne DOSY 1D. Przetworzone JDF oraz eksperymenty 3D COSY-DOSY, COSY-iDOSY i HMQC-DOSY są poza tym zakresem.
- **Wynik** zależy od prawidłowej kalibracji gradientów, poprawnego modelu **b**, jakości całek oraz adekwatności modelu zaniku. SE pojedynczego dopasowania i rozrzut **D** między sygnałami opisują różne wielkości.

## Kod i testy

| Pliki | Rola |
| --- | --- |
| `app.py`, `ui.py`, `spectrum_view.py` | Uruchamianie aplikacji, główne okno i obsługa widma |
| `jeol.py`, `core.py` | Odczyt JDF, model **b**, processing, całkowanie i dopasowanie |
| `phase_window.py`, `processing_panel.py`, `spectrum_settings.py` | Fazowanie, ustawienia processingu i wygląd widma |
| `d_statistics.py`, `statistics_window.py`, `statistics_svg.py` | Statystyka, diagnostyka ROI i wykresy |
| `exports.py`, `experiment_info.py`, `fullscreen_plot.py` | Projekty i eksport, metadane oraz widoki wykresów |
| `examples/generate_demo.py`, `test_*.py` | Syntetyczny JDF i testy obliczeń oraz interfejsu |

~~~powershell
.\.venv\Scripts\python.exe -m unittest discover -v
~~~

Testy interfejsu uruchamiają Qt w trybie offscreen. Publiczny zestaw używa danych syntetycznych. Dodatkowe testy regresji na lokalnych plikach pomiarowych są pomijane, kiedy tych plików brakuje.

## Licencja

Projekt jest udostępniony na [licencji MIT](LICENSE).
