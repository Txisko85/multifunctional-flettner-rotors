# Multifunctional Flettner rotors: reduced-order models and data

Code and data behind the article *Multifunctional Flettner rotors for wind-assisted ship propulsion: energy and
economic assessment of extension, turbine and Savonius-conversion topologies*, by Francisco Javier Asensio, Gaizka
Saldaña, Mikel González-Pérez, Oier Oñederra and Ander Ordoño, submitted to *Applied Energy*. The full reference
will be added when the article is published.

The models compare a conventional Flettner rotor with rotors that take extra functions (a vertical-axis wind
turbine that assists their spin or generates electricity, a telescopic extension, or opening into a Savonius
turbine) and with non-patented alternatives (a taller, a wider, a tiltable and a fifth rotor, and a tiltable rotor with an
extension). The case study is four 4 m × 24 m rotors on the ro-pax ferry *Galicia* on the Portsmouth–Bilbao route, at
sea and at berth, in its April–October service, with the 2026 timetable and hourly ERA5 winds for 2015–2024. A reduced-order
power balance converts every energy flow into fuel where it enters, including the wind resistance of the ship and
the resistance induced by lateral forces through a sway–yaw balance with rudder. The ship model is calibrated against the ship's EU MRV reports, and the
rotor model is checked against the full-scale figure of the Viking Grace, to which a further coefficient scenario
(A-VG) is fitted. The repository also contains the turbine
model (double-multiple streamtube) and its validation, a lifting-line model of deck-mounted rotors checked against URANS
results, a sway–yaw balance with rudder, the hourly simulation of optimal, causal and predictive switching between modes, the economic
and regulatory valuation (EU ETS, UK ETS, FuelEU Maritime, with the wind reward factor computed by the method of
MEPC.1/Circ.896) and a Sobol sensitivity analysis.

## Contents

| Path | Contents |
| --- | --- |
| `ruta/` | Route model and every analysis of the article, with their inputs and outputs (JSON), and the figure scripts |
| `ruta/datos_coef/` | Rotor coefficients digitised from published figures, and the global wind probability matrix of MEPC.1/Circ.896 (Appendix 2, Table 1) |
| `ruta/datos_mapa/` | Coastlines for the route map, from Natural Earth (public domain), clipped around the map |
| `turbina/` | Turbine model (double-multiple streamtube), its validation against wind-tunnel and field data, and the coefficient maps that the route model reads |
| `figuras/` | Scripts for the two schematic figures (Figs. 1 and 2) and the graphical abstract |
| `datos/era5_horario_ruta.csv.gz` | Hourly ERA5 wind at the 10 route points and the two ports of the Galicia, 2015–2024 (8.5 MB) |
| `datos/era5_horario_viking_grace.csv.gz` | Hourly ERA5 wind at the 10 points of the route of the Viking Grace, 2015–2024 (7.1 MB) |
| `descargar_era5_horario.py`, `descargar_era5_viking_grace.py` | Download those files again from the Open-Meteo API and write them next to the script (about 10 minutes each) |

Code comments, variable names and console output are in Spanish; this README and the figures are in English.

## Requirements

Python 3.11 or later and the packages in `requirements.txt`:

    pip install -r requirements.txt

`tzdata` provides the time zones of the timetable on Windows. Only Fig. 2 needs a browser:
`pip install playwright` and `python -m playwright install chromium`. `ruta_recorte_costa.py`, which only rebuilds the
clipped coastlines of the route map, needs `shapely`. The article's figures use the Inter typeface
(https://rsms.me/inter/); the scripts use it if it is installed and Matplotlib's default font otherwise.

## Reproducing the results

Run the scripts from their own folder (`cd turbina` or `cd ruta`). Each script overwrites the output file listed;
the repository includes those files as used in the article. The turbine scripts come first only if the turbine
maps are to be rebuilt; the route model reads `turbina/mapa_turbina.json` as included. Section, table and figure
numbers refer to the article; S numbers to its Supplementary material.

| Script (in `turbina/`) | What it does | Output | In the article |
| --- | --- | --- | --- |
| `validacion_dmst.py` | Turbine model against the wind-tunnel data of Battisti et al. (2018) and the field peak of Kjellin et al. (2011) (about 7 minutes; `--figura` redraws the figure from the JSON) | `validacion_dmst.json` | Section S2; Fig. S1 |
| `mapa_turbina.py` | Power, thrust and side-force coefficients of the case-study turbine against tip-speed ratio and wind speed (about 15 minutes) | `mapa_turbina.json` | Section 3.2; Section S2 |
| `dmst.py` | The model itself; run alone, it prints the case-study turbine's coefficients | console | Section S2 |

| Script (in `ruta/`) | What it does | Output | In the article |
| --- | --- | --- | --- |
| `coeficientes.py` | Rotor coefficient sets and scenarios A, A-VG, A+, AB, B, B′ and B″; ranges of the sampled inputs | (module) | Section 3.2; Tables 3 and 5 |
| `ruta_modelo.py` | The route model (module); run alone, it prints the reference case | console | Sections 3.1–3.7 |
| `ruta_geometria.py` | Route legs and the 10 ERA5 sampling points | `ruta.json` | Section 3.7 |
| `ruta_recorte_costa.py FOLDER` | Clips the Natural Earth coastlines around the route map (optional; `FOLDER` holds `ne_10m_land.geojson` and `ne_10m_minor_islands.geojson` from https://github.com/nvkelso/natural-earth-vector, folder `geojson/`) | `datos_mapa/ne_10m_tierra.geojson` | Fig. 4 |
| `era5_horario.py ../datos/era5_horario_ruta.csv.gz` | Rebuilds the all-hours wind histograms of `era5_ruta.json` (earlier versions) from the hourly file | `era5_ruta_desde_horario.json` | Section 3.7 |
| `era5_horario_condicionado.py ../datos/era5_horario_ruta.csv.gz` | Wind histograms conditioned on the timetable: by crossing and month, only the hours and points the ship sails through, and the berth hours at each port; the route model reads them | `era5_ruta_horario.json` | Section 3.7 |
| `ruta_calibracion_mrv.py` | Ship without rotors against its EU MRV reports at equal speed; sea margin | `calibracion_mrv.json` | Section 4.1; Table 6 |
| `validacion_viking_grace.py` | Rotor model on the route of the Viking Grace, by coefficient scenario, with sensitivities (about 2 minutes) | `validacion_viking_grace.json` | Sections 3.8 and 4.1; Tables 3 and 6; Section S6 |
| `linea_sustentadora.py` | Lifting-line model: lift factors η<sub>B</sub>, f<sub>X</sub> and f<sub>30</sub>, interaction ε, inflow to the turbine | `linea_sustentadora.json` | Section 3.3; Section S3 |
| `linea_sustentadora_validacion.py` | Lifting-line model against the URANS aspect-ratio results of De Marco et al. (2016); empirical f<sub>30</sub> and the resulting range of f<sub>X</sub> | `linea_sustentadora_validacion.json` | Section 3.3; Section S3 |
| `ruta_analisis.py` | Reference case, comparators, rules, break-even interaction ε*, scenarios and sensitivities (several minutes) | `resultados_ruta.json` | Sections 4.3 and 4.4; Tables 3, 7 and 8; Section S7 |
| `ruta_incertidumbre.py` | Sobol indices (Saltelli and Jansen estimators, N = 1,024, 18 inputs) and percentiles within scenario A, ranking probabilities in every scenario (about 80 minutes on two processes; it saves partial results and resumes if interrupted); `--muestras` repeats only the analysis from the saved samples | `incertidumbre_ruta.json`, `incertidumbre_muestras.npz` | Section 4.8; Section S8 |
| `tablas_sobol_tex.py` | LaTeX tables of the Supplementary material from `incertidumbre_ruta.json` | console | Tables S10–S13 |
| `ruta_horaria.py ../datos/era5_horario_ruta.csv.gz` | Hour-by-hour simulation of 2,087 crossings: switching costs, optimal, persistence and predictive control, year-to-year variability, berth generation | `resultados_horarios.json` | Section 4.5; Table 9; Section S5 |
| `ruta_horaria_diag.py ../datos/era5_horario_ruta.csv.gz` | Hourly simulation against the histograms (departure time, binning) | `diagnostico_horario.json` | Section 4.5 |
| `fueleu_viento.py` | FuelEU wind reward factor by the method of MEPC.1/Circ.896 (available effective power over the most favourable half of the global wind matrix) | `fueleu_viento.json` | Sections 3.9 and 4.7; Table 11; Section S9 |
| `economia.py` | Value of each option under fuel prices, EU ETS, UK ETS and FuelEU Maritime; largest justified investment; value of the FuelEU wind reward, EEXI and CII; rotor cost from EMSA | `economia.json` | Sections 3.9 and 4.7; Tables 10 and 11 |
| `ruta_momentos.py` | Base bending moments, structural limit of the extension, heeling moment and air draft | `momentos.json` | Section 4.10; Table 12 |
| `maniobra.py` | Sway–yaw balance with rudder (Kijima and Nakiri hull derivatives, MMG rudder model): the quadratic form of the induced resistance that the route model uses, the equivalent e<sub>H</sub> and the strong cases | `maniobra.json` | Sections 3.4 and 4.10; Section S4 |
| `ruta_mapa_adimensional.py` | Gain of each added function against wind-to-ship speed ratio and wind angle, at 18 and 12 kn | `mapa_adimensional.json` | Section 4.9; Fig. 9 |
| `inferencia_e.py` | Hull–rudder efficiency e<sub>H</sub> inferred from Bickert et al. (2025), for the lifting-surface sensitivity | console | Section 4.3 |

The figure scripts write `fig_*.pdf` and `fig_*.png` next to themselves:

| Script | Figure |
| --- | --- |
| `figuras/fig_topologias.py` | Fig. 1 |
| `figuras/fig_flujo.py` | Fig. 2 |
| `ruta/ruta_figura_coeficientes.py` | Fig. 3 |
| `ruta/ruta_figura_mapa.py ../datos/era5_horario_ruta.csv.gz` | Fig. 4 |
| `ruta/ruta_figura_validacion.py` (after `validacion_viking_grace.py`) | Fig. 5 |
| `ruta/ruta_figura_viento.py` | Fig. 6 |
| `ruta/ruta_figura_horaria.py` (after `ruta_horaria.py`) | Fig. 7 |
| `ruta/ruta_figura_puerto.py` | Fig. 8 |
| `ruta/ruta_mapa_adimensional.py` | Fig. 9 |
| `turbina/validacion_dmst.py` | Fig. S1 |
| `figuras/fig_resumen_grafico.py` (reads the outputs of `ruta_analisis.py` and `ruta_incertidumbre.py`) | Graphical abstract (also `fig_resumen_grafico.tif`, 600 dpi) |

Reference values for checking a run (scenario A, April–October service):

- Ship without rotors: 11,400 t of fuel a year. Four conventional rotors save 245.8 t/yr, after paying 27.0 t/yr of
  induced resistance; the telescopic extension (T3x) adds 66.7 t/yr and the retractable turbine (T4) 1.9 t/yr
  (`resultados_ruta.json`). With scenario A-VG, the four rotors save 403.3 t/yr.
- Rotor of the Viking Grace: mean thrust power 141.0 kW in scenario A; lift factor of scenario A-VG 1.25
  (`validacion_viking_grace.json`).
- Hourly simulation, 1,223 April–October crossings: the four rotors save 246.0 t/yr (`resultados_horarios.json`).
- Sobol analysis within scenario A: the four rotors save 202–298 t/yr (5–95 %, median 246 t/yr); total Sobol index of
  the extension's lift factor f<sub>X</sub> for the extension, 0.83 (`incertidumbre_ruta.json`).
- FuelEU wind reward: P<sub>wind</sub>/P<sub>prop</sub> = 0.071 for four rotors in scenario A, f<sub>wind</sub> = 0.99
  (`fueleu_viento.json`).

## Data sources and licences

The code is released under the MIT licence (`LICENSE`). The data files are released under CC BY 4.0 and include
third-party data that must be credited as set out in `DATA_LICENSE.md`: ERA5 (Copernicus Climate Change Service,
ECMWF) retrieved through Open-Meteo, the ship's EU MRV reports (EMSA), rotor coefficients digitised from
Bordogna et al. (2019), Seifert (2012) and De Marco et al. (2016), the wind-tunnel data of Battisti et al. (2018)
and the NACA 0021 section data of Sheldahl and Klimas (1981) as distributed with Sandia's CACTUS code (BSD 3-Clause
licence, `turbina/datos/LICENSE_CACTUS.txt`). Prices and regulatory parameters in `economia.py` come from the
sources cited in the article.

## Citation

Please cite the article and this repository. `CITATION.cff` gives the repository reference; the DOI of each
release is shown on its Zenodo record.

## Funding

This work was supported by the Basque Government through the GISEL research group of the University of the Basque
Country (UPV/EHU) [grant number IT1799-26].

## Use of AI

The code was written with the assistance of Claude (Anthropic). The authors reviewed it and take full responsibility
for it.
