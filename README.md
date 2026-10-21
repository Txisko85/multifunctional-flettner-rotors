# Multifunctional Flettner rotors: reduced-order models and data

Code and data behind the article *Multifunctional Flettner rotors for wind-assisted ship propulsion: a comparative
energy assessment of extension, turbine and Savonius-conversion topologies*, by F.J. Asensio, G. Saldaña,
M. González-Pérez, O. Oñederra and A. Ordoño, submitted to *Applied Energy*. The full reference will be added when
the article is published.

The models compare a conventional Flettner rotor with rotors that take extra functions: a vertical-axis wind turbine
that assists their spin or generates electricity, a telescopic extension, or opening into a Savonius turbine. The
case study is four 4 m × 24 m rotors on the ro-pax ferry *Galicia* on the Portsmouth–Bilbao route, at sea and at
berth, with the 2026 timetable and hourly ERA5 winds for 2015–2024. A reduced-order power balance converts every
energy flow into fuel where it enters, including the induced resistance of the lateral force. The turbine's power,
thrust and side-force coefficients come from a double-multiple streamtube model validated against wind-tunnel and
field data. Optimal switching between modes is compared with wind-direction rules, and coefficient uncertainty is
sampled by Monte Carlo.

## Contents

| Path | Contents |
| --- | --- |
| `ruta/` | Route model, hourly simulation, Monte Carlo analysis and figure scripts, with their inputs and outputs (JSON) |
| `ruta/datos_coef/` | Rotor coefficients digitised from published figures |
| `ruta/datos_mapa/` | Coastlines for the route map, from Natural Earth (public domain), clipped around the map |
| `turbina/` | Turbine model (double-multiple streamtube), its validation against wind-tunnel and field data, and the coefficient maps that the route model reads |
| `figuras/` | Scripts for the two schematic figures (Figs. 1 and 2) |
| `datos/era5_horario_ruta.csv.gz` | Hourly ERA5 wind at the 10 route points and the two ports, 2015–2024 (8.5 MB) |
| `descargar_era5_horario.py` | Downloads that file again from the Open-Meteo API and writes it next to the script (about 10 minutes) |

Code comments, variable names and console output are in Spanish; this README and the figures are in English.

## Requirements

Python 3.11 or later and the packages in `requirements.txt`:

    pip install -r requirements.txt

`tzdata` provides the time zones of the timetable on Windows. Only Fig. 2 needs a browser:
`pip install playwright` and `python -m playwright install chromium`. `ruta_recorte_costa.py`, which only rebuilds the
clipped coastlines of the route map, needs `shapely`. The article's figures use the Inter typeface
(https://rsms.me/inter/); the scripts use it if it is installed and Matplotlib's default font otherwise.

The scripts were checked from a clean copy of this repository with Python 3.11 and 3.13: every output file is
identical to the one included, and the figures match those of the article.

## Reproducing the results

Run the scripts from their own folder (`cd turbina` or `cd ruta`). Each script overwrites the output file listed;
the repository includes those files as used in the article. The turbine scripts come first only if the turbine
maps are to be rebuilt; the route model reads `turbina/mapa_turbina.json` as included.

| Script (in `turbina/`) | What it does | Output | In the article |
| --- | --- | --- | --- |
| `validacion_dmst.py` | Turbine model against the wind-tunnel data of Battisti et al. (2018) and the field peak of Kjellin et al. (2011) (about 7 minutes; `--figura` redraws the figure from the JSON) | `validacion_dmst.json` | Appendix B; Fig. B.1 |
| `mapa_turbina.py` | Power, thrust and side-force coefficients of the case-study turbine against tip-speed ratio and wind speed (about 15 minutes) | `mapa_turbina.json` | Section 3.2; Appendix B |
| `dmst.py` | The model itself; run alone, it prints the case-study turbine's coefficients | console | Appendix B |

| Script (in `ruta/`) | What it does | Output | In the article |
| --- | --- | --- | --- |
| `ruta_geometria.py` | Route legs and the 10 ERA5 sampling points | `ruta.json` | Section 3.6 |
| `ruta_recorte_costa.py FOLDER` | Clips the Natural Earth coastlines around the route map (optional; `FOLDER` holds `ne_10m_land.geojson` and `ne_10m_minor_islands.geojson` from https://github.com/nvkelso/natural-earth-vector, folder `geojson/`) | `datos_mapa/ne_10m_tierra.geojson` | Fig. 4 |
| `ruta_analisis.py` | Central case, step-by-step changes from the first model version, sensitivity cases | `resultados_ruta.json` | Sections 4.2, 4.3 and 4.5; Tables 5 and 6 |
| `ruta_montecarlo.py` | Monte Carlo analysis, 1,000 samples on two processes (several minutes) | `montecarlo_ruta.json` | Section 4.6; Tables 5 and 6 |
| `ruta_calibracion_mrv.py` | Ship without rotors against its EU MRV reports | console | Section 3.6 |
| `inferencia_e.py` | Hull–rudder efficiency *e*<sub>H</sub> inferred from Bickert et al. (2025) | console | Section 3.3; Table 3 |
| `ruta_momentos.py` | Bending moments at the rotor base | `momentos.json` | Section 4.7; Table 8 |
| `era5_horario.py ../datos/era5_horario_ruta.csv.gz` | Rebuilds the wind histograms of `era5_ruta.json` from the hourly file | `era5_ruta_desde_horario.json` | Section 3.6 |
| `ruta_horaria.py ../datos/era5_horario_ruta.csv.gz` | Hour-by-hour simulation of 2,087 crossings: switching cost, interannual variability, berth generation | `resultados_horarios.json` | Section 4.4; Table 7 |
| `ruta_horaria_diag.py ../datos/era5_horario_ruta.csv.gz` | Hourly simulation against the histograms (departure time, binning) | `diagnostico_horario.json` | Section 4.4 |

The figure scripts write `fig_*.pdf` and `fig_*.png` next to themselves:

| Script | Figure |
| --- | --- |
| `figuras/fig_topologias.py` | Fig. 1 |
| `figuras/fig_flujo.py` | Fig. 2 |
| `ruta/ruta_figura_coeficientes.py` | Fig. 3 |
| `ruta/ruta_figura_mapa.py ../datos/era5_horario_ruta.csv.gz` | Fig. 4 |
| `ruta/ruta_figura_viento.py` | Fig. 5 |
| `ruta/ruta_figura_horaria.py` (after `ruta_horaria.py`) | Fig. 6 |
| `ruta/ruta_figura_puerto.py` | Fig. 7 |
| `turbina/validacion_dmst.py` | Fig. B.1 |

Reference values for checking a run (central coefficients, year-round service):

- Ship without rotors: 17,236 t of fuel a year. Four conventional rotors save 332.4 t/yr, after paying
  98.3 t/yr of induced resistance (`resultados_ruta.json`).
- Hourly simulation over 2,087 crossings: the four rotors save 322.7 t/yr (`resultados_horarios.json`).
- Monte Carlo, 5–95 % interval: 113–587 t/yr (`montecarlo_ruta.json`).
- Retractable turbine (T4), relative to the conventional rotors: 5.3 t/yr; telescopic extension (T3x): 64.1 t/yr
  (`resultados_ruta.json`).

## Data sources and licences

The code is released under the MIT licence (`LICENSE`). The data files are released under CC BY 4.0 and include
third-party data that must be credited as set out in `DATA_LICENSE.md`: ERA5 (Copernicus Climate Change Service,
ECMWF) retrieved through Open-Meteo, the ship's EU MRV reports (EMSA), rotor coefficients digitised from
Bordogna et al. (2019), Seifert (2012) and De Marco et al. (2016), the wind-tunnel data of Battisti et al. (2018)
and the NACA 0021 section data of Sheldahl and Klimas (1981) as distributed with Sandia's CACTUS code (BSD 3-Clause
licence, `turbina/datos/LICENSE_CACTUS.txt`).

## Citation

Please cite the article and this repository. `CITATION.cff` gives the repository reference; the DOI of each
release is shown on its Zenodo record.

## Funding

This work was supported by the Basque Government through the GISEL research group of the University of the Basque
Country (UPV/EHU) [grant number IT1799-26].

## Use of AI

The code was written with the assistance of Claude (Anthropic). The authors reviewed it and take full responsibility
for it.
