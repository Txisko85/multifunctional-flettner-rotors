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
| `route/` | Route model and every analysis of the article, with their inputs and outputs (JSON), and the figure scripts |
| `route/coef_data/` | Rotor coefficients digitised from published figures, and the global wind probability matrix of MEPC.1/Circ.896 (Appendix 2, Table 1) |
| `route/map_data/` | Coastlines for the route map, from Natural Earth (public domain), clipped around the map |
| `turbine/` | Turbine model (double-multiple streamtube), its validation against wind-tunnel and field data, and the coefficient maps that the route model reads |
| `figures/` | Scripts for the two schematic figures (Figs. 1 and 2) and the graphical abstract |
| `data/era5_hourly_route.csv.gz` | Hourly ERA5 wind at the 10 route points and the two ports of the Galicia, 2015–2024 (8.5 MB) |
| `data/era5_hourly_viking_grace.csv.gz` | Hourly ERA5 wind at the 10 points of the route of the Viking Grace, 2015–2024 (7.1 MB) |
| `download_era5_route.py`, `download_era5_viking_grace.py` | Download those files again from the Open-Meteo API into `data/` (about 10 minutes each) |

## Requirements

Python 3.11 or later and the packages in `requirements.txt`:

    pip install -r requirements.txt

`tzdata` provides the time zones of the timetable on Windows. Only Fig. 2 needs a browser:
`pip install playwright` and `python -m playwright install chromium`. `clip_coastline.py`, which only rebuilds the
clipped coastlines of the route map, needs `shapely`. The article's figures use the Inter typeface
(https://rsms.me/inter/); the scripts use it if it is installed and Matplotlib's default font otherwise.

## Reproducing the results

Run the scripts from their own folder (`cd turbine` or `cd route`). Each script overwrites the output file listed;
the repository includes those files as used in the article. The turbine scripts come first only if the turbine
maps are to be rebuilt; the route model reads `turbine/turbine_map.json` as included. Section, table and figure
numbers refer to the article; S numbers to its Supplementary material.

| Script (in `turbine/`) | What it does | Output | In the article |
| --- | --- | --- | --- |
| `dmst_validation.py` | Turbine model against the wind-tunnel data of Battisti et al. (2018) and the field peak of Kjellin et al. (2011) (about 7 minutes; `--figure` redraws the figure from the JSON) | `dmst_validation.json`, `turbine_validation.pdf` | Section S2; Fig. S1 |
| `turbine_map.py` | Power, thrust and side-force coefficients of the case-study turbine against tip-speed ratio and wind speed (about 15 minutes) | `turbine_map.json` | Section 3.2; Section S2 |
| `dmst.py` | The model itself; run alone, it prints the case-study turbine's coefficients | console | Section S2 |

| Script (in `route/`) | What it does | Output | In the article |
| --- | --- | --- | --- |
| `coefficients.py` | Rotor coefficient sets and scenarios A, A-VG, A+, AB, B, B′ and B″; ranges of the sampled inputs | (module) | Section 3.2; Tables 3 and 5 |
| `route_model.py` | The route model (module); run alone, it prints the reference case | console | Sections 3.1–3.7 |
| `route_geometry.py` | Route legs and the 10 ERA5 sampling points | `route.json` | Section 3.7 |
| `clip_coastline.py FOLDER` | Clips the Natural Earth coastlines around the route map (optional; `FOLDER` holds `ne_10m_land.geojson` and `ne_10m_minor_islands.geojson` from https://github.com/nvkelso/natural-earth-vector, folder `geojson/`) | `map_data/ne_10m_land_clipped.geojson` | Fig. 4 |
| `era5_hourly.py ../data/era5_hourly_route.csv.gz` | Rebuilds the all-hours wind histograms of `era5_route.json` (earlier versions) from the hourly file | `era5_route_from_hourly.json` | Section 3.7 |
| `era5_conditioned.py ../data/era5_hourly_route.csv.gz` | Wind histograms conditioned on the timetable: by crossing and month, only the hours and points the ship sails through, and the berth hours at each port; the route model reads them | `era5_route_conditioned.json` | Section 3.7 |
| `mrv_calibration.py` | Ship without rotors against its EU MRV reports at equal speed; sea margin | `mrv_calibration.json` | Section 4.1; Table 6 |
| `viking_grace_check.py` | Rotor model on the route of the Viking Grace, by coefficient scenario, with sensitivities (about 2 minutes) | `viking_grace_check.json` | Sections 3.8 and 4.1; Tables 3 and 6; Section S6 |
| `lifting_line.py` | Lifting-line model: lift factors η<sub>B</sub>, f<sub>X</sub> and f<sub>30</sub>, interaction ε, inflow to the turbine | `lifting_line.json` | Section 3.3; Section S3 |
| `lifting_line_validation.py` | Lifting-line model against the URANS aspect-ratio results of De Marco et al. (2016); empirical f<sub>30</sub> and the resulting range of f<sub>X</sub> | `lifting_line_validation.json` | Section 3.3; Section S3 |
| `route_analysis.py` | Reference case, comparators, rules, break-even interaction ε*, scenarios and sensitivities (several minutes) | `route_results.json` | Sections 4.3 and 4.4; Tables 3, 7 and 8; Section S7 |
| `uncertainty.py` | Sobol indices (Saltelli and Jansen estimators, N = 1,024, 18 inputs) and percentiles within scenario A, ranking probabilities in every scenario (about 80 minutes on two processes; it saves partial results and resumes if interrupted); `--samples` repeats only the analysis from the saved samples | `uncertainty.json`, `uncertainty_samples.npz` | Section 4.8; Section S8 |
| `sobol_tables_tex.py` | LaTeX tables of the Supplementary material from `uncertainty.json` | console | Tables S10–S13 |
| `hourly_voyages.py ../data/era5_hourly_route.csv.gz` | Hour-by-hour simulation of 2,087 crossings: switching costs, optimal, persistence and predictive control, year-to-year variability, berth generation | `hourly_results.json` | Section 4.5; Table 9; Section S5 |
| `hourly_diagnostics.py ../data/era5_hourly_route.csv.gz` | Hourly simulation against the histograms (departure time, binning) | `hourly_diagnostics.json` | Section 4.5 |
| `fueleu_wind.py` | FuelEU wind reward factor by the method of MEPC.1/Circ.896 (available effective power over the most favourable half of the global wind matrix) | `fueleu_wind.json` | Sections 3.9 and 4.7; Table 11; Section S9 |
| `economics.py` | Value of each option under fuel prices, EU ETS, UK ETS and FuelEU Maritime; largest justified investment; value of the FuelEU wind reward, EEXI and CII; rotor cost from EMSA | `economics.json` | Sections 3.9 and 4.7; Tables 10 and 11 |
| `moments.py` | Base bending moments, structural limit of the extension, heeling moment and air draft | `moments.json` | Section 4.10; Table 12 |
| `sway_yaw.py` | Sway–yaw balance with rudder (Kijima and Nakiri hull derivatives, MMG rudder model): the quadratic form of the induced resistance that the route model uses, the equivalent e<sub>H</sub> and the strong cases | `sway_yaw.json` | Sections 3.4 and 4.10; Section S4 |
| `function_map.py` | Gain of each added function against wind-to-ship speed ratio and wind angle, at 18 and 12 kn | `function_map.json` | Section 4.9; Fig. 9 |
| `infer_eh.py` | Hull–rudder efficiency e<sub>H</sub> inferred from Bickert et al. (2025), for the lifting-surface sensitivity | console | Section 4.3 |

The figure scripts write a PDF and a PNG next to themselves, with the file names used in the article's LaTeX source:

| Script | Output | Figure |
| --- | --- | --- |
| `figures/fig_topologies.py` | `topologies.pdf` | Fig. 1 |
| `figures/fig_workflow.py` | `workflow.pdf` (and `workflow.svg`) | Fig. 2 |
| `route/fig_coefficients.py` | `rotor_coefficients.pdf` | Fig. 3 |
| `route/fig_route_map.py ../data/era5_hourly_route.csv.gz` | `route_map.pdf` | Fig. 4 |
| `route/fig_fullscale.py` (after `viking_grace_check.py`) | `fullscale_check.pdf` | Fig. 5 |
| `route/fig_wind.py` | `wind_roses.pdf` | Fig. 6 |
| `route/fig_hourly.py` (after `hourly_voyages.py`) | `hourly_voyages.pdf` | Fig. 7 |
| `route/fig_berth.py` | `berth_generation.pdf` | Fig. 8 |
| `route/function_map.py` | `function_map.pdf` | Fig. 9 |
| `turbine/dmst_validation.py` | `turbine_validation.pdf` | Fig. S1 |
| `figures/fig_graphical_abstract.py` (reads the outputs of `route_analysis.py` and `uncertainty.py`) | `graphical_abstract.pdf` (also `.png` and a 600 dpi `.tif`) | Graphical abstract |

Reference values for checking a run (scenario A, April–October service):

- Ship without rotors: 11,400 t of fuel a year. Four conventional rotors save 245.8 t/yr, after paying 27.0 t/yr of
  induced resistance; the telescopic extension (T3x) adds 66.7 t/yr and the retractable turbine (T4) 1.9 t/yr
  (`route_results.json`). With scenario A-VG, the four rotors save 403.3 t/yr.
- Rotor of the Viking Grace: mean thrust power 141.0 kW in scenario A; lift factor of scenario A-VG 1.25
  (`viking_grace_check.json`).
- Hourly simulation, 1,223 April–October crossings: the four rotors save 246.0 t/yr (`hourly_results.json`).
- Sobol analysis within scenario A: the four rotors save 202–298 t/yr (5–95 %, median 246 t/yr); total Sobol index of
  the extension's lift factor f<sub>X</sub> for the extension, 0.83 (`uncertainty.json`).
- FuelEU wind reward: P<sub>wind</sub>/P<sub>prop</sub> = 0.071 for four rotors in scenario A, f<sub>wind</sub> = 0.99
  (`fueleu_wind.json`).

## Data sources and licences

The code is released under the MIT licence (`LICENSE`). The data files are released under CC BY 4.0 and include
third-party data that must be credited as set out in `DATA_LICENSE.md`: ERA5 (Copernicus Climate Change Service,
ECMWF) retrieved through Open-Meteo, the ship's EU MRV reports (EMSA), rotor coefficients digitised from
Bordogna et al. (2019), Seifert (2012) and De Marco et al. (2016), the wind-tunnel data of Battisti et al. (2018)
and the NACA 0021 section data of Sheldahl and Klimas (1981) as distributed with Sandia's CACTUS code (BSD 3-Clause
licence, `turbine/data/LICENSE_CACTUS.txt`). Prices and regulatory parameters in `economics.py` come from the
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
