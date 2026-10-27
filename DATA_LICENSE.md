# Data licence and attributions

The data files in this repository (`data/`, and the JSON files in `route/`, `route/coef_data/`, `turbine/` and
`turbine/data/`) are released under the Creative Commons Attribution 4.0 International licence (CC BY 4.0,
https://creativecommons.org/licenses/by/4.0/), except the NACA 0021 section data, which keep the BSD 3-Clause
licence of their source (see Turbine data), and the coastlines in `route/map_data/`, which are in the public domain
(see Coastlines). When you reuse them, cite this repository and credit the original sources below.

## ERA5 wind

Files: `data/era5_hourly_route.csv.gz`, `data/era5_hourly_viking_grace.csv.gz`, `route/era5_route.json`,
`route/era5_route_from_hourly.json`, `route/era5_route_conditioned.json`.

ERA5 hourly reanalysis (Hersbach et al., 2020, https://doi.org/10.1002/qj.3803), © European Centre for Medium-Range
Weather Forecasts (ECMWF), produced by the Copernicus Climate Change Service and licensed under CC BY 4.0.
Retrieved through the Open-Meteo Historical Weather API with `models=era5` (weather data by Open-Meteo.com,
https://open-meteo.com/; Zippenfenig, 2023, https://doi.org/10.5281/zenodo.7970649), licensed under CC BY 4.0.

Changes: wind speed at 10 m and 100 m and wind direction at 100 m, 2015–2024, written to CSV, at 12 locations on the
route of the Galicia and 10 on the route of the Viking Grace.
`route/era5_route.json` and `route/era5_route_from_hourly.json` aggregate these series into wind histograms by route
group and season, and `route/era5_route_conditioned.json` into histograms of the hours of the timetable, by crossing
and month.

## EU MRV reports of the ship

File: `route/mrv_galicia.json`.

Source: European Maritime Safety Agency (EMSA), THETIS-MRV public CO2 emission reports for the GALICIA
(IMO 9856189), reporting periods 2020–2025, https://mrv.emsa.europa.eu/#public/emission-report (retrieved
18 October 2026). EMSA allows reproduction provided the source is acknowledged.

## Digitised rotor coefficients

The following values were read from published figures. They approximate the published data; the original
publications remain the reference and should be cited.

- `route/coef_data/bordogna2019_fig4_fig7.json`: lift and drag coefficients from Figs. 4 and 7 of
  Bordogna G, Muggiasca S, Giappino S, Belloli M, Keuning JA, Huijsmans RHM, et al. Experiments on a Flettner rotor
  at critical and supercritical Reynolds numbers. J Wind Eng Ind Aerodyn 2019;188:19–29.
  https://doi.org/10.1016/j.jweia.2019.02.006
- `route/coef_data/ackeret_seifert_fig27.json`: Ackeret's data as reproduced in Fig. 27 of Seifert J. A review of the
  Magnus effect in aeronautics. Prog Aerosp Sci 2012;55:17–45. https://doi.org/10.1016/j.paerosci.2012.07.001
- `route/coef_data/demarco2016_fig10.json` (and the values used in `route/coefficients.py`): lift and drag
  coefficients read from Fig. 10 of De Marco A, Mancini S, Pensa C, Calise G, De Luca F. Flettner rotor concept for
  marine applications: a systematic study. Int J Rotating Mach 2016;2016:3458750.
  https://doi.org/10.1155/2016/3458750

## Turbine data

- `turbine/data/NACA_0021_SK.dat`: lift and drag of the NACA 0021 section through 180° at eleven Reynolds numbers,
  from Sheldahl RE, Klimas PC. Aerodynamic characteristics of seven symmetrical airfoil sections through 180-degree
  angle of attack for use in aerodynamic analysis of vertical axis wind turbines. Sandia National Laboratories,
  SAND80-2114; 1981 (OSTI 6548367). The file is copied unchanged from the CACTUS code of Sandia National
  Laboratories (https://github.com/sandialabs/CACTUS, `test/Airfoil_Section_Data/NACA_0021.dat`, commit 4ae4faf),
  distributed under the BSD 3-Clause licence reproduced in `turbine/data/LICENSE_CACTUS.txt`. Sheldahl and Klimas
  derived the NACA 0021 values from measurements on thinner sections; they are not direct measurements.
- `turbine/data/battisti2018_table8.json`: Table 8 of Battisti L, Persico G, Dossena V, Paradiso B, Raciti
  Castelli M, Brighenti A, Benini E. Experimental benchmark data for H-shaped and troposkien VAWT architectures.
  Renew Energy 2018;125:425–444. https://doi.org/10.1016/j.renene.2018.02.098, transcribed from the accepted
  manuscript in the repository of the Politecnico di Milano. The original publication remains the reference.

## Coastlines

File: `route/map_data/ne_10m_land_clipped.geojson`.

Made with Natural Earth. Free vector and raster map data @ naturalearthdata.com (public domain). Land and minor
islands at 1:10 m, from https://github.com/nvkelso/natural-earth-vector (folder `geojson/`, commit ca96624), clipped
to 11°W–3°E and 41.5°N–52.5°N and simplified by `route/clip_coastline.py`.

## Model outputs

The other JSON files in `route/` and `turbine/` are outputs of the models in this repository.
