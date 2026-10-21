#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Recorta la costa de Natural Earth a la zona del mapa de la ruta (A1, 21-10-2026).

Natural Earth, escala 1:10 m, dominio público (https://www.naturalearthdata.com/), ficheros de la carpeta geojson/
de https://github.com/nvkelso/natural-earth-vector (commit ca96624): ne_10m_land.geojson y
ne_10m_minor_islands.geojson.
Recorta a 11°W–3°E y 41,5°N–52,5°N (más que la zona de estudio, que el mapa recorta después), simplifica con una
tolerancia de 0,002° (unos 200 m, invisible a la escala del mapa) y redondea las coordenadas a 1e-4°.
El mapa no dibuja fronteras: la política de mapas de Elsevier pide que muestren solo la zona de estudio.
Uso: python ruta_recorte_costa.py CARPETA_CON_LOS_GEOJSON
Salida: datos_mapa/ne_10m_tierra.geojson. Necesita shapely; la figura (ruta_figura_mapa.py) no lo necesita.
"""
import json
import os
import sys
from shapely.geometry import box, mapping, shape
from shapely.ops import unary_union

HERE = os.path.dirname(os.path.abspath(__file__))
BBOX = box(-11.0, 41.5, 3.0, 52.5)
TOL = 0.002


def rnd(obj):
    if isinstance(obj, (list, tuple)):
        if obj and isinstance(obj[0], (int, float)):
            return [round(float(v), 4) for v in obj]
        return [rnd(o) for o in obj]
    return obj


def load(path):
    return json.load(open(path, encoding="utf-8"))["features"]


def main():
    src = sys.argv[1]
    out_dir = os.path.join(HERE, "datos_mapa")
    os.makedirs(out_dir, exist_ok=True)
    polys = []
    for name in ("ne_10m_land.geojson", "ne_10m_minor_islands.geojson"):
        for f in load(os.path.join(src, name)):
            g = shape(f["geometry"])
            if g.intersects(BBOX):
                polys.append(g.intersection(BBOX))
    land = unary_union(polys).simplify(TOL, preserve_topology=True)
    feats = [dict(type="Feature", properties={}, geometry=dict(type=g["type"], coordinates=rnd(g["coordinates"])))
             for g in [mapping(land)]]
    note = ("Made with Natural Earth (public domain), 1:10 m land and minor islands, "
            "github.com/nvkelso/natural-earth-vector commit ca96624; clipped to 11W-3E, 41.5N-52.5N, simplified 0.002 deg")
    out = os.path.join(out_dir, "ne_10m_tierra.geojson")
    json.dump(dict(type="FeatureCollection", name="ne_10m_tierra", note=note, features=feats), open(out, "w"),
              separators=(",", ":"))
    print("ne_10m_tierra.geojson", os.path.getsize(out), "bytes")


if __name__ == "__main__":
    main()
