#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Clips the Natural Earth coastline to the area of the route map (A1, 21-10-2026).

Natural Earth, scale 1:10 m, public domain (https://www.naturalearthdata.com/), files of the geojson/ folder of
https://github.com/nvkelso/natural-earth-vector (commit ca96624): ne_10m_land.geojson and ne_10m_minor_islands.geojson.
Clips to 11°W–3°E and 41.5°N–52.5°N (more than the study area, which the map crops afterwards), simplifies with a
tolerance of 0.002° (about 200 m, invisible at the scale of the map) and rounds the coordinates to 1e-4°.
The map draws no borders: the map policy of Elsevier asks that maps show only the study area.
Usage: python clip_coastline.py FOLDER_WITH_THE_GEOJSON
Output: map_data/ne_10m_land_clipped.geojson. Needs shapely; the figure (fig_route_map.py) does not.
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
    out_dir = os.path.join(HERE, "map_data")
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
    out = os.path.join(out_dir, "ne_10m_land_clipped.geojson")
    json.dump(dict(type="FeatureCollection", name="ne_10m_land_clipped", note=note, features=feats), open(out, "w"),
              separators=(",", ":"))
    print("ne_10m_land_clipped.geojson", os.path.getsize(out), "bytes")


if __name__ == "__main__":
    main()
