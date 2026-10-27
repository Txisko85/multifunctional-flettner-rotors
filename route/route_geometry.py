"""Geometry of the Portsmouth-Bilbao route (approximate, public coastal-navigation waypoints).
Computes great-circle distances, headings and ERA5 sampling points every ~60 nm.
Output: route.json. The 10 sampling points are those of download_era5_route.py (C1-C4 in the Channel, B1-B6 in the
Bay of Biscay), their weights in nm are those of era5_route.json and the sums by group (Channel 208.4 nm, Bay of
Biscay 331.4 nm) are NM_OPEN in route_model.py."""
import math, json
R_NM = 3440.065
WP = [  # (name, lat, lon)
    ("Portsmouth", 50.790, -1.110),
    ("Nab", 50.670, -0.950),
    ("St Catherine S", 50.500, -1.300),
    ("Start Point S", 49.950, -3.500),
    ("Ushant SW", 48.650, -5.750),
    ("Bilbao approach", 43.450, -3.080),
    ("Bilbao", 43.360, -3.030),
]
def gc(lat1, lon1, lat2, lon2):
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dl = math.radians(lon2 - lon1)
    d = 2 * math.asin(math.sqrt(math.sin((p2 - p1) / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2))
    brg = math.degrees(math.atan2(math.sin(dl) * math.cos(p2), math.cos(p1) * math.sin(p2) - math.sin(p1) * math.cos(p2) * math.cos(dl))) % 360
    return d * R_NM, brg
def interp(lat1, lon1, lat2, lon2, f):
    return lat1 + f * (lat2 - lat1), lon1 + f * (lon2 - lon1)   # enough for short legs
legs, total = [], 0.0
for (n1, a1, o1), (n2, a2, o2) in zip(WP[:-1], WP[1:]):
    d, b = gc(a1, o1, a2, o2)
    legs.append({"from": n1, "to": n2, "nm": round(d, 1), "heading_south": round(b, 1)})
    total += d
print(json.dumps(legs, indent=1, ensure_ascii=False)); print("Total (nm):", round(total, 1))
# sampling points: centre of sub-legs of ~60 nm in the open-sea legs (St Catherine -> Bilbao approach)
pts = []
for (n1, a1, o1), (n2, a2, o2) in zip(WP[2:-2], WP[3:-1]):
    d, b = gc(a1, o1, a2, o2)
    n = max(1, round(d / 60))
    for k in range(n):
        f = (k + 0.5) / n
        la, lo = interp(a1, o1, a2, o2, f)
        pts.append(dict(lat=round(la, 3), lon=round(lo, 3), heading_south=round(b, 1), weight_nm=round(d / n, 2)))
print(json.dumps(pts, indent=1)); print("Points:", len(pts), " sum nm:", round(sum(p['weight_nm'] for p in pts), 1))
json.dump(dict(legs=legs, total_nm=total, points=pts, ports=[dict(name="Portsmouth", lat=50.79, lon=-1.11), dict(name="Bilbao", lat=43.36, lon=-3.03)]), open("route.json", "w"), indent=1, ensure_ascii=False)
