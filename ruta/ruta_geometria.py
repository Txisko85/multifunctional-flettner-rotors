"""Geometria de la ruta Portsmouth-Bilbao (aproximada, waypoints publicos de navegacion costera).
Calcula distancias ortodromicas, rumbos y puntos de muestreo de ERA5 cada ~60 mn.
Salida: ruta.json. Los 10 puntos de muestreo son los de descargar_era5_horario.py (C1-C4 en el Canal, B1-B6 en
el golfo), sus pesos en mn son los de era5_ruta.json y las sumas por grupo (Canal 208,4 mn, golfo 331,4 mn) son
NM_OPEN en ruta_modelo.py."""
import math, json
R_NM = 3440.065
WP = [  # (nombre, lat, lon)
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
    return lat1 + f * (lat2 - lat1), lon1 + f * (lon2 - lon1)   # suficiente para tramos cortos
legs, total = [], 0.0
for (n1, a1, o1), (n2, a2, o2) in zip(WP[:-1], WP[1:]):
    d, b = gc(a1, o1, a2, o2)
    legs.append(dict(desde=n1, hasta=n2, mn=round(d, 1), rumbo_sur=round(b, 1)))
    total += d
print(json.dumps(legs, indent=1, ensure_ascii=False)); print("Total (mn):", round(total, 1))
# puntos de muestreo: centro de subtramos de ~60 mn en los tramos de mar abierto (St Catherine -> Bilbao approach)
pts = []
for (n1, a1, o1), (n2, a2, o2) in zip(WP[2:-2], WP[3:-1]):
    d, b = gc(a1, o1, a2, o2)
    n = max(1, round(d / 60))
    for k in range(n):
        f = (k + 0.5) / n
        la, lo = interp(a1, o1, a2, o2, f)
        pts.append(dict(lat=round(la, 3), lon=round(lo, 3), rumbo_sur=round(b, 1), peso_mn=round(d / n, 2)))
print(json.dumps(pts, indent=1)); print("Puntos:", len(pts), " suma mn:", round(sum(p['peso_mn'] for p in pts), 1))
json.dump(dict(legs=legs, total_nm=total, puntos=pts, puertos=[dict(nombre="Portsmouth", lat=50.79, lon=-1.11), dict(nombre="Bilbao", lat=43.36, lon=-3.03)]), open("ruta.json", "w"), indent=1, ensure_ascii=False)
