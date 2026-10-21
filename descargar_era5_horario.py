#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Descarga del viento horario ERA5 en los 12 puntos del modelo de ruta del artículo A1
(ro-pax Galicia, Portsmouth-Bilbao): 10 puntos de mar y los dos puertos.

Fuente: Open-Meteo, Historical Weather API con models=era5 (los mismos datos que dieron era5_ruta.json).
Periodo 2015-01-01 a 2024-12-31, horas UTC. Variables: velocidad a 10 m y 100 m, dirección a 100 m.

Uso (Windows):  py descargar_era5_horario.py      (o: python descargar_era5_horario.py)
Uso (macOS/Linux): python3 descargar_era5_horario.py
Requisitos:  Python 3.8 o posterior, sin paquetes adicionales, y conexión a internet.
Duración:    unos 10 minutos. Espera 35 s entre peticiones: cada una pide 10 años, que la API gratuita
             cuenta como unas 260 llamadas (límites: 600 por minuto, 5 000 por hora, 10 000 por día).
Reanudable:  cada punto se guarda en la carpeta era5_cache; si se corta, al volver a lanzarlo sigue
             donde lo dejó.
Salida:      era5_horario_ruta.csv.gz junto al script (unos 8 MB).

Formato de la salida (CSV comprimido con gzip; las líneas que empiezan por '#' son metadatos):
  punto,h,u10,u100,dir100
  h = horas desde 2015-01-01 00:00 UTC; u10 y u100 en m/s; dir100 = procedencia a 100 m (grados);
  campo vacío = dato ausente.

Licencias: datos de Open-Meteo con licencia CC BY 4.0 (https://open-meteo.com/en/license);
ERA5 de Copernicus Climate Change Service (Hersbach et al. 2020, doi:10.1002/qj.3803).
"""
import csv
import gzip
import json
import os
import ssl
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone

PUNTOS = [  # id, latitud, longitud (mismos puntos que el modelo de ruta)
    ("C1", 50.362, -1.850), ("C2", 50.088, -2.950), ("C3", 49.625, -4.062), ("C4", 48.975, -5.188),
    ("B1", 48.217, -5.527), ("B2", 47.350, -5.082), ("B3", 46.483, -4.638), ("B4", 45.617, -4.192),
    ("B5", 44.750, -3.748), ("B6", 43.883, -3.303),
    ("PORTSMOUTH", 50.790, -1.110), ("BILBAO", 43.360, -3.030),
]
INICIO, FIN = "2015-01-01", "2024-12-31"
N_HORAS = 87672           # 3 653 días x 24 h
T0 = 1420070400           # 2015-01-01 00:00 UTC en segundos Unix
URL = ("https://archive-api.open-meteo.com/v1/archive?latitude={lat}&longitude={lon}"
       "&start_date=" + INICIO + "&end_date=" + FIN +
       "&hourly=wind_speed_10m,wind_speed_100m,wind_direction_100m"
       "&models=era5&wind_speed_unit=ms&timezone=GMT&timeformat=unixtime")
AQUI = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(AQUI, "era5_cache")
SALIDA = os.path.join(AQUI, "era5_horario_ruta.csv.gz")
PAUSA = 35  # s entre peticiones

try:  # certificados de certifi si está instalado (evita fallos SSL en algunas instalaciones)
    import certifi
    CTX = ssl.create_default_context(cafile=certifi.where())
except ImportError:
    CTX = ssl.create_default_context()


class Parar(Exception):
    pass


def pedir(url):
    req = urllib.request.Request(url, headers={"User-Agent": "A1-Flettner-UPV-EHU-research"})
    esperas_hora = 0
    for intento in range(12):
        try:
            with urllib.request.urlopen(req, timeout=300, context=CTX) as r:
                return json.loads(r.read().decode("utf-8"))
        except urllib.error.HTTPError as e:
            cuerpo = e.read().decode("utf-8", "replace")[:300]
            if e.code == 429:
                if "Daily" in cuerpo or "daily" in cuerpo:
                    raise Parar("la API ha alcanzado el límite diario de esta conexión. Vuelve a lanzar el "
                                "script mañana: seguirá donde lo dejó.")
                if "Hourly" in cuerpo or "hourly" in cuerpo:
                    esperas_hora += 1
                    if esperas_hora > 7:
                        raise Parar("la API sigue limitando tras una hora de espera. Vuelve a lanzarlo más tarde.")
                    print("    límite horario de la API: espero 10 minutos y sigo solo", flush=True)
                    time.sleep(600)
                    continue
                print("    límite por minuto de la API: espero 65 s", flush=True)
                time.sleep(65)
                continue
            if 400 <= e.code < 500:
                raise Parar(f"la API rechaza la petición (HTTP {e.code}): {cuerpo}")
            print(f"    error del servidor (HTTP {e.code}); reintento en 60 s", flush=True)
            time.sleep(60)
        except urllib.error.URLError as e:
            if isinstance(getattr(e, "reason", None), ssl.SSLCertVerificationError):
                raise Parar("error de certificados SSL. Instala certifi (py -m pip install certifi) "
                            "y vuelve a lanzar el script.")
            print(f"    sin conexión ({e.reason}); reintento en 60 s", flush=True)
            time.sleep(60)
        except (TimeoutError, ConnectionError, OSError) as e:
            print(f"    {e}; reintento en 60 s", flush=True)
            time.sleep(60)
    raise Parar("sin respuesta de la API tras varios intentos. Comprueba la conexión y vuelve a lanzarlo.")


def cache_ok(path):
    try:
        with gzip.open(path, "rt", encoding="utf-8") as f:
            j = json.load(f)
        return len(j["hourly"]["time"]) == N_HORAS
    except Exception:
        return False


def fmt(v):
    if v is None:
        return ""
    if isinstance(v, float):
        return f"{v:.2f}".rstrip("0").rstrip(".")
    return str(v)


def descargar():
    os.makedirs(CACHE, exist_ok=True)
    pendientes = [p for p in PUNTOS if not cache_ok(os.path.join(CACHE, p[0] + ".json.gz"))]
    hechos = len(PUNTOS) - len(pendientes)
    if hechos:
        print(f"Ya descargados: {hechos} de {len(PUNTOS)} puntos (carpeta era5_cache).", flush=True)
    for i, (pid, lat, lon) in enumerate(pendientes, 1):
        print(f"[{hechos + i}/{len(PUNTOS)}] {pid} ({lat}, {lon}) ...", flush=True)
        j = pedir(URL.format(lat=lat, lon=lon))
        n = len(j.get("hourly", {}).get("time", []))
        if n != N_HORAS:
            raise Parar(f"respuesta incompleta para {pid}: {n} horas en lugar de {N_HORAS}.")
        tmp = os.path.join(CACHE, pid + ".json.gz.tmp")
        with gzip.open(tmp, "wt", encoding="utf-8") as f:
            json.dump(j, f)
        os.replace(tmp, os.path.join(CACHE, pid + ".json.gz"))
        print(f"    {n} horas guardadas", flush=True)
        if i < len(pendientes):
            time.sleep(PAUSA)


def ensamblar():
    resumen = []
    with gzip.open(SALIDA, "wt", encoding="utf-8", newline="") as f:
        f.write(f"# ERA5 horario vía Open-Meteo (archive-api.open-meteo.com, models=era5), {INICIO} a {FIN}, UTC; "
                f"ensamblado el {datetime.now(timezone.utc).isoformat(timespec='seconds')}\n")
        w = csv.writer(f)
        w.writerow(["punto", "h", "u10", "u100", "dir100"])
        for pid, lat, lon in PUNTOS:
            with gzip.open(os.path.join(CACHE, pid + ".json.gz"), "rt", encoding="utf-8") as g:
                j = json.load(g)
            h = j["hourly"]
            faltan = sum(1 for a, b, c in zip(h["wind_speed_10m"], h["wind_speed_100m"], h["wind_direction_100m"])
                         if a is None or b is None or c is None)
            f.write(f"# {pid}: pedido {lat},{lon}; celda {j.get('latitude')},{j.get('longitude')}; "
                    f"elevación {j.get('elevation')}; horas {len(h['time'])}; incompletas {faltan}\n")
            for t, a, b, c in zip(h["time"], h["wind_speed_10m"], h["wind_speed_100m"], h["wind_direction_100m"]):
                w.writerow([pid, (int(t) - T0) // 3600, fmt(a), fmt(b), fmt(c)])
            resumen.append((pid, len(h["time"]), faltan))
    return resumen


def main():
    t0 = time.time()
    print("Descarga del viento horario ERA5 (2015-2024) para el artículo A1.\n", flush=True)
    descargar()
    resumen = ensamblar()
    print(f"\nListo en {(time.time() - t0) / 60:.1f} min. Fichero: {SALIDA} "
          f"({os.path.getsize(SALIDA) / 1e6:.1f} MB)")
    for pid, n, faltan in resumen:
        print(f"  {pid:10s} {n} horas, {faltan} incompletas")


if __name__ == "__main__":
    codigo = 0
    try:
        main()
    except Parar as e:
        print(f"\nNo se ha completado: {e}")
        codigo = 1
    except KeyboardInterrupt:
        print("\nInterrumpido. Al volver a lanzarlo, seguirá donde lo dejó.")
        codigo = 1
    if os.name == "nt":  # en Windows, que la ventana no se cierre sola si se abrió con doble clic
        try:
            input("\nPulsa Intro para cerrar.")
        except EOFError:
            pass
    sys.exit(codigo)
