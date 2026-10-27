#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Descarga del viento horario ERA5 en la ruta del Viking Grace (Turku-Mariehamn-Estocolmo), para el caso de
validación del rotor convencional del artículo A1 (rotor Norsepower de 24 m x 4 m instalado en 2018).

Fuente: Open-Meteo, Historical Weather API con models=era5 (el mismo producto que la ruta Portsmouth-Bilbao).
Periodo 2015-01-01 a 2024-12-31, horas UTC. Variables: velocidad a 10 m y 100 m, dirección a 100 m.
Puntos: nueve a lo largo de la ruta (archipiélago de Turku, Skiftet, Åland, mar de Åland y archipiélago de
Estocolmo) y el puerto de Mariehamn; sus coordenadas son una aproximación de la ruta, no el track real.

Uso (Windows):  py descargar_era5_viking_grace.py      (o: python descargar_era5_viking_grace.py)
Uso (macOS/Linux): python3 descargar_era5_viking_grace.py
Requisitos:  Python 3.8 o posterior, sin paquetes adicionales, y conexión a internet.
Duración:    unos 10 minutos. Pide los datos año a año (100 peticiones pequeñas, con 5 s de espera entre
             ellas, dentro de los límites de la API gratuita) y empieza con una petición de prueba de un día.
Reanudable:  cada punto y año se guarda en la carpeta era5_cache_vg; si se corta, al volver a lanzarlo sigue
             donde lo dejó.
Si la API devuelve algo que no es JSON (por ejemplo, una página de bloqueo de la red), el script muestra el
principio de la respuesta y se detiene: en ese caso, prueba desde otra conexión (por ejemplo, fuera de la red
de la universidad).
Salida:      datos/era5_horario_viking_grace.csv.gz (unos 6 MB), con el mismo formato que
             datos/era5_horario_ruta.csv.gz:
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

PUNTOS = [  # id, latitud, longitud (aproximación de la ruta Turku-Mariehamn-Estocolmo)
    ("VG1", 60.380, 22.000), ("VG2", 60.240, 21.620), ("VG3", 60.150, 21.150), ("VG4", 60.050, 20.650),
    ("VG5", 60.050, 20.200), ("MARIEHAMN", 60.095, 19.935), ("VG6", 59.950, 19.500), ("VG7", 59.800, 19.150),
    ("VG8", 59.620, 18.850), ("VG9", 59.420, 18.400),
]
INICIO, FIN = "2015-01-01", "2024-12-31"
N_HORAS = 87672           # 3 653 días x 24 h
T0 = 1420070400           # 2015-01-01 00:00 UTC en segundos Unix
URL = ("https://archive-api.open-meteo.com/v1/archive?latitude={lat}&longitude={lon}"
       "&start_date={ini}&end_date={fin}"
       "&hourly=wind_speed_10m,wind_speed_100m,wind_direction_100m"
       "&models=era5&wind_speed_unit=ms&timezone=GMT&timeformat=unixtime")
ANIOS = list(range(2015, 2025))
AQUI = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(AQUI, "era5_cache_vg")
SALIDA = os.path.join(AQUI, "datos", "era5_horario_viking_grace.csv.gz")
PAUSA = 5  # s entre peticiones (cada una pide un año: unas 26 llamadas de la API gratuita)

try:  # certificados de certifi si está instalado (evita fallos SSL en algunas instalaciones)
    import certifi
    CTX = ssl.create_default_context(cafile=certifi.where())
except ImportError:
    CTX = ssl.create_default_context()


class Parar(Exception):
    pass


def pedir(url):
    req = urllib.request.Request(url, headers={"User-Agent": "A1-Flettner-UPV-EHU-research",
                                               "Accept": "application/json", "Accept-Encoding": "gzip"})
    esperas_hora, malas = 0, 0
    for intento in range(12):
        try:
            with urllib.request.urlopen(req, timeout=300, context=CTX) as r:
                raw, estado = r.read(), r.status
                tipo = r.headers.get("Content-Type", "")
                cod = (r.headers.get("Content-Encoding") or "").lower()
            if cod == "gzip" or raw[:2] == b"\x1f\x8b":
                raw = gzip.decompress(raw)
            try:
                return json.loads(raw.decode("utf-8"))
            except (UnicodeDecodeError, ValueError):
                muestra = raw[:400].decode("utf-8", "replace").replace("\n", " ")
                print(f"    la respuesta no es JSON (HTTP {estado}, tipo '{tipo}', {len(raw)} bytes): {muestra!r}",
                      flush=True)
                malas += 1
                if malas >= 3:
                    raise Parar("la API (o la red) devuelve una respuesta que no es JSON; el principio de la "
                                "respuesta está arriba. Si es una página de la red de la universidad, prueba desde "
                                "otra conexión. Al volver a lanzarlo, seguirá donde lo dejó.")
                print("    reintento en 60 s", flush=True)
                time.sleep(60)
                continue
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


def horas_anio(y):
    return (366 if y % 4 == 0 else 365) * 24


def cache_ok(path, horas):
    try:
        with gzip.open(path, "rt", encoding="utf-8") as f:
            j = json.load(f)
        return len(j["hourly"]["time"]) == horas
    except Exception:
        return False


def ruta_cache(pid, y):
    return os.path.join(CACHE, f"{pid}_{y}.json.gz")


def fmt(v):
    if v is None:
        return ""
    if isinstance(v, float):
        return f"{v:.2f}".rstrip("0").rstrip(".")
    return str(v)


def prueba():
    """Petición de un día en el primer punto: comprueba que la API responde con JSON antes de empezar."""
    pid, lat, lon = PUNTOS[0]
    print("Prueba de conexión (un día en el primer punto)...", flush=True)
    j = pedir(URL.format(lat=lat, lon=lon, ini=INICIO, fin=INICIO))
    n = len(j.get("hourly", {}).get("time", []))
    if n != 24:
        raise Parar(f"la prueba devuelve {n} horas en lugar de 24: {str(j)[:300]}")
    print("    la API responde correctamente\n", flush=True)


def descargar():
    os.makedirs(CACHE, exist_ok=True)
    tareas = [(p, y) for p in PUNTOS for y in ANIOS]
    pendientes = [(p, y) for p, y in tareas if not cache_ok(ruta_cache(p[0], y), horas_anio(y))]
    hechas = len(tareas) - len(pendientes)
    if hechas:
        print(f"Ya descargados: {hechas} de {len(tareas)} años-punto (carpeta era5_cache_vg).", flush=True)
    if pendientes:
        prueba()
    for i, ((pid, lat, lon), y) in enumerate(pendientes, 1):
        print(f"[{hechas + i}/{len(tareas)}] {pid} ({lat}, {lon}) {y} ...", flush=True)
        j = pedir(URL.format(lat=lat, lon=lon, ini=f"{y}-01-01", fin=f"{y}-12-31"))
        n = len(j.get("hourly", {}).get("time", []))
        if n != horas_anio(y):
            raise Parar(f"respuesta incompleta para {pid} en {y}: {n} horas en lugar de {horas_anio(y)}.")
        tmp = ruta_cache(pid, y) + ".tmp"
        with gzip.open(tmp, "wt", encoding="utf-8") as f:
            json.dump(j, f)
        os.replace(tmp, ruta_cache(pid, y))
        if i < len(pendientes):
            time.sleep(PAUSA)


def ensamblar():
    resumen = []
    os.makedirs(os.path.dirname(SALIDA), exist_ok=True)
    with gzip.open(SALIDA, "wt", encoding="utf-8", newline="") as f:
        f.write(f"# ERA5 horario vía Open-Meteo (archive-api.open-meteo.com, models=era5), {INICIO} a {FIN}, UTC; "
                f"pedido año a año; ensamblado el {datetime.now(timezone.utc).isoformat(timespec='seconds')}\n")
        w = csv.writer(f)
        w.writerow(["punto", "h", "u10", "u100", "dir100"])
        for pid, lat, lon in PUNTOS:
            n, faltan, celda, filas = 0, 0, None, []
            for y in ANIOS:
                with gzip.open(ruta_cache(pid, y), "rt", encoding="utf-8") as g:
                    j = json.load(g)
                h = j["hourly"]
                celda = celda or (j.get("latitude"), j.get("longitude"), j.get("elevation"))
                faltan += sum(1 for a, b, c in zip(h["wind_speed_10m"], h["wind_speed_100m"], h["wind_direction_100m"])
                              if a is None or b is None or c is None)
                filas += [[pid, (int(t) - T0) // 3600, fmt(a), fmt(b), fmt(c)]
                          for t, a, b, c in zip(h["time"], h["wind_speed_10m"], h["wind_speed_100m"],
                                                h["wind_direction_100m"])]
                n += len(h["time"])
            if n != N_HORAS:
                raise Parar(f"{pid}: {n} horas en total en lugar de {N_HORAS}.")
            f.write(f"# {pid}: pedido {lat},{lon}; celda {celda[0]},{celda[1]}; elevación {celda[2]}; "
                    f"horas {n}; incompletas {faltan}\n")
            w.writerows(filas)
            resumen.append((pid, n, faltan))
    return resumen


def main():
    t0 = time.time()
    print("Descarga del viento horario ERA5 (2015-2024) en la ruta del Viking Grace, artículo A1.\n", flush=True)
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
