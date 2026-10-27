#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Download of the hourly ERA5 wind at the 12 points of the route model of article A1
(ro-pax Galicia, Portsmouth-Bilbao): 10 sea points and the two ports.

Source: Open-Meteo, Historical Weather API with models=era5 (the same data that gave route/era5_route.json).
Period 2015-01-01 to 2024-12-31, UTC hours. Variables: speed at 10 m and 100 m, direction at 100 m.

Usage (Windows):  py download_era5_route.py      (or: python download_era5_route.py)
Usage (macOS/Linux): python3 download_era5_route.py
Requirements: Python 3.8 or later, no additional packages, and an internet connection.
Duration:     about 10 minutes. It waits 35 s between requests: each one asks for 10 years, which the free API
              counts as about 260 calls (limits: 600 per minute, 5,000 per hour, 10,000 per day).
Resumable:    each point is saved in the folder era5_cache; if interrupted, running it again continues where it
              stopped.
Output:       data/era5_hourly_route.csv.gz (about 8 MB).

Format of the output (gzip-compressed CSV; lines starting with '#' are metadata):
  point,h,u10,u100,dir100
  h = hours since 2015-01-01 00:00 UTC; u10 and u100 in m/s; dir100 = direction the wind comes from at 100 m
  (degrees); empty field = missing value.

Licences: Open-Meteo data under CC BY 4.0 (https://open-meteo.com/en/license);
ERA5 from the Copernicus Climate Change Service (Hersbach et al. 2020, doi:10.1002/qj.3803).
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

POINTS = [  # id, latitude, longitude (the same points as the route model)
    ("C1", 50.362, -1.850), ("C2", 50.088, -2.950), ("C3", 49.625, -4.062), ("C4", 48.975, -5.188),
    ("B1", 48.217, -5.527), ("B2", 47.350, -5.082), ("B3", 46.483, -4.638), ("B4", 45.617, -4.192),
    ("B5", 44.750, -3.748), ("B6", 43.883, -3.303),
    ("PORTSMOUTH", 50.790, -1.110), ("BILBAO", 43.360, -3.030),
]
START, END = "2015-01-01", "2024-12-31"
N_HOURS = 87672           # 3,653 days x 24 h
T0 = 1420070400           # 2015-01-01 00:00 UTC in Unix seconds
URL = ("https://archive-api.open-meteo.com/v1/archive?latitude={lat}&longitude={lon}"
       "&start_date=" + START + "&end_date=" + END +
       "&hourly=wind_speed_10m,wind_speed_100m,wind_direction_100m"
       "&models=era5&wind_speed_unit=ms&timezone=GMT&timeformat=unixtime")
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "era5_cache")
OUTPUT = os.path.join(HERE, "data", "era5_hourly_route.csv.gz")
PAUSE = 35  # s between requests

try:  # certifi certificates if installed (avoids SSL failures in some installations)
    import certifi
    CTX = ssl.create_default_context(cafile=certifi.where())
except ImportError:
    CTX = ssl.create_default_context()


class Stop(Exception):
    pass


def fetch(url):
    req = urllib.request.Request(url, headers={"User-Agent": "A1-Flettner-UPV-EHU-research"})
    hourly_waits = 0
    for attempt in range(12):
        try:
            with urllib.request.urlopen(req, timeout=300, context=CTX) as r:
                return json.loads(r.read().decode("utf-8"))
        except urllib.error.HTTPError as e:
            body = e.read().decode("utf-8", "replace")[:300]
            if e.code == 429:
                if "Daily" in body or "daily" in body:
                    raise Stop("the API has reached the daily limit of this connection. Run the script again "
                               "tomorrow: it will continue where it stopped.")
                if "Hourly" in body or "hourly" in body:
                    hourly_waits += 1
                    if hourly_waits > 7:
                        raise Stop("the API keeps limiting after an hour of waiting. Run it again later.")
                    print("    hourly limit of the API: waiting 10 minutes and continuing on my own", flush=True)
                    time.sleep(600)
                    continue
                print("    per-minute limit of the API: waiting 65 s", flush=True)
                time.sleep(65)
                continue
            if 400 <= e.code < 500:
                raise Stop(f"the API rejects the request (HTTP {e.code}): {body}")
            print(f"    server error (HTTP {e.code}); retrying in 60 s", flush=True)
            time.sleep(60)
        except urllib.error.URLError as e:
            if isinstance(getattr(e, "reason", None), ssl.SSLCertVerificationError):
                raise Stop("SSL certificate error. Install certifi (py -m pip install certifi) and run the "
                           "script again.")
            print(f"    no connection ({e.reason}); retrying in 60 s", flush=True)
            time.sleep(60)
        except (TimeoutError, ConnectionError, OSError) as e:
            print(f"    {e}; retrying in 60 s", flush=True)
            time.sleep(60)
    raise Stop("no answer from the API after several attempts. Check the connection and run it again.")


def cache_ok(path):
    try:
        with gzip.open(path, "rt", encoding="utf-8") as f:
            j = json.load(f)
        return len(j["hourly"]["time"]) == N_HOURS
    except Exception:
        return False


def fmt(v):
    if v is None:
        return ""
    if isinstance(v, float):
        return f"{v:.2f}".rstrip("0").rstrip(".")
    return str(v)


def download():
    os.makedirs(CACHE, exist_ok=True)
    pending = [p for p in POINTS if not cache_ok(os.path.join(CACHE, p[0] + ".json.gz"))]
    done = len(POINTS) - len(pending)
    if done:
        print(f"Already downloaded: {done} of {len(POINTS)} points (folder era5_cache).", flush=True)
    for i, (pid, lat, lon) in enumerate(pending, 1):
        print(f"[{done + i}/{len(POINTS)}] {pid} ({lat}, {lon}) ...", flush=True)
        j = fetch(URL.format(lat=lat, lon=lon))
        n = len(j.get("hourly", {}).get("time", []))
        if n != N_HOURS:
            raise Stop(f"incomplete answer for {pid}: {n} hours instead of {N_HOURS}.")
        tmp = os.path.join(CACHE, pid + ".json.gz.tmp")
        with gzip.open(tmp, "wt", encoding="utf-8") as f:
            json.dump(j, f)
        os.replace(tmp, os.path.join(CACHE, pid + ".json.gz"))
        print(f"    {n} hours saved", flush=True)
        if i < len(pending):
            time.sleep(PAUSE)


def assemble():
    summary = []
    os.makedirs(os.path.dirname(OUTPUT), exist_ok=True)
    with gzip.open(OUTPUT, "wt", encoding="utf-8", newline="") as f:
        f.write(f"# ERA5 hourly via Open-Meteo (archive-api.open-meteo.com, models=era5), {START} to {END}, UTC; "
                f"assembled {datetime.now(timezone.utc).isoformat(timespec='seconds')}\n")
        w = csv.writer(f)
        w.writerow(["point", "h", "u10", "u100", "dir100"])
        for pid, lat, lon in POINTS:
            with gzip.open(os.path.join(CACHE, pid + ".json.gz"), "rt", encoding="utf-8") as g:
                j = json.load(g)
            h = j["hourly"]
            missing = sum(1 for a, b, c in zip(h["wind_speed_10m"], h["wind_speed_100m"], h["wind_direction_100m"])
                          if a is None or b is None or c is None)
            f.write(f"# {pid}: requested {lat},{lon}; cell {j.get('latitude')},{j.get('longitude')}; "
                    f"elevation {j.get('elevation')}; hours {len(h['time'])}; incomplete {missing}\n")
            for t, a, b, c in zip(h["time"], h["wind_speed_10m"], h["wind_speed_100m"], h["wind_direction_100m"]):
                w.writerow([pid, (int(t) - T0) // 3600, fmt(a), fmt(b), fmt(c)])
            summary.append((pid, len(h["time"]), missing))
    return summary


def main():
    t0 = time.time()
    print("Download of the hourly ERA5 wind (2015-2024) for article A1.\n", flush=True)
    download()
    summary = assemble()
    print(f"\nDone in {(time.time() - t0) / 60:.1f} min. File: {OUTPUT} "
          f"({os.path.getsize(OUTPUT) / 1e6:.1f} MB)")
    for pid, n, missing in summary:
        print(f"  {pid:10s} {n} hours, {missing} incomplete")


if __name__ == "__main__":
    code = 0
    try:
        main()
    except Stop as e:
        print(f"\nNot completed: {e}")
        code = 1
    except KeyboardInterrupt:
        print("\nInterrupted. Running it again will continue where it stopped.")
        code = 1
    if os.name == "nt":  # on Windows, keep the window open if it was started with a double click
        try:
            input("\nPress Enter to close.")
        except EOFError:
            pass
    sys.exit(code)
