#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Download of the hourly ERA5 wind on the route of the Viking Grace (Turku-Mariehamn-Stockholm), for the validation
case of the conventional rotor of article A1 (Norsepower rotor of 24 m x 4 m installed in 2018).

Source: Open-Meteo, Historical Weather API with models=era5 (the same product as the Portsmouth-Bilbao route).
Period 2015-01-01 to 2024-12-31, UTC hours. Variables: speed at 10 m and 100 m, direction at 100 m.
Points: nine along the route (Turku archipelago, Skiftet, Åland, Sea of Åland and Stockholm archipelago) and the
port of Mariehamn; their coordinates approximate the route, not the actual track.

Usage (Windows):  py download_era5_viking_grace.py      (or: python download_era5_viking_grace.py)
Usage (macOS/Linux): python3 download_era5_viking_grace.py
Requirements: Python 3.8 or later, no additional packages, and an internet connection.
Duration:     about 10 minutes. It asks for the data year by year (100 small requests, with 5 s of waiting between
              them, within the limits of the free API) and starts with a one-day test request.
Resumable:    each point and year is saved in the folder era5_cache_vg; if interrupted, running it again continues
              where it stopped.
If the API returns something that is not JSON (for instance, a blocking page of the network), the script shows the
start of the answer and stops: in that case, try from another connection (for instance, outside the university
network).
Output:       data/era5_hourly_viking_grace.csv.gz (about 7 MB), with the same format as
              data/era5_hourly_route.csv.gz:
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

POINTS = [  # id, latitude, longitude (approximation of the Turku-Mariehamn-Stockholm route)
    ("VG1", 60.380, 22.000), ("VG2", 60.240, 21.620), ("VG3", 60.150, 21.150), ("VG4", 60.050, 20.650),
    ("VG5", 60.050, 20.200), ("MARIEHAMN", 60.095, 19.935), ("VG6", 59.950, 19.500), ("VG7", 59.800, 19.150),
    ("VG8", 59.620, 18.850), ("VG9", 59.420, 18.400),
]
START, END = "2015-01-01", "2024-12-31"
N_HOURS = 87672           # 3,653 days x 24 h
T0 = 1420070400           # 2015-01-01 00:00 UTC in Unix seconds
URL = ("https://archive-api.open-meteo.com/v1/archive?latitude={lat}&longitude={lon}"
       "&start_date={ini}&end_date={fin}"
       "&hourly=wind_speed_10m,wind_speed_100m,wind_direction_100m"
       "&models=era5&wind_speed_unit=ms&timezone=GMT&timeformat=unixtime")
YEARS = list(range(2015, 2025))
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "era5_cache_vg")
OUTPUT = os.path.join(HERE, "data", "era5_hourly_viking_grace.csv.gz")
PAUSE = 5  # s between requests (each one asks for one year: about 26 calls of the free API)

try:  # certifi certificates if installed (avoids SSL failures in some installations)
    import certifi
    CTX = ssl.create_default_context(cafile=certifi.where())
except ImportError:
    CTX = ssl.create_default_context()


class Stop(Exception):
    pass


def fetch(url):
    req = urllib.request.Request(url, headers={"User-Agent": "A1-Flettner-UPV-EHU-research",
                                               "Accept": "application/json", "Accept-Encoding": "gzip"})
    hourly_waits, bad = 0, 0
    for attempt in range(12):
        try:
            with urllib.request.urlopen(req, timeout=300, context=CTX) as r:
                raw, status = r.read(), r.status
                ctype = r.headers.get("Content-Type", "")
                enc = (r.headers.get("Content-Encoding") or "").lower()
            if enc == "gzip" or raw[:2] == b"\x1f\x8b":
                raw = gzip.decompress(raw)
            try:
                return json.loads(raw.decode("utf-8"))
            except (UnicodeDecodeError, ValueError):
                sample = raw[:400].decode("utf-8", "replace").replace("\n", " ")
                print(f"    the answer is not JSON (HTTP {status}, type '{ctype}', {len(raw)} bytes): {sample!r}",
                      flush=True)
                bad += 1
                if bad >= 3:
                    raise Stop("the API (or the network) returns an answer that is not JSON; the start of the "
                               "answer is above. If it is a page of the university network, try from another "
                               "connection. Running it again will continue where it stopped.")
                print("    retrying in 60 s", flush=True)
                time.sleep(60)
                continue
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


def year_hours(y):
    return (366 if y % 4 == 0 else 365) * 24


def cache_ok(path, hours):
    try:
        with gzip.open(path, "rt", encoding="utf-8") as f:
            j = json.load(f)
        return len(j["hourly"]["time"]) == hours
    except Exception:
        return False


def cache_path(pid, y):
    return os.path.join(CACHE, f"{pid}_{y}.json.gz")


def fmt(v):
    if v is None:
        return ""
    if isinstance(v, float):
        return f"{v:.2f}".rstrip("0").rstrip(".")
    return str(v)


def test():
    """One-day request at the first point: checks that the API answers with JSON before starting."""
    pid, lat, lon = POINTS[0]
    print("Connection test (one day at the first point)...", flush=True)
    j = fetch(URL.format(lat=lat, lon=lon, ini=START, fin=START))
    n = len(j.get("hourly", {}).get("time", []))
    if n != 24:
        raise Stop(f"the test returns {n} hours instead of 24: {str(j)[:300]}")
    print("    the API answers correctly\n", flush=True)


def download():
    os.makedirs(CACHE, exist_ok=True)
    tasks = [(p, y) for p in POINTS for y in YEARS]
    pending = [(p, y) for p, y in tasks if not cache_ok(cache_path(p[0], y), year_hours(y))]
    done = len(tasks) - len(pending)
    if done:
        print(f"Already downloaded: {done} of {len(tasks)} point-years (folder era5_cache_vg).", flush=True)
    if pending:
        test()
    for i, ((pid, lat, lon), y) in enumerate(pending, 1):
        print(f"[{done + i}/{len(tasks)}] {pid} ({lat}, {lon}) {y} ...", flush=True)
        j = fetch(URL.format(lat=lat, lon=lon, ini=f"{y}-01-01", fin=f"{y}-12-31"))
        n = len(j.get("hourly", {}).get("time", []))
        if n != year_hours(y):
            raise Stop(f"incomplete answer for {pid} in {y}: {n} hours instead of {year_hours(y)}.")
        tmp = cache_path(pid, y) + ".tmp"
        with gzip.open(tmp, "wt", encoding="utf-8") as f:
            json.dump(j, f)
        os.replace(tmp, cache_path(pid, y))
        if i < len(pending):
            time.sleep(PAUSE)


def assemble():
    summary = []
    os.makedirs(os.path.dirname(OUTPUT), exist_ok=True)
    with gzip.open(OUTPUT, "wt", encoding="utf-8", newline="") as f:
        f.write(f"# ERA5 hourly via Open-Meteo (archive-api.open-meteo.com, models=era5), {START} to {END}, UTC; "
                f"requested year by year; assembled {datetime.now(timezone.utc).isoformat(timespec='seconds')}\n")
        w = csv.writer(f)
        w.writerow(["point", "h", "u10", "u100", "dir100"])
        for pid, lat, lon in POINTS:
            n, missing, cell, rows = 0, 0, None, []
            for y in YEARS:
                with gzip.open(cache_path(pid, y), "rt", encoding="utf-8") as g:
                    j = json.load(g)
                h = j["hourly"]
                cell = cell or (j.get("latitude"), j.get("longitude"), j.get("elevation"))
                missing += sum(1 for a, b, c in zip(h["wind_speed_10m"], h["wind_speed_100m"], h["wind_direction_100m"])
                               if a is None or b is None or c is None)
                rows += [[pid, (int(t) - T0) // 3600, fmt(a), fmt(b), fmt(c)]
                         for t, a, b, c in zip(h["time"], h["wind_speed_10m"], h["wind_speed_100m"],
                                               h["wind_direction_100m"])]
                n += len(h["time"])
            if n != N_HOURS:
                raise Stop(f"{pid}: {n} hours in total instead of {N_HOURS}.")
            f.write(f"# {pid}: requested {lat},{lon}; cell {cell[0]},{cell[1]}; elevation {cell[2]}; "
                    f"hours {n}; incomplete {missing}\n")
            w.writerows(rows)
            summary.append((pid, n, missing))
    return summary


def main():
    t0 = time.time()
    print("Download of the hourly ERA5 wind (2015-2024) on the route of the Viking Grace, article A1.\n", flush=True)
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
