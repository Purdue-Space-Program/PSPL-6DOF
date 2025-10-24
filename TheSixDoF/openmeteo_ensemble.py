#!/usr/bin/env python3
"""
Open-Meteo Ensemble API fetcher — single date *or* date range.

Usage examples
--------------
# Single calendar date (local time if timezone=auto):
python openmeteo_ensemble.py --lat 28.5 --lon -80.6 --date 2025-11-18 \
  --hourly temperature_2m,wind_speed_10m,wind_gusts_10m,wind_speed_100m,wind_direction_100m,cloud_cover,precipitation,cape

# Date range (inclusive ISO dates, up to 35 days ahead of today):
python openmeteo_ensemble.py --lat 28.5 --lon -80.6 \
  --start-date 2025-11-10 --end-date 2025-11-20 \
  --hourly wind_speed_10m,wind_gusts_10m,wind_speed_100m,wind_direction_100m,pressure_msl,cloud_cover,precipitation,cape

Notes
-----
• This endpoint supports absolute time windows via &start_date=&end_date= (ISO YYYY-MM-DD) and also the relative
  &forecast_days= / &past_days= alternatives. We use absolute dates for clarity. The forecast horizon is
  available **up to 35 days ahead** (model-dependent). Beyond that, the API won’t return future data. 
• If you request a date/window partly beyond a model’s horizon, that model will return no values there.
• Timezone: 'auto' returns hours starting at 00:00 *local* time for the coordinates.
• Ensemble model selection: kept configurable via --models (default 'gfs_seamless'); your MATLAB caller locks this down.

Docs: https://open-meteo.com/en/docs/ensemble-api
"""

import argparse
import datetime as dt
import json
import sys
from typing import List, Optional
from urllib.parse import urlencode
from urllib.request import urlopen, Request
from urllib.error import HTTPError

API_URL = "https://ensemble-api.open-meteo.com/v1/ensemble"

# Rocket-friendly defaults (you can override with --hourly)
DEFAULT_HOURLY = [
    "wind_speed_10m","wind_gusts_10m","wind_direction_10m",
    "wind_speed_100m","wind_direction_100m",
    "temperature_2m","relative_humidity_2m","pressure_msl",
    "cloud_cover","visibility","precipitation","rain",
    "cape","freezing_level_height"
]

def _hourly_to_list(hourly_str: Optional[str]) -> List[str]:
    if not hourly_str:
        return DEFAULT_HOURLY[:]
    return [p.strip() for p in hourly_str.replace(";", ",").split(",") if p.strip()]

def _parse_date(s: str) -> dt.date:
    try:
        return dt.date.fromisoformat(s)
    except Exception as e:
        raise SystemExit(f"Invalid date '{s}' (use YYYY-MM-DD). Error: {e}")

def _validate_window(start: dt.date, end: dt.date) -> None:
    if end < start:
        raise SystemExit("end-date must be >= start-date.")
    today = dt.date.today()
    if start < today:
        raise SystemExit("Requested start-date is in the past. Use Past/Previous APIs for historical data.")
    if (end - today).days > 35:
        raise SystemExit("Window goes beyond 35 days ahead (API limit). Choose earlier dates.")

def _build_params(lat, lon, models, hourly_list, timezone, start_date, end_date) -> List[tuple]:
    # Return a list of (key, value) tuples so we can repeat &hourly= many times
    params = [
        ("latitude", lat),
        ("longitude", lon),
        ("models", models),
        ("timezone", timezone),
        ("timeformat", "iso8601"),
        ("start_date", start_date),
        ("end_date", end_date),
    ]
    for v in hourly_list:
        params.append(("hourly", v))  # repeated params to avoid CSV parsing BS
    return params

def fetch_json(params_kv: List[tuple], timeout: int = 60) -> dict:
    url = f"{API_URL}?{urlencode(params_kv, doseq=True)}"
    req = Request(url, headers={"User-Agent": "openmeteo-ensemble-client/1.2"})
    try:
        with urlopen(req, timeout=timeout) as r:
            body = r.read().decode("utf-8")
        return json.loads(body)
    except HTTPError as e:
        try:
            err_body = e.read().decode("utf-8")
        except Exception:
            err_body = ""
        msg = f"HTTP {e.code}"
        if err_body:
            try:
                msg += f": {json.dumps(json.loads(err_body), ensure_ascii=False)}"
            except Exception:
                msg += f": {err_body[:400]}"
        raise SystemExit(msg)
    except Exception as e:
        raise SystemExit(f"Request failed: {e}")

def main(argv: Optional[List[str]] = None) -> None:
    ap = argparse.ArgumentParser(description="Open-Meteo Ensemble API (single date or date range)")
    ap.add_argument("--lat", "--latitude", dest="lat", type=float, required=True)
    ap.add_argument("--lon", "--longitude", dest="lon", type=float, required=True)
    grp = ap.add_mutually_exclusive_group(required=True)
    grp.add_argument("--date", type=str, help="Single date (YYYY-MM-DD).")
    grp.add_argument("--start-date", dest="start_date", type=str, help="Start date (YYYY-MM-DD) — use with --end-date")
    ap.add_argument("--end-date", dest="end_date", type=str, help="End date (YYYY-MM-DD) — required with --start-date")
    ap.add_argument("--hourly", type=str, default=",".join(DEFAULT_HOURLY),
                    help="Comma- or semicolon-separated hourly variables (see docs).")
    ap.add_argument("--models", type=str, default="gfs_seamless",
                    help="Comma-separated model codes (default: gfs_seamless).")
    ap.add_argument("--timezone", type=str, default="auto", help="e.g., 'auto' or 'America/New_York'.")
    args = ap.parse_args(argv)

    # Resolve dates
    if args.date:
        d = _parse_date(args.date)
        start, end = d, d
    else:
        if not args.start_date or not args.end_date:
            raise SystemExit("Both --start-date and --end-date are required for a range.")
        start = _parse_date(args.start_date)
        end = _parse_date(args.end_date)

    _validate_window(start, end)
    hourly_list = _hourly_to_list(args.hourly)
    params_kv = _build_params(
        lat=args.lat, lon=args.lon, models=args.models,
        hourly_list=hourly_list, timezone=args.timezone,
        start_date=start.isoformat(), end_date=end.isoformat()
    )
    data = fetch_json(params_kv)
    json.dump(data, sys.stdout, ensure_ascii=False)
    sys.stdout.write("\n")

if __name__ == "__main__":
    main()