#!/usr/bin/env python3
"""
Open-Meteo Ensemble API fetcher — single date *or* date range.

- Past-only windows -> Historical Weather API (/v1/archive)
- Today/Future-only windows -> Ensemble API (/v1/ensemble)
- Mixed windows (return is kinda weird, try to avoid) -> split return {"archive": {...}, "ensemble": {...}}

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
• This supports absolute time windows via &start_date=&end_date= (ISO YYYY-MM-DD) and also the relative
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
from typing import List, Optional, Tuple, Dict, Any
from urllib.parse import urlencode
from urllib.request import urlopen, Request
from urllib.error import HTTPError

ENSEMBLE_URL = "https://ensemble-api.open-meteo.com/v1/ensemble"
ARCHIVE_URL  = "https://archive-api.open-meteo.com/v1/archive"

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

def _validate_future_window(start: dt.date, end: dt.date, today: dt.date) -> None:
    if start < today:
        raise SystemExit("Internal error: future window starts before today.")
    if end < start:
        raise SystemExit("end-date must be >= start-date.")
    if (end - today).days > 35:
        raise SystemExit("Future window extends beyond 35 days ahead (Ensemble API limit).")

def _build_params_common(lat: float, lon: float, timezone: str) -> List[Tuple[str, Any]]:
    return [("latitude", lat), ("longitude", lon), ("timezone", timezone), ("timeformat", "iso8601")]

def _build_params_hourly(params: List[Tuple[str, Any]], hourly_list: List[str]) -> None:
    # Use repeated &hourly=... to avoid CSV enum parsing issues
    for v in hourly_list:
        params.append(("hourly", v))

def _fetch(url: str, params_kv: List[Tuple[str, Any]], timeout: int = 60) -> Dict[str, Any]:
    req_url = f"{url}?{urlencode(params_kv, doseq=True)}"
    req = Request(req_url, headers={"User-Agent": "openmeteo-router/1.0"})
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

def _fetch_ensemble(lat, lon, hourly_list, timezone, start, end, models) -> Dict[str, Any]:
    params = _build_params_common(lat, lon, timezone)
    params += [("models", models), ("start_date", start.isoformat()), ("end_date", end.isoformat())]
    _build_params_hourly(params, hourly_list)
    return _fetch(ENSEMBLE_URL, params)

def _fetch_archive(lat, lon, hourly_list, timezone, start, end) -> Dict[str, Any]:
    params = _build_params_common(lat, lon, timezone)
    params += [("start_date", start.isoformat()), ("end_date", end.isoformat())]
    _build_params_hourly(params, hourly_list)
    return _fetch(ARCHIVE_URL, params)

def main(argv: Optional[List[str]] = None) -> None:
    ap = argparse.ArgumentParser(description="Open-Meteo router: Ensemble (future) + Archive (past)")
    ap.add_argument("--lat", "--latitude", dest="lat", type=float, required=True)
    ap.add_argument("--lon", "--longitude", dest="lon", type=float, required=True)

    grp = ap.add_mutually_exclusive_group(required=True)
    grp.add_argument("--date", type=str, help="Single date (YYYY-MM-DD).")
    grp.add_argument("--start-date", dest="start_date", type=str, help="Start date (YYYY-MM-DD) — use with --end-date")
    ap.add_argument("--end-date", dest="end_date", type=str, help="End date (YYYY-MM-DD) — required with --start-date")

    ap.add_argument("--hourly", type=str, default=",".join(DEFAULT_HOURLY),
                    help="Comma-/semicolon-separated hourly variables (exact API names).")
    ap.add_argument("--models", type=str, default="gfs_seamless",
                    help="Ensemble models for FUTURE part (ignored for archive).")
    ap.add_argument("--timezone", type=str, default="auto", help="e.g., 'auto' or 'America/New_York'.")
    args = ap.parse_args(argv)

    # Resolve dates
    if args.date:
        start = end = _parse_date(args.date)
    else:
        if not args.start_date or not args.end_date:
            raise SystemExit("Both --start-date and --end-date are required for a range.")
        start, end = _parse_date(args.start_date), _parse_date(args.end_date)
        if end < start:
            raise SystemExit("end-date must be >= start-date.")

    today = dt.date.today()
    hourly_list = _hourly_to_list(args.hourly)

    # Case 1: all past
    if end < today:
        data = _fetch_archive(args.lat, args.lon, hourly_list, args.timezone, start, end)
        json.dump(data, sys.stdout, ensure_ascii=False); sys.stdout.write("\n")
        return

    # Case 2: all future/today
    if start >= today:
        _validate_future_window(start, end, today)
        data = _fetch_ensemble(args.lat, args.lon, hourly_list, args.timezone, start, end, args.models)
        json.dump(data, sys.stdout, ensure_ascii=False); sys.stdout.write("\n")
        return

    # Case 3: mixed window (split)
    past_end = min(end, today - dt.timedelta(days=1))
    future_start = max(start, today)
    _validate_future_window(future_start, end, today)

    archive_part  = _fetch_archive(args.lat, args.lon, hourly_list, args.timezone, start, past_end)
    ensemble_part = _fetch_ensemble(args.lat, args.lon, hourly_list, args.timezone, future_start, end, args.models)

    out = {
        "note": "Window crosses past/future boundary; results are split by API.",
        "archive": archive_part,
        "ensemble": ensemble_part
    }
    json.dump(out, sys.stdout, ensure_ascii=False); sys.stdout.write("\n")

if __name__ == "__main__":
    main()