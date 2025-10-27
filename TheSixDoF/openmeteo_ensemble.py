#!/usr/bin/env python3
"""
Open-Meteo Ensemble API fetcher — single date *or* date range.

- Past-only windows -> Historical Weather API (/v1/archive)
- Today/Future-only windows -> Ensemble API (/v1/ensemble)
- Mixed windows (return is kinda weird, try to avoid) -> split return {"archive": {...}, "ensemble": {...}}

USAGE EXAMPLE (what MATLAB will run for you):
python fetch_openmeteo_ensemble.py \
    --lat 28.5 --lon -80.6 \
    --start-date 2025-11-10 --end-date 2025-11-12 \
    --models gfs_seamless \
    --timezone auto \
    --hourly wind_speed_10m wind_speed_100m wind_gusts_10m wind_direction_10m surface_pressure cloud_cover precipitation \
    --wind-speed-unit ms --precipitation-unit mm --temperature-unit celsius
NOTES:
- "ensemble-api" can forecast out to ~35 days ahead with many ensemble members
  (e.g. GFS ensemble ~31 members). Each member is a plausible future. We'll grab
  all members. The hourly arrays will look like [N_hours x N_members].
- "archive-api" is past weather from reanalysis and usually returns a single
  deterministic value for each hour (so effectively N_members = 1).

Docs: https://open-meteo.com/en/docs/ensemble-api
"""

import argparse
import datetime as dt
import json
import sys
import requests


def call_api(url, params):
    """Perform GET, raise if HTTP error, return parsed JSON."""
    r = requests.get(url, params=params, timeout=30)
    r.raise_for_status()
    return r.json()


def main():
    parser = argparse.ArgumentParser(
        description="Fetch Open-Meteo ensemble/historical data for MATLAB."
    )
    parser.add_argument("--lat", type=float, required=True, help="Latitude (deg)")
    parser.add_argument("--lon", type=float, required=True, help="Longitude (deg)")

    parser.add_argument(
        "--start-date",
        dest="start_date",
        required=True,
        help="YYYY-MM-DD (local calendar date start)"
    )
    parser.add_argument(
        "--end-date",
        dest="end_date",
        required=False,
        help="YYYY-MM-DD (local calendar date end, inclusive). " +
             "If omitted we use start_date."
    )

    # Ensemble model name to use for future forecasts.
    # We'll default to 'gfs_seamless' which is a long-range GFS-style ensemble.
    parser.add_argument(
        "--models",
        type=str,
        default="gfs_seamless",
        help="Comma-separated ensemble model list for future forecast portion."
    )

    # Hourly variables. We accept multiple tokens after --hourly.
    parser.add_argument(
        "--hourly",
        nargs="+",
        default=[
            "wind_speed_10m",
            "wind_speed_100m",
            "wind_gusts_10m",
            "wind_direction_10m",
            "surface_pressure",
            "cloud_cover",
            "precipitation",
        ],
        help="One or more Open-Meteo hourly variable names (see docs)."
    )

    parser.add_argument(
        "--timezone",
        type=str,
        default="auto",
        help="Open-Meteo timezone param. 'auto' = infer local tz from lat/lon."
    )
    parser.add_argument(
        "--wind-speed-unit",
        dest="wind_speed_unit",
        type=str,
        default="ms",
        help="ms, kmh, mph, kn"
    )
    parser.add_argument(
        "--precipitation-unit",
        dest="precipitation_unit",
        type=str,
        default="mm",
        help="mm or inch"
    )
    parser.add_argument(
        "--temperature-unit",
        dest="temperature_unit",
        type=str,
        default="celsius",
        help="celsius or fahrenheit"
    )

    args = parser.parse_args()

    if args.end_date is None:
        args.end_date = args.start_date

    # Parse and sanity-check dates
    try:
        start_dt = dt.date.fromisoformat(args.start_date)
        end_dt = dt.date.fromisoformat(args.end_date)
    except ValueError as e:
        print(json.dumps({"error": True, "reason": f"Invalid start/end date: {e}"}))
        sys.exit(1)

    # Swap if user flipped them
    if end_dt < start_dt:
        start_dt, end_dt = end_dt, start_dt

    # "today" in UTC. We treat anything strictly before today as "past"
    # (historical reanalysis) and anything from today forward as "future"
    # (ensemble forecast).
    today_utc = dt.datetime.utcnow().date()

    # We'll build a result dict. Could contain:
    #   {"archive": {...}, "ensemble": {...}}
    # or just one of them, or an error.
    results = {}

    # -------- Past chunk (strictly before today) --------
    if start_dt < today_utc:
        past_end_dt = min(end_dt, today_utc - dt.timedelta(days=1))

        archive_params = {
            "latitude": args.lat,
            "longitude": args.lon,
            "start_date": start_dt.isoformat(),
            "end_date": past_end_dt.isoformat(),
            # IMPORTANT: For Open-Meteo, if you pass a list for hourly,
            # requests will encode ?hourly=var1&hourly=var2&...
            "hourly": args.hourly,
            "timezone": args.timezone,
            "timeformat": "iso8601",
            "temperature_unit": args.temperature_unit,
            "wind_speed_unit": args.wind_speed_unit,
            "precipitation_unit": args.precipitation_unit,
        }

        try:
            results["archive"] = call_api(
                "https://archive-api.open-meteo.com/v1/archive",
                archive_params
            )
        except requests.HTTPError as e:
            print(json.dumps({
                "error": True,
                "source": "archive",
                "status": e.response.status_code,
                "reason": e.response.text
            }))
            sys.exit(1)

    # -------- Future/forecast chunk (today and after) --------
    if end_dt >= today_utc:
        future_start_dt = max(start_dt, today_utc)

        ensemble_params = {
            "latitude": args.lat,
            "longitude": args.lon,
            "start_date": future_start_dt.isoformat(),
            "end_date": end_dt.isoformat(),
            "models": args.models,
            "hourly": args.hourly,
            "timezone": args.timezone,
            "timeformat": "iso8601",
            "temperature_unit": args.temperature_unit,
            "wind_speed_unit": args.wind_speed_unit,
            "precipitation_unit": args.precipitation_unit,
        }

        try:
            results["ensemble"] = call_api(
                "https://ensemble-api.open-meteo.com/v1/ensemble",
                ensemble_params
            )
        except requests.HTTPError as e:
            print(json.dumps({
                "error": True,
                "source": "ensemble",
                "status": e.response.status_code,
                "reason": e.response.text
            }))
            sys.exit(1)

    # If we only got one side, unwrap it so MATLAB gets a clean struct.
    if "archive" in results and "ensemble" in results:
        final_payload = results
    elif "archive" in results:
        final_payload = results["archive"]
    elif "ensemble" in results:
        final_payload = results["ensemble"]
    else:
        final_payload = {
            "error": True,
            "reason": "No data returned for given date range."
        }

    # Print one-line JSON for MATLAB's jsondecode
    print(json.dumps(final_payload))


if __name__ == "__main__":
    main()