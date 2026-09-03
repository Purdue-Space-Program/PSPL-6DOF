import requests
import requests_cache
import numpy as np
import pandas as pd
import time
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

def _make_session():
    session = requests_cache.CachedSession('.cache', expire_after=3600)
    retry_policy = Retry(total=10, backoff_factor=0.5,
                         status_forcelist=[429, 500, 502, 503, 504])
    session.mount('https://', HTTPAdapter(max_retries=retry_policy))
    return session

def _get(url, params, max_attempts=5):
    session = _make_session()
    for attempt in range(max_attempts):
        try:
            resp = session.get(url, params=params, timeout=30)
            resp.raise_for_status()
            return resp.json()
        except Exception as e:
            if attempt < max_attempts - 1:
                time.sleep(2 ** attempt)
            else:
                raise

def getCurrentWeather(lat=39.7392, lon=-104.9847):
    url = "https://api.open-meteo.com/v1/forecast"
    params = {
        "latitude": lat,
        "longitude": lon,
        "models": "best_match",
        "current": [
            "temperature_2m",
            "surface_pressure",
            "relative_humidity_2m",
            "wind_speed_10m",
            "wind_direction_10m",
            "wind_gusts_10m",
        ],
        "timezone": "auto",
        "forecast_days": 1,
        "timeformat": "unixtime",
        "wind_speed_unit": "ms",
    }

    data = _get(url, params)
    current = data.get("current", {})

    return {
        "latitude":             data.get("latitude"),
        "longitude":            data.get("longitude"),
        "elevation":            data.get("elevation"),
        "timezone":             data.get("timezone"),
        "timezone_abbrev":      data.get("timezone_abbreviation"),
        "utc_offset_seconds":   float(data.get("utc_offset_seconds", 0)),
        "time":                 current.get("time"),
        "temperature_2m":       current.get("temperature_2m"),
        "surface_pressure":     current.get("surface_pressure"),
        "relative_humidity_2m": current.get("relative_humidity_2m"),
        "wind_speed_10m":       current.get("wind_speed_10m"),
        "wind_direction_10m":   current.get("wind_direction_10m"),
        "wind_gusts_10m":       current.get("wind_gusts_10m"),
    }


def getHourlyWeather(lat, lon):
    url = "https://api.open-meteo.com/v1/forecast"
    pressure_levels = [1000,975,950,925,900,850,800,700,600,500,400,300,250,200,150,100,70,50,30]
    hourly_vars = []
    for p in pressure_levels:
        hourly_vars += [
            f"temperature_{p}hPa",
            f"wind_speed_{p}hPa",
            f"wind_direction_{p}hPa",
            f"geopotential_height_{p}hPa",
        ]
    hourly_vars += ["temperature_2m", "surface_pressure", "wind_gusts_10m"]

    params = {
        "latitude": lat,
        "longitude": lon,
        "hourly": hourly_vars,
        "models": "best_match",
        "timezone": "auto",
        "forecast_days": 1,
        "timeformat": "unixtime",
        "wind_speed_unit": "ms",
    }

    data = _get(url, params)
    hourly = data.get("hourly", {})
    utc_offset = float(data.get("utc_offset_seconds", 0))

    times = hourly.get("time", [])
    hourly_numpy = {"date": np.array(times, dtype=float).tolist(),
                    "timeOffset": utc_offset}

    for p in pressure_levels:
        for var in ["temperature", "wind_speed", "wind_direction", "geopotential_height"]:
            key = f"{var}_{p}hPa"
            val = hourly.get(key, [])
            hourly_numpy[key] = np.array(val, dtype=float).tolist()

    hourly_numpy["temperature_2m"]    = np.array(hourly.get("temperature_2m", []),    dtype=float).tolist()
    hourly_numpy["surface_pressure"]  = np.array(hourly.get("surface_pressure", []),  dtype=float).tolist()
    hourly_numpy["wind_gusts_10m"]    = np.array(hourly.get("wind_gusts_10m", []),    dtype=float).tolist()

    return hourly_numpy
