import openmeteo_requests
import requests_cache
from retry_requests import retry
import numpy as np

FORECAST_URL    = "https://api.open-meteo.com/v1/forecast"
PRESSURE_LEVELS = [1000, 925, 850, 700, 500, 300, 200, 100, 50]

cache_session = requests_cache.CachedSession('.omcache', expire_after=3600)
retry_session = retry(cache_session, retries=5, backoff_factor=0.2)
openmeteo     = openmeteo_requests.Client(session=retry_session)


def getHourlyWeather(lat, lon):
    hourly_vars = []
    for p in PRESSURE_LEVELS:
        hourly_vars += [
            f"temperature_{p}hPa",
            f"wind_speed_{p}hPa",
            f"wind_direction_{p}hPa",
            f"geopotential_height_{p}hPa",
        ]
    hourly_vars += ["temperature_2m", "surface_pressure", "wind_gusts_10m"]

    params = {
        "latitude":        lat,
        "longitude":       lon,
        "hourly":          hourly_vars,
        "forecast_days":   1,
        "timezone":        "auto",
        "timeformat":      "unixtime",
        "wind_speed_unit": "ms",
    }

    responses = openmeteo.weather_api(FORECAST_URL, params=params)
    response  = responses[0]

    hourly     = response.Hourly()
    utc_offset = float(response.UtcOffsetSeconds())
    times      = np.arange(hourly.Time(), hourly.TimeEnd(), hourly.Interval())

    result = {
        "date":       times.tolist(),
        "timeOffset": utc_offset,
    }

    # Variables come back in the same order as requested
    for i, key in enumerate(hourly_vars):
        result[key] = hourly.Variables(i).ValuesAsNumpy().tolist()

    return result


def getCurrentWeather(lat, lon):
    current_vars = [
        "temperature_2m", "surface_pressure", "relative_humidity_2m",
        "wind_speed_10m", "wind_direction_10m", "wind_gusts_10m",
    ]

    params = {
        "latitude":        lat,
        "longitude":       lon,
        "current":         current_vars,
        "forecast_days":   1,
        "timezone":        "auto",
        "timeformat":      "unixtime",
        "wind_speed_unit": "ms",
    }

    responses = openmeteo.weather_api(FORECAST_URL, params=params)
    response  = responses[0]
    current   = response.Current()

    return {
        "latitude":             response.Latitude(),
        "longitude":            response.Longitude(),
        "elevation":            response.Elevation(),
        "timezone":             str(response.Timezone()),
        "utc_offset_seconds":   float(response.UtcOffsetSeconds()),
        "time":                 current.Time(),
        "temperature_2m":       current.Variables(0).Value(),
        "surface_pressure":     current.Variables(1).Value(),
        "relative_humidity_2m": current.Variables(2).Value(),
        "wind_speed_10m":       current.Variables(3).Value(),
        "wind_direction_10m":   current.Variables(4).Value(),
        "wind_gusts_10m":       current.Variables(5).Value(),
    }
