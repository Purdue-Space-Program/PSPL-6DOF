import openmeteo_requests
import requests_cache
from retry_requests import retry

def getCurrentWeather(lat=39.7392, lon=-104.9847):
    """Fetch current weather data from Open-Meteo and return a dict
       (automatically converts to MATLAB struct)."""

    # Setup Open-Meteo API client with caching and retry
    cache_session = requests_cache.CachedSession('.cache', expire_after=3600)
    retry_session = retry(cache_session, retries=5, backoff_factor=0.2)
    openmeteo = openmeteo_requests.Client(session=retry_session)

    # Define query parameters
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
        "wind_speed_unit": "ms",
    }

    # Make API request
    responses = openmeteo.weather_api(url, params=params)
    response = responses[0]
    current = response.Current()

    # Return results as a Python dict (MATLAB converts to struct)
    data = {
        "latitude": response.Latitude(),
        "longitude": response.Longitude(),
        "elevation": response.Elevation(),
        "timezone": response.Timezone(),
        "timezone_abbrev": response.TimezoneAbbreviation(),
        "utc_offset_seconds": response.UtcOffsetSeconds(),
        "time": current.Time(),
        "temperature_2m": current.Variables(0).Value(),
        "surface_pressure": current.Variables(1).Value(),
        "relative_humidity_2m": current.Variables(2).Value(),
        "wind_speed_10m": current.Variables(3).Value(),
        "wind_direction_10m": current.Variables(4).Value(),
        "wind_gusts_10m": current.Variables(5).Value(),
    }

    return data
