import requests
import json

def fetch_windy_data(api_key, lat, lon, model, parameters, levels=None):
    url = "https://api.windy.com/api/point-forecast/v2"
    headers = {
        'Content-Type': 'application/json'
    }
    body = {
        "lat": lat,
        "lon": lon,
        "model": model,
        "parameters": parameters,
        "levels": levels,
        "key": api_key
    }
    
    response = requests.post(url, headers=headers, data=json.dumps(body))
    
    if response.status_code == 200:
        data = response.json()
        return data
    else:
        raise Exception(f"Error {response.status_code}: {response.text}")

# Example usage
api_key = "HcwZQ8f1w4Qu2sSifMc6qt1VnyeRjsgm"
lat = 35.009
lon = -115.473
model = "gfs"
parameters = ["wind", "windGust", "gh", "pressure", "temp"]
levels = ["surface", "1000h","950h","925h","900h","850h","800h", "300h"]

try:
    data = fetch_windy_data(api_key, lat, lon, model, parameters, levels)
    print(json.dumps(data, indent=4))
except Exception as e:
    print(e)