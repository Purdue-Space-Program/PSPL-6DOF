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
        
        # Convert the JSON to a Python dictionary that is easy to handle in MATLAB
        matlab_data = convert_to_matlab_friendly_format(data)
        return matlab_data
    else:
        raise Exception(f"Error {response.status_code}: {response.text}")

def convert_to_matlab_friendly_format(data):
    """
    Converts JSON data into a dictionary-like structure that MATLAB can interpret.
    """
    matlab_data = {}
    
    if isinstance(data, dict):
        for key, value in data.items():
            if isinstance(value, dict):
                matlab_data[key] = convert_to_matlab_friendly_format(value)
            elif isinstance(value, list):
                matlab_data[key] = value
            else:
                matlab_data[key] = value
    elif isinstance(data, list):
        matlab_data = data
    
    return matlab_data

def main(api_key, lat, lon, model, parameters, levels):
    matlab_friendly_data = fetch_windy_data(api_key, lat, lon, model, parameters, levels)
    
    # Print the data as JSON for MATLAB to capture and process
    print(json.dumps(matlab_friendly_data))  # MATLAB will capture this output

if __name__ == "__main__":
    api_key = "HcwZQ8f1w4Qu2sSifMc6qt1VnyeRjsgm"
    lat = 40.439
    lon = -86.9668
    model = "namConus"
    parameters = ["wind", "windGust", "gh", "pressure", "temp"]
    levels = ["surface", "1000h","950h","925h","900h","850h","800h","700h","600h","500h", "300h"]
    
    main(api_key, lat, lon, model, parameters, levels)
