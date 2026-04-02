import json

import numpy as np

def slant_rate(OBJs, oi, data_ext_dict, data_ext_prev_dict, t, dT,
              min_elevation_deg, type, metadata=None):
    # Plugin that applies a rate equal to the maximum one scaled by the rate (zenit_altitude/distance)^2, where distance is calculated from the delay and zenit_altitude is the altitude at which the rate is zenith_rate. 
    # Additional atmosferic loss can be applied at zenith. This is done by multiplying the rate by 10**(-zenit_atm_loss_db * (relative_air_mass-1) / 10), where zenit_atm_loss_db is the atmospheric loss in dB at zenith.
    # This models is valid for satellite-ground and satellite-user links for low SNR regime where ln(1+x) ~ x, and elevation angle greather than 10 degrees. 
    if metadata is None:
        return None
    try:
        with open(metadata, "r") as f:
            metadata_dict = json.load(f)
    except Exception as e:
        print(f"Error loading metadata from {metadata}: {e}")
        return None
    zenith_rate = metadata_dict.get("zenith_rate_mbps", 50) # max rate in Mbit/s when the elevation angle is 90 degrees
    zenith_altitude = metadata_dict.get("zenith_altitude_km", 1200) # altitude in km at which the rate is zenith_rate.
    zenit_atm_loss_db = metadata_dict.get("zenith_atm_loss_db", 0) # atmospheric loss in dB at zenith.
    delay_data = data_ext_dict.get("delay")
    elevation_angle_data = data_ext_dict.get("angle")
    rate_updated = np.zeros(data_ext_dict["rate"].shape[0])  # initialize the array of updated rates with zeros.
    for i in range(data_ext_dict["rate"].shape[0]):
        if delay_data[oi, i] == 0: # if there is no link, skip
            continue
        distance_km = delay_data[oi, i] * 3e5 # distance in km (assuming speed of light in km/s)
        elevation_angle_deg = elevation_angle_data[oi, i]
        FSPL_factor = (zenith_altitude / distance_km)**2
        relative_air_mass = 1.0 * distance_km / zenith_altitude
        #relative_air_mass = 1/np.sin(np.radians(elevation_angle_deg)) # relative air mass based on elevation angle with respect to zenith
        atm_loss_factor = 10**(-zenit_atm_loss_db * (relative_air_mass-1) / 10) # atmospheric loss factor based on zenith atmospheric loss and relative air mass. Note: atmospheric loss increases as exp(-alpha * relative_air_mass), where alpha is the atmospheric loss coefficient in dB and relative_air_mass is the relative air mass. (Lambert-Beer law)
        
        rate_updated[i] = zenith_rate * FSPL_factor * atm_loss_factor 
    return rate_updated
