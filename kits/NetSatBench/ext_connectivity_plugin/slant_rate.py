import json

import numpy as np

def slant_rate(OBJs, oi, data_ext_dict, data_ext_prev_dict, t, dT,
              min_elevation_deg, type, metadata=None):
    # Plugin that applies a rate equal to the zenith one scaled by free-space and atmosferic losses using Shannon capacity
  
    if metadata is None:
        return None
    try:
        with open(metadata, "r") as f:
            metadata_dict = json.load(f)
    except Exception as e:
        print(f"Error loading metadata from {metadata}: {e}")
        return None
    zenith_rate_mbps = metadata_dict.get("zenith_rate_mbps", 50) # max rate in Mbit/s when the elevation angle is 90 degrees
    zenith_altitude_km = metadata_dict.get("zenith_altitude_km", 1200) # altitude in km at which the rate is zenith_rate.
    zenit_atm_loss_db = metadata_dict.get("zenith_atm_loss_db", 0.5) # atmospheric loss in dB at zenith ITU-R (P.676 e P.840) .
    zenith_SNR_db = metadata_dict.get("zenith_SNR_db", 12) # SNR in dB at zenith.
    delay_data = data_ext_dict.get("delay")
    elevation_angle_data = data_ext_dict.get("angle")
    rate_updated = np.zeros(data_ext_dict["rate"].shape[0])  # initialize the array of updated rates with zeros.
    for i in range(data_ext_dict["rate"].shape[0]):
        if delay_data[oi, i] == 0: # if there is no link, skip
            continue
        slant_distance_km = delay_data[oi, i] * 3e5 # distance in km (assuming speed of light in km/s)
        elevation_angle_deg = elevation_angle_data[oi, i]
        FSPL_factor = (zenith_altitude_km / slant_distance_km)**2
        relative_air_mass = 1.0 * slant_distance_km / zenith_altitude_km
        atm_PL_factor = 10**(-zenit_atm_loss_db * (relative_air_mass-1) / 10) # atmospheric loss factor based on zenith atmospheric loss and relative air mass. Note: atmospheric loss increases as exp(-gamma * relative_air_mass), where gamma is the atmospheric loss coefficient in dB and relative_air_mass is the relative air mass. (Lambert-Beer law)
        A_tot = FSPL_factor * atm_PL_factor
        zenith_SNR = 10**(zenith_SNR_db / 10) # convert zenith SNR from dB to linear scale
        rate_updated[i] = zenith_rate_mbps * np.log2(1+zenith_SNR * A_tot) /  np.log2(1+zenith_SNR)
    return rate_updated
