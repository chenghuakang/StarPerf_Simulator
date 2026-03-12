import numpy as np

def norm_rate(OBJs, oi, data_ext_dict, data_ext_prev_dict, t, dT,
              min_elevation_deg, type, metadata=None):
    # plugin that applies a rate equal to a normal distribution with mean metadata["mean"] and std metadata["std"]
    # the different rate is applied only for new links (i.e., links that were not present in the previous snapshot)
    if metadata is None:
        return None
    mean = metadata.get("mean", 50)
    std = metadata.get("std", 10)
    delay_data = data_ext_dict.get("delay")
    delay_data_prev = data_ext_prev_dict.get("delay")
    linked_sats = np.where(delay_data[oi, :] != 0)[0]
    linked_sats_prev = np.where(delay_data_prev[oi, :] != 0)[0]
    rate_updated = np.zeros(data_ext_dict["rate"].shape[0])  # initialize the array of updated rates with zeros.
    for i in range(data_ext_dict["rate"].shape[0]):
        if i not in linked_sats_prev and i in linked_sats:
            rate = np.random.normal(mean, std)
            rate = max(0, rate)  # ensure non-negative rate
            rate = min(rate, 1000)  # cap rate to a maximum value (e.g., 1000 Mbit/s)
            rate_updated[i] = rate
    return rate_updated
