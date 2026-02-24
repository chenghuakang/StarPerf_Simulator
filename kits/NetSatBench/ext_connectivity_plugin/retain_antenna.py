# #!/usr/bin/env python3
import numpy as np

def retain_antenna(OBJs, oi, data_ext_dict, data_ext_prev_dict, t, dT,
                min_elevation_deg, type, metadata=None):
    """
    Antenna management plugin with link-retention policy.

    Enforces antenna constraints for user and ground station objects and
    updates the connectivity (del_ext[oi, :]) accordingly to minimize handhovers.

    Policy
    ------
    - Each object has `antenna_count` antennas.
    - Keep at most `antenna_count` active links

    Link selection (when more candidates than allowed)
    --------------------------------------------------
    1) Prefer already active links (minimize unnecessary handovers).
    2) Prefer links with increasing elevation angle (angle rising).
    3) Prefer links with higher elevation angle.

    Returns
    -------
    numpy.ndarray or None
        Updated row del_ext[oi, :]. Returning None means no change.
    """


    # Only apply to user / ground station
    if type not in ("gs", "user"):
        return None

    obj = OBJs[oi]
    antenna_count = int(getattr(obj, "antenna_count", 0))

    # If no antennas, drop all links (or return None: choose what your framework expects)
    if antenna_count <= 0:
        # If you prefer "no change" instead, replace with: return None
        delay = data_ext_dict.get("delay")
        if isinstance(delay, np.ndarray) and delay.ndim == 2:
            out = delay[oi, :].copy()
            out[:] = 0
            return out
        return None


    delay_raw = data_ext_dict.get("delay")
    angle_raw = data_ext_dict.get("angle")

    if not isinstance(delay_raw, np.ndarray) or delay_raw.ndim != 2:
        return None  # cannot operate
    if not isinstance(angle_raw, np.ndarray) or angle_raw.ndim != 2:
        return None  # cannot operate

    # Work on copies so inputs are never mutated
    delay_data = delay_raw.copy()
    angle_data = angle_raw.copy()

    # Previous snapshot (copy as well to avoid any chance of shared-memory surprises)
    delay_prev_raw = data_ext_prev_dict.get("delay")
    angle_prev_raw = data_ext_prev_dict.get("angle")

    if isinstance(delay_prev_raw, np.ndarray) and delay_prev_raw.shape == delay_data.shape:
        delay_data_prev = delay_prev_raw.copy()
    else:
        delay_data_prev = delay_data.copy()

    if isinstance(angle_prev_raw, np.ndarray) and angle_prev_raw.shape == angle_data.shape:
        angle_data_prev = angle_prev_raw.copy()
    else:
        angle_data_prev = angle_data.copy()

    # --- Core logic  ---
    linked_sats = np.where(delay_data[oi, :] != 0)[0]
    linked_sats_prev = np.where(delay_data_prev[oi, :] != 0)[0]

    link_sat_old = np.intersect1d(linked_sats, linked_sats_prev)
    linked_sat_old_rising = link_sat_old[angle_data[oi, link_sat_old] > angle_data_prev[oi, link_sat_old]]
    linked_sat_old_rising = linked_sat_old_rising[np.argsort(angle_data[oi, linked_sat_old_rising])]

    linked_sat_old_setting = link_sat_old[angle_data[oi, link_sat_old] <= angle_data_prev[oi, link_sat_old]]
    linked_sat_old_setting = linked_sat_old_setting[np.argsort(-angle_data[oi, linked_sat_old_setting])]

    linked_sat_new = np.setdiff1d(linked_sats, linked_sats_prev)
    linked_sat_new_rising = linked_sat_new[angle_data[oi, linked_sat_new] > angle_data_prev[oi, linked_sat_new]]
    linked_sat_new_rising = linked_sat_new_rising[np.argsort(angle_data[oi, linked_sat_new_rising])]

    linked_sat_new_setting = linked_sat_new[angle_data[oi, linked_sat_new] <= angle_data_prev[oi, linked_sat_new]]
    linked_sat_new_setting = linked_sat_new_setting[np.argsort(-angle_data[oi, linked_sat_new_setting])]

    linked_sats_sorted = np.concatenate(
        (linked_sat_old_rising, linked_sat_old_setting, linked_sat_new_rising, linked_sat_new_setting)
    )

    linked_sats_updated = linked_sats_sorted[:antenna_count]
    linked_sat_to_delete = np.setdiff1d(linked_sats, linked_sats_updated)

    delay_data[oi, linked_sat_to_delete] = 0

    # Return only the row (copy is already detached from original inputs)
    return delay_data[oi, :]

    
    

    
    
