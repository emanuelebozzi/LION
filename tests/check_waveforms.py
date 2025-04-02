# %%
#check waveforms 

import obspy 

import matplotlib.pyplot as plt

from obspy import read

from obspy import Stream

import copy

import numpy as np

import pandas as pd


# %% [markdown]
# ##STATIONS

# %%
path_st = '/home/emanuele/data/emanuele/loki-das/cefalonia/Waveforms_Emanuele/noa2024prefy.mseed'
path_save_st = '/home/emanuele/data/emanuele/loki-das/cefalonia/events/prefy/stations/noa2024prefy_filt.mseed'
path =  '/home/emanuele/data/emanuele/loki-das/cefalonia/Waveforms_Emanuele/noa2024prefy_das.mseed'
save_path = '/home/emanuele/data/emanuele/loki-das/cefalonia/events/prefy/fiber/noa2024prefy_das_sel.mseed'
save_path_subsampling = '/home/emanuele/data/emanuele/loki-das/cefalonia/events/prefy/hybrid/noa2024prefy_das_sel_subsampling.mseed'

#path_st = '/home/emanuele/data/emanuele/loki-das/cefalonia/Waveforms_Emanuele/noa2024qalpq.mseed'
#path_save_st = '/home/emanuele/data/emanuele/loki-das/cefalonia/events/qalpq/stations/noa2024qalpq_filt.mseed'
#path =  '/home/emanuele/data/emanuele/loki-das/cefalonia/Waveforms_Emanuele/noa2024qalpq_das.mseed'
#save_path = '/home/emanuele/data/emanuele/loki-das/cefalonia/events/qalpq/fiber/noa2024qalpq_das_sel.mseed'
#save_path_subsampling = '/home/emanuele/data/emanuele/loki-das/cefalonia/events/qalpq/hybrid/noa2024qalpq_das_sel_subsampling.mseed'

#path_st = '/home/emanuele/data/emanuele/loki-das/cefalonia/Waveforms_Emanuele/noa2024pqvcv.mseed'
#path_save_st = '/home/emanuele/data/emanuele/loki-das/cefalonia/events/pqvcv/stations/noa2024pqvcv_filt.mseed'
#path =  '/home/emanuele/data/emanuele/loki-das/cefalonia/Waveforms_Emanuele/noa2024pqvcv_das.mseed'
#save_path = '/home/emanuele/data/emanuele/loki-das/cefalonia/events/pqvcv/fiber/noa2024pqvcv_das_sel.mseed'
#save_path_subsampling = '/home/emanuele/data/emanuele/loki-das/cefalonia/events/pqvcv/hybrid/noa2024pqvcv_das_sel_subsampling.mseed' 


#path_st = '/home/emanuele/data/emanuele/loki-das/cefalonia/Waveforms_Emanuele/noa2024rgjis.mseed'
#path_save_st = '/home/emanuele/data/emanuele/loki-das/cefalonia/events/rgjis/stations/noa2024rgjis_filt.mseed'
#path =  '/home/emanuele/data/emanuele/loki-das/cefalonia/Waveforms_Emanuele/noa2024rgjis_das.mseed'
#save_path = '/home/emanuele/data/emanuele/loki-das/cefalonia/events/rgjis/fiber/noa2024rgjis_das_sel.mseed'
#save_path_subsampling = '/home/emanuele/data/emanuele/loki-das/cefalonia/events/rgjis/hybrid/noa2024rgjis_das_sel_subsampling.mseed'

#path_st = '/home/emanuele/data/emanuele/loki-das/cefalonia/Waveforms_Emanuele/noa2024prlri.mseed'
#path_save_st = '/home/emanuele/data/emanuele/loki-das/cefalonia/events/prlri/stations/noa2024prlri_filt.mseed'
#path =  '/home/emanuele/data/emanuele/loki-das/cefalonia/events/noa2024prlri_das.mseed'
#save_path = '/home/emanuele/data/emanuele/loki-das/cefalonia/events/prlri/fiber/noa2024prlri_das_sel.mseed'
#save_path_subsampling = '/home/emanuele/data/emanuele/loki-das/cefalonia/events/prlri/hybrid/noa2024prlri_das_sel_subsampling.mseed'

#path_st = '/home/emanuele/data/emanuele/loki-das/cefalonia/Waveforms_Emanuele/noa2024qmbey.mseed'
#path_save_st = '/home/emanuele/data/emanuele/loki-das/cefalonia/events/qmbey/stations/noa2024qmbey_filt.mseed'
#path =  '/home/emanuele/data/emanuele/loki-das/cefalonia/Waveforms_Emanuele/noa2024qmbey_das.mseed'
#save_path = '/home/emanuele/data/emanuele/loki-das/cefalonia/events/qmbey/fiber/noa2024qmbey_das_sel.mseed'
#save_path_subsampling = '/home/emanuele/data/emanuele/loki-das/cefalonia/events/qmbey/hybrid/noa2024qmbey_das_sel_subsampling.mseed'


#path_st = '/home/emanuele/data/emanuele/loki-das/cefalonia/Waveforms_Emanuele/noa2024qsmcp.mseed'
#path_save_st = '/home/emanuele/data/emanuele/loki-das/cefalonia/events/qsmcp/stations/noa2024qsmcp_filt.mseed'
#path =  '/home/emanuele/data/emanuele/loki-das/cefalonia/Waveforms_Emanuele/noa2024qsmcp_das.mseed'
#save_path = '/home/emanuele/data/emanuele/loki-das/cefalonia/events/qsmcp/fiber/noa2024qsmcp_das_sel.mseed'
#save_path_subsampling = '/home/emanuele/data/emanuele/loki-das/cefalonia/events/qsmcp/hybrid/noa2024qsmcp_das_sel_subsampling.mseed'

st = read(path_st)

print(st)

sta_filt = st.copy()

sta_filt = sta_filt.detrend("demean")  
sta_filt = sta_filt.taper(0.05, type='cosine')
sta_filt = sta_filt.filter("bandpass", freqmin=2.0, freqmax = 15)  
sta_filt.normalize()

stations_sf = sta_filt[0].stats.sampling_rate 

print(stations_sf)

# Plot normalized traces
sta_filt.plot()


sta_filt.write(path_save_st, format="MSEED")

sta_filt.plot()


# %% [markdown]
# ##FIBER

# %%
st = read(path)

print(st)

# %%
#import traveltimes 

# %%
#stack traces (stack every gauge length) and associate the stacked trace to the middle

# %%
#remove loops (just use the id of the channel which have a location)


# Load CSV file
file_path = "/home/emanuele/data/emanuele/loki-das/cefalonia/geometry/das_channels_edited.csv"
df = pd.read_csv(file_path)

# Function to compute azimuth (bearing) between two lat/lon points
def calculate_azimuth(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    delta_lon = lon2 - lon1
    x = np.sin(delta_lon) * np.cos(lat2)
    y = np.cos(lat1) * np.sin(lat2) - (np.sin(lat1) * np.cos(lat2) * np.cos(delta_lon))
    azimuth = np.degrees(np.arctan2(x, y))
    return (azimuth + 360) % 360  # Normalize to 0-360 degrees

# Compute azimuth for each consecutive pair of points
azimuths = np.full(len(df), np.nan)  # Initialize with NaN
for i in range(1, len(df)):
    azimuths[i] = calculate_azimuth(df.loc[i-1, "latitude"], df.loc[i-1, "longitude"],
                                    df.loc[i, "latitude"], df.loc[i, "longitude"])

# Apply moving average smoothing
window_size = 20  # Adjust for desired smoothing
smoothed_azimuths = np.convolve(azimuths, np.ones(window_size)/window_size, mode='same')

# Plot azimuths

plt.plot(smoothed_azimuths, label="Smoothed azimuth")

# %%
import numpy as np
import matplotlib.pyplot as plt

# Define signal and noise windows (in seconds)
noise_start, noise_end = 0, 10  # Example: first 10 seconds as noise
signal_start, signal_end = 10, 20  # Example: 20-30s as signal

snr_db_sum = []  # List to store SNR values
snr_dict = {}  # Dictionary to store SNR values with trace ID as the key
processed_traces = []  # List to store processed traces

# Compute SNR for each trace
for idx, tr in enumerate(st):
    tr_new = tr.copy()  # Copy the original trace

    # Detrend, taper, filter, and normalize the trace
    tr_new = tr_new.detrend("demean")  
    tr_new = tr_new.taper(0.05, type='cosine')
    tr_new = tr_new.filter("bandpass", freqmin=2, freqmax=15)  
    tr_new.normalize()

    # Extract noise and signal windows
    noise_window = tr_new.slice(starttime=tr_new.stats.starttime + noise_start,
                                endtime=tr_new.stats.starttime + noise_end).data
    signal_window = tr_new.slice(starttime=tr_new.stats.starttime + signal_start,
                                 endtime=tr_new.stats.starttime + signal_end).data
    
    # Compute RMS values
    rms_noise = np.sqrt(np.mean(noise_window ** 2))
    rms_signal = np.sqrt(np.mean(signal_window ** 2))

    # Compute SNR (Ratio and dB)
    if rms_noise > 0:
        snr = rms_signal / rms_noise
        snr_db = 20 * np.log10(snr)
    else:
        snr = np.inf  # Avoid division by zero
        snr_db = np.inf

    # Store the SNR in both the list and dictionary with the trace ID (index)
    snr_db_sum.append(snr_db)
    snr_dict[idx] = snr_db  # Store SNR with trace ID as key

    # Add the processed trace to the new list (only change the data, not other attributes)
    tr_new.data = tr_new.data  # Replace the data with processed data (already done in the pipeline)
    processed_traces.append(tr_new)  # Append the processed trace to the list

# Convert the list to a numpy array for further use
snr_db_sum = np.array(snr_db_sum)

# Create a new Stream from the processed traces
from obspy import Stream
processed_stream = Stream(traces=processed_traces)

# Plot the SNR values
plt.figure(figsize=(10, 6))
plt.plot(snr_db_sum)
plt.xlabel("Trace Index")
plt.ylabel("SNR (dB)")
plt.title("SNR for Each Trace")
plt.show()

# Optionally, you can print the SNR dictionary to see trace IDs and their SNR values
# print(snr_dict)




# %%
import numpy as np
import matplotlib.pyplot as plt
from obspy import Stream

# Calculate the 80th percentile of the SNR values
snr_percentile_80 = np.percentile(snr_db_sum, 90)

# Divide into 10 subarrays with constraints
num_subarrays = 10
max_channels = len(azimuths) // num_subarrays  # Limit channels per group

sorted_indices = np.argsort(smoothed_azimuths)  # Sort by azimuth value
subarrays = np.array_split(sorted_indices, num_subarrays)  # Split into 10 parts

# Visualization of azimuth clustering
plt.figure(figsize=(10, 5))
colors = plt.cm.viridis(np.linspace(0, 1, num_subarrays))

for i, indices in enumerate(subarrays):
    plt.scatter(indices, smoothed_azimuths[indices], color=colors[i], label=f"Group {i+1}")

plt.xlabel("Channel Index")
plt.ylabel("Smoothed Azimuth")
plt.title("Azimuth Clustering into 10 Subarrays")

# Select the most azimuthally distinct and medium-high SNR point from each subarray
selected_indices = []
for subset in subarrays:
    if len(subset) > 0:
        # Calculate azimuthal range for the subset
        az_range = smoothed_azimuths[subset].max() - smoothed_azimuths[subset].min()
        
        # Get SNR values for the subset
        snr_values = snr_db_sum[subset]
        
        # Filter out the high SNR values (above the 80th percentile)
        valid_indices = subset[snr_values <= snr_percentile_80]
        
        # If there are valid indices, calculate the combined score for the remaining points
        if len(valid_indices) > 0:
            # Calculate the SNR for the valid points
            valid_snr_range = snr_db_sum[valid_indices]
            
            # Select the point that maximizes both azimuthal range and medium-high SNR
            combined_score = (smoothed_azimuths[valid_indices].max() - smoothed_azimuths[valid_indices].min()) + valid_snr_range
            selected_idx = valid_indices[np.argmax(combined_score)]
            selected_indices.append(selected_idx)

# Plot the selected points on the lat/lon trajectory
plt.figure(figsize=(8, 6))
plt.plot(df["longitude"], df["latitude"], color='gray', alpha=0.5, label="Trajectory")
plt.scatter(df.loc[selected_indices, "longitude"], df.loc[selected_indices, "latitude"], c='red', marker='o', label="Selected Points")

plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.title("Selected Points for Maximum Azimuthal Coverage and Medium-High SNR")
plt.legend()
plt.show()

# Now, create a new Stream for the selected processed traces
selected_traces = [processed_traces[idx] for idx in selected_indices]  # Extract the selected processed traces

# Convert selected traces into a Stream object
stream = Stream(selected_traces)

# Plot the selected traces using Obspy's Stream.plot()
stream.plot()


# %%
#Within those sections identify max SNR and the associated channel

# %%
print(str(selected_indices))

# %%
selected_stations = list(map(str, selected_indices))

selected_stations_hybrid = selected_stations

#selected_stations = ["0223", "0500", "0800", "1000", "1200", "1400", "1663", "1723", "2000", "2487", "3039", "4000", "5000", "6000", "7150", "7185"]  # Replace with your station IDs
#selected_stations = ["0223", "0500", "1200", "3039", "4000", "5000", "7150", "7185", "7350"]  # Replace with your station IDs
#selected_stations_hybrid = ["0223",  "1200", "3039", "5000",  "7185"]  # Replace with your station IDs


filtered_stream = Stream()
subsampled_filtered_stream = Stream()

g = st.copy()

for a in selected_stations: 
    
    b = g.select(station=str(a))
    filtered_stream.append(b[0])

ch = ["HHE", "HHN", "HHZ"]

e = st.copy()

for a in selected_stations_hybrid: 
    
    b_new = e.select(station=str(a))
    c_new = b_new.resample(stations_sf)

    trace_N = copy.deepcopy(c_new[0]) 
    trace_E = copy.deepcopy(c_new[0]) 
    trace_N.stats.channel = "HHE"
    trace_E.stats.channel = "HHN"

    subsampled_filtered_stream.append(c_new[0])
    subsampled_filtered_stream += trace_N
    subsampled_filtered_stream += trace_E


filtered_stream.write(save_path, format="MSEED")

#sta_filt_only_z = sta_filt.copy()

#sta_filt_only_z = sta_filt_only_z.select(channel="*Z")



subsampled_filtered_stream_hydrid = subsampled_filtered_stream + sta_filt

print(subsampled_filtered_stream)


subsampled_filtered_stream_hydrid.write(save_path_subsampling, format="MSEED")





# %%
das = read(save_path_subsampling)

print(das)

das_filt = das.copy()

das = das.detrend("demean")  
das = das.taper(0.05, type='cosine')
das_filt = das.filter("bandpass", freqmin=2, freqmax = 15)  
das_filt.normalize()

# Plot normalized traces
das_filt.plot()

das_filt.write(save_path_subsampling, format="MSEED")

# %%



