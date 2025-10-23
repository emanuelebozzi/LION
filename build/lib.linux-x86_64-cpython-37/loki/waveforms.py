import os
from obspy.core import read
from datetime import datetime


class Waveforms:

    # 1. __init__ remains the same, but 'comps' will be handled in load_waveforms
    def __init__(self, event_path, extension_sta='*', comps=None, freq=None):
        if not os.path.isdir(event_path):
            raise ValueError('Error: data path does not exist')
        try:
            self.load_waveforms(event_path, extension_sta, comps, freq)
        except Exception as e:
            # Added e for better error reporting
            raise WaveformLoadingError('Error: data not read for the event: %s. Details: %s' %(event_path, e))
        self.station_list()

    def station_list(self):
        data_stalist=[]
        for comp in (self.stream).keys():
            for sta in (self.stream[comp]).keys():
                if sta not in data_stalist:
                    data_stalist.append(sta)
        self.data_stations=set(data_stalist)

        #print('data_stations', self.data_stations)
#---------------------------------------------------------------------------------------------------
    # 2. Key changes are here: dynamic component loading and optional filtering
    def load_waveforms(self, event_path, extension_sta, comps, freq):
        files=os.path.join(event_path,extension_sta)
        #print('files', files)
        
        # Read all traces from files
        traces=read(files)
        
        if freq:
            traces.detrend('demean')
            traces.detrend('linear')
            if len(freq) == 1:
                traces.filter("highpass", freq=freq[0])
            elif len(freq) == 2:
                traces.filter("bandpass", freqmin=freq[0], freqmax=freq[1])
        
        self.stream={}
        
        # Iterate over all read traces
        for tr in traces:
            # Dynamically determine the component (e.g., 'Z', 'N', 'E' from 'BHZ')
            current_comp = tr.stats.channel[-1] 
            
            # Check if filtering is required (i.e., 'comps' was provided)
            # If 'comps' is None, it reads ALL components.
            # If 'comps' is provided, it only processes components in the list.
            if comps is not None and current_comp not in comps:
                continue # Skip this trace
            
            # Initialize the dictionary for this component if it doesn't exist
            if current_comp not in self.stream:
                self.stream[current_comp] = {}
            
            # Store the data
            dtime=datetime.strptime(str(tr.stats.starttime),"%Y-%m-%dT%H:%M:%S.%fZ")
            self.stream[current_comp][tr.stats.station]=[dtime, tr.stats.delta, tr.data]


class WaveformLoadingError(Exception):
    pass