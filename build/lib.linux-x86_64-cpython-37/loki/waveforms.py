import os
from obspy.core import read
from datetime import datetime
import h5py
from scipy.signal import detrend, butter, filtfilt
from obspy import Stream, Trace, UTCDateTime


'''
This class is used to read the station and DAS waveforms. 
They are kept separate to add eventually the coherece matrices 

#first load the available waveforms for the stations, then create the station list 
#then again for DAS 
'''

class Waveforms:

    def __init__(self, tobj, data_path, event_path, extension_sta='*', extension_das='CANDAS2_2023-01-07_10-48-10.h5', comps=['E','N','Z'], freq=None):

        #seismometers 

        if not os.path.isdir(event_path):
            raise ValueError('Error: data path does not exist')
        try:
            self.load_sta_waveforms(event_path, extension_sta, comps, freq)
        except:
            raise WaveformLoadingError('Error: data (station) not read for the event: %s' %(event_path))
        self.station_list()  

        #DAS 

        try:
            self.load_das_waveforms(data_path, extension_das, tobj, freq)


        except:
            raise WaveformLoadingError('Error: data (DAS) not read for the event: %s' %(event_path))
        self.channel_list()  

    #calle methods (in order)

    def load_sta_waveforms(self, event_path, extension_sta, comps, freq):

 
        files=os.path.join(event_path,extension_sta)

        print('Loading station event: ' + str(files))
        traces_sta=read(files)
        
        if freq:
            traces_sta.detrend('demean')
            traces_sta.detrend('linear')
            if len(freq) == 1:
                traces_sta.filter("highpass", freq=freq[0])
            elif len(freq) == 2:
                traces_sta.filter("bandpass", freqmin=freq[0], freqmax=freq[1])

        #create the stream for stations

        self.stream_sta={}
        for comp in comps:
            self.stream_sta[comp]={}
            for tr in traces_sta:
                if tr.stats.channel[-1]==comp:
                    dtime=datetime.strptime(str(tr.stats.starttime),"%Y-%m-%dT%H:%M:%S.%fZ")
                    self.stream_sta[comp][tr.stats.station]=[dtime, tr.stats.delta, tr.data]



    def load_das_waveforms(self, data_path, extension_das,tobj, freq):
        id_das_stations = tobj.id_das_stations
        delta_das =  tobj.delta_das
        #read the events in .h5 format 
        files=os.path.join(data_path,extension_das)
        
        print('Loading DAS event: ' + str(files))
        #traces_das=read(files)      

        # Open the .h5 file in read mode
        with h5py.File(files, 'r') as file:
             
            #this is valid for the Gran Canaria dataset 
            dataset_name = 'data' 
            dataset = file[dataset_name]
            data = dataset[:]  
            #print("Data:", data)

        num_channel = len(data[:,1])


#utilizzare scipy

        if freq:
            for i in range(0, num_channel): 

                data[i,:] = detrend(data[i,:], type='constant')
                data[i,:] = detrend(data[i,:], type='linear')
                b, a = butter(N=4, Wn=freq[0] if len(freq) == 1 else [freq[0], freq[1]], btype='high' if len(freq) == 1 else 'band')
                data[i, :] = filtfilt(b, a, data[i, :])

        #create the stream for DAS

        self.stream_das = Stream() 

        for i in range(0, len(id_das_stations)-6):
            
            #print(data.shape)
            #print(i)
            #print(len(id_das_stations))
            #print(id_das_stations[i])

            trace = Trace(data=data[id_das_stations[i], :])  # Use `row` directly
            trace.stats.delta = delta_das
            #starttime = list(self.stream_sta.values())[0].stats.starttime
            #trace.stats.starttime = starttime
            trace.stats.starttime = UTCDateTime("2025-01-01T00:00:00")
            trace.stats.station = str(id_das_stations[i])  # Unique station name for each trace
            trace.stats.channel = "E"  # Same channel for all traces
            
            # Add the trace to the Stream
            self.stream_das.append(trace)

            #print('the DAS stream is: ', self.stream_das)

    def station_list(self):
        data_stalist=[]
        for comp in (self.stream_sta).keys():
            for sta in (self.stream_sta[comp]).keys():
                if sta not in data_stalist:
                    data_stalist.append(sta)
        self.data_stations=set(data_stalist)

    def channel_list(self):
        data_chlist = []
        for trace in self.stream_das:
            station = trace.stats.station  # Get the station name from trace metadata
            if station not in data_chlist:
                data_chlist.append(station)
        self.data_channels = set(data_chlist)



# forme d'onda da tenere separate, per pesare la DAS e sismometri ugualmente 
#poi sommare le matrici di coerenza 

class WaveformLoadingError(Exception):
    pass
