import numpy as num
import DET_STALTA     # C
import LOC_STALTA     # C
from datetime import datetime
from obspy import UTCDateTime

class Stacktraces:

    def __init__(self, tobj, wobj, **inputs):
        #check input objects tobj=traveltime_object wobj=Raw_waveform_object
        self.check_sampling_rate_sta(wobj)
        self.check_sampling_rate_das(wobj)
        self.check_starting_time_sta(wobj)
        self.check_starting_time_das(wobj)


        #loki_input to change
        
        if ('vfunc' in inputs) or ('hfunc' in inputs):
            # need to calculate characteristic functions from input data
            vfunc=inputs['vfunc']
            hfunc=inputs['hfunc']
            epsilon=inputs['epsilon']
            derivative=inputs['derivative']
            self.loki_input(wobj, tobj, derivative)
            self.characteristic_function(vfunc, hfunc, epsilon)
        else:
            # no need to calculate characteristic function from input data
            # directly input data of P and S component
            if 'normthrd' in inputs:
                normalize = inputs['normthrd']
            else:
                normalize = False
            self.loki_input(wobj, tobj, derivative=False, direct_input=True, normalize=normalize)
            
            # compute array element wise power over the input probabilities if needed
            if 'ppower' in inputs:
                self.obs_dataV_sta = self.obs_dataV_sta**inputs['ppower']
                self.obs_dataH_sta  = self.obs_dataH_sta**inputs['ppower']
                self.obs_dataV_das = self.obs_dataV_das**inputs['ppower']
                self.obs_dataH_das  = self.obs_dataH_das**inputs['ppower']

    def check_sampling_rate_sta(self,wobj):
        intsamp=1E6
        deltas_sta=[]
        for comp in (wobj.stream_sta).keys():
            for sta in (wobj.stream_sta[comp]).keys():
                deltas_sta.append(wobj.stream_sta[comp][sta][1])
        deltas_sta=num.array(deltas_sta)
        ideltas_sta=num.unique((deltas_sta*intsamp).astype(int))
        if num.size(ideltas_sta)==1:
            self.deltat_sta=deltas_sta[0]
        else:
            raise ValueError('(STA) Error!! All trace must have the same sampling rate')


    def check_sampling_rate_das(self,wobj):
        intsamp=1E6
        deltas_das = []
        for trace in wobj.stream_das:  # Iterate over traces in the Stream
            deltas_das.append(trace.stats.delta)  # Sampling interval in seconds

        deltas_das = num.array(deltas_das)
        ideltas_das = num.unique((deltas_das * intsamp).astype(int))
        if num.size(ideltas_das) == 1:
            self.deltat_das = deltas_das[0]
        else:
            raise ValueError('Error (DAS)!! All traces must have the same sampling rate')


    def check_starting_time_sta(self,wobj):
        dtimes_sta=[]
        self.ns_sta=0
        for comp in (wobj.stream_sta).keys():
            for sta in (wobj.stream_sta[comp]).keys():
                dtimes_sta.append(wobj.stream_sta[comp][sta][0])
                if self.ns_sta<num.size(wobj.stream_sta[comp][sta][2]):
                    self.ns_sta=num.size(wobj.stream_sta[comp][sta][2])
        self.dtime_max_sta=max(dtimes_sta)

        self.evid_sta=(self.dtime_max_sta).isoformat()

    def check_starting_time_das(self, wobj):
        dtimes_das = []
        self.ns_das = 0

        # Iterate over traces in the Stream

        for trace in wobj.stream_das:

            # Assuming each trace has a stats.starttime attribute and data array
            #dtimes_das.append(trace.stats.starttime)  # Collect starting times

            das_start_time = trace.stats.starttime.strftime("%Y-%m-%dT%H:%M:%S.%f")
            datetime_obj = datetime.strptime(das_start_time, "%Y-%m-%dT%H:%M:%S.%f")

                # Convert to the desired format
            formatted_date = datetime_obj.strftime("%Y-%m-%d %H:%M:%S")
            dtimes_das.append(formatted_date)
            current_ns = len(trace.data)  # Assuming trace.data is array-like
            if self.ns_das < current_ns:
                self.ns_das = current_ns

        # Determine the maximum starting time
        self.dtime_max_das = max(dtimes_das)
        self.dtime_max_das = datetime.strptime(self.dtime_max_das, "%Y-%m-%d %H:%M:%S")
  
        self.evid_das = self.dtime_max_das.isoformat()


# before changing this check the attributed of wobj 

    def loki_input(self, wobj, tobj, derivative, direct_input=False, normalize=True):
        if direct_input:
            # directly use input data as characteristic function

            self.obs_dataV_sta = self.select_data('P', wobj, wobj.data_stations, derivative, normalize)
            self.obs_dataH_sta = self.select_data('S', wobj, wobj.data_stations, derivative, normalize)


            self.obs_dataV_das = self.select_data('P', wobj, wobj.data_channels, derivative, normalize)
            self.obs_dataH_das = self.select_data('S', wobj, wobj.data_channels, derivative, normalize)
        else:
            # normal input, input 1- or 3-component data for calculating characteristic
            # function later
            self.comp_sta=tuple((wobj.stream_sta).keys())
            if len(self.comp_sta)==3:
                self.xtr_sta=self.select_data_sta(self.comp_sta[0], wobj, wobj.data_stations, derivative, normalize)
                self.ytr_sta=self.select_data_sta(self.comp_sta[1], wobj, wobj.data_stations, derivative, normalize)
                self.ztr_sta=self.select_data_sta(self.comp_sta[2], wobj, wobj.data_stations, derivative, normalize)
            elif len(self.comp)==1:
                self.ztr_sta=self.select_data_sta(self.comp_sta[0], wobj, wobj.data_stations, derivative, normalize)
            else:
                raise ValueError('Traces (STA) must have 1 or 3 components!')
            

            self.comp_das=tuple('E')
            #if len(self.comp)==3:
            #    self.xtr_das=self.select_data(self.comp[0], wobj, wobj.data_channels, derivative, normalize)
            #    self.ytr_das=self.select_data(self.comp[1], wobj, wobj.data_channels, derivative, normalize)
            #    self.ztr_das=self.select_data(self.comp[2], wobj, wobj.data_channels, derivative, normalize)
            #elif len(self.comp)==1:
            self.xtr_das=self.select_data_das(self.comp_das, wobj, wobj.data_channels, derivative, normalize)
            #else:
            #    raise ValueError('Traces (DAS) must have 1 or 3 components!')



    def select_data_sta(self, comp, wobj, db_stations, derivative, normalize):

    #stations 

        self.stations=tuple(wobj.data_stations & db_stations)  # find stations that are in common
        self.nstation=num.size(self.stations)
        tr_sta=num.zeros([self.nstation,self.ns_sta])
        stream=wobj.stream_sta


        for i, sta in enumerate(self.stations):
            found = False  # Flag to check if the station is found in any component
            for comp in stream.keys():  # Iterate over components ('E', 'N', 'Z')
                if sta in stream[comp]:  # Check if the station exists in the current component
                    #print('self.dtime_max_sta:', self.dtime_max_sta, stream[comp][sta][0])
                    # Get station data
                    nstr = num.size(stream[comp][sta][2])  # Length of the trace array
                    idt = int((self.dtime_max_sta - stream[comp][sta][0]).total_seconds() / self.deltat_sta)  # Time offset index
                    tr_sta[i, 0:nstr-idt] = stream[comp][sta][2][idt:]  # Assign trace data
                    
                    found = True  # Mark as found
                    break  # Exit loop once the station is found

            if not found:
                print(f"Warning: Station {sta} not found in any component")
            
            if derivative:
                # calculate derivatives of input data
                tr_sta[i,1:self.ns_sta]=((tr_sta[i,1:]-tr_sta[i,0:self.ns_sta-1])/self.deltat_sta)
                tr_sta[i,0]=0.
            
            if isinstance(normalize, float):
                # normalize data only if the absolute maxima is largar than a 
                # certain input threshold
                trmax = num.max(num.abs(tr_sta[i,:]))
                if trmax >= normalize:
                    tr_sta[i,:]=tr_sta[i,:]/trmax
            elif normalize:
                # normalize data by the absolute data maxima
                if num.max(num.abs(tr_sta[i,:])) > 0:
                    tr_sta[i,:]=tr_sta[i,:]/num.max(num.abs(tr_sta[i,:]))
        
        return tr_sta

    def select_data_das(self, comp_das, wobj, db_stations, derivative, normalize): 
        self.channels = tuple(wobj.data_channels & db_stations)  # Find stations that are in common
        self.nchannel = num.size(self.channels)
        tr_das = num.zeros([self.nchannel, self.ns_das])  # Initialize output array
        stream = wobj.stream_das  # Stream object containing traces

        for i, sta in enumerate(self.channels):

            # Convert station index to integer
            sta_idx = int(sta)

            # Validate that the station index is within the stream bounds
            if sta_idx < 0 or sta_idx >= len(stream):
                print(f"Invalid station index: {sta_idx}. Skipping.")
                continue

            # Get the trace for the station
            trace = stream[sta_idx]

            # Access the trace data properly
            trace_data = trace.data

            # Get trace length
            nstr = len(trace_data)

            # Calculate the index offset (idt)
            idt = int((UTCDateTime(self.dtime_max_das) - trace.stats.starttime) / self.deltat_das)

            # Validate slice boundaries
            if idt < 0 or idt >= nstr:
                print(f"Invalid slice index {idt} for trace {sta_idx}. Skipping.")
                continue

            # Assign the data to the output array
            tr_das[i, 0:nstr-idt] = trace_data[idt:]

            
            if derivative:
                # calculate derivatives of input data
                tr_das[i,1:self.ns_das]=((tr_das[i,1:]-tr_das[i,0:self.ns_das-1])/self.deltat_das)
                tr_das[i,0]=0.
            
            if isinstance(normalize, float):
                # normalize data only if the absolute maxima is largar than a 
                # certain input threshold
                trmax = num.max(num.abs(tr_das[i,:]))
                if trmax >= normalize:
                    tr_das[i,:]=tr_das[i,:]/trmax
            elif normalize:
                # normalize data by the absolute data maxima
                if num.max(num.abs(tr_das[i,:])) > 0:
                    tr_das[i,:]=tr_das[i,:]/num.max(num.abs(tr_das[i,:]))
                    

        return tr_das

    def analytic_signal_sta(self, trace):
        tracef=num.fft.fft(trace)
        nsta,nfreq=num.shape(tracef)
        freqs=num.fft.fftfreq(nfreq,self.deltat_sta)
        traceh=-1j*num.sign(freqs).T*tracef
        trace=trace+1j*num.fft.ifft(traceh).real
        return trace

    def analytic_signal_das(self, trace):
        tracef=num.fft.fft(trace)
        nsta,nfreq=num.shape(tracef)
        freqs=num.fft.fftfreq(nfreq,self.deltat_das)
        traceh=-1j*num.sign(freqs).T*tracef
        trace=trace+1j*num.fft.ifft(traceh).real
        return trace
    

    def time_extractor_sta(self, tp, ts):
        nxyz= num.size(tp[self.stations[0]])
        tp_mod_sta=num.zeros([nxyz,self.nstation])
        ts_mod_sta=num.zeros([nxyz,self.nstation])
        for i,sta in enumerate(self.stations):
            tp_mod_sta[:,i]=tp[sta]
            ts_mod_sta[:,i]=ts[sta]
        return (tp_mod_sta, ts_mod_sta)
    
    def time_extractor_das(self, tp, ts):
        nxyz= num.size(tp[self.channels[0]])
        tp_mod_das=num.zeros([nxyz,self.nstation])
        ts_mod_das=num.zeros([nxyz,self.nstation])
        for i,sta in enumerate(self.channels):
            tp_mod_das[:,i]=tp[sta]
            ts_mod_das[:,i]=ts[sta]
        return (tp_mod_das, ts_mod_das)


    def time_extractor_sta(self, tp, ts):
        nxyz= num.size(tp[self.stations[0]])
        tp_mod=num.zeros([nxyz,self.nstation])
        ts_mod=num.zeros([nxyz,self.nstation])
        for i,sta in enumerate(self.stations):
            tp_mod[:,i]=tp[sta]
            ts_mod[:,i]=ts[sta]
        return (tp_mod, ts_mod)
    

    def time_extractor_das(self, tp, ts):
        nxyz= num.size(tp[self.channels[0]])
        tp_mod=num.zeros([nxyz,self.nchannel])
        ts_mod=num.zeros([nxyz,self.nchannel])
        for i,sta in enumerate(self.channels):
            tp_mod[:,i]=tp[sta]
            ts_mod[:,i]=ts[sta]
        return (tp_mod, ts_mod)

    def characteristic_function(self, vfunc='erg', hfunc='pca', epsilon=0.001):
        if vfunc=='erg' and hfunc=='pca':
            self.cfunc_erg_sta(True)
            self.cfunc_erg_das(True)
            self.cfunc_pca_sta(epsilon)
            self.cfunc_pca_das(epsilon)
        elif vfunc=='pca' and hfunc=='pca':
            self.cfunc_pcafull_sta(epsilon)
            self.cfunc_pcafull_das(epsilon)
        elif vfunc=='erg' and hfunc=='erg':
            self.cfunc_erg_sta(False)
            self.cfunc_erg_das(False)
        elif vfunc=='erg' and hfunc=='null':
            self.cfunc_erg_sta(True)
            self.cfunc_erg_das(True)
        else:
            print('wrong characterstic functions, energy used as default')
            self.cfunc_erg_sta(False)
            self.cfunc_erg_das(False)

    def cfunc_erg_sta(self, ergz):
        if ergz:
            obs_dataV_sta=(self.ztr_sta**2)
            for i in range(self.nstation):
                if num.max(obs_dataV_sta[i,:]) > 0:
                    obs_dataV_sta[i,:]=(obs_dataV_sta[i,:]/num.max(obs_dataV_sta[i,:]))
            self.obs_dataV_sta=obs_dataV_sta
        else:
            obs_dataV_sta=(self.ztr_sta**2)
            obs_dataH_sta=(self.xtr_sta**2)+(self.ytr_sta**2)
            for i in range(self.nstation):
                if abs(num.max(obs_dataH_sta[i,:])) > 0:
                    obs_dataH_sta[i,:]=(obs_dataH_sta[i,:]/num.max(obs_dataH_sta[i,:]))
                if abs(num.max(obs_dataV_sta[i,:])) > 0:
                    obs_dataV_sta[i,:]=(obs_dataV_sta[i,:]/num.max(obs_dataV_sta[i,:]))
            self.obs_dataH_sta=obs_dataH_sta
            self.obs_dataV_sta=obs_dataV_sta

    def cfunc_erg_das(self, ergz):
        if ergz:
            obs_dataV_das=(self.xtr_das**2)
            for i in range(self.nchannel):
                if num.max(obs_dataV_das[i,:]) > 0:
                    obs_dataV_das[i,:]=(obs_dataV_das[i,:]/num.max(obs_dataV_das[i,:]))
            self.obs_dataV_das=obs_dataV_das
        else:
            obs_dataV_das=(self.xtr_das**2)
            obs_dataH_das=(self.xtr_das**2)+(self.xtr_das**2)
            for i in range(self.nchannel):
                if abs(num.max(obs_dataH_das[i,:])) > 0:
                    obs_dataH_das[i,:]=(obs_dataH_das[i,:]/num.max(obs_dataH_das[i,:]))
                if abs(num.max(obs_dataV_das[i,:])) > 0:
                    obs_dataV_das[i,:]=(obs_dataV_das[i,:]/num.max(obs_dataV_das[i,:]))
            self.obs_dataH_das=obs_dataH_das
            self.obs_dataV_das=obs_dataV_das

    def cfunc_pcafull_sta(self, epsilon):
        obs_dataH_sta=num.zeros([self.nstation,self.ns_sta]); obs_dataV_sta=num.zeros([self.nstation,self.ns_sta])
        obs_dataH1_sta=self.analytic_signal(self.xtr_sta); obs_dataH2_sta=self.analytic_signal(self.ytr_sta); obs_dataH3_sta=self.analytic_signal(self.ztr_sta)
        obs_dataH1C_sta=num.conjugate(obs_dataH1_sta); obs_dataH2C_sta=num.conjugate(obs_dataH2_sta); obs_dataH3C_sta=num.conjugate(obs_dataH3_sta)
        xx_sta=obs_dataH1_sta*obs_dataH1C_sta; xy_sta=obs_dataH1_sta*obs_dataH2C_sta; xz_sta=obs_dataH1_sta*obs_dataH3C_sta
        yx_sta=obs_dataH2_sta*obs_dataH1C_sta; yy_sta=obs_dataH2_sta*obs_dataH2C_sta; yz_sta=obs_dataH2_sta*obs_dataH2C_sta
        zx_sta=obs_dataH3_sta*obs_dataH1C_sta; zy_sta=obs_dataH3_sta*obs_dataH2C_sta; zz_sta=obs_dataH3_sta*obs_dataH3C_sta
        for i in range(self.nstation):
            for j in range(self.ns_sta):
                cov3d_sta=num.array([[xx_sta[i,j], xy_sta[i,j], xz_sta[i,j]],[yx_sta[i,j], yy_sta[i,j], yz_sta[i,j]],[zx_sta[i,j],zy_sta[i,j],zz_sta[i,j]]])
                cov2d_sta=num.array([[xx_sta[i,j], xy_sta[i,j]],[yx_sta[i,j], yy_sta[i,j]]])
                U2d_sta, s2d_sta, V2d_sta = num.linalg.svd(cov2d_sta, full_matrices=True)
                U3d_sta, s3d_sta, V3d_sta = num.linalg.svd(cov3d_sta, full_matrices=True)
                obs_dataV_sta[i,j]=(s3d_sta[0]**2)*(num.abs(V3d_sta[0][2]))
                obs_dataH_sta[i,j]=(s2d_sta[0]**2)*(1-num.abs(V3d_sta[0][2]))
                
            if abs(num.max(obs_dataH_sta[i,:])) > 0:
                obs_dataH_sta[i,:]=(obs_dataH_sta[i,:]/num.max(obs_dataH_sta[i,:]))+epsilon
            if abs(num.max(obs_dataV_sta[i,:])) > 0:
                obs_dataV_sta[i,:]=(obs_dataV_sta[i,:]/num.max(obs_dataV_sta[i,:]))
        self.obs_dataH_sta=obs_dataH_sta
        self.obs_dataV_sta=obs_dataV_sta


    def cfunc_pcafull_das(self, epsilon):
        obs_dataH_das=num.zeros([self.nstation,self.ns_das]); obs_dataV_das=num.zeros([self.nstation,self.ns_das])
        obs_dataH1_das=self.analytic_signal(self.xtr_das); obs_dataH2_das=self.analytic_signal(self.ytr_das); obs_dataH3_das=self.analytic_signal(self.ztr_das)
        obs_dataH1C_das=num.conjugate(obs_dataH1_das); obs_dataH2C_das=num.conjugate(obs_dataH2_das); obs_dataH3C_das=num.conjugate(obs_dataH3_das)
        xx_das=obs_dataH1_das*obs_dataH1C_das; xy_das=obs_dataH1_das*obs_dataH2C_das; xz_das=obs_dataH1_das*obs_dataH3C_das
        yx_das=obs_dataH2_das*obs_dataH1C_das; yy_das=obs_dataH2_das*obs_dataH2C_das; yz_das=obs_dataH2_das*obs_dataH2C_das
        zx_das=obs_dataH3_das*obs_dataH1C_das; zy_das=obs_dataH3_das*obs_dataH2C_das; zz_das=obs_dataH3_das*obs_dataH3C_das
        for i in range(self.nstation):
            for j in range(self.ns_das):
                cov3d_das=num.array([[xx_das[i,j], xy_das[i,j], xz_das[i,j]],[yx_das[i,j], yy_das[i,j], yz_das[i,j]],[zx_das[i,j],zy_das[i,j],zz_das[i,j]]])
                cov2d_das=num.array([[xx_das[i,j], xy_das[i,j]],[yx_das[i,j], yy_das[i,j]]])
                U2d_das, s2d_das, V2d_das = num.linalg.svd(cov2d_das, full_matrices=True)
                U3d_das, s3d_das, V3d_das = num.linalg.svd(cov3d_das, full_matrices=True)
                obs_dataV_das[i,j]=(s3d_das[0]**2)*(num.abs(V3d_das[0][2]))
                obs_dataH_das[i,j]=(s2d_das[0]**2)*(1-num.abs(V3d_das[0][2]))
                
            if abs(num.max(obs_dataH_das[i,:])) > 0:
                obs_dataH_das[i,:]=(obs_dataH_das[i,:]/num.max(obs_dataH_das[i,:]))+epsilon
            if abs(num.max(obs_dataV_das[i,:])) > 0:
                obs_dataV_das[i,:]=(obs_dataV_das[i,:]/num.max(obs_dataV_das[i,:]))
        self.obs_dataH_das=obs_dataH_das
        self.obs_dataV_das=obs_dataV_das


    def cfunc_pca_sta(self, epsilon):
        obs_dataH_sta=num.zeros([self.nstation,self.ns_sta])
        obs_dataH1_sta=self.analytic_signal_sta(self.xtr_sta)
        obs_dataH2_sta=self.analytic_signal_sta(self.ytr_sta)
        obs_dataH1C_sta=num.conjugate(obs_dataH1_sta)
        obs_dataH2C_sta=num.conjugate(obs_dataH2_sta)
        xx_sta=obs_dataH1_sta*obs_dataH1C_sta; xy_sta=obs_dataH1_sta*obs_dataH2C_sta
        yx_sta=obs_dataH2_sta*obs_dataH1C_sta; yy_sta=obs_dataH2_sta*obs_dataH2C_sta
        for i in range(self.nstation):
            for j in range(self.ns_sta):
                cov_sta=num.array([[xx_sta[i,j], xy_sta[i,j]],[yx_sta[i,j], yy_sta[i,j]]])
                U_sta, s_sta, V_sta = num.linalg.svd(cov_sta, full_matrices=True)
                obs_dataH_sta[i,j]=(s_sta[0]**2)
                
            if abs(num.max(obs_dataH_sta[i,:])) > 0:
                obs_dataH_sta[i,:]=(obs_dataH_sta[i,:]/num.max(obs_dataH_sta[i,:]))+epsilon
        self.obs_dataH_sta=obs_dataH_sta
    

    def cfunc_pca_das(self, epsilon):
        obs_dataH_das=num.zeros([self.nstation,self.ns_das])
        obs_dataH1_das=self.analytic_signal_das(self.xtr_das)
        obs_dataH2_das=self.analytic_signal_das(self.xtr_das)
        obs_dataH1C_das=num.conjugate(obs_dataH1_das)
        obs_dataH2C_das=num.conjugate(obs_dataH2_das)
        xx_das=obs_dataH1_das*obs_dataH1C_das; xy_das=obs_dataH1_das*obs_dataH2C_das
        yx_das=obs_dataH2_das*obs_dataH1C_das; yy_das=obs_dataH2_das*obs_dataH2C_das
        for i in range(self.nstation):
            for j in range(self.ns_das):
                cov_das=num.array([[xx_das[i,j], xy_das[i,j]],[yx_das[i,j], yy_das[i,j]]])
                U_das, s_das, V_das = num.linalg.svd(cov_das, full_matrices=True)
                obs_dataH_das[i,j]=(s_das[0]**2)
                
            if abs(num.max(obs_dataH_das[i,:])) > 0:
                obs_dataH_das[i,:]=(obs_dataH_das[i,:]/num.max(obs_dataH_das[i,:]))+epsilon
        self.obs_dataH_das=obs_dataH_das


    # def pca_win(xtr, ytr, ztr, iwin):
    #     for i in range(self.nstation):
    #          for j in range(iwin,self.ns-iwin):
    #              X=self.xtr[i,j-iwin:j+iwin]-num.mean(self.xtr[i,j-iwin:j+iwin])
    #              Y=self.ytr[i,j-iwin:j+iwin]-num.mean(self.ytr[i,j-iwin:j+iwin])
    #              Z=self.ztr[i,j-iwin:j+iwin]-num.mean(self.ztr[i,j-iwin:j+iwin])
    #              cov=num.vstack((X,Y,Z))
    #              C=num.dot(cov,cov.T)
    #              U, s, V = num.linalg.svd(C, full_matrices=True)
    #              obs_dataH[i,j]=1-((s[1]+s[2])/(2*s[0])))
    #              azm.append(num.arctan2(V[0][1],V[0][0])*(180/num.pi))
    #              inc.append(num.arccos(num.abs(V[0][2]))*(180/num.pi))
    #           obs_dataH[i,j]
    #     dol=num.array(dol); azm=num.array(azm); inc=num.array(inc)
    #     azm[azm<0]=360+azm[azm<0]
    # return dol,azm,inc


    def loc_stalta_sta(self, nshort_p, nshort_s, slrat, norm=1):
        tshort_p=nshort_p*self.deltat_sta; tshort_s=nshort_s*self.deltat_sta
        tlong_p=tshort_p*slrat; tlong_s=tshort_s*slrat
        ks_p=self.deltat_sta/tshort_p; kl_p=self.deltat_sta/tlong_p;
        ks_s=self.deltat_sta/tshort_s; kl_s=self.deltat_sta/tlong_s;
        obs_dataP=LOC_STALTA.recursive_stalta(tshort_p, tlong_p, self.deltat_sta, self.obs_dataV_sta, kl_p, ks_p, norm)
        obs_dataS=LOC_STALTA.recursive_stalta(tshort_s, tlong_s, self.deltat_sta, self.obs_dataH_sta, kl_s, ks_s, norm)
        return obs_dataP, obs_dataS

    def loc_stalta_das(self, nshort_p, nshort_s, slrat, norm=1):
        tshort_p=nshort_p*self.deltat_das; tshort_s=nshort_s*self.deltat_das
        tlong_p=tshort_p*slrat; tlong_s=tshort_s*slrat
        ks_p=self.deltat_das/tshort_p; kl_p=self.deltat_das/tlong_p;
        ks_s=self.deltat_das/tshort_s; kl_s=self.deltat_das/tlong_s;
        obs_dataP=LOC_STALTA.recursive_stalta(tshort_p, tlong_p, self.deltat_das, self.obs_dataV_das, kl_p, ks_p, norm)
        obs_dataS=LOC_STALTA.recursive_stalta(tshort_s, tlong_s, self.deltat_das, self.obs_dataH_das, kl_s, ks_s, norm)
        return obs_dataP, obs_dataS

    def det_stalta_sta(self, nshort_p, nshort_s, slrat, staltap0, staltas0, thres=0.):
        tshort_p=nshort_p*self.deltat_sta; tshort_s=nshort_s*self.deltat_sta
        tlong_p=tshort_p*slrat; tlong_s=tshort_s*slrat
        obs_dataP=num.zeros([self.nstation,self.ns_sta])
        obs_dataS=num.zeros([self.nstation,self.ns_sta])
        for i in range(self.nstation):
            obs_dataP[i,:],stltp0=DET_STALTA.recursive_stalta(staltap0[i][0], staltap0[i][1], tshort_p, tlong_p, self.deltat_sta, self.obs_dataV_sta[i,:], thres)
            obs_dataS[i,:],stlts0=DET_STALTA.recursive_stalta(staltas0[i][0], staltas0[i][1], tshort_s, tlong_s, self.deltat_sta, self.obs_dataH_sta[i,:], thres)
            staltap0[i][0]=stltp0[0]; staltap0[i][1]=stltp0[1]
            staltas0[i][0]=stlts0[0]; staltas0[i][1]=stlts0[1]
        return obs_dataP, obs_dataS, staltap0, staltas0
    
    def det_stalta_das(self, nshort_p, nshort_s, slrat, staltap0, staltas0, thres=0.):
        tshort_p=nshort_p*self.deltat; tshort_s=nshort_s*self.deltat
        tlong_p=tshort_p*slrat; tlong_s=tshort_s*slrat
        obs_dataP=num.zeros([self.nstation,self.ns_das])
        obs_dataS=num.zeros([self.nstation,self.ns_das])
        for i in range(self.nstation):
            obs_dataP[i,:],stltp0=DET_STALTA.recursive_stalta(staltap0[i][0], staltap0[i][1], tshort_p, tlong_p, self.deltat_das, self.obs_dataV_das[i,:], thres)
            obs_dataS[i,:],stlts0=DET_STALTA.recursive_stalta(staltas0[i][0], staltas0[i][1], tshort_s, tlong_s, self.deltat_das, self.obs_dataH_das[i,:], thres)
            staltap0[i][0]=stltp0[0]; staltap0[i][1]=stltp0[1]
            staltas0[i][0]=stlts0[0]; staltas0[i][1]=stlts0[1]
        return obs_dataP, obs_dataS, staltap0, staltas0


'''

import numpy as num
import DET_STALTA     # C
import LOC_STALTA     # C
from obspy import UTCDateTime

class Stacktraces:

    def __init__(self, tobj, wobj, **inputs):
        # Initialize time-related parameters
        
        self.dtime_max_sta = None  # Initialize to None before calculating
        self.dtime_max_das = None  # Initialize to None before calculating

        self.deltat = tobj.dtt

        # Initialize the station and DAS samples (ns)
        self.ns_sta = None
        self.ns_das = None

        # Process station data if available
        if hasattr(wobj, 'stream_sta') and wobj.stream_sta:
            first_component = next(iter(wobj.stream_sta.keys()))
            first_station = next(iter(wobj.stream_sta[first_component].keys()))
            self.ns_sta = len(wobj.stream_sta[first_component][first_station][2])
        else:
            raise ValueError("No station waveform data found in wobj")

        # Process DAS data if available
        if hasattr(wobj, 'stream_das') and wobj.stream_das:
            self.ns_das = len(wobj.stream_das[0].data)
        else:
            raise ValueError("No DAS waveform data found in wobj")

        if not self.ns_sta and not self.ns_das:
            raise ValueError("Cannot determine 'ns' from the waveform data. Neither 'sta' nor 'das' data found.")

        # Call separate input handling methods for sta and das
        self.loki_input_sta(wobj, tobj, inputs.get('derivative', False), 
                            inputs.get('direct_input', False), 
                            inputs.get('normalize', True))
        
        self.loki_input_das(wobj, tobj, inputs.get('derivative', False), 
                            inputs.get('direct_input', False), 
                            inputs.get('normalize', True))

        # Check sampling rates for both sta and das
        self.check_sampling_rate_sta(wobj)
        self.check_sampling_rate_das(wobj)

        # Ensure dtime_max_sta is set correctly by calling check_starting_time
        self.check_starting_time(wobj)  # Make sure this is called to set `dtime_max_sta`
        
        # Process characteristic functions if provided
        if ('vfunc' in inputs) or ('hfunc' in inputs):
            vfunc = inputs['vfunc']
            hfunc = inputs['hfunc']
            epsilon = inputs['epsilon']
            derivative = inputs['derivative']
            
            # Reinitialize data if characteristic functions need to be calculated
            self.loki_input_sta(wobj, tobj, derivative)
            self.loki_input_das(wobj, tobj, derivative)
            self.characteristic_function(vfunc, hfunc, epsilon)
        else:
            # Directly input data of P and S components
            normalize = inputs.get('normthrd', False)
            self.loki_input_sta(wobj, tobj, derivative=False, direct_input=True, normalize=normalize)
            self.loki_input_das(wobj, tobj, derivative=False, direct_input=True, normalize=normalize)
            
            # Compute power if needed for both `sta` and `das`
            if 'ppower' in inputs:
                self.obs_dataV = self.obs_dataV ** inputs['ppower']  # Vertical
                self.obs_dataH = self.obs_dataH ** inputs['ppower']  # Horizontal

            # Additional processing for DAS
            if hasattr(self, 'obs_data_das'):
                self.obs_data_das = self.apply_hilbert_transform(self.obs_data_das)


    def check_sampling_rate_sta(self, wobj):
        intsamp = 1E6
        deltas = []
        if hasattr(wobj, 'stream_sta'):
            for comp in wobj.stream_sta.keys():
                for sta in wobj.stream_sta[comp].keys():
                    deltas.append(wobj.stream_sta[comp][sta][1])  # Sampling rate
        deltas = num.array(deltas)
        ideltas = num.unique((deltas * intsamp).astype(int))
        if num.size(ideltas) == 1:
            self.deltat = deltas[0]
        else:
            raise ValueError('Error (stations)!! All traces must have the same sampling rate')

    def check_sampling_rate_das(self, wobj):
        intsamp = 1E6
        deltas = []
        if hasattr(wobj, 'stream_das'):
            for trace in wobj.stream_das:
                deltas.append(trace.stats.delta)  # Sampling interval (delta)
        deltas = num.array(deltas)
        ideltas = num.unique((deltas * intsamp).astype(int))
        if num.size(ideltas) == 1:
            self.deltat = deltas[0]
        else:
            raise ValueError('Error (DAS) !! All traces must have the same sampling rate')

    def check_starting_time(self, wobj):
        # Initialize an empty list to store time values for DAS and STA separately
        dtimes_das = []
        dtimes_sta = []

        # Check for DAS data and extract the times
        for comp in wobj.stream_das.keys():
            for sta in wobj.stream_das[comp].keys():
                dtimes_das.append(wobj.stream_das[comp][sta][0])  # Time is at index 0
                # Update the size for DAS if needed
                if self.ns_das < num.size(wobj.stream_das[comp][sta][2]):
                    self.ns_das = num.size(wobj.stream_das[comp][sta][2])

        # Check for STA data and extract the times
        for comp in wobj.stream_sta.keys():
            for sta in wobj.stream_sta[comp].keys():
                dtimes_sta.append(wobj.stream_sta[comp][sta][0])  # Time is at index 0
                # Update the size for STA if needed
                if self.ns_sta < num.size(wobj.stream_sta[comp][sta][2]):
                    self.ns_sta = num.size(wobj.stream_sta[comp][sta][2])

        # Calculate dtime_max for DAS and STA independently
        self.dtime_max_das = max(dtimes_das) if dtimes_das else None
        self.dtime_max_sta = max(dtimes_sta) if dtimes_sta else None

        # Set evid using ISO format of the max times for DAS and STA
        if self.dtime_max_das:
            self.evid_das = self.dtime_max_das.isoformat()
        else:
            self.evid_das = None

        if self.dtime_max_sta:
            self.evid_sta = self.dtime_max_sta.isoformat()
        else:
            self.evid_sta = None

        # Ensure dtime_max_sta is initialized
        if self.dtime_max_sta is None:
            self.dtime_max_sta = UTCDateTime(0)  # Or set to a default value if needed


    def loki_input_sta(self, wobj, tobj, derivative, direct_input=False, normalize=True):
        """
        Process station data only.
        """
        if direct_input:
            # Handle direct input for station data
            self.obs_dataV_sta = self.select_data('P', wobj, tobj.db_stations, derivative, normalize, 'station')
            self.obs_dataH_sta = self.select_data('S', wobj, tobj.db_stations, derivative, normalize, 'station')
        else:
            # Handle stream data for stations
            station_components = []

            if hasattr(wobj, 'stream_sta') and wobj.stream_sta:
                station_components = list(wobj.stream_sta.keys())  # Only consider station components

            self.comp_sta = tuple(set(station_components))  # Only store station components

            if len(self.comp_sta) == 3:
                self.xtr_sta = self.select_data(self.comp_sta[0], wobj, tobj.db_stations, derivative, normalize, 'station')
                self.ytr_sta = self.select_data(self.comp_sta[1], wobj, tobj.db_stations, derivative, normalize, 'station')
                self.ztr_sta = self.select_data(self.comp_sta[2], wobj, tobj.db_stations, derivative, normalize, 'station')
            elif len(self.comp_sta) == 1:
                self.ztr_sta = self.select_data(self.comp_sta[0], wobj, tobj.db_stations, derivative, normalize, 'station')
            else:
                raise ValueError('Station traces must have 3 components!')

    def loki_input_das(self, wobj, tobj, derivative, direct_input=False, normalize=True):
        """
        Process DAS data only.
        """
        if direct_input:
            # Handle direct input for DAS data
            self.obs_dataV_das = self.select_data('P', wobj, tobj.db_stations, derivative, normalize, 'das')
            self.obs_dataH_das = self.select_data('S', wobj, tobj.db_stations, derivative, normalize, 'das')
        else:
            # Handle stream data for DAS
            das_components = []

            if hasattr(wobj, 'stream_das') and wobj.stream_das:
                das_components = list(set(trace.stats.channel for trace in wobj.stream_das))  # Only consider DAS components

            self.comp_das = tuple(set(das_components))  # Only store DAS components

            if len(self.comp_das) == 3:
                self.xtr_das = self.select_data(self.comp_das[0], wobj, tobj.db_stations, derivative, normalize, 'das')
                self.ytr_das = self.select_data(self.comp_das[1], wobj, tobj.db_stations, derivative, normalize, 'das')
                self.ztr_das = self.select_data(self.comp_das[2], wobj, tobj.db_stations, derivative, normalize, 'das')
            elif len(self.comp_das) == 1:
                self.ztr_das = self.select_data(self.comp_das[0], wobj, tobj.db_stations, derivative, normalize, 'das')
            else:
                raise ValueError('DAS traces must have 3 components!')
            

    def select_data(self, comp, wobj, db_stations, derivative, normalize, data_type):
    # Ensure dtime_max is not None before proceeding
        if self.dtime_max_sta is None and self.dtime_max_das is None:
            raise ValueError("Both dtime_max_sta and dtime_max_das are None, unable to proceed.")

        # Choose the appropriate max time based on data type
        if data_type == 'das':
            dtime_max = self.dtime_max_das
        else:  # For 'sta'
            dtime_max = self.dtime_max_sta

        # Proceed if dtime_max is valid
        if dtime_max is not None:
            # Compute time index safely
            idt = int((dtime_max - trace.stats.starttime).total_seconds() / self.deltat)
        else:
            raise ValueError(f"{data_type} dtime_max is None, unable to compute time index.")
        # Process station data
        if 'station' in data_type:
            self.stations_sta = tuple(wobj.data_stations & db_stations)  # Find common stations between wobj.data_stations and db_stations
            stream_sta = wobj.stream_sta[comp]  # Get the stream of station data for the specified component
            ns_sta = self.ns_sta  # Use the number of samples specific to the station data
            dtime_max = self.dtime_max_sta  # Use the station maximum time for station data
        else:
            self.stations_sta = []  # No station data to process

        # Process DAS data
        if 'das' in data_type:
            self.stations_das = tuple(wobj.data_channels)  # Use DAS stations from the stream_das
            stream_das = wobj.stream_das  # Get the stream of DAS data
            ns_das = self.ns_das  # Use the number of samples specific to the DAS data
            dtime_max = self.dtime_max_das  # Use the DAS maximum time for DAS data
        else:
            self.stations_das = []  # No DAS data to process

        # Initialize the result array for both station and DAS data
        tr_sta = num.zeros([len(self.stations_sta), ns_sta]) if self.stations_sta else None
        tr_das = num.zeros([len(self.stations_das), ns_das]) if self.stations_das else None

        # Process station data if present
        if tr_sta is not None:
            for i, sta in enumerate(self.stations_sta):
                nstr = len(stream_sta[sta][2])  # Number of samples for this station
                idt = int((dtime_max - stream_sta[sta][0]).total_seconds() / self.deltat)  # Time index
                tr_sta[i, 0:nstr - idt] = stream_sta[sta][2][idt:]  # Store data, adjusting for time offset

                # Apply derivative if specified
                if derivative:
                    tr_sta[i, 1:ns_sta] = (tr_sta[i, 1:] - tr_sta[i, 0:ns_sta-1]) / self.deltat
                    tr_sta[i, 0] = 0.  # Initial point is 0 after derivative calculation

                # Normalize the data
                if isinstance(normalize, float):
                    # Normalize based on a threshold value
                    trmax = num.max(num.abs(tr_sta[i, :]))
                    if trmax >= normalize:
                        tr_sta[i, :] = tr_sta[i, :] / trmax
                elif normalize:
                    # Normalize data by the maximum value
                    if num.max(num.abs(tr_sta[i, :])) > 0:
                        tr_sta[i, :] = tr_sta[i, :] / num.max(num.abs(tr_sta[i, :]))

        # Process DAS data if present
        if tr_das is not None:
            for i, sta in enumerate(self.stations_das):
                trace = [tr for tr in stream_das if tr.stats.station == sta][0]  # Find the corresponding DAS trace by station name
                nstr = len(trace.data)  # Number of samples for this DAS trace
                idt = int((dtime_max - trace.stats.starttime).total_seconds() / self.deltat)  # Time index
                tr_das[i, 0:nstr - idt] = trace.data[idt:]  # Store the DAS data, adjusting for time offset

                # Apply derivative if specified
                if derivative:
                    tr_das[i, 1:ns_das] = (tr_das[i, 1:] - tr_das[i, 0:ns_das-1]) / self.deltat
                    tr_das[i, 0] = 0.  # Initial point is 0 after derivative calculation

                # Normalize the data
                if isinstance(normalize, float):
                    # Normalize based on a threshold value
                    trmax = num.max(num.abs(tr_das[i, :]))
                    if trmax >= normalize:
                        tr_das[i, :] = tr_das[i, :] / trmax
                elif normalize:
                    # Normalize data by the maximum value
                    if num.max(num.abs(tr_das[i, :])) > 0:
                        tr_das[i, :] = tr_das[i, :] / num.max(num.abs(tr_das[i, :]))

        # Return processed station and DAS data separately
        return tr_sta, tr_das


    def analytic_signal(self, trace):
        tracef=num.fft.fft(trace)
        nsta,nfreq=num.shape(tracef)
        freqs=num.fft.fftfreq(nfreq,self.deltat)
        traceh=-1j*num.sign(freqs).T*tracef
        trace=trace+1j*num.fft.ifft(traceh).real
        return trace


    def time_extractor(self, tp, ts):
        nxyz= num.size(tp[self.stations[0]])
        tp_mod=num.zeros([nxyz,self.nstation])
        ts_mod=num.zeros([nxyz,self.nstation])
        for i,sta in enumerate(self.stations):
            tp_mod[:,i]=tp[sta]
            ts_mod[:,i]=ts[sta]
        return (tp_mod, ts_mod)


    def characteristic_function(self, vfunc='erg', hfunc='pca', epsilon=0.001):
        if vfunc=='erg' and hfunc=='pca':
            self.cfunc_erg(True)
            self.cfunc_pca(epsilon)
        elif vfunc=='pca' and hfunc=='pca':
            self.cfunc_pcafull(epsilon)
        elif vfunc=='erg' and hfunc=='erg':
            self.cfunc_erg(False)
        elif vfunc=='erg' and hfunc=='null':
            self.cfunc_erg(True)
        else:
            print('wrong characterstic functions, energy used as default')
            self.cfunc_erg(False)

    def cfunc_erg(self, ergz):
        if ergz:
            obs_dataV=(self.ztr**2)
            for i in range(self.nstation):
                if num.max(obs_dataV[i,:]) > 0:
                    obs_dataV[i,:]=(obs_dataV[i,:]/num.max(obs_dataV[i,:]))
            self.obs_dataV=obs_dataV
        else:
            obs_dataV=(self.ztr**2)
            obs_dataH=(self.xtr**2)+(self.ytr**2)
            for i in range(self.nstation):
                if abs(num.max(obs_dataH[i,:])) > 0:
                    obs_dataH[i,:]=(obs_dataH[i,:]/num.max(obs_dataH[i,:]))
                if abs(num.max(obs_dataV[i,:])) > 0:
                    obs_dataV[i,:]=(obs_dataV[i,:]/num.max(obs_dataV[i,:]))
            self.obs_dataH=obs_dataH
            self.obs_dataV=obs_dataV

    def cfunc_pcafull(self, epsilon):
        obs_dataH=num.zeros([self.nstation,self.ns]); obs_dataV=num.zeros([self.nstation,self.ns])
        obs_dataH1=self.analytic_signal(self.xtr); obs_dataH2=self.analytic_signal(self.ytr); obs_dataH3=self.analytic_signal(self.ztr)
        obs_dataH1C=num.conjugate(obs_dataH1); obs_dataH2C=num.conjugate(obs_dataH2); obs_dataH3C=num.conjugate(obs_dataH3)
        xx=obs_dataH1*obs_dataH1C; xy=obs_dataH1*obs_dataH2C; xz=obs_dataH1*obs_dataH3C
        yx=obs_dataH2*obs_dataH1C; yy=obs_dataH2*obs_dataH2C; yz=obs_dataH2*obs_dataH2C
        zx=obs_dataH3*obs_dataH1C; zy=obs_dataH3*obs_dataH2C; zz=obs_dataH3*obs_dataH3C
        for i in range(self.nstation):
            for j in range(self.ns):
                cov3d=num.array([[xx[i,j], xy[i,j], xz[i,j]],[yx[i,j], yy[i,j], yz[i,j]],[zx[i,j],zy[i,j],zz[i,j]]])
                cov2d=num.array([[xx[i,j], xy[i,j]],[yx[i,j], yy[i,j]]])
                U2d, s2d, V2d = num.linalg.svd(cov2d, full_matrices=True)
                U3d, s3d, V3d = num.linalg.svd(cov3d, full_matrices=True)
                obs_dataV[i,j]=(s3d[0]**2)*(num.abs(V3d[0][2]))
                obs_dataH[i,j]=(s2d[0]**2)*(1-num.abs(V3d[0][2]))
                
            if abs(num.max(obs_dataH[i,:])) > 0:
                obs_dataH[i,:]=(obs_dataH[i,:]/num.max(obs_dataH[i,:]))+epsilon
            if abs(num.max(obs_dataV[i,:])) > 0:
                obs_dataV[i,:]=(obs_dataV[i,:]/num.max(obs_dataV[i,:]))
        self.obs_dataH=obs_dataH
        self.obs_dataV=obs_dataV


    def cfunc_pca(self, epsilon):
        obs_dataH=num.zeros([self.nstation,self.ns])
        obs_dataH1=self.analytic_signal(self.xtr)
        obs_dataH2=self.analytic_signal(self.ytr)
        obs_dataH1C=num.conjugate(obs_dataH1)
        obs_dataH2C=num.conjugate(obs_dataH2)
        xx=obs_dataH1*obs_dataH1C; xy=obs_dataH1*obs_dataH2C
        yx=obs_dataH2*obs_dataH1C; yy=obs_dataH2*obs_dataH2C
        for i in range(self.nstation):
            for j in range(self.ns):
                cov=num.array([[xx[i,j], xy[i,j]],[yx[i,j], yy[i,j]]])
                U, s, V = num.linalg.svd(cov, full_matrices=True)
                obs_dataH[i,j]=(s[0]**2)
                
            if abs(num.max(obs_dataH[i,:])) > 0:
                obs_dataH[i,:]=(obs_dataH[i,:]/num.max(obs_dataH[i,:]))+epsilon
        self.obs_dataH=obs_dataH


    # def pca_win(xtr, ytr, ztr, iwin):
    #     for i in range(self.nstation):
    #          for j in range(iwin,self.ns-iwin):
    #              X=self.xtr[i,j-iwin:j+iwin]-num.mean(self.xtr[i,j-iwin:j+iwin])
    #              Y=self.ytr[i,j-iwin:j+iwin]-num.mean(self.ytr[i,j-iwin:j+iwin])
    #              Z=self.ztr[i,j-iwin:j+iwin]-num.mean(self.ztr[i,j-iwin:j+iwin])
    #              cov=num.vstack((X,Y,Z))
    #              C=num.dot(cov,cov.T)
    #              U, s, V = num.linalg.svd(C, full_matrices=True)
    #              obs_dataH[i,j]=1-((s[1]+s[2])/(2*s[0])))
    #              azm.append(num.arctan2(V[0][1],V[0][0])*(180/num.pi))
    #              inc.append(num.arccos(num.abs(V[0][2]))*(180/num.pi))
    #           obs_dataH[i,j]
    #     dol=num.array(dol); azm=num.array(azm); inc=num.array(inc)
    #     azm[azm<0]=360+azm[azm<0]
    # return dol,azm,inc


    def loc_stalta(self, nshort_p, nshort_s, slrat, norm=1):
        tshort_p=nshort_p*self.deltat; tshort_s=nshort_s*self.deltat
        tlong_p=tshort_p*slrat; tlong_s=tshort_s*slrat
        ks_p=self.deltat/tshort_p; kl_p=self.deltat/tlong_p;
        ks_s=self.deltat/tshort_s; kl_s=self.deltat/tlong_s;
        obs_dataP=LOC_STALTA.recursive_stalta(tshort_p, tlong_p, self.deltat, self.obs_dataV, kl_p, ks_p, norm)
        obs_dataS=LOC_STALTA.recursive_stalta(tshort_s, tlong_s, self.deltat, self.obs_dataH, kl_s, ks_s, norm)
        return obs_dataP, obs_dataS


    def det_stalta(self, nshort_p, nshort_s, slrat, staltap0, staltas0, thres=0.):
        tshort_p=nshort_p*self.deltat; tshort_s=nshort_s*self.deltat
        tlong_p=tshort_p*slrat; tlong_s=tshort_s*slrat
        obs_dataP=num.zeros([self.nstation,self.ns])
        obs_dataS=num.zeros([self.nstation,self.ns])
        for i in range(self.nstation):
            obs_dataP[i,:],stltp0=DET_STALTA.recursive_stalta(staltap0[i][0], staltap0[i][1], tshort_p, tlong_p, self.deltat, self.obs_dataV[i,:], thres)
            obs_dataS[i,:],stlts0=DET_STALTA.recursive_stalta(staltas0[i][0], staltas0[i][1], tshort_s, tlong_s, self.deltat, self.obs_dataH[i,:], thres)
            staltap0[i][0]=stltp0[0]; staltap0[i][1]=stltp0[1]
            staltas0[i][0]=stlts0[0]; staltas0[i][1]=stlts0[1]
        return obs_dataP, obs_dataS, staltap0, staltas0

'''