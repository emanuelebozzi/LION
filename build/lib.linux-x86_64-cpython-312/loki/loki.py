# %%
import os
import math
import numpy as num
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import datetime
import copy
import gc
from loki import ioformatting
from loki import traveltimes
from loki import waveforms
from loki import stacktraces
from loki import latlon2cart
import tt_processing                       # C
import location_t0                         # C  for multiplying the P- and S-stacking values using this
#import location_t0_plus                   # C  for adding the P- and S-stacking values using this


class Loki:
    """docstring for Loki"""

    def __init__(self, data_path, output_path, db_path, hdr_filename, geometry_filename, mode='locator'):
        self.data_path = data_path
        self.output_path = output_path
        self.db_path = db_path
        self.hdr_filename = hdr_filename
        self.geometry_filename = geometry_filename
        if mode == 'locator':
            self.data_tree, self.events = self.location_data_struct(self.data_path, self.output_path)
        elif mode == 'detector':
            self.data_tree, self.events = self.detection_data_struct(self.data_path, self.output_path)
        else:
            raise ValueError('mode must be "detector" or "locator"')

    def location_data_struct(self, data_path, output_path):
        events=[]
        data_tree=[]
        for root, dirs, files in os.walk(data_path):
            if not dirs:
                data_tree.append(root)
        data_tree.sort()
        events = [idtree.split('/')[-1] for idtree in data_tree]
        if not os.path.isdir(output_path):
            os.mkdir(output_path)
        return data_tree, events

    def detection_data_struct(self, data_path, output_path):
        events = []
        data_tree = []
        for root, dirs, files in os.walk(data_path):
            if not dirs:
                data_tree.append(root)
        data_tree.sort()
        events = [idtree.split('/')[-1] for idtree in data_tree]
        if not os.path.isdir(output_path):
            os.mkdir(output_path)
        return data_tree, events

    def location(self, comp=['E', 'N', 'Z'], precision='single', **inputs):
        if 'tshortp_min' in inputs:
            # need to calculate STA/LTA for stacking
            STALTA = True
            tshortp_min = inputs['tshortp_min']
            tshortp_max = inputs['tshortp_max']
            tshorts_min = inputs['tshorts_min']
            tshorts_max = inputs['tshorts_max']
            slrat = inputs['slrat']
            ntrial = inputs['ntrial']
            #extension_sta = inputs['extension_sta']
            #extension_das = inputs['extension_das']
            #delta_das = inputs['delta_das']
            tshortp = num.linspace(tshortp_min, tshortp_max, ntrial)
            tshorts = num.linspace(tshorts_min, tshorts_max, ntrial)
        else:
            # no need to calculate STA/LTA ratio for stacking
            STALTA = False
            ntrial = 1
            
        npr = inputs['npr']
        model = inputs['model']
        if 'freq' in inputs:
            freq = inputs['freq']
        else:
            freq = None  # no filtering
        if 'opsf' in inputs:
            opsf = inputs['opsf']
        else:
            opsf = False
        
        #for each source I have P and S traveltimes for 2D velocity models

        # Synthetic Traveltimes and metadata generated with NonLinLoc and 2D 

        tobj = traveltimes.Traveltimes(self.db_path, self.hdr_filename, self.geometry_filename)

        print('The traveltime object is:', tobj)

        attributes = [name for name in dir(tobj) if not callable(getattr(tobj, name)) and not name.startswith("__")]
        methods = [name for name in dir(tobj) if callable(getattr(tobj, name)) and not name.startswith("__")]

        print("Attributes tobj:", attributes)
        print("Methods tobj:", methods)

        #load the traveltimes

        tp = tobj.load_traveltimes('P', model, precision) 
        ts = tobj.load_traveltimes('S', model, precision)

        
        #for each event 

        for event_path in self.data_tree:

            data_path = self.data_path

            #Reading the observed wavefields (stations and DAS)

            wobj = waveforms.Waveforms(tobj=tobj, data_path = data_path, event_path=event_path, extension_sta="*", extension_das='CANDAS2_2023-01-07_10-48-10.h5', freq=None)

            print('The waveforms object is:', wobj)

            attributes = [name for name in dir(wobj) if not callable(getattr(wobj, name)) and not name.startswith("__")]
            methods = [name for name in dir(wobj) if callable(getattr(wobj, name)) and not name.startswith("__")]

            print("Attributes wobj:", attributes)
            print("Methods wobj:", methods)
            print('The stream DAS is:', wobj.stream_das)
            print('The stream station is:', wobj.stream_sta)

            #object of the class stacktraces  

            sobj = stacktraces.Stacktraces(tobj, wobj, **inputs)

            print('The stacktraces object is:', sobj)

            attributes = [name for name in dir(sobj) if not callable(getattr(sobj, name)) and not name.startswith("__")]
            methods = [name for name in dir(sobj) if callable(getattr(sobj, name)) and not name.startswith("__")]

            print("Attributes sobj:", attributes)
            print("Methods sobj:", methods)


            event = event_path.split('/')[-1]

            print('Processing to the event folder: ', event_path, event)
            if os.path.isdir(self.output_path+'/'+event):
                print('directory already exists')
                #continue
            else:
                os.mkdir(self.output_path+'/'+event)

   
            tpxz=tp['HM01'].reshape(tobj.nxz, 1)
            tsxz=tp['HM01'].reshape(tobj.nxz, 1)


            tp_modse = num.ascontiguousarray(tpxz)
            ts_modse = num.ascontiguousarray(tsxz)

            print('tp_modse prima di tt_processing', tp_modse, tp_modse.shape)

            ########################################

            tp_mod_sta, ts_mod_sta = tt_processing.tt_f2i(sobj.deltat_sta, tp_modse, ts_modse, npr)  # traveltime table in time sample, for each imaging point traveltimes have substracted the minimal P traveltime
            tp_mod_das, ts_mod_das = tt_processing.tt_f2i(sobj.deltat_das, tp_modse, ts_modse, npr)  # traveltime table in time sample, for each imaging point traveltimes have substracted the minimal P traveltime


# %%


            cmax_pre = -1.0
            for i in range(ntrial):
                if STALTA:
                    # need to calculate STA/LTA from the characteristic funtion
                    # then stack the STA/LTA for imaging
                    nshort_p_sta = int(tshortp[i]//sobj.deltat_sta)
                    nshort_s_sta = int(tshorts[i]//sobj.deltat_sta)
                    nshort_p_das = int(tshortp[i]//sobj.deltat_das)
                    nshort_s_das = int(tshorts[i]//sobj.deltat_das)
                    obs_dataP_sta, obs_dataS_sta = sobj.loc_stalta_sta(nshort_p_sta, nshort_s_sta, slrat, norm=1)
                    obs_dataP_das, obs_dataS_das = sobj.loc_stalta_das(nshort_p_das, nshort_s_das, slrat, norm=1)

                else:
                    # no need to calculate STA/LTA 
                    # directly stack the characteristic function for imaging
                    obs_dataP_sta = sobj.obs_dataV_sta  # vertical -> P
                    obs_dataS_sta = sobj.obs_dataH_sta  # horizontal -> S
                    obs_dataP_das = sobj.obs_dataV_das  # vertical -> P
                    obs_dataS_das = sobj.obs_dataH_das  # horizontal -> S

                if opsf:
                    # output the characteristic functions for stacking
                    datainfo = {}
                    datainfo['dt_sta'] = sobj.deltat_sta
                    datainfo['starttime_sta'] = sobj.dtime_max_sta
                    datainfo['dt_das'] = sobj.deltat_das
                    datainfo['starttime_das'] = sobj.dtime_max_das
                    for ista, sta in enumerate(sobj.stations):
                        datainfo['station_name'] = sta
                        datainfo['channel_name'] = 'CFP'  # note maximum three characters, the last one must be 'P'
                        ioformatting.vector2trace(datainfo, obs_dataP_sta[ista,:], self.output_path+'/'+event+'/cf/trial{}'.format(i))
                        datainfo['channel_name'] = 'CFS'  # note maximum three characters, the last one must be 'S'
                        ioformatting.vector2trace(datainfo, obs_dataS_sta[ista,:], self.output_path+'/'+event+'/cf/trial{}'.format(i))
 
                    for ista, sta in enumerate(sobj.channels):
                        datainfo['station_name'] = sta
                        datainfo['channel_name'] = 'CFP'  # note maximum three characters, the last one must be 'P'
                        ioformatting.vector2trace(datainfo, obs_dataP_das[ista,:], self.output_path+'/'+event+'/cf/trial{}'.format(i))
                        datainfo['channel_name'] = 'CFS'  # note maximum three characters, the last one must be 'S'
                        ioformatting.vector2trace(datainfo, obs_dataS_das[ista,:], self.output_path+'/'+event+'/cf/trial{}'.format(i))

                ######## modify 3D>>2D ##############


                print('input (STA) before locator', tp_mod_sta[0,0], ts_mod_sta[0,0], obs_dataP_sta[0,0], obs_dataS_sta[0,0], obs_dataP_sta.shape)

                print('input (DAS) before locator', tp_mod_das[0,0], ts_mod_das[0,0], obs_dataP_das[0,0], obs_dataS_das[0,0], obs_dataP_das.shape)
                
                ############à new 




                def validate_input_array(arr, name):
                    """ Validates that an array contains valid numeric values (no NaN, Inf, negative values). """
                    if isinstance(arr, num.ndarray):
                        # If it's a NumPy array, check for NaN or Inf values
                        if num.any(num.isnan(arr)) or num.any(num.isinf(arr)):
                            print(f"Error: {name} contains invalid value (NaN or Inf).")
                            return False
                        if num.any(arr < 0):
                            print(f"Error: {name} contains negative value.")
                            return False
                    elif isinstance(arr, list):
                        # If it's a Python list, iterate through it and validate
                        for row in arr:
                            if not isinstance(row, list):
                                print(f"Error: {name} contains non-list elements.")
                                return False
                            for val in row:
                                if math.isnan(val) or math.isinf(val):
                                    print(f"Error: {name} contains invalid value (NaN or Inf): {val}")
                                    return False
                                if val < 0:
                                    print(f"Error: {name} contains negative value: {val}")
                                    return False
                    else:
                        print(f"Error: {name} is neither a NumPy array nor a list.")
                        return False
                    return True

                def validate_npr(npr):
                    """ Validates that `npr` is a valid number. """
                    if not isinstance(npr, (int, float)):
                        print("Error: npr should be a number.")
                        return False
                    if npr <= 0:
                        print(f"Error: npr should be positive, but got {npr}.")
                        return False
                    return True


                # Validate the inputs
                if validate_input_array(tp_mod_sta, "tp_mod_sta") and \
                validate_input_array(ts_mod_sta, "ts_mod_sta") and \
                validate_input_array(obs_dataP_sta, "obs_dataP_sta") and \
                validate_input_array(obs_dataS_sta, "obs_dataS_sta") and \
                validate_npr(npr):

                    print('good, i am locating now')
                    # Proceed with the computation if all validations pass
                    iloctime_sta, corrmatrix_sta = location_t0.stacking(tp_mod_sta, ts_mod_sta, obs_dataP_sta, obs_dataS_sta, npr)
                    iloctime_das, corrmatrix_das = location_t0.stacking(tp_mod_das, ts_mod_das, obs_dataP_das[0:100, :], obs_dataS_das[0:100, :], npr)  # iloptime_das[0]: grid index; iloptime_das[1]: time index

                
                else:
                    print("Error: One or more inputs are invalid. Computation skipped.")

                #iloctime_sta, corrmatrix_sta = location_t0.stacking(tp_mod_sta, ts_mod_sta, obs_dataP_sta, obs_dataS_sta, npr)  # iloctime[0]: the grid index of the maximum stacking point; iloctime[1]: the time ndex at the maximum stacking point

                print('output location stations:', iloctime_sta, corrmatrix_sta[0,0], corrmatrix_sta.shape)
                print('output location DAS channels:', iloctime_das, corrmatrix_das[0,0], corrmatrix_das.shape)
                # 1. Compute evtpmin_sta for all stations in tp_modse_sta

                # corrmatrix is the stacking matrix, in 1D format but can be 
                # reformat to 3D format, each point saves the maximum stacking 
                # value during this calculation time period


                evtpmin_sta = {}
                for station_data in tp_mod_sta:
                    print(tp_mod_sta.shape)
                    print(station_data)
                    for station_key, station_array in station_data.items():
                        # Compute the minimum for each station's array
                        evtpmin_sta[station_key] = num.amin(station_array)

                # 2. Compute evtpmin_das for all stations in tp_modse_das
               # evtpmin_das = {}
               # for station_data in tp_mod_das:
               #     for station_key, station_array in station_data.items():
               #         # Compute the minimum for each DAS station's array
               #         evtpmin_das[station_key] = num.amin(station_array)

                # 3. Calculate event origin time for each station in tp_modse_sta (using evtpmin_sta)
#                for station_key_sta, evtpmin_sta_value in evtpmin_sta.items():
#                    if evtpmin_sta_value is not None:
#                        # Cast the numpy.float32 to a native Python float
#                        evtpmin_sta_value = float(evtpmin_sta_value)
#                        
#                        event_t0_sta = sobj.dtime_max_sta + datetime.timedelta(seconds=iloctime_sta[1]*sobj.deltat_sta) - datetime.timedelta(seconds=evtpmin_sta_value)  # event origin time for sta
#                        event_t0s_sta = event_t0_sta.isoformat()
#                        print(f"Event origin time for station {station_key_sta}: {event_t0s_sta}")
#                    else:
#                        print(f"Station {station_key_sta} has no minimum value in evtpmin_sta")
#
#                # 4. Process the stacking function for DAS (using evtpmin_das for tp_modse_das)
#                # Call the stacking function to get iloptime_das and corrmatrix_das


 #               print('tp_mod_das, ts_mod_das', tp_mod, ts_mod)

                


  #              # Now calculate event origin time for each DAS station
   #             for station_key_das, evtpmin_das_value in evtpmin_das.items():
    #                if evtpmin_das_value is not None:
     #                   # Cast the numpy.float32 to a native Python float
      #                  evtpmin_das_value = float(evtpmin_das_value)
       #                 
        #                # Assuming iloptime_das[0] is the grid index and iloptime_das[1] is the time index
         #               event_t0_das = sobj.dtime_max_das + datetime.timedelta(seconds=iloptime_das[1]*sobj.deltat_das) - datetime.timedelta(seconds=evtpmin_das_value)  # event origin time for das
          #              event_t0s_das = event_t0_das.isoformat()
           #             print(f"Event origin time for DAS station {station_key_das}: {event_t0s_das}")
            #        else:
             #           print(f"Station {station_key_das} has no minimum value in evtpmin_das")




                cmax_sta = num.max(corrmatrix_sta)
                cmax_das = num.max(corrmatrix_das)


                corrmatrix_sta = num.reshape(corrmatrix_sta,(tobj.nx,tobj.nx,tobj.nz))
                corrmatrix_das = num.reshape(corrmatrix_das,(tobj.nx,tobj.nx,tobj.nz))

                corrmatrix = corrmatrix_sta + corrmatrix_das


                (ixloc_sta, iyloc_sta, izloc_sta) = num.unravel_index(iloctime_sta[0],(tobj.nx,tobj.nx,tobj.nz))
                xloc_sta = tobj.x[ixloc_sta]
                yloc_sta = tobj.y[iyloc_sta]
                zloc_sta = tobj.z[izloc_sta]


                (ixloc_das, iyloc_das, izloc_das) = num.unravel_index(iloctime_das[0],(tobj.nx,tobj.nx,tobj.nz))
                xloc_das = tobj.x[ixloc_das]
                yloc_das = tobj.y[iyloc_das]
                zloc_das = tobj.z[izloc_das]
                
                # output the current location result
                if ntrial > 1:
                    cmfilename = self.output_path+'/'+event+'/'+event
                else:
                    cmfilename = self.output_path+'/'+event+'/'+event_t0s_sta
                out_file = open(cmfilename+'.loc', 'a')
                if STALTA:
                    out_file.write(str(i)+' '+str(xloc_sta)+' '+str(yloc_sta)+' '+str(zloc_sta)+' '+str(cmax_sta)+' '+str(nshort_p_sta)+' '+str(nshort_s_sta)+' '+str(slrat)+'\n')
                else:
                    out_file.write(str(i)+' '+str(xloc_sta)+' '+str(yloc_sta)+' '+str(zloc_sta)+' '+str(cmax_sta)+'\n')
                out_file.close()

                
                # output the current location result
                if ntrial > 1:
                    cmfilename = self.output_path+'/'+event+'/'+event
                else:
                    cmfilename = self.output_path+'/'+event+'/'+event_t0s_das
                out_file = open(cmfilename+'.loc', 'a')
                if STALTA:
                    out_file.write(str(i)+' '+str(xloc_das)+' '+str(yloc_das)+' '+str(zloc_das)+' '+str(cmax_das)+' '+str(nshort_p_das)+' '+str(nshort_s_das)+' '+str(slrat)+'\n')
                else:
                    out_file.write(str(i)+' '+str(xloc_das)+' '+str(yloc_das)+' '+str(zloc_das)+' '+str(cmax_das)+'\n')
                out_file.close()  

                # save the stacked coherence matrix
                num.save(self.output_path+'/'+event+'/'+'corrmatrix_trial_'+str(i),corrmatrix)
                
                # plot migration profiles
                self.coherence_plot(self.output_path+'/'+event, corrmatrix, tobj.x, tobj.y, tobj.z, i)
            
                # output theoretical P- and S-wave arrivaltimes
                fname_sta = cmfilename + '_trial_sta{}.phs'.format(i)
                self.write_phasetime(sobj.stations, event_t0_sta, tp_modse, ts_modse, iloctime_sta[0], fname_sta)

                # output theoretical P- and S-wave arrivaltimes
                fname_das = cmfilename + '_trial_das{}.phs'.format(i)
                self.write_phasetime(sobj.stations, event_t0_das, tp_modse, ts_modse, iloctime_das[0], fname_das)

                if cmax_sta > cmax_pre:
                    event_t0s_final_sta = copy.deepcopy(event_t0s_sta)
                    cmax_pre = copy.deepcopy(cmax_sta)

                if cmax_das > cmax_pre:
                    event_t0s_final_das = copy.deepcopy(event_t0s_das)
                    cmax_pre = copy.deepcopy(cmax_das)
            
            event_t0s_final = event_t0s_final_sta + event_t0s_final_das

            print('Total correlation matrix:', corrmatrix)

            self.catalogue_creation(event, event_t0s_final, tobj.lat0, tobj.lon0, ntrial, corrmatrix)
        print('Location process completed!!!')
        gc.collect()

    def catalogue_creation(self, event, event_t0s, lat0, lon0, ntrial, corrmatrix, refell=23):
        latref=lat0; lonref=lon0; eleref=0.
        origin=latlon2cart.Coordinates(latref,lonref,eleref)
        if (ntrial > 1):
            ev_file = self.output_path+'/'+event+'/'+event+'.loc'
            data = num.loadtxt(ev_file)
            w = num.sum(data[:, 4])
            xb = ((num.dot(data[:, 1], data[:, 4])/w)*1000)
            yb = ((num.dot(data[:, 2], data[:, 4])/w)*1000)
            late,lone,elev=origin.cart2geo(xb,yb,eleref)
            zb = num.dot(data[:, 3], data[:, 4])/w  # depth in km
            cb = num.mean(data[:, 4])  # the mean coherence over the ntrial realizations
            cmax = num.max(data[:, 4])  # the maximum coherence over the ntrial realizations
            merr = num.vstack((data[:, 1], data[:, 2], data[:, 3]))
            err = num.cov(merr)
            errmax = num.sqrt(num.max(num.linalg.eigvals(err)))
        else:
            ev_file = self.output_path+'/'+event+'/'+event_t0s+'.loc'
            data = num.loadtxt(ev_file)
            xb=data[1]*1000
            yb=data[2]*1000
            zb = data[3]  # depth in km
            late,lone,elev=origin.cart2geo(xb,yb,eleref) # latitude, longitude
            cmax = data[4]  # the maximum coherence over the 3D corrmatrix
            
            # nomalize corrmatrix first, let minimal->1, maximum->2
            n1 = 1.0  # minimal limit
            n2 = 2.0  # maximum limit
            dmax = num.amax(corrmatrix, axis=None, keepdims=True)
            dmin = num.amin(corrmatrix, axis=None, keepdims=True)
            k = (n2-n1)/(dmax-dmin)
            b = (dmax*n1-dmin*n2)/(dmax-dmin)
            corrmatrix = k*corrmatrix + b
            
            errmax = num.std(corrmatrix, axis=None)  # the coherence standard deviation over the 3D corrmatrix
            cb = num.median(corrmatrix, axis=None)  # the median coherence over the 3D corrmatrix
            
        f = open(self.output_path+'/'+'catalogue', 'a')
        f.write(event_t0s+'    '+str(late)+'   '+str(lone)+'   '+str(zb)+'   '+str(errmax)+'   '+str(cb)+'   '+str(cmax)+'\n')
        f.close()

    def coherence_plot(self, event_path, corrmatrix, xax, yax, zax, itrial, normalization=False):
        nx, ny, nz = num.shape(corrmatrix)
        CXY = num.zeros([ny, nx])
        for i in range(ny):
            for j in range(nx):
                CXY[i,j]=num.max(corrmatrix[j,i,:])

        CXZ = num.zeros([nz, nx])
        for i in range(nz):
            for j in range(nx):
                CXZ[i, j] = num.max(corrmatrix[j,:,i])

        CYZ = num.zeros([nz, ny])
        for i in range(nz):
            for j in range(ny):
                CYZ[i, j] = num.max(corrmatrix[:, j, i])

        if normalization:
            nrm = Normalize(vmin=0., vmax=1.)
        else:
            nrm = None


        xticks=num.min(xax)+num.arange(6)*(num.max(xax)-num.min(xax))/5
        yticks=num.min(yax)+num.arange(6)*(num.max(yax)-num.min(yax))/5
        zticks=num.min(zax)+num.arange(6)*(num.max(zax)-num.min(zax))/5

        fig, axs = plt.subplots(1,3, figsize=(15, 7.5))
        fig.suptitle('Coherence matrices trial '+str(itrial), fontsize=14, fontweight='bold')
        cmap = plt.cm.get_cmap('viridis', 100)
    
        ax1 = axs[0]
        cs1=ax1.contourf(xax, yax, CXY, 20, cmap=cmap, norm=nrm)
        ax1.set_xticks(xticks)
        ax1.set_yticks(yticks)
        ax1.set_xlabel('X (km)')
        ax1.set_ylabel('Y (km)')
        ax1.set_aspect('auto')
        
        ax2 = axs[1]
        cs2=ax2.contourf(yax, zax, CYZ, 20, cmap=cmap, norm=nrm)
        ax2.set_xticks(yticks)
        ax2.set_yticks(zticks)
        ax2.set_xlabel('Y (km)')
        ax2.set_ylabel('Z (km)')
        ax2.set_aspect('auto')
        ax2.invert_yaxis()

        ax3 = axs[2]
        cs3=ax3.contourf(xax, zax, CXZ, 20, cmap=cmap, norm=nrm)
        ax3.set_xticks(xticks)
        ax3.set_yticks(zticks)
        ax3.set_xlabel('X (km)')
        ax3.set_ylabel('Z (km)')
        ax3.set_aspect('auto')
        ax3.invert_yaxis()
        

        cbar=plt.colorbar(cs1, ax=axs, orientation='horizontal', shrink=0.6)
        cbar.set_label('Coherence')

        plt.savefig(event_path+'/'+'Coherence_matrix_'+str(itrial)+'.eps')
        plt.close("all")
        
    
    def write_phasetime(self, stations, event_t0, tp_modse, ts_modse, grididx, fname):
        """
        Calculate the theoretical arrival-times of P- and S-phases for the located
        event and output to a text file.

        Parameters
        ----------
        stations : list of str
            station names.
        event_t0 : datetime
            event origin time.
        tp_modse : numpy array, shape: n_stations*n_grids
            P-wave traveltime table in second.
        ts_modse : numpy array, shape: n_stations*n_grids
            S-wave traveltime table in second.
        grididx : int
            grid index where the seismic event is located.
        fname : str
            output filename including path.

        Returns
        -------
        None.

        """
        
        ofile = open(fname, 'a')
        ofile.write('# station    P_arrivaltime    S_arrivaltime \n')
        
        for ii, sta in enumerate(stations):
            # loop over each station to output the theoretical arrival-times for the P- and S-phases
            tp_tavt = event_t0 + datetime.timedelta(seconds=tp_modse[grididx,ii])  # P_arrival-time = event_origin_time + P_traveltime
            ts_tavt = event_t0 + datetime.timedelta(seconds=ts_modse[grididx,ii])  # S_arrival-time = event_origin_time + S_traveltime
            ofile.write(sta+' '+tp_tavt.isoformat()+' '+ts_tavt.isoformat()+'\n')
            ofile.flush()
            
        ofile.close()        
        
        return None

