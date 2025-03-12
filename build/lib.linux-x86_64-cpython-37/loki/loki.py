# %%
import os
from obspy.core import read
import math
import numpy as num
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.pyplot as plt
import datetime
import copy
import gc
from loki import ioformatting
from loki import traveltimes
from loki import waveforms
from loki import stacktraces
from loki import latlon2cart
from loki import location_t0_py
import tt_processing                       # C
import location_t0                         # C  for multiplying the P- and S-stacking values using this
#import location_t0_plus                   # C  for adding the P- and S-stacking values using this


class Loki:
    """docstring for Loki"""

    def __init__(self, data_path, output_path, db_path, hdr_filename, geometry_filename_fiber, geometry_filename_stat, mode='locator'):
        self.data_path = data_path
        self.output_path = output_path
        self.db_path = db_path
        self.hdr_filename = hdr_filename
        self.geometry_filename_stat = geometry_filename_stat
        self.geometry_filename_fiber = geometry_filename_fiber
        
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
        data_tree = []
        
        for root, dirs, files in os.walk(data_path):
            if "stations" in dirs and "fiber" in dirs:  # Ensure both folders exist
                data_tree.append({
                    "event_path": root,
                    "stations": os.path.join(root, "stations"),
                    "fiber": os.path.join(root, "fiber")
                })
        
        data_tree.sort(key=lambda x: x["event_path"])  # Sort based on event path

        events = [os.path.basename(entry["event_path"]) for entry in data_tree]

        if not os.path.isdir(output_path):
            os.mkdir(output_path)

        return data_tree, events

    def location(self, comp=['E', 'N', 'Z'], precision='single', **inputs):
        if 'tshortp_min_sta' in inputs:
            # need to calculate STA/LTA for stacking
            STALTA = True
            tshortp_min = inputs['tshortp_min_sta']
            tshortp_max = inputs['tshortp_max_sta']
            tshorts_min = inputs['tshorts_min_sta']
            tshorts_max = inputs['tshorts_max_sta']
            tshortp_min_fiber = inputs['tshortp_min_fiber']
            tshortp_max_fiber = inputs['tshortp_max_fiber']
            tshorts_min_fiber = inputs['tshorts_min_fiber']
            tshorts_max_fiber = inputs['tshorts_max_fiber']
            slrat = inputs['slrat']
            ntrial = inputs['ntrial']
            #extension_sta = inputs['extension_sta']
            #extension_das = inputs['extension_das']
            #delta_das = inputs['delta_das']
            tshortp = num.linspace(tshortp_min, tshortp_max, ntrial)
            tshorts = num.linspace(tshorts_min, tshorts_max, ntrial)
            tshortp_fiber = num.linspace(tshortp_min_fiber, tshortp_max_fiber, ntrial)
            tshorts_fiber = num.linspace(tshorts_min_fiber, tshorts_max_fiber, ntrial)
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


        #print(self.geometry_filename_fiber)

        tobj = traveltimes.Traveltimes(self.db_path, self.hdr_filename, self.geometry_filename_fiber, self.geometry_filename_stat)


        tp = tobj.load_traveltimes('P', model, precision) 


        #print('dimension of the tt', tp['HM00'].shape)

        ts = tobj.load_traveltimes('S', model, precision)


        #load id of stations, channels and their location 
    
        tobj.load_station_info()
        #tobj.load_channel_info()


        #for each event, locate with stations, fiber, save the two, then merge the correlation 

        print('here!!')
        
        for event_path in self.data_tree:

            #print('This is the data tree:', self.data_tree)

            #print('I am reading the data in this path:', event_path)


            self.subdata_path = event_path
        
            #for folder in ["fiber", "stations"]:   #access stations and fiber independently
            #    self.folder_path = os.path.join(event_path, folder)

                
            #    if os.path.exists(self.folder_path):  # Ensure the folder exists


            #print('I am reading data from this directory: ', self.subdata_path)

            #Reading the observed wavefields (stations and DAS)

            #wobj = waveforms.Waveforms(tobj=tobj, data_path = data_path, event_path=event_path, extension_sta="*", extension_das='CANDAS2_2023-01-07_10-48-10.h5', freq=None)
            
            last_folder = os.path.basename(self.subdata_path)  # Get the last folder name

            print(last_folder)


            if last_folder == "hybrid":
            
                label = "hybrid"
                
                print(f"Rete ibrida!")

                st = read(os.path.join(event_path,"*"))

                components = set(tr.stats.channel for tr in st)

                if len(components) == 1: 
                    comps = ['Z']

                if len(components) == 3: 
                    comps = ['E','N','Z']

                print('number of components in hybrid:', comps)
   
                wobj = waveforms.Waveforms(self.subdata_path, extension_sta="*", comps=comps, freq=None)

                sobj = stacktraces.Stacktraces(tobj, wobj, **inputs)


            if last_folder == "stations":
            
                label = "station"
                
                print(f"Stazioni!")
                    
                wobj = waveforms.Waveforms(self.subdata_path, extension_sta="*", comps=['E','N','Z'], freq=None)

                sobj = stacktraces.Stacktraces(tobj, wobj, **inputs)

            if last_folder == 'fiber':
            
                print(f"Fibre!")
                    
                label = "fibre"
                        
                wobj = waveforms.Waveforms(self.subdata_path, extension_sta="*", comps=['Z'], freq=None)

                sobj = stacktraces.Stacktraces(tobj, wobj, **inputs)

            print('current sobj.deltat', sobj.deltat)
            event = event_path.split('/')[-1]

            print('Processing to the event folder: ', event_path, event)
            if os.path.isdir(self.output_path+'/'+event):
                print('directory already exists')
                #continue
            else:
                os.mkdir(self.output_path+'/'+event)

    
            tpxz=tp['HM00'].reshape(tobj.nxz, 1)
            tsxz=ts['HM00'].reshape(tobj.nxz, 1)

            print('tpxz:',  tpxz.shape, tpxz.dtype)
            print('tsxz:',  tsxz.shape, tsxz.dtype)


            tpxz = num.asarray(tpxz, dtype=num.float64)
            tsxz = num.asarray(tsxz, dtype=num.float64)


            tp_modse = num.ascontiguousarray(tpxz)
            ts_modse = num.ascontiguousarray(tsxz)


            tp_mod_sta, ts_mod_sta = tt_processing.tt_f2i(sobj.deltat, tp_modse, ts_modse, npr)  # traveltime table in time sample


            cmax_pre = -1.0
            for i in range(ntrial):
                if STALTA:

                    # need to calculate STA/LTA from the characteristic funtion
                    # then stack the STA/LTA for imaging

                    if last_folder == 'fiber':

                        nshort_p_sta = int(tshortp_fiber[i]//sobj.deltat)
                        nshort_s_sta = int(tshorts_fiber[i]//sobj.deltat)

                    else:


                        nshort_p_sta = int(tshortp[i]//sobj.deltat)
                        nshort_s_sta = int(tshorts[i]//sobj.deltat)


                    obs_dataP_sta, obs_dataS_sta = sobj.loc_stalta(nshort_p_sta, nshort_s_sta, slrat, norm=1)

                else:

                    # no need to calculate STA/LTA 
                    # directly stack the characteristic function for imaging
                    obs_dataP_sta = sobj.obs_dataV_sta  # vertical -> P
                    obs_dataS_sta = sobj.obs_dataH_sta  # horizontal -> S


                if opsf:

                    datainfo = {}
                    datainfo['dt'] = sobj.deltat
                    datainfo['starttime'] = sobj.dtime_max
                    for ista, sta in enumerate(sobj.stations):
                        print('sta:', sta)
                        print('ista:', ista)
                        datainfo['station_name'] = sta
                        datainfo['channel_name'] = 'CFP'  # note maximum three characters, the last one must be 'P'
                        ioformatting.vector2trace(datainfo, obs_dataP_sta[ista,:], self.output_path+'/'+event+'/cf/trial{}'.format(i))
                        datainfo['channel_name'] = 'CFS'  # note maximum three characters, the last one must be 'S'
                        ioformatting.vector2trace(datainfo, obs_dataS_sta[ista,:], self.output_path+'/'+event+'/cf/trial{}'.format(i))


                print("tp_mod_sta shape:", tp_mod_sta.shape)
                print("ts_mod_sta shape:", ts_mod_sta.shape)
                print("obs_dataP_sta shape:", obs_dataP_sta.shape)
                print("obs_dataS_sta shape:", obs_dataS_sta.shape)
                print("npr:", npr)
                
                
                #python stack 

                #iloc, itime, corrmatrix = location_t0.stacking(itp, its, stalta_p, stalta_s, nproc)
                #stacking = location_t0_py.WaveformStacking(tobj, sobj, npr, tp_mod_sta, ts_mod_sta, obs_dataP_sta[:,:], obs_dataS_sta[:,:], obs_dataP_das[0:2,:], obs_dataS_das[0:2,:])
                #iloc_sta, iloc_ch, iloc, itime, corrmatrix_sta, corrmatrix_ch, corrmatrix = stacking.locate_event()
                
                x_stations_tot = num.array(tobj.lon_stations, dtype=float) #change name, otherwise misleading (x,y,z)
                y_stations_tot = num.array(tobj.lat_stations, dtype=float) 
                z_stations_tot = num.array(tobj.depth_stations, dtype=float)


                tp_mod_sta = tp_mod_sta.reshape(tobj.nx,tobj.nz)
                ts_mod_sta = ts_mod_sta.reshape(tobj.nx,tobj.nz)

                nsta = len(obs_dataP_sta[:, 1])
                
                x_stations =[]
                y_stations =[]
                z_stations =[]


            #LOOP ON THE STATIONS 
                for j in range(nsta):

                    
                    current_sta = sobj.stations[j]

                    #print(current_sta)

                    # Access the tuple (lon, lat, depth) from the dictionary
                    lon, lat, depth = tobj.stations_coordinates.get(current_sta, (None, None, None))

                    #print(tobj.stations_coordinates)
                    #print(lon,lat,depth)

                    x_stations.append(lon)
                    y_stations.append(lat)
                    z_stations.append(depth)
 
                    #if current_sta[2] == '0':
                    #    current_sta = num.int(current_sta[3]) -1
                    #else:
                    #    current_sta = num.int(current_sta[2:4]) -1

                    #x_stations.append(x_stations_tot[current_sta])
                    #y_stations.append(y_stations_tot[current_sta])
                    #z_stations.append(z_stations_tot[current_sta])

                x_stations = num.array(x_stations, dtype=float)
                y_stations = num.array(y_stations, dtype=float)
                z_stations = num.array(z_stations, dtype=float)



                # Ensure inputs are contiguous and have the correct types
                tp_mod_sta = num.ascontiguousarray(tp_mod_sta, dtype=num.int32)
                ts_mod_sta = num.ascontiguousarray(ts_mod_sta, dtype=num.int32)
                x_stations = num.ascontiguousarray(x_stations, dtype=num.float64)
                y_stations = num.ascontiguousarray(y_stations, dtype=num.float64)
                z_stations = num.ascontiguousarray(z_stations, dtype=num.float64)
                tobj_x = num.ascontiguousarray(tobj.x, dtype=num.float64)
                tobj_y = num.ascontiguousarray(tobj.y, dtype=num.float64)
                tobj_z = num.ascontiguousarray(tobj.z, dtype=num.float64)
                obs_dataP_sta = num.ascontiguousarray(obs_dataP_sta, dtype=num.float64)
                obs_dataS_sta = num.ascontiguousarray(obs_dataS_sta, dtype=num.float64)
                
                print('I am checking for memory contiguity')
                
                print(tp_mod_sta.flags['C_CONTIGUOUS'])  # True if C-contiguous)
                print(ts_mod_sta.flags['C_CONTIGUOUS'])  # True if C-contiguous)
                print(x_stations.flags['C_CONTIGUOUS'])  # True if C-contiguous)
                print(y_stations.flags['C_CONTIGUOUS'])  # True if C-contiguous)
                print(z_stations.flags['C_CONTIGUOUS'])  # True if C-contiguous)
                print(tobj_x.flags['C_CONTIGUOUS'])  # True if C-contiguous)
                print(tobj_y.flags['C_CONTIGUOUS'])  # True if C-contiguous)
                print(tobj_z.flags['C_CONTIGUOUS'])  # True if C-contiguous)
                print(obs_dataP_sta.flags['C_CONTIGUOUS'])  # True if C-contiguous)
                print(obs_dataS_sta.flags['C_CONTIGUOUS'])  # True if C-contiguous)

                iloctime, corrmatrix = location_t0.stacking(tp_mod_sta, ts_mod_sta, x_stations, y_stations, z_stations, tobj_x, tobj_y, tobj_z, obs_dataP_sta, obs_dataS_sta, npr)

                # Step 2: Save the 3D array
                num.save(os.path.dirname(event_path).rsplit("/", 1)[0] + "/" + os.path.dirname(event_path).rstrip("/").split("/")[-1]  + '_' + label + "_coherence_matrix.npy", corrmatrix)

                
            print('Now creating hybrid coherence map')
                

            # Construct the path to the coherence_fibre file
            coherence_fibre_path = os.path.dirname(event_path).rsplit("/", 1)[0] + "/" + os.path.dirname(event_path).rstrip("/").split("/")[-1] + "_fibre_coherence_matrix.npy"
            coherence_stations_path = os.path.dirname(event_path).rsplit("/", 1)[0] + "/" + os.path.dirname(event_path).rstrip("/").split("/")[-1] + "_station_coherence_matrix.npy"
            coherence_full_hybrid_path = os.path.dirname(event_path).rsplit("/", 1)[0] + "/" + os.path.dirname(event_path).rstrip("/").split("/")[-1] + "_hybrid_coherence_matrix.npy"


            # Check if the file exists
            if os.path.exists(coherence_full_hybrid_path):
                # Load the file if it exists
                coherence_full_hybrid = num.load(os.path.dirname(event_path).rsplit("/", 1)[0] + "/" + os.path.dirname(event_path).rstrip("/").split("/")[-1] + "_hybrid_coherence_matrix.npy")
                #coherence_hybrid = (coherence_stations + coherence_fibre)/2
                #num.save(os.path.dirname(event_path).rsplit("/", 1)[0] + "/" + os.path.dirname(event_path).rstrip("/").split("/")[-1] + "_hybrid_coherence_matrix.npy", coherence_hybrid)
            else:
                print(f"fully hybrid")

            # Check if the file exists
            if os.path.exists(coherence_fibre_path):
                # Load the file if it exists
                coherence_fibre = num.load(os.path.dirname(event_path).rsplit("/", 1)[0] + "/" + os.path.dirname(event_path).rstrip("/").split("/")[-1] + "_fibre_coherence_matrix.npy")
                #coherence_hybrid = (coherence_stations + coherence_fibre)/2
                #num.save(os.path.dirname(event_path).rsplit("/", 1)[0] + "/" + os.path.dirname(event_path).rstrip("/").split("/")[-1] + "_hybrid_coherence_matrix.npy", coherence_hybrid)
            else:
                print(f"no fiber for this event")

            # Check if the file exists
            if os.path.exists(coherence_stations_path):
                # Load the file if it exists
                coherence_stations = num.load(os.path.dirname(event_path).rsplit("/", 1)[0] + "/" + os.path.dirname(event_path).rstrip("/").split("/")[-1] + "_station_coherence_matrix.npy")

            else:
                print(f"no stations for this event")

            if os.path.exists(coherence_fibre_path) & os.path.exists(coherence_stations_path):

                print(f"stations and fiber for this event!")

                coherence_stations_norm = coherence_stations/max(coherence_stations)
                coherence_fibre_norm = coherence_fibre/max(coherence_fibre)
                
                coherence_hybrid = (coherence_stations_norm + coherence_fibre_norm)/2
                num.save(os.path.dirname(event_path).rsplit("/", 1)[0] + "/" + os.path.dirname(event_path).rstrip("/").split("/")[-1] + "_hybrid_only_sum_coherence_matrix.npy", coherence_hybrid)

'''                
                iloc0, iloc1, iloc2, itime = iloctime
                print(f"Best location index: {iloc0}")
                print(f"Best location indices (tt): ({iloc1}, {iloc2})")
                print(f"Best time index: {itime}")

                tp_modse_2d = num.reshape(tp_modse,(tobj.nx,tobj.nz))
                ts_modse_2d = num.reshape(ts_modse,(tobj.nx,tobj.nz))
                evtpmin = num.amin(tp_modse_2d[iloctime[1],iloctime[2]])
                event_t0 = sobj.dtime_max + datetime.timedelta(seconds=iloctime[3]*sobj.deltat) - datetime.timedelta(seconds=evtpmin)  # event origin time
                event_t0s = (event_t0).isoformat()
                # corrmatrix is the stacking matrix, in 1D format but can be 
                # reformat to 3D format, each point saves the maximum stacking 
                # value during this calculation time period
                cmax = num.max(corrmatrix)
                corrmatrix = num.reshape(corrmatrix,(tobj.nx,tobj.nx,tobj.nz))
                (ixloc, iyloc, izloc) = num.unravel_index(iloctime[0],(tobj.nx,tobj.nx,tobj.nz))
                xloc = tobj.x[ixloc]
                yloc = tobj.y[iyloc]
                zloc = tobj.z[izloc]
                
                # output the current location result
                if ntrial > 1:
                    cmfilename = self.output_path+'/'+event+'/'+event
                else:
                    cmfilename = self.output_path+'/'+event+'/'+event_t0s
                out_file = open(cmfilename+'.loc', 'a')
                if STALTA:
                    out_file.write(str(i)+' '+str(xloc)+' '+str(yloc)+' '+str(zloc)+' '+str(cmax)+' '+str(nshort_p_sta)+' '+str(nshort_s_sta)+' '+str(slrat)+'\n')
                else:
                    out_file.write(str(i)+' '+str(xloc)+' '+str(yloc)+' '+str(zloc)+' '+str(cmax)+'\n')
                out_file.close()
                
                # save the stacked coherence matrix
                num.save(self.output_path+'/'+event+'/'+'corrmatrix_trial_'+str(i),corrmatrix)
                
                # plot migration profiles
                self.coherence_plot(self.output_path+'/'+event, corrmatrix, tobj.x, tobj.y, tobj.z, i)
            
                # output theoretical P- and S-wave arrivaltimes
                fname = cmfilename + '_trial{}.phs'.format(i)
                self.write_phasetime(sobj.stations, iloctime, event_t0, tp_modse_2d,ts_modse_2d, fname)

                if cmax > cmax_pre:
                    event_t0s_final = copy.deepcopy(event_t0s)
                    cmax_pre = copy.deepcopy(cmax)
            


            self.catalogue_creation(event, event_t0s_final, tobj.lat0, tobj.lon0, ntrial, corrmatrix)
        print('Location process completed!!!')
        gc.collect()

#methods
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
        
    
    def write_phasetime(self, stations, iloctime, event_t0, tp_modse, ts_modse, fname):
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
            tp_tavt = event_t0 + datetime.timedelta(seconds=tp_modse[iloctime[1],iloctime[2]])  # P_arrival-time = event_origin_time + P_traveltime
            ts_tavt = event_t0 + datetime.timedelta(seconds=ts_modse[iloctime[1],iloctime[2]])  # S_arrival-time = event_origin_time + S_traveltime
            ofile.write(sta+' '+tp_tavt.isoformat()+' '+ts_tavt.isoformat()+'\n')
            ofile.flush()
            
        ofile.close()        
        
        return None
'''

# %%
