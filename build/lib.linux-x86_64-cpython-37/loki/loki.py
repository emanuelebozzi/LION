# %%
import os
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
#import location_t0                         # C  for multiplying the P- and S-stacking values using this
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


        print(self.geometry_filename_fiber)

        tobj = traveltimes.Traveltimes(self.db_path, self.hdr_filename, self.geometry_filename_fiber, self.geometry_filename_stat)


        tp = tobj.load_traveltimes('P', model, precision) 
        ts = tobj.load_traveltimes('S', model, precision)

        #load id of stations, channels and their location 
    
        tobj.load_station_info()
        tobj.load_channel_info()


        #for each event 

        for event_path in self.data_tree:

            data_path = self.data_path

            #Reading the observed wavefields (stations and DAS)

            wobj = waveforms.Waveforms(tobj=tobj, data_path = data_path, event_path=event_path, extension_sta="*", extension_das='CANDAS2_2023-01-07_10-48-10.h5', freq=None)


            #object of the class stacktraces  

            sobj = stacktraces.Stacktraces(tobj, wobj, **inputs)


            event = event_path.split('/')[-1]

            print('Processing to the event folder: ', event_path, event)
            if os.path.isdir(self.output_path+'/'+event):
                print('directory already exists')
                #continue
            else:
                os.mkdir(self.output_path+'/'+event)

   
            tpxz=tp['HM01'].reshape(tobj.nxz, 1)
            tsxz=ts['HM01'].reshape(tobj.nxz, 1)

            print('tpxz:',  tpxz.shape, tpxz.dtype)
            print('tsxz:',  tsxz.shape, tsxz.dtype)


            tpxz = num.asarray(tpxz, dtype=num.float64)
            tsxz = num.asarray(tsxz, dtype=num.float64)


            tp_modse = num.ascontiguousarray(tpxz)
            ts_modse = num.ascontiguousarray(tsxz)

            print('sobj.deltat_sta', sobj.deltat_sta)
            print('sobj.deltat_das', sobj.deltat_das)
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

                ######## pre-location plots ##############

                # Save the observed P-wave data for STA (Station)
                # Assuming obs_dataP_sta is the data with shape (6591, 8000)
                plt.figure(figsize=(10, 10))
                plt.plot(tp_modse, '*')
                plt.title('tp mod sta')


                # Save the figure as a PNG file
                plt.tight_layout()
                plt.savefig('tp_modse.png')  # Save the plot to a file
                plt.close()  # Close the plot to free up memory


                plt.figure(figsize=(10, 10))
                plt.plot(tp_mod_sta, '*')
                plt.title('tp mod sta')


                # Save the figure as a PNG file
                plt.tight_layout()
                plt.savefig('tp_mod_sta.png')  # Save the plot to a file
                plt.close()  # Close the plot to free up memory

                # Save the observed P-wave data for STA (Station)
                # Assuming obs_dataP_sta is the data with shape (6591, 8000)
                plt.figure(figsize=(12, 8))
                plt.imshow(obs_dataP_sta, aspect='auto', cmap='seismic', origin='lower', interpolation='none')
                plt.colorbar(label='Amplitude')
                plt.title('Observed P-Wave Data (STA) - Space vs Time')
                plt.xlabel('Time Samples')
                plt.ylabel('Station Index')

                # Save the figure as a PNG file
                plt.tight_layout()
                plt.savefig('observed_p_wave_sta.png')  # Save the plot to a file
                plt.close()  # Close the plot to free up memory

                # Save the observed P-wave data for DAS (Distributed Acoustic Sensing)
                # Assuming obs_dataP_das is the data with shape (6591, 8000)
                plt.figure(figsize=(12, 8))
                plt.imshow(obs_dataP_das, aspect='auto', cmap='seismic', origin='lower', interpolation='none')
                plt.colorbar(label='Amplitude')
                plt.title('Observed P-Wave Data (DAS) - Space vs Time')
                plt.xlabel('Time Samples')
                plt.ylabel('DAS Channel Index')

                # Save the figure as a PNG file
                plt.tight_layout()
                plt.savefig('observed_p_wave_das.png')  # Save the plot to a file
                plt.close()  # Close the plot to free up memory

                ############ 

                print("tp_mod_sta shape:", tp_mod_sta.shape)
                print("ts_mod_sta shape:", ts_mod_sta.shape)
                print("obs_dataP_sta shape:", obs_dataP_sta.shape)
                print("obs_dataS_sta shape:", obs_dataS_sta.shape)
                print("npr:", npr)

                #iloc, itime, corrmatrix = location_t0.stacking(itp, its, stalta_p, stalta_s, nproc)

                stacking = location_t0_py.WaveformStacking(tobj, npr, tp_mod_sta, ts_mod_sta, obs_dataP_sta[0:2,:], obs_dataS_sta[0:2,:], obs_dataP_das[0:2,:], obs_dataS_das[0:2,:])
                iloc_sta, iloc_ch, iloc, itime, corrmatrix_sta, corrmatrix_ch, corrmatrix = stacking.locate_event()
                 
                #save 



           
                # Step 2: Save the 3D array
                num.save("array_3d.npy", corrmatrix)

                # Step 3: Extract slices
                slice_xy = corrmatrix[:, :, 25]  # Middle slice along Z-axis
                slice_xz = corrmatrix[:, 25, :]  # Middle slice along Y-axis
                slice_yz = corrmatrix[25, :, :]  # Middle slice along X-axis

                # Step 4: Plot the slices in subplots
                fig, axes = plt.subplots(1, 3, figsize=(12, 4))

                axes[0].imshow(slice_xy, cmap="gray")
                axes[0].set_title("XY Plane (Z=25)")

                axes[1].imshow(slice_xz, cmap="gray")
                axes[1].set_title("XZ Plane (Y=25)")

                axes[2].imshow(slice_yz, cmap="gray")
                axes[2].set_title("YZ Plane (X=25)")

                for ax in axes:
                    ax.axis("off")  # Hide axes for better visualization

                # Step 5: Save the figure
                plt.savefig("slices_subplot.png", dpi=300, bbox_inches="tight")
                plt.show()


                ####


'''

                #stations 
                evtpmin_sta = num.amin(tp_modse[iloc_sta[0],:])
                event_t0_sta = sobj.dtime_max + datetime.timedelta(seconds=iloc_sta[1]*sobj.deltat_sta) - datetime.timedelta(seconds=evtpmin_sta)  # event origin time
                event_t0s_sta = (event_t0_sta).isoformat()
                # corrmatrix is the stacking matrix, in 1D format but can be 
                # reformat to 3D format, each point saves the maximum stacking 
                # value during this calculation time period
                cmax_sta = num.max(corrmatrix_sta)
                #corrmatrix = num.reshape(corrmatrix,(tobj.nx,tobj.ny,tobj.nz))
                (ixloc_sta, iyloc_sta, izloc_sta) = num.unravel_index(iloc_sta[0],(tobj.nx,tobj.nx,tobj.nz))
                xloc_sta = tobj.x[ixloc_sta]
                yloc_sta = tobj.y[iyloc_sta]
                zloc_sta = tobj.z[izloc_sta]
                
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
                
                # save the stacked coherence matrix
                num.save(self.output_path+'/'+event+'/'+'corrmatrix_trial_'+str(i),corrmatrix_sta)
                
                # plot migration profiles
                self.coherence_plot(self.output_path+'/'+event, corrmatrix_sta, tobj.x, tobj.y, tobj.z, i)
            
                # output theoretical P- and S-wave arrivaltimes
                fname_sta = cmfilename + '_trial{}.phs'.format(i)
                self.write_phasetime(sobj.stations, event_t0_sta, tp_modse, ts_modse, iloc_sta[0], fname_sta)
                
                if cmax_sta > cmax_pre:
                    event_t0s_final_sta = copy.deepcopy(event_t0s_sta)
                    cmax_pre = copy.deepcopy(cmax_sta)


                #fiber  
                evtpmin_ch = num.amin(tp_modse[iloc_ch[0],:])
                event_t0_ch = sobj.dtime_max + datetime.timedelta(seconds=iloc_sta[1]*sobj.deltat_das) - datetime.timedelta(seconds=evtpmin_ch)  # event origin time
                event_t0s_ch = (event_t0_ch).isoformat()
                # corrmatrix is the stacking matrix, in 1D format but can be 
                # reformat to 3D format, each point saves the maximum stacking 
                # value during this calculation time period
                cmax_ch = num.max(corrmatrix_ch)
                #corrmatrix = num.reshape(corrmatrix,(tobj.nx,tobj.ny,tobj.nz))
                (ixloc_ch, iyloc_ch, izloc_ch) = num.unravel_index(iloc_ch[0],(tobj.nx,tobj.nx,tobj.nz))
                xloc_ch = tobj.x[ixloc_ch]
                yloc_ch = tobj.y[iyloc_ch]
                zloc_ch = tobj.z[izloc_ch]
                
                # output the current location result
                if ntrial > 1:
                    cmfilename = self.output_path+'/'+event+'/'+event
                else:
                    cmfilename = self.output_path+'/'+event+'/'+event_t0s_ch
                out_file = open(cmfilename+'.loc', 'a')
                if STALTA:
                    out_file.write(str(i)+' '+str(xloc_ch)+' '+str(yloc_ch)+' '+str(zloc_ch)+' '+str(cmax_ch)+' '+str(nshort_p_sta)+' '+str(nshort_s_sta)+' '+str(slrat)+'\n')
                else:
                    out_file.write(str(i)+' '+str(xloc_ch)+' '+str(yloc_ch)+' '+str(zloc_ch)+' '+str(cmax_sta)+'\n')
                out_file.close()
                
                # save the stacked coherence matrix
                num.save(self.output_path+'/'+event+'/'+'corrmatrix_trial_'+str(i),corrmatrix_sta)
                
                # plot migration profiles
                self.coherence_plot(self.output_path+'/'+event, corrmatrix_ch, tobj.x, tobj.y, tobj.z, i)
            
                # output theoretical P- and S-wave arrivaltimes
                fname_ch = cmfilename + '_trial{}.phs'.format(i)
                self.write_phasetime(sobj.stations, event_t0_ch, tp_modse, ts_modse, iloc_ch[0], fname_ch)
                
                if cmax_ch > cmax_pre:
                    event_t0s_final_ch = copy.deepcopy(event_t0s_ch)
                    cmax_pre = copy.deepcopy(cmax_ch)

                #totale
                evtpmin = num.amin(tp_modse[iloc[0],:])
                event_t0 = sobj.dtime_max + datetime.timedelta(seconds=iloc[1]*sobj.deltat_sta) - datetime.timedelta(seconds=evtpmin)  # event origin time
                event_t0s = (event_t0).isoformat()
                # corrmatrix is the stacking matrix, in 1D format but can be 
                # reformat to 3D format, each point saves the maximum stacking 
                # value during this calculation time period
                cmax = num.max(corrmatrix)
                #corrmatrix = num.reshape(corrmatrix,(tobj.nx,tobj.ny,tobj.nz))
                (ixloc, iyloc, izloc) = num.unravel_index(iloc[0],(tobj.nx,tobj.nx,tobj.nz))
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
                self.write_phasetime(sobj.stations, event_t0, tp_modse, ts_modse, iloc[0], fname)
                
                if cmax_sta > cmax_pre:
                    event_t0s_final_tot = copy.deepcopy(event_t0s)
                    cmax_pre = copy.deepcopy(cmax)


            
            self.catalogue_creation(event, event_t0s_final_sta, tobj.lat0, tobj.lon0, ntrial, corrmatrix_sta)
            print('Location process completed (STA)!!!')
            self.catalogue_creation(event, event_t0s_final_ch, tobj.lat0, tobj.lon0, ntrial, corrmatrix_ch)
            print('Location process completed (FIBER)!!!')
            self.catalogue_creation(event, event_t0s_final_tot, tobj.lat0, tobj.lon0, ntrial, corrmatrix)
            print('Location process completed (HYBRID)!!!')


            
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

'''