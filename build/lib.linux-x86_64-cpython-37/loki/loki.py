# %%
import os
import glob
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
import location_t0_sta                     
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
            tshortp = num.linspace(tshortp_min, tshortp_max, ntrial)
            tshorts = num.linspace(tshorts_min, tshorts_max, ntrial)
            tshortp_fiber = num.linspace(tshortp_min_fiber, tshortp_max_fiber, ntrial)
            tshorts_fiber = num.linspace(tshorts_min_fiber, tshorts_max_fiber, ntrial)
        else:
            STALTA = False
            ntrial = 1

        npr = inputs['npr']
        model = inputs['model']
        freq = inputs.get('freq', None)
        opsf = inputs.get('opsf', False)

        tobj = traveltimes.Traveltimes(self.db_path, self.hdr_filename, self.geometry_filename_fiber, self.geometry_filename_stat)
        tp = tobj.load_traveltimes('P', model, precision)
        ts = tobj.load_traveltimes('S', model, precision)
        tobj.load_station_info()

        #tobj.x = tobj.x[0:len(tobj.x)-1]. #if nz different from nx,ny
        #tobj.y = tobj.y[0:len(tobj.y)-1]



        for event_path in self.data_tree:
            self.subdata_path = event_path    #this is stations, diber
            rel_path = os.path.relpath(event_path, self.data_path)  # e.g., 'closest_event' or 'ww2_bomb'
            output_dir = os.path.join(self.output_path, rel_path)
            os.makedirs(output_dir, exist_ok=True)
            last_folder = os.path.basename(self.subdata_path)

            print(f"Processing event: {event_path}")
            print(f"Processing network type: {last_folder}")

            if "hybrid" in last_folder:
                label = "hybrid"
                st = read(os.path.join(event_path, "*"))
                # get last-letter channel codes present (E,N,Z,...)
                components = set(tr.stats.channel[-1] for tr in st)
                print("components in folder:", components)
                # preserve preferred ordering E, N, Z
                pref = ['E', 'N', 'Z']
                comps = [c for c in pref if c in components]
                # if there are other component letters (rare), append them
                other = sorted([c for c in components if c not in pref])
                comps += other
                print("final components used:", comps)

            elif last_folder == "hybrid_strain_to_vel":
                label = "hybrid_strain_to_vel"
                st = read(os.path.join(event_path, "*"))
                components = set(tr.stats.channel[-1] for tr in st)
                print("components in folder:", components)
                pref = ['E', 'N', 'Z']
                comps = [c for c in pref if c in components]
                other = sorted([c for c in components if c not in pref])
                comps += other
                

            elif last_folder == "stations":
                label = "stations"
                comps = ['E', 'N', 'Z']

            elif "fiber" in last_folder:
                label = "fiber"
                comps = ['E']

            elif last_folder == "hybrid_gilbert":
                label = "hybrid_gilbert"
                st = read(os.path.join(event_path, "*"))
                components = set(tr.stats.channel[-1] for tr in st)
                print("components in folder:", components)
                pref = ['E', 'N', 'Z']
                comps = [c for c in pref if c in components]
                other = sorted([c for c in components if c not in pref])
                comps += other

            else:
                continue

            print('waveforms class')

            wobj = waveforms.Waveforms(self.subdata_path, extension_sta="*", comps=comps, freq=None)

            print('stacktraces class')

            sobj = stacktraces.Stacktraces(tobj, wobj, **inputs)
            event = event_path.split('/')[-1]
            #print('event', event)
            #output_dir = os.path.join(self.output_path, event)
            #os.makedirs(output_dir, exist_ok=True)

            print(tobj.nxz)

            tp_modse = num.ascontiguousarray(tp['HM00'].reshape(tobj.nxz, 1), dtype=num.float64)
            ts_modse = num.ascontiguousarray(ts['HM00'].reshape(tobj.nxz, 1), dtype=num.float64)
            tp_mod_sta, ts_mod_sta = tt_processing.tt_f2i(sobj.deltat, tp_modse, ts_modse, npr)

            cmax_pre = -1.0
            event_t0s_final = None
            corrmatrix_best = None


            # Clean previous output files only in this event's output folder
            for f in glob.glob(os.path.join(output_dir, '*.loc')):
                os.remove(f)
            for f in glob.glob(os.path.join(output_dir, '*.npy')):
                os.remove(f)

            catalogue_path = os.path.join(output_dir, '*.catalogue')
            if os.path.exists(catalogue_path):
                os.remove(catalogue_path)

            png_path = os.path.join(output_dir, '*.png')
            if os.path.exists(png_path):
                os.remove(png_path)

            for i in range(ntrial):
                if STALTA:

                    print('hey!')
                    if label == 'fiber':
                        nshort_p_sta = int(tshortp_fiber[i] // sobj.deltat)
                        nshort_s_sta = int(tshorts_fiber[i] // sobj.deltat)
                    else:
                        nshort_p_sta = int(tshortp[i] // sobj.deltat)
                        nshort_s_sta = int(tshorts[i] // sobj.deltat)
                        
                    obs_dataP_sta, obs_dataS_sta = sobj.loc_stalta(nshort_p_sta, nshort_s_sta, slrat, norm=1)
                else:
                    obs_dataP_sta = sobj.obs_dataV_sta
                    obs_dataS_sta = sobj.obs_dataH_sta

                x_stations, y_stations, z_stations = [], [], []
                for sta in sobj.stations:
                    lon, lat, depth = tobj.stations_coordinates.get(sta, (None, None, None))
                    x_stations.append(lon)
                    y_stations.append(lat)
                    z_stations.append(depth)

                
                import numpy as np
                import matplotlib.pyplot as plt

                # ----------------------------------------
                # ZERO CHECKING
                # ----------------------------------------

                print("🔍 Checking for zeros in obs_dataP_sta ...")

                # Does the entire array contain at least one zero?
                has_any_zero = np.any(obs_dataS_sta == 0)
                print("Does the array contain at least one zero? ->", has_any_zero)

                # Count zeros per row
                zero_counts = np.sum(obs_dataS_sta == 0, axis=1)

                for i, count in enumerate(zero_counts):
                    print(f"Row {i}: {count} zeros")

                # Detect fully zero rows
                fully_zero_rows = np.where(zero_counts == obs_dataP_sta.shape[1])[0]
                print("Rows entirely zero:", fully_zero_rows.tolist())

                # ----------------------------------------
                # PLOTTING SECTION
                # ----------------------------------------

                fig, axes = plt.subplots(
                    obs_dataS_sta.shape[0],
                    1,
                    figsize=(10, 2 * obs_dataS_sta.shape[0]),
                    sharex=True
                )


                if obs_dataS_sta.shape[0] <100:  
                    for i in range(obs_dataS_sta.shape[0]):
                        axes[i].plot(obs_dataS_sta[i])
                        axes[i].set_ylabel(f"Row {i}")

                        # Mark subplot in red title if it contains zeros
                        if zero_counts[i] > 0:
                            axes[i].set_title(f"{zero_counts[i]} zeros", color="red", fontsize=8)

                    axes[-1].set_xlabel("Index")
                    plt.tight_layout()

                    # ----------------------------------------
                    # SAVE TO PDF
                    # ----------------------------------------

                    pdf_path = "/home/emanuele/data/emanuele/loki-das/bedretto/output_m05/obs_dataP_sta_plots.pdf"
                    print("Saving figure to:", pdf_path)

                    plt.savefig(pdf_path)
                    plt.close()

                    print("✅ PDF saved successfully.")


                #print(z_stations)

                tp_mod_sta = num.ascontiguousarray(tp_mod_sta.reshape(tobj.nx, tobj.nz), dtype=num.int32)
                ts_mod_sta = num.ascontiguousarray(ts_mod_sta.reshape(tobj.nx, tobj.nz), dtype=num.int32)
                x_stations = num.ascontiguousarray(x_stations, dtype=num.float64)
                y_stations = num.ascontiguousarray(y_stations, dtype=num.float64)
                z_stations = num.ascontiguousarray(z_stations, dtype=num.float64)
                obs_dataP_sta = num.ascontiguousarray(obs_dataP_sta, dtype=num.float64)
                obs_dataS_sta = num.ascontiguousarray(obs_dataS_sta, dtype=num.float64)

                print('aaaa',obs_dataP_sta.shape)




                print('bbbb', len(tp_mod_sta))

                print


                print(tobj.nx, tobj.nz)
                print('Starting stacking')



                if last_folder == "stations":  

                    print('now only stations and thus stacking P*S')                                      

                    iloctime, corrmatrix = location_t0.stacking(tp_mod_sta, ts_mod_sta,
                                                                x_stations, y_stations, z_stations,
                                                                tobj.x, tobj.y, tobj.z,
                                                                obs_dataP_sta, obs_dataS_sta, npr)

                elif last_folder == "hybrid":  

                    print('now hybrid thus stacking P * S')                                      

                    iloctime, corrmatrix = location_t0.stacking(tp_mod_sta, ts_mod_sta,
                                                                x_stations, y_stations, z_stations,
                                                                tobj.x, tobj.y, tobj.z,
                                                                obs_dataP_sta, obs_dataS_sta, npr)                    

                elif last_folder == "hybrid_gilbert":  

                    print('now hybrid_gilbert thus stacking P * S')                                      

                    iloctime, corrmatrix = location_t0.stacking(tp_mod_sta, ts_mod_sta,
                                                                x_stations, y_stations, z_stations,
                                                                tobj.x, tobj.y, tobj.z,
                                                                obs_dataP_sta, obs_dataS_sta, npr)    

                else:

                    iloctime, corrmatrix = location_t0.stacking(tp_mod_sta, ts_mod_sta,
                                                                x_stations, y_stations, z_stations,
                                                                tobj.x, tobj.y, tobj.z,
                                                                obs_dataP_sta, obs_dataS_sta, npr)



                evtpmin = num.amin(tp_modse)
                event_t0 = sobj.dtime_max + datetime.timedelta(seconds=iloctime[3] * sobj.deltat) - datetime.timedelta(seconds=evtpmin)
                event_t0s = event_t0.isoformat()
                cmax = num.max(corrmatrix)

                # Use consistent .loc filenames: event name only, no label appended
                cmfilename = f"{output_dir}/{event}" if ntrial > 1 else f"{output_dir}/{event}"

                #print('cmfilename', cmfilename)

                mode = 'a' if i > 0 else 'w'
                with open(cmfilename + '.loc', mode) as out_file:
                    if STALTA:
                        out_file.write(f"{i} {x_stations[0]} {y_stations[0]} {z_stations[0]} {cmax} {nshort_p_sta} {nshort_s_sta} {slrat}\n")
                    else:
                        out_file.write(f"{i} {x_stations[0]} {y_stations[0]} {z_stations[0]} {cmax}\n")

                num.save(f"{output_dir}/corrmatrix_trial_{i}_{label}.npy", corrmatrix)

                if cmax > cmax_pre:
                    event_t0s_final = copy.deepcopy(event_t0s)
                    corrmatrix_best = copy.deepcopy(corrmatrix)
                    cmax_pre = cmax
                    best_trial_index = i  # <--- Add this line

            if event_t0s_final and corrmatrix_best is not None:
                print(tobj.x, tobj.y, tobj.z)
                try:
                    self.catalogue_creation(cmfilename, event, event_t0s_final, tobj.lat0, tobj.lon0, ntrial, corrmatrix_best, label)
                    # ---- Coherence Plot ----
                    self.coherence_plot(
                        cmfilename, event_path,  # event_path
                        corrmatrix_best,tobj,                        # coherence matrix of best trial
                        tobj.x, tobj.y, tobj.z,                 # grid axes
                        best_trial_index,                       # trial number
                        normalization=False                      # optional
                    )
                except Exception as e:
                    print(f"Catalogue creation or coherence plot failed for event {event} - {label}: {e}")

#methods
    def catalogue_creation(self, cmfilename, event, event_t0s, lat0, lon0, ntrial, corrmatrix, refell=23):
        latref = lat0
        lonref = lon0
        eleref = 1
        origin = latlon2cart.Coordinates(latref, lonref, eleref)

        if (ntrial > 1):
            ev_file = cmfilename + '.loc'
            data = num.loadtxt(ev_file)
            w = num.sum(data[:, 4])
            xb = ((num.dot(data[:, 1], data[:, 4]) / w) * 1000)
            yb = ((num.dot(data[:, 2], data[:, 4]) / w) * 1000)
            print(f"xb: {xb}, yb: {yb}, eleref: {eleref}")
            try:
                result = origin.cart2geo(xb, yb, eleref)
                if isinstance(result, (tuple, list)) and len(result) == 3:
                    late, lone, elev = result
                else:
                    raise ValueError(f"Invalid cart2geo output: {result}")
            except Exception as e:
                print(f"[ERROR] cart2geo failed for event {event} ({event_t0s}): {e}")
                late, lone, elev = 0.0, 0.0, 0.0  # Or consider skipping this catalogue entry

            zb = num.dot(data[:, 3], data[:, 4]) / w  # depth in km
            cb = num.mean(data[:, 4])  # mean coherence
            cmax = num.max(data[:, 4])  # max coherence
            merr = num.vstack((data[:, 1], data[:, 2], data[:, 3]))
            err = num.cov(merr)
            errmax = num.sqrt(num.max(num.linalg.eigvals(err)))
        else:
            ev_file = cmfilename + '.loc'
            data = num.loadtxt(ev_file)
            if data.ndim > 1:
               data = data[0]  # Get first row
            xb = data[1] * 1000
            yb = data[2] * 1000
            zb = data[3]  # depth in km
            late, lone, elev = origin.cart2geo(xb, yb, eleref)
            #print('late, lone, elev', late, lone, elev)
            cmax = data[4]

            # Normalize corrmatrix: min -> 1, max -> 2
            n1 = 1.0
            n2 = 2.0
            dmax = num.amax(corrmatrix)
            dmin = num.amin(corrmatrix)
            k = (n2 - n1) / (dmax - dmin)
            b = (dmax * n1 - dmin * n2) / (dmax - dmin)
            corrmatrix = k * corrmatrix + b

            errmax = num.std(corrmatrix)
            cb = num.median(corrmatrix)

        # Ensure output per subfolder
        catalogue_path = cmfilename + '.catalogue'
        with open(catalogue_path, 'a') as f:
            f.write(f"{event_t0s}    {late}   {lone}   {zb}   {errmax}   {cb}   {cmax}\n")


    def coherence_plot(self, cmfilename, event_path, corrmatrix, tobj, xax, yax, zax, itrial, normalization=False):
        
        nx = tobj.nx
        ny = tobj.nx
        nz = tobj.nz

         

        corrmatrix = corrmatrix.reshape((nx, ny, nz))
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

        plt.savefig(cmfilename + 'Coherence_matrix_'+str(itrial)+'.png', dpi =300)
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


# %%
