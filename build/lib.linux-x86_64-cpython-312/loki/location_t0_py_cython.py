import numpy as np
import logging
import time
from concurrent.futures import ProcessPoolExecutor
from numba import jit, prange
from scipy.spatial import cKDTree
from waveform_stacking import compute_sensor_stack  # Import the Cythonized function

class WaveformStacking:
    def __init__(self, tobj, nproc, ttp, tts, obsp_sta, obss_sta, obsp_ch, obss_ch):
        required_attrs = ["nx", "nz", "dx", "dz", "lon_stations", "lat_stations", 
                          "depth_stations", "lon_channels", "lat_channels", "depth_channels"]
        for attr in required_attrs:
            if not hasattr(tobj, attr):
                raise AttributeError(f"Missing required attribute: {attr} in tobj")

        self.nproc = nproc
        self.nx = tobj.nx
        self.nz = tobj.nz
        self.dx = tobj.dx
        self.dz = tobj.dz
        self.ttp = ttp
        self.tts = tts
        self.obsp_sta = obsp_sta
        self.obss_sta = obss_sta
        self.obsp_ch = obsp_ch
        self.obss_ch = obss_ch

        # Adapt location of sensors to the dx-dz of the domain
        self.lon_stations = np.array(tobj.lon_stations, dtype=float) * 0.0001
        self.lat_stations = np.array(tobj.lat_stations, dtype=float) * 0.0001
        self.depth_stations = np.array(tobj.depth_stations, dtype=float) * 0.0001
        self.lon_channels = np.array(tobj.lon_channels, dtype=float) * 0.0001
        self.lat_channels = np.array(tobj.lat_channels, dtype=float) * 0.0001
        self.depth_channels = np.array(tobj.depth_channels, dtype=float) * 0.0001

        # Prepare KD-tree for nearest neighbor search
        self.station_tree = cKDTree(np.column_stack([self.lon_stations, self.lat_stations, self.depth_stations]))
        self.channel_tree = cKDTree(np.column_stack([self.lon_channels, self.lat_channels, self.depth_channels]))

    def location_domain(self):
        # Define the 3D grid location domain, starting from the 2D traveltime table 
        extx_sub = np.abs(np.minimum(self.lon_stations.min(), self.lon_channels.min()) - 
                           np.maximum(self.lon_stations.max(), self.lon_channels.max()))
        exty_sub = np.abs(np.minimum(self.lat_stations.min(), self.lat_channels.min()) - 
                           np.maximum(self.lat_stations.max(), self.lat_channels.max()))
        extz_sub = np.abs(np.minimum(self.depth_stations.min(), self.depth_channels.min()) - 
                           np.maximum(self.depth_stations.max(), self.depth_channels.max()))

        extx_tt = (self.nx * self.dx)
        extz_tt = (self.nz * self.dz)

        extx = extx_tt / np.sqrt(2)
        exty = extx_tt / np.sqrt(2)
        extz = extz_tt

        diff_subx = extx - extx_sub
        diff_suby = exty - exty_sub
        diff_subz = extz - extz_sub

        print('diff', extx_sub, exty_sub, extz_sub, extx_tt, extz_tt, diff_subx, diff_suby, diff_subz)

        self.lon_stations_rel = (self.lon_stations + diff_subx / 2) - self.lon_stations.min()
        self.lat_stations_rel = self.lat_stations + diff_suby / 2 - self.lat_stations.min()
        self.depth_stations_rel = self.depth_stations + diff_subz / 2 - self.depth_stations.min()
        self.lon_channels_rel = self.lon_channels + diff_subx / 2 - self.lon_channels.min()
        self.lat_channels_rel = self.lat_channels + diff_suby / 2 - self.lat_channels.min()
        self.depth_channels_rel = self.depth_channels + diff_subz / 2 - self.depth_channels.min()

    def get_closest_travel_time(self, horizontal_distance, depth_value, tt_table):
        tt_2d = tt_table.reshape(self.nx, self.nz)
        horiz_dist = np.linspace(0, self.nx * self.dx, self.nx)
        depth = np.linspace(0, self.nz * self.dz, self.nz)

        closest_x_idx = np.argmin(np.abs(horiz_dist - horizontal_distance))
        closest_z_idx = np.argmin(np.abs(depth - depth_value))

        return closest_x_idx, closest_z_idx, tt_2d[closest_x_idx, closest_z_idx]

    def stacking(self, lon_sensors, lat_sensors, depth_sensors, itp, its, stalta_p, stalta_s):
        logging.info("Stacking function called.")
        start_time = time.time()  # Start timer

        nxyz = self.nx * self.nx * self.nz
        nsta = len(stalta_p[:, 1])
        nsamples = len(stalta_p[1, :])

        corrmatrix = np.zeros(nxyz)
        corrmax = -1.0
        iloc, itime = 0, 0

        progress_step = nxyz // 300000  # 2% progress intervals

        # Precompute the sensor distances to avoid repeated calculation
        sensor_distances = np.sqrt(np.array(lon_sensors)**2 + np.array(lat_sensors)**2)

        with ProcessPoolExecutor(max_workers=self.nproc) as executor:
            futures = []
            for i in prange(0, nxyz):  # Parallelize using prange
                if i % progress_step == 0:
                    print(f"Processing: {100 * i / nxyz:.6f}% completed")

                for j in prange(nsta):  # Parallelize inner loop
                    current_horizontal_distance = sensor_distances[j]
                    current_depth = depth_sensors[j]

                    _, _, ttp_val = self.get_closest_travel_time(current_horizontal_distance, current_depth, itp)
                    _, _, tts_val = self.get_closest_travel_time(current_horizontal_distance, current_depth, its)

                    ttp_val = int(ttp_val)
                    tts_val = int(tts_val)

                    futures.append(executor.submit(compute_sensor_stack, j, stalta_p, stalta_s, nsamples, ttp_val, tts_val))

            # Collect the results and update corrmatrix
            for future in futures:
                stk0p, stk0s = future.result()
                if stk0p * stk0s > corrmax:
                    corrmax = stk0p * stk0s
                    iloc = i  # Save the best matching index

        end_time = time.time()
        logging.info(f"Stacking function completed in {end_time - start_time:.2f} seconds.")

        return (iloc, itime), corrmatrix

    def locate_event(self):
        self.location_domain()
        print("Location domain set up.")

        # corrmatrix stations-fiber
        iloc_sta, corrmatrix_sta = self.stacking(self.lon_stations_rel, self.lat_stations_rel, self.depth_stations_rel, self.ttp, self.tts, self.obsp_sta, self.obss_sta)
        iloc_ch, corrmatrix_ch = self.stacking(self.lon_channels_rel, self.lat_channels_rel, self.depth_channels_rel, self.ttp, self.tts, self.obsp_ch, self.obss_ch)

        iloc = (iloc_sta[0] + iloc_ch[0]) / 2
        itime = (iloc_sta[1] + iloc_ch[1]) / 2

        # Hybrid corrmatrix
        corrmatrix = (corrmatrix_sta + corrmatrix_ch) / 2

        print("Event located successfully!")

        return iloc_sta, iloc_ch, iloc, itime, corrmatrix_sta.reshape(self.nx, self.nx, self.nz), corrmatrix_ch.reshape(self.nx, self.nx, self.nz), corrmatrix.reshape((self.nx, self.nx, self.nz))
