import numpy as num

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

        # Convert lists to NumPy arrays
        self.lon_stations = num.array(tobj.lon_stations, dtype=float)
        self.lat_stations = num.array(tobj.lat_stations, dtype=float)
        self.depth_stations = num.array(tobj.depth_stations, dtype=float)
        self.lon_channels = num.array(tobj.lon_channels, dtype=float)
        self.lat_channels = num.array(tobj.lat_channels, dtype=float)
        self.depth_channels = num.array(tobj.depth_channels, dtype=float)

    def location_domain(self):
        extx_sub = num.abs(num.minimum(self.lon_stations.min(), self.lon_channels.min()) - 
                           num.maximum(self.lon_stations.max(), self.lon_channels.max()))
        exty_sub = num.abs(num.minimum(self.lat_stations.min(), self.lat_channels.min()) - 
                           num.maximum(self.lat_stations.max(), self.lat_channels.max()))
        extz_sub = num.abs(num.minimum(self.depth_stations.min(), self.depth_channels.min()) - 
                           num.maximum(self.depth_stations.max(), self.depth_channels.max()))

        extx_tt = self.nx * self.dx
        extz_tt = self.nz * self.dz

        extx = extx_tt / num.sqrt(2)
        exty = extx_tt / num.sqrt(2)
        extz = extz_tt

        diff_subx = extx - extx_sub
        diff_suby = exty - exty_sub
        diff_subz = extz - extz_sub

        self.lon_stations_rel = self.lon_stations + diff_subx / 2 - self.lon_stations.min()
        self.lat_stations_rel = self.lat_stations + diff_suby / 2 - self.lat_stations.min()
        self.depth_stations_rel = self.depth_stations + diff_subz / 2 - self.depth_stations.min()
        self.lon_channels_rel = self.lon_channels + diff_subx / 2 - self.lon_channels.min()
        self.lat_channels_rel = self.lat_channels + diff_suby / 2 - self.lat_channels.min()
        self.depth_channels_rel = self.depth_channels + diff_subz / 2 - self.depth_channels.min()

    def get_closest_travel_time(self, horizontal_distance, depth_value, tt_table):
        tt_2d = tt_table.reshape(self.nx, self.nz)
        horiz_dist = num.linspace(0, self.nx * self.dx, self.nx)
        depth = num.linspace(0, self.nz * self.dz, self.nz)
        closest_x_idx = num.argmin(num.abs(horiz_dist - horizontal_distance))
        closest_z_idx = num.argmin(num.abs(depth - depth_value))
        travel_time = tt_2d[closest_x_idx, closest_z_idx]
        return closest_x_idx, closest_z_idx, travel_time

    def stacking(self, lon_sensors, lat_sensors, depth_sensors, itp, its, stalta_p, stalta_s):
        nxyz = self.nx * self.nx * self.nz
        nsta = len(stalta_p[:, 1])
        nsamples = len(stalta_p[1, :])
        
        corrmatrix = num.zeros(nxyz)
        corrmax = -1.0
        iloc, itime = 0, 0

        progress_step = nxyz // 300000  # 2% progress intervals

        for i in range(0, nxyz, 1000):  # Perform computations every 1000 steps
            if i % progress_step == 0:
                print(f"Processing: {100 * i / nxyz:.6f}% completed")

            stkmax = -1.0
            kmax = 0
            for k in range(nsamples):
                stk0p = 0.0
                stk0s = 0.0
                for j in range(nsta):
                    current_horizontal_distance = num.sqrt(lon_sensors[j]**2 + lat_sensors[j]**2)
                    current_depth = depth_sensors[j]
                    _, _, ttp_val = self.get_closest_travel_time(current_horizontal_distance, current_depth, itp)
                    _, _, tts_val = self.get_closest_travel_time(current_horizontal_distance, current_depth, its)
                    
                    ip = int(ttp_val + k)
                    is_ = int(tts_val + k)

                    if (ip < nsamples) & (is_ < nsamples):
                        stk0p += stalta_p[j, ip]
                        stk0s += stalta_s[j, is_]

                if stk0p * stk0s > stkmax:
                    stkmax = stk0p * stk0s
                    kmax = k

            corrmatrix[i] = num.sqrt(stkmax) / nsta
            if corrmatrix[i] > corrmax:
                corrmax = corrmatrix[i]
                iloc = i
                itime = kmax

        return (iloc, itime), corrmatrix

    def locate_event(self):
        self.location_domain()
        print("Location domain set up.")

        iloc_sta, corrmatrix_sta = self.stacking(self.lon_stations_rel, self.lat_stations_rel, self.depth_stations_rel, self.ttp, self.tts, self.obsp_sta, self.obss_sta)
        iloc_ch, corrmatrix_ch = self.stacking(self.lon_channels_rel, self.lat_channels_rel, self.depth_channels_rel, self.ttp, self.tts, self.obsp_ch, self.obss_ch)

        iloc = (iloc_sta[0] + iloc_ch[0]) / 2
        itime = (iloc_sta[1] + iloc_ch[1]) / 2
        corrmatrix = (corrmatrix_sta + corrmatrix_ch) / 2

        print("Event located successfully!")

        return iloc, itime, corrmatrix_sta.reshape(self.nx, self.nx, self.nz), corrmatrix_ch.reshape(self.nx, self.nx, self.nz), corrmatrix.reshape((self.nx, self.nx, self.nz))
