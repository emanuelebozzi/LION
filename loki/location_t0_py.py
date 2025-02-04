import numpy as num

class WaveformStacking:
    def __init__(self, tobj, nproc, nx, nz, dx, dz):
        self.tobj = tobj
        self.nproc = nproc
        self.nx = nx
        self.nz = nz
        self.dx = dx
        self.dz = dz

    def location_domain(self):
        tobj = self.tobj
        extx_sub = num.abs(num.min([tobj.lon_stations.min(), tobj.lon_channels.min()]) - 
                          num.max([tobj.lon_stations.max(), tobj.lon_channels.max()]))
        exty_sub = num.abs(num.min([tobj.lat_stations.min(), tobj.lat_channels.min()]) - 
                          num.max([tobj.lat_stations.max(), tobj.lat_channels.max()]))
        extz_sub = num.abs(num.min([tobj.depth_stations.min(), tobj.depth_channels.min()]) - 
                          num.max([tobj.depth_stations.max(), tobj.depth_channels.max()]))

        extx_tt = self.nx * self.dx
        extz_tt = self.nz * self.dz

        extx = extx_tt / num.sqrt(2)
        exty = extx_tt / num.sqrt(2)
        extz = extz_tt

        diff_subx = extx - extx_sub
        diff_suby = exty - exty_sub
        diff_subz = extz - extz_sub

        self.lon_stations_rel = tobj.lon_stations + diff_subx / 2 - tobj.lon_stations.min()
        self.lat_stations_rel = tobj.lat_stations + diff_suby / 2 - tobj.lat_stations.min()
        self.depth_stations_rel = tobj.depth_stations + diff_subz / 2 - tobj.depth_stations.min()
        self.lon_channels_rel = tobj.lon_channels + diff_subx / 2 - tobj.lon_channels.min()
        self.lat_channels_rel = tobj.lat_channels + diff_suby / 2 - tobj.lat_channels.min()
        self.depth_channels_rel = tobj.depth_channels + diff_subz / 2 - tobj.depth_channels.min()

    def get_closest_travel_time(self, horizontal_distance, depth_value, tt_table):
        tt_2d = tt_table.reshape(self.nx, self.nz)
        horiz_dist = num.linspace(0, self.nx * self.dx, self.nx)
        depth = num.linspace(0, self.nz * self.dz, self.nz)
        closest_x_idx = num.argmin(num.abs(horiz_dist - horizontal_distance))
        closest_z_idx = num.argmin(num.abs(depth - depth_value))
        travel_time = tt_2d[closest_x_idx, closest_z_idx]
        return closest_x_idx, closest_z_idx, travel_time

    def stacking(self, nx, nz, itp, its, stalta_p, stalta_s, nproc):
        nxyz = nx * nx * nz
        nsta = len(self.lon_stations_rel)
        nsamples = len(itp)

        corrmatrix = num.zeros(nxyz)
        corrmax = -1.0
        iloc, itime = 0, 0

        for i in range(nxyz):
            stkmax = -1.0
            kmax = 0
            for k in range(nsamples):
                stk0p = 0.0
                stk0s = 0.0
                for j in range(nsta):
                    current_horizontal_distance = num.sqrt(self.lon_stations_rel[j]**2 + self.lat_stations_rel[j]**2)
                    current_depth = self.depth_stations_rel[j]
                    _, _, ttp_val = self.get_closest_travel_time(current_horizontal_distance, current_depth, itp)
                    _, _, tts_val = self.get_closest_travel_time(current_horizontal_distance, current_depth, its)

                    ip = int(ttp_val + k)
                    is_ = int(tts_val + k)
                    if is_ < nsamples:
                        stk0p += stalta_p[j, ip]
                        stk0s += stalta_s[j, is_]
                    else:
                        stk0p += 0.0
                        stk0s += 0.0

                if stk0p * stk0s > stkmax:
                    stkmax = stk0p * stk0s
                    kmax = k

            corrmatrix[i] = num.sqrt(stkmax) / nsta
            if corrmatrix[i] > corrmax:
                corrmax = corrmatrix[i]
                iloc = i
                itime = kmax

        return (iloc, itime), corrmatrix

    def locate_event(self, ttp, tts, obsp_sta, obss_sta, obsp_ch, obss_ch):
        self.location_domain()

        iloc_sta, corrmatrix_sta = self.stacking(self.nx, self.nz, ttp, tts, obsp_sta, obss_sta, self.nproc)
        iloc_ch, corrmatrix_ch = self.stacking(self.nx, self.nz, ttp, tts, obsp_ch, obss_ch, self.nproc)

        iloc = num.mean([iloc_sta[0], iloc_ch[0]])
        itime = num.mean([iloc_sta[1], iloc_ch[1]])
        corrmatrix = num.mean([corrmatrix_sta, corrmatrix_ch], axis=0)

        return iloc, itime, corrmatrix.reshape((self.nx, self.nx, self.nz))

# Example usage
# nx, nz, dx, dz = 50, 50, 1.0, 1.0
# tobj = YourObjectWithAttributes()  # Replace with actual object
# travel_time_flat = num.random.random(nx * nz)
# stacking = WaveformStacking(tobj, nproc=4, nx=nx, nz=nz, dx=dx, dz=dz)
# iloc, itime, corrmatrix = stacking.locate_event(travel_time_flat, travel_time_flat, obsp_sta, obss_sta, obsp_ch, obss_ch)
# print(f"Event located at index {iloc} with time {itime}")
