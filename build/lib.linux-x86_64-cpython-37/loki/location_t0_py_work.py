import numpy as num
import logging 
import time 

#logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s")


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

        #adapt location of sensors to the dx-dz of the domain (to be fixed sistematically)
        self.lon_stations = num.array(tobj.lon_stations, dtype=float) * 0.0001
        self.lat_stations = num.array(tobj.lat_stations, dtype=float) *0.0001
        self.depth_stations = num.array(tobj.depth_stations, dtype=float)*0.0001
        self.lon_channels = num.array(tobj.lon_channels, dtype=float)*0.0001
        self.lat_channels = num.array(tobj.lat_channels, dtype=float)*0.0001
        self.depth_channels = num.array(tobj.depth_channels, dtype=float)*0.0001


    def location_domain(self):

        '''
        This function define the 3D grid location domain, starting from the 2D traveltime table 
        spanning the diagonal of the domain
        '''

        #define the current extension of the sensors on x,y,z (necessary for relative location)
        extx_sub = num.abs(num.minimum(self.lon_stations.min(), self.lon_channels.min()) - 
                           num.maximum(self.lon_stations.max(), self.lon_channels.max()))
        exty_sub = num.abs(num.minimum(self.lat_stations.min(), self.lat_channels.min()) - 
                           num.maximum(self.lat_stations.max(), self.lat_channels.max()))
        extz_sub = num.abs(num.minimum(self.depth_stations.min(), self.depth_channels.min()) - 
                           num.maximum(self.depth_stations.max(), self.depth_channels.max()))

        #define the current extension of the traveltime domain 
        extx_tt = (self.nx * self.dx)
        extz_tt = (self.nz * self.dz)
        
        #define the definitive dimension of the location domain 
        extx = extx_tt / num.sqrt(2)
        exty = extx_tt / num.sqrt(2)
        extz = extz_tt

        diff_subx = extx - extx_sub
        diff_suby = exty - exty_sub
        diff_subz = extz - extz_sub


        print('diff',extx_sub, exty_sub, extz_sub, extx_tt, extz_tt, diff_subx, diff_suby, diff_subz)
        
        #define all the relative location of the sensors (to the 0,0,0)
        #position the stations at the center of the investigated domain 

        self.lon_stations_rel = (self.lon_stations + diff_subx / 2) - self.lon_stations.min()
        self.lat_stations_rel = self.lat_stations + diff_suby / 2 - self.lat_stations.min()
        self.depth_stations_rel = self.depth_stations + diff_subz / 2 - self.depth_stations.min()
        self.lon_channels_rel = self.lon_channels + diff_subx / 2 - self.lon_channels.min()
        self.lat_channels_rel = self.lat_channels + diff_suby / 2 - self.lat_channels.min()
        self.depth_channels_rel = self.depth_channels + diff_subz / 2 - self.depth_channels.min()



    def get_closest_travel_time(self, horizontal_distance, depth_value, tt_table):
        start_time = time.time()  # Start timing
        tt_2d = tt_table.reshape(self.nx, self.nz)

        # Log the input values
        logging.debug(f"Computing travel time for Distance={horizontal_distance:.4f}, Depth={depth_value:.4f}")

        horiz_dist = num.linspace(0, self.nx * self.dx, self.nx)
        depth = num.linspace(0, self.nz * self.dz, self.nz)

        closest_x_idx = num.argmin(num.abs(horiz_dist - horizontal_distance))
        closest_z_idx = num.argmin(num.abs(depth - depth_value))

        travel_time = tt_2d[closest_x_idx, closest_z_idx]

        end_time = time.time()  # End timing
        elapsed_time = end_time - start_time

        # Log execution time
        logging.debug(f"Travel time lookup took {elapsed_time:.6f} seconds")

        return closest_x_idx, closest_z_idx, travel_time
    

    def stacking(self, lon_sensors, lat_sensors, depth_sensors, itp, its, stalta_p, stalta_s):
        """Function stacking energy along predicted travel times."""
        
        logging.info("Stacking function called.")
        start_time = time.time()  # Start timer

        nxyz = self.nx * self.nx * self.nz
        nsta = len(stalta_p[:, 1])
        nsamples = len(stalta_p[1, :])

        corrmatrix = num.zeros(nxyz)
        corrmax = -1.0
        iloc, itime = 0, 0

        progress_step = nxyz // 300000  # 2% progress intervals

        # **Monitor outer loop**
        outer_start = time.time()

        # Precompute the sensor distances to avoid repeated calculation
        sensor_distances = num.sqrt(num.array(lon_sensors)**2 + num.array(lat_sensors)**2)

        for i in range(0, nxyz):  # Limit loop for debugging
            loop_start = time.time()  # Start timer for outer loop iteration
            if i % progress_step == 0:
                print(f"Processing: {100 * i / nxyz:.6f}% completed")

            stkmax = -1.0
            kmax = 0

            # **Monitor sensor loop**
            sensor_loop_start = time.time()

            for j in range(nsta):
                # Precompute travel times for each sensor
                current_horizontal_distance = sensor_distances[j]
                current_depth = depth_sensors[j]

                ttp_start = time.time()  # Start timing the travel time lookup
                _, _, ttp_val = self.get_closest_travel_time(current_horizontal_distance, current_depth, itp)
                ttp_end = time.time()  # End timing

                tts_start = time.time()  # Start timing the travel time lookup
                _, _, tts_val = self.get_closest_travel_time(current_horizontal_distance, current_depth, its)
                tts_end = time.time()  # End timing

                ttp_val = int(ttp_val)
                tts_val = int(tts_val)

                # **Monitor time loop**
                time_loop_start = time.time()
                stk0p = 0.0
                stk0s = 0.0

                # Use numpy to avoid repeated range checks
                ip_range = range(ttp_val, min(ttp_val + nsamples, nsamples))
                is_range = range(tts_val, min(tts_val + nsamples, nsamples))

                for ip, is_ in zip(ip_range, is_range):
                    stk0p += stalta_p[j, ip]  # assuming j is defined and valid
                    stk0s += stalta_s[j, is_]  # assuming j is defined and valid

                time_loop_end = time.time()
                if stk0p * stk0s > stkmax:
                    stkmax = stk0p * stk0s
                    kmax = ip  # Save the best matching index

                logging.debug(f"Time loop for sensor {j+1} and iteration {i} took {time_loop_end - time_loop_start:.6f} seconds")

            sensor_loop_end = time.time()
            logging.debug(f"Sensor loop for outer loop iteration {i} took {sensor_loop_end - sensor_loop_start:.6f} seconds")

            corrmatrix[i] = num.sqrt(stkmax) / nsta

            if corrmatrix[i] > corrmax:
                corrmax = corrmatrix[i]
                iloc = i
                itime = kmax

            loop_end = time.time()  # End timer for outer loop iteration
            logging.debug(f"Outer loop iteration {i} took {loop_end - loop_start:.6f} seconds")

        outer_end = time.time()
        logging.info(f"Outer loop duration: {outer_end - outer_start:.2f} sec")

        end_time = time.time()
        logging.info(f"Stacking function completed in {end_time - start_time:.2f} seconds.")

        return (iloc, itime), corrmatrix



    def locate_event(self):

        '''
        Main code
        '''

        self.location_domain()
        print("Location domain set up.")

        #corrmatrix stations-fiber 

        iloc_sta, corrmatrix_sta = self.stacking(self.lon_stations_rel, self.lat_stations_rel, self.depth_stations_rel, self.ttp, self.tts, self.obsp_sta, self.obss_sta)
        iloc_ch, corrmatrix_ch = self.stacking(self.lon_channels_rel, self.lat_channels_rel, self.depth_channels_rel, self.ttp, self.tts, self.obsp_ch, self.obss_ch)
        

        iloc = (iloc_sta[0] + iloc_ch[0]) / 2
        itime = (iloc_sta[1] + iloc_ch[1]) / 2

        #hybrid corrmatrix 
        corrmatrix = (corrmatrix_sta + corrmatrix_ch) / 2

        print("Event located successfully!")

        return iloc_sta, iloc_ch, iloc, itime, corrmatrix_sta.reshape(self.nx, self.nx, self.nz), corrmatrix_ch.reshape(self.nx, self.nx, self.nz), corrmatrix.reshape((self.nx, self.nx, self.nz))
