import numpy as num
import logging 
import time 
from pathlib import Path

#logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s")


class WaveformStacking:

    #inizia la classe 

    def __init__(self, tobj, sobj, nproc, ttp, tts, obsp_sta, obss_sta, obsp_ch, obss_ch):
        required_attrs = ["nx", "nz", "dx", "dz", "lon_stations", "lat_stations", 
                          "depth_stations", "lon_channels", "lat_channels", "depth_channels"]
        for attr in required_attrs:
            if not hasattr(tobj, attr):
                raise AttributeError(f"Missing required attribute: {attr} in tobj")


        #attributi 

        self.stations = sobj.stations  #id_stations 
        #self.channels = sobj.channels  #id_stations 
        self.deltat = sobj.deltat      #sampling rate 

        self.nproc = nproc
        #traveltime table dimension
        self.x0 = tobj.x0   #orgin of the traveltime table (0,0,0)
        self.y0 = tobj.y0
        self.z0 = tobj.z0

        self.nx = tobj.nx #dimension of the traveltime table
        self.nz = tobj.nz 

        self.dx = tobj.dx #spacing of the traveltime table
        self.dz = tobj.dz
        #domain (from the assumption that the domain is defined from the dimension of the traveltime table, representing its diagonal)
        self.x = tobj.x
        self.y = tobj.y
        self.z = tobj.z
        #traveltime table
        self.ttp = ttp
        self.tts = tts
        #observed waveform stations
        self.obsp_sta = obsp_sta
        self.obss_sta = obss_sta

        #observed waveform fiber 
        #self.obsp_ch = obsp_ch
        #self.obss_ch = obss_ch

        #adapt location of sensors to the dx-dz of the domain (to be fixed sistematically)
        
        self.lon_stations = num.array(tobj.lon_stations, dtype=float) #change name, otherwise misleading (x,y,z)
        self.lat_stations = num.array(tobj.lat_stations, dtype=float) 
        self.depth_stations = num.array(tobj.depth_stations, dtype=float)
        #self.lon_channels = num.array(tobj.lon_channels, dtype=float)
        #self.lat_channels = num.array(tobj.lat_channels, dtype=float)
        #self.depth_channels = num.array(tobj.depth_channels, dtype=float)


    def location_domain(self):  #prende tutti gli attributi iniziali della classe 

        '''
        This function define the 3D grid location domain, starting from the 2D traveltime table 
        spanning the diagonal of the domain
        '''
        #number of points in the 3D location grid 
        nxyz = self.nx * self.nx * self.nz

        #define the current extension of the sensors on x,y,z (necessary for relative location)

        #extx_sub = num.abs(self.lon_stations.min() - self.lon_stations.max())
        #exty_sub = num.abs(self.lat_stations.min() - self.lat_stations.max())

        #extx_sub = extx_sub 
        #exty_sub = extx_sub

        #define the current extension of the traveltime domain 
        #extx_tt = (self.nx * self.dx)
        #extz_tt = (self.nz * self.dz)
        
        #define the definitive dimension of the location domain 
        #extx = num.max(self.x)  #extx_tt / num.sqrt(2)
        #exty = num.max(self.x)  #extx_tt / num.sqrt(2)
        #extz = num.max(self.z)  #extz_tt

        # Generate grid coordinates
        x_vals = self.x
        y_vals = self.y
        z_vals = self.z

        print('this is the dimension on x,y,z of the location domain', x_vals, y_vals, z_vals )
        

        # Create a 3D mesh grid
        self.X, self.Y, self.Z = num.meshgrid(x_vals, y_vals, z_vals, indexing='ij')

        #Alternative formulation avoiding meshgrid 
        #for i in range(nx):
        #    for j in range(ny):
        #        for k in range(nz):
        #            x=i*dx
        #            y=j*dy
        #            z=k*dz


        # Flatten to get 1D arrays of grid points
        self.lon_source = self.X.ravel()
        self.lat_source = self.Y.ravel()
        self.depth_source = self.Z.ravel()

        #print(len(x_vals), len(y_vals), len(z_vals), len(self.lon_source))
        #print(self.lon_source)

        #for k in range(0,nxyz):
        #    print(k)
        #    for x in range(0,self.nx):
        #        for y in range(0, self.nx):
        #            for z in range(0,self.nz):
        #                self.lon_source[k] = x*self.nx
        #                self.lat_source[k] = y*self.nx
        #                self.depth_source[k] = z*self.nz


        #difference between the domain spanned by the sensors and the location domain 
        #self.diff_subx = extx - extx_sub
        #self.diff_suby = exty - exty_sub


        #define all the relative location of the sensors (to the 0,0,0)
        #position the stations at the center of the investigated domain 

        self.lon_stations_rel = self.lon_stations 
        self.lat_stations_rel = self.lat_stations 
        self.depth_stations_rel = self.depth_stations
        #self.lon_channels_rel = self.lon_channels 
        #self.lat_channels_rel = self.lat_channels
        #self.depth_channels_rel = self.depth_channels




    def get_closest_travel_time(self, horizontal_distance, depth_value, tt_table):
        

        start_time = time.time()  # Start timing to monitor the procedure 

        tt_2d = tt_table.reshape(self.nx, self.nz)  #reshape the traveltimes using nx and nz dimenisons 

        # Log the input values
        logging.debug(f"Computing travel time for Distance={horizontal_distance:.4f}, Depth={depth_value:.4f}")

        #horizontal and depth distances  of the 2d traveltime table
        horiz_dist = num.arange(0, (self.nx * self.dx) -self.dx, (self.dx))
        depth = num.arange(0, (self.nz * self.dz) -self.dz, self.dz)

        #identify where the current depth, horizontal distance are most symilar to the one the traveltime

        a = num.abs(horiz_dist - horizontal_distance) 

        closest_x_idx = num.argmin(a)  

        diffx = horiz_dist[closest_x_idx] - horizontal_distance


        c = num.abs(depth - depth_value)  

        closest_z_idx = num.argmin(c)  

        ratiox = diffx/self.dx  #ratio between the distance of the closest point and the distance of the point of interest
        dtt = tt_2d[closest_x_idx, closest_z_idx] - tt_2d[closest_x_idx-1, closest_z_idx]  #dtt between the closest point and the point of interest
        diff_tt = dtt*ratiox #difference in travel time between the closest point and the point of interest
        
        #interpolation 

        travel_time = tt_2d[closest_x_idx, closest_z_idx] - diff_tt #if diff_tt is positive, the travel time is less than the closest point, otherwise is greater
        
        #no interpolation
        
        #travel_time = tt_2d[closest_x_idx, closest_z_idx]  #if diff_tt is positive, the travel time is less than the closest point, otherwise is greater

        end_time = time.time()  # End timing
        elapsed_time = end_time - start_time

        # Log execution time
        logging.debug(f"Travel time lookup took {elapsed_time:.6f} seconds")

        return horiz_dist, depth, closest_x_idx, closest_z_idx, travel_time
    

    def stacking(self, lon_sensors, lat_sensors, depth_sensors, itp, its, stalta_p, stalta_s, sensor = 'stations'):

        """Function stacking energy along predicted travel times."""
        
        logging.info("Stacking function called.")
        start_time = time.time()  # Start timer

        #define the location domain 
        nxyz = self.nx * self.nx * self.nz
        nsta = len(stalta_p[:, 1])
        nsamples = len(stalta_p[1, :])
        
        #correlation matrix 3D for the location
        corrmatrix = num.zeros(nxyz)

        corrmax = -1.0
        iloc, itime = 0, 0

        progress_step = nxyz//100 #10000 #nxyz // 300000  # 2% progress intervals

        # **Monitor outer loop**
        outer_start = time.time()

        #OUTHER LOOP ON THE GRID 
        #nxyz = 10
        for i in range(0, nxyz):  # Limit loop for debugging

            stk0p = num.zeros((nsta, nsamples))
            stk0s = num.zeros((nsta, nsamples))
            
            loop_start = time.time()  # Start timer for outer loop iteration
            if i % progress_step == 0:
                print(f"Processing: {100 * i / nxyz:.6f}% completed")

            stkmax = -1.0
            kmax = 0

            # **Monitor sensor loop**
            sensor_loop_start = time.time()
            
            
            #LOOP ON THE STATIONS 
            for j in range(nsta):

                #this is highly dependent on the format of the stations (it works now but not)

                if sensor == 'stations':
                    current_sta = self.stations[j]

                if sensor == 'fiber':
                    current_sta = self.channels[j]

                if current_sta[2] == '0':
                    current_sta = num.int(current_sta[3]) -1
                else:
                    current_sta = num.int(current_sta[2:4]) -1

                
         
                current_horizontal_distance = num.sqrt(num.abs(lon_sensors[current_sta]-self.lon_source[i])**2 + num.abs(lat_sensors[current_sta]-self.lat_source[i])**2)
                current_depth = num.abs(depth_sensors[current_sta] - self.depth_source[i])

                #extract the traveltime 

                ttp_start = time.time()  # Start timing the travel time lookup (TP)
                horiz, dep, _, _, ttp_val = self.get_closest_travel_time(current_horizontal_distance, current_depth, itp)
                ttp_end = time.time()  # End timing

                
                tts_start = time.time()  # Start timing the travel time lookup (TS)
                horiz2, dep2, _, _, tts_val = self.get_closest_travel_time(current_horizontal_distance, current_depth, its)
                tts_end = time.time()  # End timing


                ttp_val = int(ttp_val)
                tts_val = int(tts_val)

                # **Monitor time loop**
                time_loop_start = time.time()

            
                ip_range = range(ttp_val, ttp_val + nsamples)
                is_range = range(tts_val, tts_val + nsamples)
                
                #stack for each station and time  

                a = 0 #initialize a counter for samples  
                
                for ip, is_ in zip(ip_range, is_range):

                    if ip < nsamples and is_ < nsamples:

                        stk0p[j,a] = stalta_p[j, ip]  
                        stk0s[j,a] = stalta_s[j, is_] 

                        a = a +1
                    else: 

                        stk0p[j,a] = 0 
                        stk0s[j,a] = 0
                        a = a +1

            #stack the contribution from all the stations
            stk0p_sta = num.sum(stk0p, axis=0)
            stk0s_sta = num.sum(stk0s, axis=0)


            for k in range(0, nsamples):

                if stk0p_sta[k]*stk0s_sta[k] > stkmax:
                    stkmax = stk0p_sta[k]*stk0s_sta[k]
                    kmax = k 

            #normalize the correlation matrix
                    
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

        #corrmatrix stations

        iloc_sta, corrmatrix_sta = self.stacking(self.lon_stations_rel, self.lat_stations_rel, self.depth_stations_rel, self.ttp, self.tts, self.obsp_sta, self.obss_sta, sensor = 'stations')
        
        itime_sta = (iloc_sta[1] + iloc_sta[1]) / 2

        corrmatrix_sta = corrmatrix_sta 

        #corrmatrix fiber

        #iloc_ch, corrmatrix_ch = self.stacking(self.lon_channels_rel, self.lat_channels_rel, self.depth_channels_rel, self.ttp, self.tts, self.obsp_ch, self.obss_ch, sensor = 'fiber')
       
        #itime_ch = (iloc_ch[1] + iloc_ch[1]) / 2

        #corrmatrix_ch = corrmatrix_ch  

        
        #hybrid network results 

        #iloc_hybrid = (iloc_sta[0] + iloc_ch[0]) / 2
        #itime_hybrid = (iloc_sta[1] + iloc_ch[1]) / 2
        #corrmatrix_hybrid = (corrmatrix_sta + corrmatrix_ch) / 2


        print("Event located successfully!")

        #return iloc_sta, iloc_ch, itime_sta, iloc_ch, iloc_hybrid, itime_hybrid, corrmatrix_sta.reshape(self.nx, self.nx, self.nz), corrmatrix_ch.reshape(self.nx, self.nx, self.nz), corrmatrix_hybrid.reshape((self.nx, self.nx, self.nz))
        
        return iloc_sta, iloc_sta, iloc_sta, itime_sta, corrmatrix_sta.reshape(self.nx, self.nx, self.nz), corrmatrix_sta.reshape(self.nx, self.nx, self.nz), corrmatrix_sta.reshape(self.nx, self.nx, self.nz)


