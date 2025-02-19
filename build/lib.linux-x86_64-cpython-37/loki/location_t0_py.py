import numpy as num
import logging 
import time 

#logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s")


class WaveformStacking:
    def __init__(self, tobj, sobj, nproc, ttp, tts, obsp_sta, obss_sta, obsp_ch, obss_ch):
        required_attrs = ["nx", "nz", "dx", "dz", "lon_stations", "lat_stations", 
                          "depth_stations", "lon_channels", "lat_channels", "depth_channels"]
        for attr in required_attrs:
            if not hasattr(tobj, attr):
                raise AttributeError(f"Missing required attribute: {attr} in tobj")
            
        self.stations = sobj.stations

        self.nproc = nproc
        #traveltime table dimension
        self.x0 = tobj.x0   #orgin of the traveltime table
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
        self.obsp_ch = obsp_ch
        self.obss_ch = obss_ch

        #adapt location of sensors to the dx-dz of the domain (to be fixed sistematically)
        
        self.lon_stations = num.array(tobj.lon_stations, dtype=float) 
        self.lat_stations = num.array(tobj.lat_stations, dtype=float) 
        self.depth_stations = num.array(tobj.depth_stations, dtype=float)
        self.lon_channels = num.array(tobj.lon_channels, dtype=float)*0.00001
        self.lat_channels = num.array(tobj.lat_channels, dtype=float)*0.00001
        self.depth_channels = num.array(tobj.depth_channels, dtype=float)*0.00001


    def location_domain(self):

        '''
        This function define the 3D grid location domain, starting from the 2D traveltime table 
        spanning the diagonal of the domain
        '''
        #number of points in the 3D location grid 
        nxyz = self.nx * self.nx * self.nz

        print(nxyz)

        #define the current extension of the sensors on x,y,z (necessary for relative location)
        #extx_sub = num.abs(num.minimum(self.lon_stations.min(), self.lon_channels.min()) - 
        #                   num.maximum(self.lon_stations.max(), self.lon_channels.max()))
        #exty_sub = num.abs(num.minimum(self.lat_stations.min(), self.lat_channels.min()) - 
        #                   num.maximum(self.lat_stations.max(), self.lat_channels.max()))
        
        #f
        extx_sub = num.abs(self.lon_stations.min() - self.lon_stations.max())
        exty_sub = num.abs(self.lat_stations.min() - self.lat_stations.max())
        #extz_sub = num.abs(num.minimum(self.depth_stations.min(), self.depth_channels.min()) - 
                          # num.maximum(self.depth_stations.max(), self.depth_channels.max()))

        extx_sub = extx_sub 
        exty_sub = extx_sub

        #define the current extension of the traveltime domain 
        extx_tt = (self.nx * self.dx)
        extz_tt = (self.nz * self.dz)
        
        #define the definitive dimension of the location domain 
        extx = num.max(self.x)  #extx_tt / num.sqrt(2)
        exty = num.max(self.x)  #extx_tt / num.sqrt(2)
        extz = num.max(self.z)  #extz_tt

        #loc_model = nxyz.reshape(self.nx, self.nx, self.nz)

        print('i am here 1')

        # Generate grid coordinates
        x_vals = self.x
        y_vals = self.y
        z_vals = self.z

        print('this is the dimension on x,y,z of the location domain', x_vals, y_vals, z_vals )
        
        
        # Create meshgrid for coordinates
        #X, Y, Z = num.meshgrid(x_vals, y_vals, z_vals, indexing='ij')

        # Create a 3D mesh grid
        self.X, self.Y, self.Z = num.meshgrid(x_vals, y_vals, z_vals, indexing='ij')

        #for i in range(nx):
        #    for j in range(ny):
        #        for k in range(nz):
        #            x=i*dx
        #            y=j*dy
        #            z=k*dz




        #

        
        # Flatten to get 1D arrays of grid points
        self.lon_source = self.X.ravel()
        self.lat_source = self.Y.ravel()
        self.depth_source = self.Z.ravel()


        # Flatten the grids to create 1D arrays for sources
        #self.lon_source = ((X * self.nx/num.sqrt(2)).flatten())*self.dx/num.sqrt(2)*self.dx/num.sqrt(2)*self.dx/num.sqrt(2)
        #self.lat_source = ((Y * self.nx/num.sqrt(2)).flatten())*self.dx/num.sqrt(2)*self.dx/num.sqrt(2)*self.dx/num.sqrt(2)
        #self.depth_source = ((Z * self.nz/num.sqrt(2)).flatten())*self.dz/num.sqrt(2)*self.dz/num.sqrt(2)*self.dz/num.sqrt(2)

        #print(len(x_vals), len(y_vals), len(z_vals), len(self.lon_source))
        #print(self.lon_source)

        print('i am here 2')


        #for k in range(0,nxyz):
        #    print(k)
        #    for x in range(0,self.nx):
        #        for y in range(0, self.nx):
        #            for z in range(0,self.nz):
        #                self.lon_source[k] = x*self.nx
        #                self.lat_source[k] = y*self.nx
        #                self.depth_source[k] = z*self.nz


        #difference between the domain spanned by the sensors and the location domain 
        self.diff_subx = extx - extx_sub
        self.diff_suby = exty - exty_sub
        #diff_subz = extz - extz_sub


        #print('diff',extx_sub, exty_sub, extx_tt, extz_tt, self.diff_subx, self.diff_suby)
        
        #define all the relative location of the sensors (to the 0,0,0)
        #position the stations at the center of the investigated domain 

        self.lon_stations_rel = (self.lon_stations) # - self.x0 )
        self.lat_stations_rel = (self.lat_stations) # - self.y0 )
        self.depth_stations_rel = (self.depth_stations) # - self.z0
        self.lon_channels_rel = (self.lon_channels) # -  self.x0)
        self.lat_channels_rel = (self.lat_channels) # - self.y0) 
        self.depth_channels_rel = (self.depth_channels) # - self.z0


        print('lon stations rel', self.lon_stations_rel)
        print('lat stations rel', self.lat_stations_rel)



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

        c = num.abs(depth - depth_value)  
        closest_z_idx = num.argmin(c)  
        
        #closest_x_idx = num.argmin(num.abs(horiz_dist - horizontal_distance))
        #closest_z_idx = num.argmin(num.abs(depth - depth_value))
        
        #extract the traveltime associated to that horizontal distance-depth 
        travel_time = tt_2d[closest_x_idx, closest_z_idx]

        end_time = time.time()  # End timing
        elapsed_time = end_time - start_time

        # Log execution time
        logging.debug(f"Travel time lookup took {elapsed_time:.6f} seconds")

        return horiz_dist, depth, closest_x_idx, closest_z_idx, travel_time
    

    def stacking(self, lon_sensors, lat_sensors, depth_sensors, itp, its, stalta_p, stalta_s):

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

        # Precompute the sensor distances to avoid repeated calculation

        #sensor_distances = num.sqrt(num.array(lon_sensors)**2 + num.array(lat_sensors)**2) #distance from 0,0,0
        
        #this must be from all the possible 3d points
        
        #print('len(sensor_distances), ', len(sensor_distances))

        #print('sensor distances', sensor_distances[0:100])


        #OUTHER LOOP ON THE GRID 
        #nxyz = 1
        for i in range(0, nxyz):  # Limit loop for debugging


            stk0p = num.zeros((nsta, nsamples))
            stk0s = num.zeros((nsta, nsamples))
            

            loop_start = time.time()  # Start timer for outer loop iteration
            if i % progress_step == 0:
                print(f"Processing: {100 * i / nxyz:.6f}% completed")
                #
                # print(f"Processing: {10 * i / nxyz:.6f}% completed")

            stkmax = -1.0
            kmax = 0

            # **Monitor sensor loop**
            sensor_loop_start = time.time()
            
            
                        #LOOP ON THE STATIONS 
            for j in range(nsta):

                current_sta = self.stations[j]

                #print(current_sta)

                if current_sta[2] == '0':
                    current_sta = num.int(current_sta[3]) -1
                else:
                    current_sta = num.int(current_sta[2:4]) -1

                #print(current_sta)

                #current_sta = num.int(current_sta[3:4]) - 1

                #print('lon_sensors[current_sta]', lon_sensors[current_sta])
                #print('lon_source[i]', self.lon_source[i])

                #print(num.int(current_sta[3:4]))


                a = 0
                #print('j', j)
                # current horizontal distance and depth 
                
                #compute euclidean distance 

                
            
                #print('current sensor', lon_sensors[j], lat_sensors[j], depth_sensors[j])

                #print('current source', self.lon_source[i], self.lat_source[i], self.depth_source[i])            

                current_horizontal_distance = num.sqrt(num.abs(lon_sensors[current_sta]-self.lon_source[i])**2 + num.abs(lat_sensors[current_sta]-self.lat_source[i])**2)

                #current_horizontal_distance = sensor_distances[j]
                #print('current_horizontal_distance',current_horizontal_distance)
                current_depth = num.abs(depth_sensors[current_sta] - self.depth_source[i])
                #print('current_depth',current_depth)

                #extract the traveltime 

                ttp_start = time.time()  # Start timing the travel time lookup (TP)
                horiz, dep, _, _, ttp_val = self.get_closest_travel_time(current_horizontal_distance, current_depth, itp)
                ttp_end = time.time()  # End timing

                
                tts_start = time.time()  # Start timing the travel time lookup (TS)
                horiz2, dep2, _, _, tts_val = self.get_closest_travel_time(current_horizontal_distance, current_depth, its)
                tts_end = time.time()  # End timing

                #CHECK THIS 
                ttp_val = int(ttp_val)
                 
                #print('horiz, dep', horiz[0:100], dep[0:100])
                #print('ttp_val', ttp_val)
                tts_val = int(tts_val)

                # **Monitor time loop**
                time_loop_start = time.time()


                # Use numpy to avoid repeated range checks
                #ip_range = range(ttp_val, min(ttp_val + nsamples, nsamples))
                #is_range = range(tts_val, min(tts_val + nsamples, nsamples))
                
                ip_range = range(ttp_val, ttp_val + nsamples)
                is_range = range(tts_val, tts_val + nsamples)
                

                #stack for each station and time  
                

                for ip, is_ in zip(ip_range, is_range):
                    if is_ < nsamples:

                        #print(j)

                        #print('ip', ip, 'is_', is_, 'j',j)   
                        stk0p[j,a] = stalta_p[j, ip]  
                        stk0s[j,a] = stalta_s[j, is_] 
                        #print('a', a, 'stalta_p[0, ip]', stalta_p[j, ip], 'stalta_s[j, is_]', stalta_s[j, is_])
                        a = a +1
                    else: 
                    
                    #stk0p=0 + stk0p
                    #stk0s=0 + stk0s
                        stk0p[j,a] = 0 
                        stk0s[j,a] = 0
                        a = a +1
                #print('stack p', stk0p)
                

                #time_loop_end = time.time()

            #print(stalta_p.shape)
            #print(stalta_p[2,0:10])

            #print('stack0 p', stk0p)
            #print('stack0 s', stk0s)


            stk0p_sta = num.sum(stk0p, axis=0)
            stk0s_sta = num.sum(stk0s, axis=0)

            #print('stack p', stk0p_sta)
            #print('stack s', stk0s_sta)

            for k in range(0, nsamples):

                if stk0p_sta[k]*stk0s_sta[k] > stkmax:
                    stkmax = stk0p_sta[k]*stk0s_sta[k]
                    kmax = k  # Save the best matching index
                    #print('stack', stkmax)

            corrmatrix[i] = num.sqrt(stkmax) / nsta
            #print('current corr matrix', corrmatrix[i])

            if corrmatrix[i] > corrmax:
                corrmax = corrmatrix[i]
                #print('current corr matrix', corrmatrix[i])
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
        #iloc_sta, corrmatrix_sta = self.stacking(self.lon_stations_rel, self.lat_stations_rel, self.depth_stations_rel, self.ttp, self.tts, self.obsp_sta, self.obss_sta)

        #iloc_ch, corrmatrix_ch = self.stacking(self.lon_channels_rel, self.lat_channels_rel, self.depth_channels_rel, self.ttp, self.tts, self.obsp_ch, self.obss_ch)
        

        #iloc = (iloc_sta[0] + iloc_ch[0]) / 2
        #itime = (iloc_sta[1] + iloc_ch[1]) / 2
        itime = (iloc_sta[1] + iloc_sta[1]) / 2

        #hybrid corrmatrix 
        #corrmatrix = (corrmatrix_sta + corrmatrix_ch) / 2
        corrmatrix = corrmatrix_sta 

        print("Event located successfully!")

        #return iloc_sta, iloc_ch, iloc, itime, corrmatrix_sta.reshape(self.nx, self.nx, self.nz), corrmatrix_ch.reshape(self.nx, self.nx, self.nz), corrmatrix.reshape((self.nx, self.nx, self.nz))
        return iloc_sta, iloc_sta, iloc_sta, itime, corrmatrix_sta.reshape(self.nx, self.nx, self.nz), corrmatrix.reshape(self.nx, self.nx, self.nz), corrmatrix.reshape(self.nx, self.nx, self.nz)


