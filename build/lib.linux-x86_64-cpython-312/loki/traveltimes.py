# %%

#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.
#       Author: Francesco Grigoli
import os
import sys
import numpy as num
import matplotlib.pyplot as plt
from pyproj import Proj, transform
import utm
from loki import latlon2cart
#import loki.latlon2cart as ll2c

class Traveltimes:

    def __init__(self, db_path, hdr_filename, geometry_filename_fiber, geometry_filename_stat):
        if not os.path.isdir(db_path):
            print('Error: data or database path do not exist')
            sys.exit()
        self.db_path = db_path
        if not os.path.isfile(db_path+'/'+hdr_filename):
            print('Error: header file does not exist')
            print('this is the path:', db_path+'/'+hdr_filename)
            sys.exit()
        self.refsta = None
        self.hdr_filename = hdr_filename
        self.load_header()
        self.geometry_filename_stat = geometry_filename_stat
        self.geometry_filename_fiber = geometry_filename_fiber
        self.load_station_info()
        self.load_channel_info()


#modify to read only two components (ny out, modify header)





    def convert_to_utm(lat, lon, ref_lat=51.64, ref_lon=7.72):
        ref_east, ref_north, _, _ = utm.from_latlon(ref_lat, ref_lon)
        east, north, _, _ = utm.from_latlon(lat, lon)
        return (east - ref_east) / 1000, (north - ref_north) / 1000

    def load_header(self):

        #this method loads information on the big 2D traveltime grid and 3D location grid 

                # Define WGS84 (lat/lon) and UTM projection (example for zone 33N)
        self.wgs84 = Proj(proj="latlong", datum="WGS84")
        self.utm33n = Proj(proj="utm", zone=32, datum="WGS84")
        #print(self.wgs84)

   
        f = open(os.path.join(self.db_path, self.hdr_filename))
        lines = f.readlines()  #read header info 
        #here info on the big traveltime 2D grid 
        self.nx, self.nz = [ int(x)   for x in lines[0].split()]  #number of points on the grid
        self.y0, self.x0, self.z0 = [ float(x) for x in lines[1].split()]  #starting point of the grid 
        self.dx, self.dz = [ float(x) for x in lines[2].split()]  #grid spacing
        self.lat0, self.lon0 = [ float(x) for x in lines[3].split()] #lat lon starting point 
        self.nttx, self.nttz=[ int(x)   for x in lines[4].split()] #traveltime points along distance and depth 
        self.dttx, self.dttz=[ float(x)   for x in lines[5].split()] #spacing in travetimes 
        self.dtt=float(lines[6])
        toks = lines[7].split()
        if len(toks) > 1:
            self.ref_station_coordinates = [eval(toks[1]) - self.lon0, eval(toks[2]) - self.lat0, eval(toks[3])]
        else:
            self.ref_station_coordinates = None  
        self.refsta = toks[0] if toks else None  

        #here info on the location 3D grid. The domain is a cube, so the grid is the same in x and y
        #the dimension of the traveltime table is the one driving the dimension of the grid
        #the traveltime table is the diagonal of the grid  (num.sqrt(2)), thus x,y,z are recontructed using the square root of 2

        #self.x =  num.arange(0, (self.nx * self.dx)/num.sqrt(2) , (self.dx/num.sqrt(2)))  #define the grid search based on the 2D traveltime grid
        #self.y =  num.arange(0, (self.nx * self.dx)/num.sqrt(2), (self.dx/num.sqrt(2) ))  #define the grid search based on the 2D traveltime grid
        #self.z =  num.arange(0, (self.nz * self.dz), self.dz) #define the grid search based on the 2D traveltime grid
        
        self.x =  num.arange(0, (self.nx * self.dx) , (self.dx))  #define the grid search based on the 2D traveltime grid
        self.y =  num.arange(0, (self.nx * self.dx), (self.dx))  #define the grid search based on the 2D traveltime grid
        self.z =  num.arange(0, (self.nz * self.dz), self.dz) #define the grid search based on the 2D traveltime grid
        

        self.nxyz=self.nx*self.nx*self.nz #number of points 
        self.nxz=self.nx*self.nz 
        self.delta_das = 0.01  #

        #origin=latlon2cart.Coordinates(self.x0, self.y0,self.z0)


        
        


        #self.x0, self.y0 = transform(self.wgs84, self.utm33n, self.x0, self.y0)
        #self.x0 = self.x0*1e-3
        #self.y0 = self.y0*1e-3


    def load_station_info(self): 
        
        #read information on the location grid and the stations 

        origin=latlon2cart.Coordinates(self.lat0, self.lon0,self.z0)

        #print(origin.X0)

        #self.x0 = origin.x0
        #self.y0 = origin.y0
        #self.z0 = 0

        self.stations_coordinates={}
        self.db_stations = []
        self.lon_stations = []
        self.lat_stations = []
        self.depth_stations = []
        f = open(os.path.join(self.db_path, self.geometry_filename_stat))
        lines = f.readlines()  #read header info
        for line in lines:
            # Assuming the file is space or tab-delimited. Adjust delimiter as needed.
            columns = line.strip().split()  # Removes any leading/trailing whitespace and splits by spaces

            # Check if the line has at least 3 columns to avoid errors
            if len(columns) >= 4:
                self.db_stations.append(str(columns[0]))
                lon_degr = float(columns[2])
                lat_degr = float(columns[1])
                depth = float(columns[3])
                #print('a')
                #print(lon_degr,lat_degr,depth)
                #late,lone,elev=origin.cart2geo(lon_degr,lat_degr,depth)
                lone,late,elev = origin.geo2cart(lat_degr, lon_degr,ele=0, relative=True, geo2enu=False)

                
                #print(lone,late,elev)
                self.lon_stations.append(lone)
                self.lat_stations.append(late)
                self.depth_stations.append(elev)

                #lon_utm, lat_utm = transform(self.wgs84, self.utm33n, lon_degr, lat_degr)
                #self.lon_stations.append(lon_utm*1e-3 - self.lon0)
                #self.lat_stations.append(lat_utm*1e-3 - self.lat0)
                #self.depth_stations.append(float(columns[3]))

                self.stations_coordinates[str(columns[0])] = (lone, late, elev)
                

    def load_channel_info(self): 


        
        #read information on the location grid and the stations 
        self.channels_coordinates={}
        self.db_channels = []
        self.lon_channels = []
        self.lat_channels = []
        self.depth_channels = []

        a = open(os.path.join(self.db_path, self.geometry_filename_fiber))
        
        lines = a.readlines()  #read header info
        for line in lines:
            # Assuming the file is space or tab-delimited. Adjust delimiter as needed.
            columns = line.strip().split()  # Removes any leading/trailing whitespace and splits by spaces

            # Check if the line has at least 3 columns to avoid errors
            if len(columns) >= 4:

                self.db_channels.append(str(columns[0]))
                self.lon_channels.append(float(columns[1]))
                self.lat_channels.append(float(columns[2]))
                self.depth_channels.append(float(columns[3]))
                self.channels_coordinates[str(columns[0])] = (self.lon_stations, self.lat_stations, self.depth_stations)
        

    def load_traveltimes(self, phase, label='layer', precision='single'):

        '''
        This function now reads the big 2d traveltime table, not a traveltime table for each station
        
        '''

        t={}
        sta = self.refsta
        fn = os.path.join(self.db_path, '%(label)s.%(phase)s.%(station)s.time.buf' %{"label":label,"phase":phase, "station":sta} )
        if (precision=='single'):
            t[sta]= num.fromfile(fn, dtype=num.float32)
        elif (precision=='double'):
            t[sta]= num.fromfile(fn, dtype=num.float64)
        elif (precision=='single_int'):
            t[sta]= num.fromfile(fn, dtype=num.int16)
        elif (precision=='double_int'):
            t[sta]= num.fromfile(fn, dtype=num.int32)
        else:
            print('Error: precision must be set to "single_int" or "double_int"!!!')
            sys.exit()
        return t

######### updated up to here #########


#modify 

    def ttdb_reduce(self,tt,l_lim,u_lim,zlim=[]):
        latref=self.lat0; lonref=self.lon0; eleref=0.
        #origin=ll2c.Coordinates(latref,lonref,eleref)
        x_l,y_l,u_l = origin.geo2cart(l_lim[0],l_lim[1],eleref)
        x_u,y_u,u_u = origin.geo2cart(u_lim[0],u_lim[1],eleref)
        x_l=x_l/1000.; y_l=y_l/1000.
        x_u=x_u/1000.; y_u=y_u/1000.
        nx_ini=int((x_l-self.x0)/self.dx); ny_ini=int((y_l-self.y0)/self.dy)
        nx_fin=int((x_u-self.x0)/self.dx); ny_fin=int((y_u-self.y0)/self.dy)
        self.x0=x_l-self.x0; self.y0=y_l-self.y0
        if zlim:
           nz_ini=int((zlim[0]-self.z0)/self.dz)
           nz_fin=int((zlim[1]-self.z0)/self.dz)
           self.z0=zlim[0]
        else:
           nz_ini=0; nz_fin=self.nz
        nx_new=nx_fin-nx_ini; ny_new=ny_fin-ny_ini; nz_new=nz_fin-nz_ini
        tt_new={}
        for key in tt.keys():
            t_arr=tt[key].reshape(self.nx,self.ny,self.nz)
            tt_new[key]=t_arr[nx_ini:nx_fin,ny_ini:ny_fin,nz_ini:nz_fin].reshape(nx_new*ny_new*nz_new)
        self.nx=nx_new; self.ny=ny_new; self.nz=nz_new
        self.x0=0.; self.y0=0.
        self.x = self.x0+(num.arange(0,self.nx)*self.dx)
        self.y = self.y0+(num.arange(0,self.ny)*self.dy)
        self.z = self.z0+(num.arange(0,self.nz)*self.dz)
        self.nxyz=self.nx*self.ny*self.nz
        self.lat0=l_lim[0]
        self.lon0=l_lim[1]

        return tt_new

    def interpolation(self, tt, dx_new, dy_new, dz_new):
        from scipy.interpolate import RegularGridInterpolator
        xr=self.dx/dx_new; yr=self.dy/dy_new; zr=self.dz/dz_new;
        x_old=num.arange(self.nx)*self.dx; y_old=num.arange(self.ny)*self.dy; z_old=num.arange(self.nz)*self.dz
        x_new=num.arange(1,(self.nx-1)*xr)*dx_new; y_new=num.arange(1,(self.ny-1)*yr)*dy_new; z_new=num.arange(1,(self.nz-1)*zr)*dz_new
        nx_new=num.size(x_new); ny_new=num.size(y_new); nz_new=num.size(z_new)
        grid=[]
        for i in range(nx_new):
              for j in range(ny_new):
                  for k in range(nz_new):
                        grid.append([x_new[i],y_new[j],z_new[k]])
        grid=num.array(grid)
        t_interp={}
        for key in tt.keys():
            print('interpolation station: '+key+'\n')
            t_arr=tt[key].reshape(self.nx,self.ny,self.nz)
            interpolation = RegularGridInterpolator((x_old, y_old, z_old), t_arr)
            t_int=interpolation(grid)
            t_interp[key]=t_int.reshape(nx_new*ny_new*nz_new)
        self.nx=nx_new; self.ny=ny_new; self.nz=nz_new
        self.dx=dx_new; self.dy=dy_new; self.dz=dz_new
        self.x = self.x0+x_new; self.y = self.y0+y_new; self.z = self.z0+z_new
        self.nxyz=self.nx*self.ny*self.nz

        return t_interp
    
#modify 

    def save_ttdb(self,tt,phase,label):

        f=open(self.db_path+'/'+label+'.header.hdr','w')
        f.write('%d %d %d' %(self.nx,self.ny,self.nz)+'\n'+
                '%f %f %f' %(self.x0, self.y0, self.z0)+'\n'+
                '%f %f %f' %(self.dx, self.dy, self.dz)+'\n'+
                '%f %f' %(self.lat0, self.lon0)+'\n')

        for key in tt.keys():
            fn = os.path.join(self.db_path, '%(label)s.%(phase)s.%(station)s.time.buf' %{"label":label, "phase":phase, "station":key})
            (tt[key].astype(dtype=num.float32)).tofile(fn)
            f.write(key+'\n')

        f.close()

        return None

    def event_indexes(self,evlat,evlon,evdepth):
        km=1000.
        latref=self.lat0; lonref=self.lon0; eleref=0.
        origin=ll2c.Coordinates(latref,lonref,eleref)
        xi,yi,ui = origin.geo2cart(evlat,evlon,eleref)
        zi=evdepth*km
        x0=self.x0*km; y0=self.y0*km; z0=0.0+self.z0
        d_spac=self.dx*km
        ix=int(num.round((xi-x0)/d_spac)); iy=int(num.round((yi-y0)/d_spac)); iz=int(num.round((zi-z0)/d_spac))
        return ix,iy,iz

    def ttdb_generator(self, velocity, phase='P'):
        # Generate traveltimes for an homogenaous medium
        if self.stations_coordinates is not None:
            for sta in self.stations_coordinates:
                print('Starting to calculate traveltimes for the station : ' + sta +' \n')
                xsta,ysta,zsta=self.stations_coordinates[sta][0:3]
                fname='homo.%(phase)s.%(station)s.time.buf' %{"phase":phase, "station":sta}
                fout=open(self.db_path+'/'+fname,'wb')
                tt=num.zeros(self.nxyz)
                print(self.nxyz,self.nx,self.ny,self.nz)
                for k in range(self.nxyz):
                    ix=k//(self.ny*self.nz)
                    iy=k//self.nz-(ix*self.ny)
                    iz=k-(iy*self.nz)-(ix*self.ny*self.nz)
                    dist=num.sqrt((self.x[ix]-xsta)**2+(self.y[iy]-ysta)**2+(self.z[iz]-zsta)**2)
                    tt[k]=dist/velocity
                tt.tofile(fout)
                fout.close()
                print('Traveltimes computation for the station : ' + sta + ' completed! \n')
        return None

    def apply_master_event_correction(self, phase, dt, label='layer', precision='single'):
        t={}
        for sta in dt.keys():
            try:
              fn = os.path.join(self.db_path, '%(label)s.%(phase)s.%(station)s.time.buf' %{"label":label,"phase":phase, "station":sta} )
              if (precision=='single'):
                  tt= num.fromfile(fn, dtype=num.float32)+dt[sta]
                  tt[tt<0.]=0.
                  t[sta]=tt
              elif (precision=='double'):
                  t[sta]= num.fromfile(fn, dtype=num.float64)+dt[sta]
              else:
                  print('Error: precision must be set to "single" or "double"!!!')
                  sys.exit()
            except:
                print('Error: reading file for station' + sta)
        return t



if __name__=='__main__':
    db_path='/home/emanuele/data/emanuele/loki-das/Traveltimes'
    tt0=Traveltimes(db_path, 'header_long.hdr', 'station_das_ign.tmp')
    tp0=tt0.load_traveltimes('P', 'das')
    #print(tt0.nxyz)
    tpxy=tp0['HM01'].reshape(tt0.nx, tt0.nz)
    #print(tpxy)
    plt.figure()
    plt.imshow(tpxy)
    plt.show(block=True)




# #    tp4=tp3[k]
# #    print(num.shape(tp4),tt3.nx*tt3.ny*tt3.nz)
# #    plt.figure()
# #    plt.imshow(tp2[:,:,0])
# #    plt.figure()
# #    plt.imshow(tp3[:,:,0])
# #    plt.show()
# tt1=traveltimes(db_path, 'header.hdr')
# tp0=tt0.load_traveltimes('P','layer')
# tp1=tt1.load_traveltimes('P','layer')
# tp_red0=tt0.ttdb_reduce(tp0,[35.5,129.00],[36.5,130.0])
# tp_red1=tt1.ttdb_reduce(tp1,[35.5,129.00],[36.5,130.0])
# tp_int1=tt1.interpolation(tp_red1,0.5,0.5,0.5)
# tt3=traveltimes(db_path, 'interpolated.header.hdr')
# tp3=tt3.load_traveltimes('P','interpolated')
# #tt1.save_ttdb(tp_int1,'P','interpolated')
#
# #400 350 65
# #-250.000000 -150.000000 -5.000000
# #31.88 31.94  130.84 130.92
# #for k in tp0.keys():
# #    #tp1=tp[k].reshape(400, 350, 65)
# #    tp2=tp_red0[k].reshape(tt0.nx, tt0.ny, tt0.nz)
# #    tp4=tp3[k]
# #    print(num.shape(tp4),tt3.nx*tt3.ny*tt3.nz)
# #    plt.figure()
# #    plt.imshow(tp2[:,:,0])
# #    plt.figure()
# #    plt.imshow(tp3[:,:,0])
# #    plt.show()

# %%
