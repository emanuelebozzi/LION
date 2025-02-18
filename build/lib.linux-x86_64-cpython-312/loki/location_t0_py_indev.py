import ctypes
import numpy as np
import logging
import time

# Assuming the C library is called 'stacking_lib.so' on Linux or 'stacking_lib.dll' on Windows
class WaveformStacking:
    def __init__(self, tobj, nproc, ttp, tts, obsp_sta, obss_sta, obsp_ch, obss_ch):
        # Initialization of required attributes
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

        # Adapt sensor locations as before
        self.lon_stations = np.array(tobj.lon_stations, dtype=float) * 0.0001
        self.lat_stations = np.array(tobj.lat_stations, dtype=float) * 0.0001
        self.depth_stations = np.array(tobj.depth_stations, dtype=float) * 0.0001
        self.lon_channels = np.array(tobj.lon_channels, dtype=float) * 0.0001
        self.lat_channels = np.array(tobj.lat_channels, dtype=float) * 0.0001
        self.depth_channels = np.array(tobj.depth_channels, dtype=float) * 0.0001




        # Load the compiled C library (adjust path to the shared library accordingly)
        self.lib = ctypes.CDLL('/home/emanuele/LOKI-DAS/loki/stacking_lib.so')  # Replace with your actual path

        # Define argument and return types for C function
        self.lib.stacking.argtypes = [
            ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
            ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,
            ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)
        ]
        self.lib.stacking.restype = None

    def stacking(self, lon_sensors, lat_sensors, depth_sensors, itp, its, stalta_p, stalta_s):
        logging.info("Stacking function called.")
        start_time = time.time()

        nxyz = self.nx * self.nx * self.nz
        nsta = stalta_p.shape[0]
        nsamples = stalta_p.shape[1]

        # Ensure all arrays are NumPy float64 arrays and C-contiguous
        lon_sensors = np.ascontiguousarray(lon_sensors, dtype=np.float64)
        lat_sensors = np.ascontiguousarray(lat_sensors, dtype=np.float64)
        depth_sensors = np.ascontiguousarray(depth_sensors, dtype=np.float64)
        itp = np.ascontiguousarray(itp, dtype=np.float64)
        its = np.ascontiguousarray(its, dtype=np.float64)
        stalta_p = np.ascontiguousarray(stalta_p, dtype=np.float64).ravel()
        stalta_s = np.ascontiguousarray(stalta_s, dtype=np.float64).ravel()

        # Prepare the result arrays
        corrmatrix = np.zeros(nxyz, dtype=np.float64)
        iloc = ctypes.c_int(0)
        itime = ctypes.c_int(0)

        # Call the C function
        self.lib.stacking(
            lon_sensors.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            lat_sensors.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            depth_sensors.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            itp.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            its.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            stalta_p.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            stalta_s.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            self.nx, self.nz, nsta, nsamples,
            corrmatrix.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ctypes.byref(iloc), ctypes.byref(itime)
        )

        end_time = time.time()
        logging.info(f"Stacking function completed in {end_time - start_time:.2f} seconds.")

        return (iloc.value, itime.value), corrmatrix


    def location_domain(self):
        '''Define the 3D grid location domain'''
        extx_sub = np.abs(np.minimum(self.lon_stations.min(), self.lon_channels.min()) - 
                          np.maximum(self.lon_stations.max(), self.lon_channels.max()))
        exty_sub = np.abs(np.minimum(self.lat_stations.min(), self.lat_channels.min()) - 
                          np.maximum(self.lat_stations.max(), self.lat_channels.max()))
        extz_sub = np.abs(np.minimum(self.depth_stations.min(), self.depth_channels.min()) - 
                          np.maximum(self.depth_stations.max(), self.depth_channels.max()))
        
        extx_tt = self.nx * self.dx
        extz_tt = self.nz * self.dz
        
        extx = extx_tt / np.sqrt(2)
        exty = extx_tt / np.sqrt(2)
        extz = extz_tt

        diff_subx = extx - extx_sub
        diff_suby = exty - exty_sub
        diff_subz = extz - extz_sub

        self.lon_stations_rel = (self.lon_stations + diff_subx / 2) - self.lon_stations.min()
        self.lat_stations_rel = self.lat_stations + diff_suby / 2 - self.lat_stations.min()
        self.depth_stations_rel = self.depth_stations + diff_subz / 2 - self.depth_stations.min()
        self.lon_channels_rel = self.lon_channels + diff_subx / 2 - self.lon_channels.min()
        self.lat_channels_rel = self.lat_channels + diff_suby / 2 - self.lat_channels.min()
        self.depth_channels_rel = self.depth_channels + diff_subz / 2 - self.depth_channels.min()


    def locate_event(self):
        '''Main code to locate the event'''
        self.location_domain()
        print("Location domain set up.")

        iloc_sta, corrmatrix_sta = self.stacking(self.lon_stations_rel, self.lat_stations_rel, self.depth_stations_rel, self.ttp, self.tts, self.obsp_sta, self.obss_sta)
        iloc_ch, corrmatrix_ch = self.stacking(self.lon_channels_rel, self.lat_channels_rel, self.depth_channels_rel, self.ttp, self.tts, self.obsp_ch, self.obss_ch)

        iloc = (iloc_sta[0] + iloc_ch[0]) / 2
        itime = (iloc_sta[1] + iloc_ch[1]) / 2

        # Hybrid corrmatrix
        corrmatrix = (corrmatrix_sta + corrmatrix_ch) / 2

        print("Event located successfully!")

        return iloc_sta, iloc_ch, iloc, itime, corrmatrix_sta.reshape(self.nx, self.nx, self.nz), corrmatrix_ch.reshape(self.nx, self.nx, self.nz), corrmatrix.reshape((self.nx, self.nx, self.nz))


# Example usage
if __name__ == "__main__":
    # Create your 'tobj' object and other necessary inputs
    # Initialize the WaveformStacking class and call the methods
    pass
