# waveform_stacking.pyx

import numpy as np
cimport numpy as np

# Cythonized function to compute the stack
cpdef tuple compute_sensor_stack(int j, np.ndarray stalta_p, np.ndarray stalta_s, int nsamples, int ttp_val, int tts_val):
    cdef double stk0p = 0.0
    cdef double stk0s = 0.0
    cdef int ip, is_
    
    # Using numpy's slicing to avoid loops in Python
    for ip in range(ttp_val, min(ttp_val + nsamples, nsamples)):
        for is_ in range(tts_val, min(tts_val + nsamples, nsamples)):
            stk0p += stalta_p[j, ip]
            stk0s += stalta_s[j, is_]

    return stk0p, stk0s
