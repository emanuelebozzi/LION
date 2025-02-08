#include <stdio.h>
#include <math.h>

void stacking(double *lon_sensors, double *lat_sensors, double *depth_sensors,
              double *itp, double *its, double *stalta_p[nsta][nsamples], double *stalta_s[nsta][nsamples],
              int nx, int nz, int nsta, int nsamples,
              double *corrmatrix, int *iloc, int *itime){

    int nxyz = nx * nx * nz;
    double corrmax = -1.0;
    *iloc = 0;
    *itime = 0;

    for (int i = 0; i < nxyz; i++) {
        double stkmax = -1.0;
        int kmax = 0;

        for (int j = 0; j < nsta; j++) {
            int ttp_val = (int)itp[j];  // Example; you'd actually compute this
            int tts_val = (int)its[j];  // Example; you'd actually compute this

            double stk0p = 0.0;
            double stk0s = 0.0;

            // Loop over samples
            for (int ip = ttp_val; ip < (ttp_val + nsamples); ip++) {
                for (int is_ = tts_val; is_ < (tts_val + nsamples); is_++) {
                    stk0p += stalta_p[j * nsamples + ip];
                    stk0s += stalta_s[j * nsamples + is_];
                }
            }

            // Store the maximum correlation for this sensor
            if (stk0p * stk0s > stkmax) {
                stkmax = stk0p * stk0s;
                kmax = ttp_val;  // Save the correct best matching index (use ttp_val as it corresponds to the start of the range)
            }
        }

        corrmatrix[i] = sqrt(stkmax) / nsta;

        if (corrmatrix[i] > corrmax) {
            corrmax = corrmatrix[i];
            *iloc = i;
            *itime = kmax;
        }
    }
}
