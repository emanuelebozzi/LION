#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#ifndef max
    #define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif
#ifndef min
    #define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

/* Prototypes */
int stacking(long int nrs, long int nzs, long int nsta, long int nx, long int ny, long int nz, long int nsamples, long int nxyz, int itp[nrs][nzs], int its[nrs][nzs], double stax[nsta], double stay[nsta], double staz[nsta], double x[nx], double y[ny], double z[nz], double stackf_p[nsta][nsamples], double stackf_s[nsta][nsamples], double corrmatrix[nxyz], long int iloc[3],long int *itime, long int nproc);

/* Python wrapper of the C function stacking */
static char module_docstring[] = "Module for computing of the location";
static char stacking_docstring[] = "location through waveform stacking";


/* wrapper */

static PyObject *py_stacking(PyObject *self, PyObject *args){
   PyArrayObject *itp, *its, *stax, *stay, *staz, *x, *y, *z, *stackf_p, *stackf_s, *corrmatrix;
   long int nrs, nzs, nsta, nx, ny, nz, nsamples, nxyz, nproc;
   long int iloc[3], itime;
   npy_intp dims[1];
   /* checking the format of the arguments */

   if(!PyArg_ParseTuple(args, "OOOOOOOOOOi", &itp, &its, &stax, &stay, &staz, &x, &y, &z, &stackf_p, &stackf_s, &nproc)){
      PyErr_SetString(PyExc_RuntimeError, "Invalid arguments for the C function stacking");
      return NULL; 
   }

   /* Checking the contiguity of the arrays */

   if(!PyArray_Check(stackf_p) || !PyArray_ISCONTIGUOUS(stackf_p)){
      PyErr_SetString(PyExc_RuntimeError, "stackf_p is not a contiguous array");
      return NULL; 
   }

   if(!PyArray_Check(stackf_s) || !PyArray_ISCONTIGUOUS(stackf_s)){
      PyErr_SetString(PyExc_RuntimeError, "stackf_s is not a contiguous array");
      return NULL; 
   }

   if(!PyArray_Check(stax) || !PyArray_ISCONTIGUOUS(stax)){
      PyErr_SetString(PyExc_RuntimeError, "stax is not a contiguous array");
      return NULL; 
   }

   if(!PyArray_Check(stay) || !PyArray_ISCONTIGUOUS(stay)){
      PyErr_SetString(PyExc_RuntimeError, "stay is not a contiguous array");
      return NULL; 
   }

   if(!PyArray_Check(staz) || !PyArray_ISCONTIGUOUS(staz)){
      PyErr_SetString(PyExc_RuntimeError, "staz is not a contiguous array");
      return NULL; 
   }

   if(!PyArray_Check(itp) || !PyArray_ISCONTIGUOUS(itp)){
      PyErr_SetString(PyExc_RuntimeError, "tp is not a contiguous array");
      return NULL; 
   }

   if(!PyArray_Check(its) || !PyArray_ISCONTIGUOUS(its)){
      PyErr_SetString(PyExc_RuntimeError, "ts is not a contiguous array");
      return NULL; 
   }



   /* Checking that obs_data and stalta are the same type of array and with the same dimensions */

   if((PyArray_NDIM(stackf_p) != 2)){
      PyErr_SetString(PyExc_RuntimeError, "stackf_p is not a 2D array");
      return NULL; 
   }

   if((PyArray_NDIM(stackf_s) != 2)){
      PyErr_SetString(PyExc_RuntimeError, "stackf_s is not a 2D array");
      return NULL; 
   }

   if((PyArray_NDIM(x) != 1)){
      PyErr_SetString(PyExc_RuntimeError, "x is not a 1D array");
      return NULL; 
   }

   if((PyArray_NDIM(y) != 1)){
      PyErr_SetString(PyExc_RuntimeError, "y is not a 1D array");
      return NULL; 
   }

   if((PyArray_NDIM(z) != 1)){
      PyErr_SetString(PyExc_RuntimeError, "z is not a 1D array");
      return NULL; 
   }

   if((PyArray_NDIM(stax) != 1)){
      PyErr_SetString(PyExc_RuntimeError, "stax is not a 1D array");
      return NULL; 
   }

   if((PyArray_NDIM(stay) != 1)){
      PyErr_SetString(PyExc_RuntimeError, "stay is not a 1D array");
      return NULL; 
   }

   if((PyArray_NDIM(staz) != 1)){
      PyErr_SetString(PyExc_RuntimeError, "staz is not a 1D array");
      return NULL; 
   }

   if((PyArray_NDIM(itp) != 2)){
      PyErr_SetString(PyExc_RuntimeError, "tp is not a 2D array");
      return NULL; 
   }

   if((PyArray_NDIM(its) != 2)){
      PyErr_SetString(PyExc_RuntimeError, "ts is not a 2D array");
      return NULL; 
   }

    /* find the dimension of obs_data and stalta */
    nsta = (long int)PyArray_DIM(stackf_p, 0);
    nsamples = (long int)PyArray_DIM(stackf_p, 1);
    nx = (long int)PyArray_DIM(x, 0);
    ny = (long int)PyArray_DIM(y, 0);
    nz = (long int)PyArray_DIM(z, 0);
    nrs = (int)PyArray_DIM(itp, 0);
    nzs = (int)PyArray_DIM(itp, 1);
    nxyz = nx * ny * nz;
    dims[0] = nxyz;
    corrmatrix = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);

    /* call stacking */
    if (0 != stacking(nrs, nzs, nsta, nx, ny, nz, nsamples, nxyz, PyArray_DATA(itp), PyArray_DATA(its),
                 PyArray_DATA(stax), PyArray_DATA(stay), PyArray_DATA(staz),
                 (double *)PyArray_DATA(x), (double *)PyArray_DATA(y), (double *)PyArray_DATA(z),
                 (double (*)[nsamples])PyArray_DATA(stackf_p), (double (*)[nsamples])PyArray_DATA(stackf_s),
                 (double *)PyArray_DATA(corrmatrix), &iloc, &itime, nproc)) {
        PyErr_SetString(PyExc_RuntimeError, "Running stacking failed."); return NULL;
      }

      PyObject *iloctime = Py_BuildValue("(iiii)", iloc[0], iloc[1], iloc[2], itime);;
      /*Py_DECREF(&iloc);*/
      /*Py_DECREF(&itime);*/
       
      PyObject *cohermat=Py_BuildValue("O",corrmatrix);
      Py_DECREF(corrmatrix);
  
      PyObject *locres=Py_BuildValue("OO",iloctime, cohermat);
      Py_DECREF(iloctime);
      Py_DECREF(cohermat);
  
     return locres; }

/* module specifications and inizialization*/

static PyMethodDef module_methods[]={
   /* {method_name, Cfunction, argument_types, docstring} */
      {"stacking", py_stacking, METH_VARARGS, stacking_docstring},
      {NULL, NULL, 0, NULL}
  };
  
  static struct PyModuleDef modlocation_t0 = {
         PyModuleDef_HEAD_INIT,
         "location_t0",
         module_docstring,
         -1,
         module_methods
  };
  
  PyMODINIT_FUNC PyInit_location_t0(void){
      PyObject *m;
      m = PyModule_Create(&modlocation_t0);
      if (m==NULL)
         return NULL;
      import_array();
      return m;
  };
  

int stacking(long int nrs, long int nzs, long int nsta, long int nx, long int ny, long int nz, long int nsamples, long int nxyz,
             int itp[nrs][nzs], int its[nrs][nzs],
             double stax[nsta], double stay[nsta], double staz[nsta],
             double x[nx], double y[ny], double z[nz],
             double stackf_p[nsta][nsamples], double stackf_s[nsta][nsamples],
             double corrmatrix[nxyz], long int iloc[3], long int *itime, long int nproc) {

    long int iter, i, j, k;
    int ix, iy, iz, ip, is;
    int rdist_ind, zdist_ind, kmax;
    long int tp[nsta], ts[nsta];
    double xdist, ydist, rdist, zdist, stk0p, stk0s, stkmax;

    double dx = x[1] - x[0];   /* traveltime table dx */
    double dz = z[1] - z[0];   /* traveltime table dz */

    iter = 0;
    double corrmax = -1.0;

    omp_set_num_threads((int)nproc);

    #pragma omp parallel for shared(corrmax, iloc, itime) private(ix, iy, iz, stkmax, kmax, rdist_ind, zdist_ind, xdist, ydist, zdist, rdist, ip, is, stk0p, stk0s, k, j, tp, ts)
    for (i = 0; i < nxyz; i++) {
        /* compute indices */
        ix = i / (ny * nz);
        iy = (i / nz) % ny;
        iz = i % nz;

        /* initialize */
        stkmax = -1.0;
        kmax = 0;

        for (j = 0; j < nsta; j++) {
            /* compute horizontal and vertical distances */
            xdist = (x[ix] - stax[j]) * (x[ix] - stax[j]);
            ydist = (y[iy] - stay[j]) * (y[iy] - stay[j]);
            zdist = z[iz] - staz[j];
            rdist = sqrt(xdist + ydist);

            /* compute floating r-index and integer indices */
            double r_idx = rdist / dx;
            rdist_ind = (int)floor(r_idx);
            zdist_ind = (int)floor(zdist / dz);
            int rdist_ind_next = rdist_ind + 1;

            /* clamp horizontal indices */
            if (rdist_ind < 0) {
                rdist_ind = 0;
                rdist_ind_next = 1;
                r_idx = 0.0;
            }
            if (rdist_ind >= (int)(nx - 1)) {
                rdist_ind = (int)(nx - 2);
                rdist_ind_next = (int)(nx - 1);
                r_idx = (double)(nx - 2);
            }

            /* clamp vertical index */
            if (zdist_ind < 0) {
                zdist_ind = 0;
            }
            if (zdist_ind >= (int)nz) {
                zdist_ind = (int)(nz - 1);
            }

            double w = r_idx - rdist_ind;
            if (w < 0.0) w = 0.0;
            if (w > 1.0) w = 1.0;

            /* safe interpolation of traveltimes (use ints but clamp) */
            double tpi = (1.0 - w) * (double)itp[rdist_ind][zdist_ind] + w * (double)itp[rdist_ind_next][zdist_ind];
            double tsi = (1.0 - w) * (double)its[rdist_ind][zdist_ind] + w * (double)its[rdist_ind_next][zdist_ind];

            /* clamp and store as long ints */
            long int tpi_i = (long int) (tpi + 0.5);
            long int tsi_i = (long int) (tsi + 0.5);

            if (tpi_i < 0) tpi_i = 0;
            if (tpi_i > nsamples - 1) tpi_i = nsamples - 1;
            if (tsi_i < 0) tsi_i = 0;
            if (tsi_i > nsamples - 1) tsi_i = nsamples - 1;

            tp[j] = tpi_i;
            ts[j] = tsi_i;
        } /* end stations loop */

        /* stacking across samples: ensure index safety when accessing stackf arrays */
        for (k = 0; k < nsamples; k++) {
            stk0p = 0.0;
            stk0s = 0.0;
            for (j = 0; j < nsta; j++) {
                ip = (int)(tp[j] + k);
                is = (int)(ts[j] + k);

                if (ip >= 0 && ip < nsamples) {
                    stk0p += stackf_p[j][ip];
                }
                /* else add 0 */
                if (is >= 0 && is < nsamples) {
                    stk0s += stackf_s[j][is];
                }
            }
            if ((stk0p + stk0s) > stkmax) {
                stkmax = stk0p + stk0s;
                kmax = (int)k;
            }
        } /* end samples loop */

        corrmatrix[i] = stkmax / ((double) nsta);

        /* update global maximum safely */
        #pragma omp critical
        {
            if (corrmatrix[i] > corrmax) {
                corrmax = corrmatrix[i];
                iloc[0] = i;
                iloc[1] = rdist_ind;
                iloc[2] = zdist_ind;
                *itime = kmax;
            }
        }

    } /* end parallel for */

    printf("\n ------ Event located ------ \n");
    return 0;
}
