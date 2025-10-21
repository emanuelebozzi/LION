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
  
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

typedef struct {
    int ix, iy, iz;
    int itime;
    double score;
} Sample;

double evaluate_point(int ix, int iy, int iz, int nsta, int nsamples, double x[], double y[], double z[], 
                      double stax[], double stay[], double staz[], int nrs, int nzs, 
                      int itp[nrs][nzs], int its[nrs][nzs],
                      double stackf_p[nsta][nsamples], double stackf_s[nsta][nsamples],
                      int *best_time)
{
    int tp[nsta], ts[nsta];
    double dx = x[1]-x[0];
    double dz = z[1]-z[0];
    double stkmax = -1.0;

    for (int j=0;j<nsta;j++){
        double rdist = sqrt(pow(x[ix]-stax[j],2) + pow(y[iy]-stay[j],2));
        double zdist = z[iz];
        int rdist_ind = floor(rdist/dx);
        int rdist_ind_next = fmin(rdist_ind+1, nrs-1);
        rdist_ind = fmin(rdist_ind, nrs-2);
        int zdist_ind = fmin((int)floor(zdist/dz), nzs-1);
        double w = rdist/dx - rdist_ind;
        double tpi = (1-w)*itp[rdist_ind][zdist_ind] + w*itp[rdist_ind_next][zdist_ind];
        double tsi = (1-w)*its[rdist_ind][zdist_ind] + w*its[rdist_ind_next][zdist_ind];
        tp[j] = (int)(tpi + 0.5);
        ts[j] = (int)(tsi + 0.5);
    }

    int best_k = 0;
    for (int k=0;k<nsamples;k++){
        double stk0p=0, stk0s=0;
        for (int j=0;j<nsta;j++){
            int ip = tp[j]+k;
            int is = ts[j]+k;
            if (ip<nsamples) stk0p += stackf_p[j][ip];
            if (is<nsamples) stk0s += stackf_s[j][is];
        }
        if (stk0p+stk0s>stkmax){
            stkmax = stk0p+stk0s;
            best_k = k;
        }
    }
    *best_time = best_k;
    return stkmax / nsta;
}

int stacking(long int nrs, long int nzs, long int nsta, long int nx, long int ny, long int nz, 
             long int nsamples, long int nxyz, int itp[nrs][nzs], int its[nrs][nzs], 
             double stax[nsta], double stay[nsta], double staz[nsta], double x[nx], double y[ny], double z[nz], 
             double stackf_p[nsta][nsamples], double stackf_s[nsta][nsamples], 
             double corrmatrix[nxyz], long int iloc[3], long int *itime, long int nproc)
{
    int initial_samples = 100;      // few random starting points
    int max_total_samples = 5000;  // total samples
    int iter = 0;
    srand(time(NULL));
    omp_set_num_threads(nproc);

    Sample *samples = malloc(max_total_samples*sizeof(Sample));
    int ns = 0;

    double corrmax = -1.0;

    // Step 1: initial random points
    #pragma omp parallel for
    for (int s=0;s<initial_samples;s++){
        int ix = rand()%nx;
        int iy = rand()%ny;
        int iz = rand()%nz;
        int best_time;
        double score = evaluate_point(ix, iy, iz, nsta, nsamples, x, y, z, stax, stay, staz,
                                      nrs, nzs, itp, its, stackf_p, stackf_s, &best_time);

        #pragma omp critical
        {
            samples[ns].ix = ix;
            samples[ns].iy = iy;
            samples[ns].iz = iz;
            samples[ns].itime = best_time;
            samples[ns].score = score;
            ns++;

            if (score > corrmax){
                corrmax = score;
                iloc[0] = ix; iloc[1] = iy; iloc[2] = iz;
                *itime = best_time;
            }
            corrmatrix[ix*ny*nz + iy*nz + iz] = score;
        }
    }

    // Step 2: iterative neighborhood sampling
    while (ns < max_total_samples){
        // pick a sample weighted by score for neighborhood
        int idx = rand() % ns;
        Sample parent = samples[idx];

        // define neighborhood size (can shrink with iteration)
        int nh = fmax(1, nx/10);  // example neighborhood
        int ix = fmin(nx-1, fmax(0, parent.ix + (rand()%(2*nh+1)-nh)));
        int iy = fmin(ny-1, fmax(0, parent.iy + (rand()%(2*nh+1)-nh)));
        int iz = fmin(nz-1, fmax(0, parent.iz + (rand()%(2*nh+1)-nh)));

        int best_time;
        double score = evaluate_point(ix, iy, iz, nsta, nsamples, x, y, z, stax, stay, staz,
                                      nrs, nzs, itp, its, stackf_p, stackf_s, &best_time);

        samples[ns].ix = ix;
        samples[ns].iy = iy;
        samples[ns].iz = iz;
        samples[ns].itime = best_time;
        samples[ns].score = score;
        ns++;

        corrmatrix[ix*ny*nz + iy*nz + iz] = score;

        if (score > corrmax){
            corrmax = score;
            iloc[0] = ix; iloc[1] = iy; iloc[2] = iz;
            *itime = best_time;
        }
    }

    free(samples);
    printf("\n ------ Event located ------ \n");
    return 0;
}
