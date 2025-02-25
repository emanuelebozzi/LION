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
int stacking(long int nrs, long int nzs, long int nsta, long int nx, long int ny, long int nz, long int nsamples, long int nxyz, int itp[nrs][nzs], int its[nrs][nzs], double stax[nsta], double stay[nsta], double staz[nsta], double x[nx], double y[ny], double z[nz], double stackf_p[nsta][nsamples], double stackf_s[nsta][nsamples], double corrmatrix[nxyz], int nproc);
/* Python wrapper of the C function stacking */

static char module_docstring[]="Module for computing of the location";
static char stacking_docstring[]="location throug waveform stacking";


/* wrapper */

static PyObject *py_stacking(PyObject *self, PyObject *args){
   PyArrayObject *itp, *its, *stax, *stay, *staz, *x, *y, *z, *stackf_p, *stackf_s, *corrmatrix;
   long int nrs, nzs, nsta, nx, ny, nz, nsamples, nxyz, nproc;
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
   nsta=(long int) PyArray_DIM(stackf_p, 0); 
   nsamples=(long int) PyArray_DIM(stackf_p, 1);
   nx=(long int) PyArray_DIM(x, 0);
   ny=(long int) PyArray_DIM(y, 0);
   nz=(long int) PyArray_DIM(z, 0);
   nrs=(long int) PyArray_DIM(itp, 0);
   nzs=(long int) PyArray_DIM(itp, 1);
   nxyz= nx*ny*nz; 
   dims[0] = nxyz; 
   corrmatrix=(PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_DOUBLE);

   /*call stacking */

   if (stacking(nrs, nzs, nsta, nx, ny, nz, nsamples, nxyz,
      (int (*) [nrs][nzs])PyArray_DATA(itp), (int (*) [nrs][nzs])PyArray_DATA(its),
      (double (*) [nsta])PyArray_DATA(stax), (double (*) [nsta])PyArray_DATA(stay), (double (*) [nsta])PyArray_DATA(staz),
      (double *)PyArray_DATA(x), (double *)PyArray_DATA(y), (double *)PyArray_DATA(z),
      (double (*)[nsta][nsamples])PyArray_DATA(stackf_p), (double (*)[nsta][nsamples])PyArray_DATA(stackf_s),
      (double (*)[nxyz])PyArray_DATA(corrmatrix), nproc) != 0) {
         PyErr_SetString(PyExc_RuntimeError, "Running stacking failed.");
         Py_DECREF(corrmatrix);
         return NULL;
      }

   
   Py_INCREF(corrmatrix);
   return PyArray_Return(corrmatrix);

}


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
  


int stacking(long int nrs, long int nzs, long int nsta, long int nx, long int ny, long int nz, long int nsamples, long int nxyz, int itp[nrs][nzs], int its[nrs][nzs], double stax[nsta], double stay[nsta], double staz[nsta], double x[nx], double y[ny], double z[nz], double stackf_p[nsta][nsamples], double stackf_s[nsta][nsamples], double corrmatrix[nxyz], int nproc){

    long int iter, i, j, k; 
    int ip, is, rdist_ind, zdist_ind, *itime, *iloc, kmax;
    int tp[nsta], ts[nsta];
    double xdist, ydist, rdist, zdist, stk0p, stk0s, stkmax;

    int dx = x[1] - x[0];   /*traveltime table dx*/
    int dz = z[1] - z[0];   /*traveltime table dz*/

    omp_set_num_threads(nproc);
    
    printf(" Location process complete at : %3d %%",0);
    
    iter=0;
    double corrmax = -1;
    
    #pragma omp parallel for shared(iter) private(ip,is,stkmax,stk0p,stk0s,k,j) 
    for(i=0;i<nxyz;i++){
       printf("\b\b\b\b\b%3ld %%", (100*iter++)/(nxyz-2));
       int ix = i / 100;   /*this is the index of the x coordinate*/
       int iy = (i / 10) % 10; /*this is the index of the y coordinate*/
       int iz = i % 10; /*this is the index of the z coordinate*/

       stkmax=-1.0;
       kmax= 0;

       /*For each grid point compute the distance to each station and extract the associated traveltime*/      

       for(j=0;j<nsta;j++){
         xdist=pow(2,(x[ix]-stax[j]));
         ydist=pow(2,(y[iy]-stay[j]));
         zdist=sqrt(pow(2,(z[iz]-staz[j])));
         rdist=sqrt(xdist+ydist); /*this is the distance in meters*/
         rdist_ind = round((double)rdist / dx);  /*round the nearest index*/
         zdist_ind = round((double)zdist / dz);
         tp[j] = itp[rdist_ind][zdist_ind];
         ts[j] = its[rdist_ind][zdist_ind];
      }
       for(k=0;k<nsamples;k++){
         stk0p=0.;
         stk0s=0.;
         for(j=0;j<nsta;j++){
            ip=tp[j] + k;
            is=ts[j] + k;
               if (is<nsamples){
                  stk0p=stackf_p[j][ip] + stk0p;
                  stk0s=stackf_s[j][is] + stk0s;
               }
               else {
                  stk0p=0. + stk0p;
                  stk0s=0. + stk0s;
               }
         }
              if (stk0p*stk0s>stkmax){
                  stkmax=stk0p*stk0s;
                   kmax=k;
              }
     }
     corrmatrix[i]=sqrt(stkmax)/((float) nsta);
     #pragma omp critical
     if (corrmatrix[i]>corrmax){
        corrmax=corrmatrix[i];
        *iloc=i;
        *itime=kmax;
     }

  }
  printf("\n ------ Event located ------ \n");
  return 0;
}