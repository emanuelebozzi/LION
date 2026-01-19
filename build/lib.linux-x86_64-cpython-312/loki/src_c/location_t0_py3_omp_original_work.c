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
int stacking(long int nxyz, long int nsta, long int nsamples, int itp[nxyz][nsta], int its[nxyz][nsta], double stalta_p[nsta][nsamples], double stalta_s[nsta][nsamples], double corrmatrix[nxyz],long int *iloc,long int *itime, int nproc);
/* Python wrapper of the C function stacking */

static char module_docstring[]="Module for computing of the location";
static char stacking_docstring[]="location throug waveform stacking";


/* wrapper */

static PyObject *py_stacking(PyObject *self, PyObject *args){
   PyArrayObject *itp, *its, *stalta_p, *stalta_s, *corrmatrix;
   long int nxyz, nsamples, nsta, nproc;
	 long int iloc, itime;
   npy_intp dims[1];
   /* checking the format of the arguments */

      PyErr_SetString(PyExc_RuntimeError, "Invalid arguments for the C function stacking");
      return NULL; ng(PyExc_RuntimeError, "Invalid arguments for the C function stacking");
   }  return NULL;
   }
   /* Checking the contiguity of the arrays */
   /* Checking the contiguity of the arrays */
   if(!PyArray_Check(stackf_p) || !PyArray_ISCONTIGUOUS(stackf_p)){
      PyErr_SetString(PyExc_RuntimeError, "stackf_p is not a contiguous array");
      return NULL; ng(PyExc_RuntimeError, "stalta_p is not a contiguous array");
   }  return NULL;
   }
   if(!PyArray_Check(stackf_s) || !PyArray_ISCONTIGUOUS(stackf_s)){
      PyErr_SetString(PyExc_RuntimeError, "stackf_s is not a contiguous array");
      return NULL; ng(PyExc_RuntimeError, "stalta_s is not a contiguous array");
   }  return NULL;
   }
   if(!PyArray_Check(stax) || !PyArray_ISCONTIGUOUS(stax)){
      PyErr_SetString(PyExc_RuntimeError, "stax is not a contiguous array");
      return NULL; ng(PyExc_RuntimeError, "tp is not a contiguous array");
   }  return NULL;
   }
   if(!PyArray_Check(stay) || !PyArray_ISCONTIGUOUS(stay)){
      PyErr_SetString(PyExc_RuntimeError, "stay is not a contiguous array");
      return NULL; ng(PyExc_RuntimeError, "ts is not a contiguous array");
   }  return NULL;
   }
   if(!PyArray_Check(staz) || !PyArray_ISCONTIGUOUS(staz)){
      PyErr_SetString(PyExc_RuntimeError, "staz is not a contiguous array");
      return NULL; 
   }   /* Checking that obs_data and stalta are the same type of array and with the same dimensions */

   if(!PyArray_Check(itp) || !PyArray_ISCONTIGUOUS(itp)){
      PyErr_SetString(PyExc_RuntimeError, "tp is not a contiguous array");ng(PyExc_RuntimeError, "stalta_p is not a 2D array");
      return NULL;   return NULL;
   }   }

   if(!PyArray_Check(its) || !PyArray_ISCONTIGUOUS(its)){
      PyErr_SetString(PyExc_RuntimeError, "ts is not a contiguous array");ng(PyExc_RuntimeError, "stalta_s is not a 2D array");
      return NULL;   return NULL;
   }   }


ng(PyExc_RuntimeError, "tp is not a 2D array");
   /* Checking that obs_data and stalta are the same type of array and with the same dimensions */  return NULL;
   }
   if((PyArray_NDIM(stackf_p) != 2)){
      PyErr_SetString(PyExc_RuntimeError, "stackf_p is not a 2D array");
      return NULL; ng(PyExc_RuntimeError, "ts is not a 2D array");
   }  return NULL;
   }
   if((PyArray_NDIM(stackf_s) != 2)){
      PyErr_SetString(PyExc_RuntimeError, "stackf_s is not a 2D array");ta */
      return NULL; 
   }ta_p, 1);
PyArray_DIM(itp, 0);
   if((PyArray_NDIM(x) != 1)){
      PyErr_SetString(PyExc_RuntimeError, "x is not a 1D array");   corrmatrix=(PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_DOUBLE);
      return NULL; 
   }
DATA(its), PyArray_DATA(stalta_p), PyArray_DATA(stalta_s), PyArray_DATA(corrmatrix), &iloc, &itime, nproc)) {
   if((PyArray_NDIM(y) != 1)){ing(PyExc_RuntimeError, "running stacking failed.");
      PyErr_SetString(PyExc_RuntimeError, "y is not a 1D array");  return NULL;
      return NULL;    }
   }
BuildValue("(i,i)", iloc, itime);
   if((PyArray_NDIM(z) != 1)){
      PyErr_SetString(PyExc_RuntimeError, "z is not a 1D array");	 /*Py_DECREF(&itime);*/
      return NULL;     
   }	 PyObject *cohermat=Py_BuildValue("O",corrmatrix);
	 Py_DECREF(corrmatrix);
   if((PyArray_NDIM(stax) != 1)){
      PyErr_SetString(PyExc_RuntimeError, "stax is not a 1D array");e, cohermat);
      return NULL; 	 Py_DECREF(iloctime);
   }

   if((PyArray_NDIM(stay) != 1)){
      PyErr_SetString(PyExc_RuntimeError, "stay is not a 1D array");
      return NULL; 
   }

   if((PyArray_NDIM(staz) != 1)){
      PyErr_SetString(PyExc_RuntimeError, "staz is not a 1D array");ations and inizialization*/
      return NULL; 
   }ethodDef module_methods[]={
nction, argument_types, docstring} */
   if((PyArray_NDIM(itp) != 2)){  {"stacking", py_stacking, METH_VARARGS, stacking_docstring},
      PyErr_SetString(PyExc_RuntimeError, "tp is not a 2D array");    {NULL, NULL, 0, NULL}
      return NULL; 
   }
t0 = {
   if((PyArray_NDIM(its) != 2)){ef_HEAD_INIT,
      PyErr_SetString(PyExc_RuntimeError, "ts is not a 2D array");",
      return NULL; ring,
   }
     module_methods
   /* find the dimension of obs_data and stalta */};
   nsta=(long int) PyArray_DIM(stackf_p, 0); 
   nsamples=(long int) PyArray_DIM(stackf_p, 1);PyMODINIT_FUNC PyInit_location_t0(void){
   nx=(long int) PyArray_DIM(x, 0);
   ny=(long int) PyArray_DIM(y, 0);    m = PyModule_Create(&modlocation_t0);
   nz=(long int) PyArray_DIM(z, 0);
   nrs=(long int) PyArray_DIM(itp, 0);ULL;
   nzs=(long int) PyArray_DIM(itp, 1);
   nxyz= nx*ny*nz; m;
   dims[0] = nxyz; };
   corrmatrix=(PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_DOUBLE);

   /*call stacking */
], int its[nxyz][nsta], double stalta_p[nsta][nsamples], double stalta_s[nsta][nsamples], double corrmatrix[nxyz], long int *iloc, long int *itime, int nproc){
   if (stacking(nrs, nzs, nsta, nx, ny, nz, nsamples, nxyz,
      (int *)PyArray_DATA(itp), (int *)PyArray_DATA(its),
      (double *)PyArray_DATA(stax), (double *)PyArray_DATA(stay), (double *)PyArray_DATA(staz),
      (double *)PyArray_DATA(x), (double *)PyArray_DATA(y), (double *)PyArray_DATA(z),, corrmax;
      (double *)PyArray_DATA(stackf_p), (double *)PyArray_DATA(stackf_s),
      (double *)PyArray_DATA(corrmatrix), nproc) != 0) {
         PyErr_SetString(PyExc_RuntimeError, "Running stacking failed.");
         Py_DECREF(corrmatrix);
         return NULL;omplete at : %3d %%",0);
      }d(iter,corrmax,iloc,itime) private(ip,is,stkmax,stk0p,stk0s,kmax,k,j)

   xyz-2));
   Py_INCREF(corrmatrix);;
   return PyArray_Return(corrmatrix);

}
.;
or(j=0;j<nsta;j++){
/* module specifications and initialization*/
      is=its[i][j] + k;
static PyMethodDef module_methods[]={
 /* {method_name, Cfunction, argument_types, docstring} */               stk0p=stalta_p[j][ip] + stk0p;
    {"stacking", py_stacking, METH_VARARGS, stacking_docstring},s;
    {NULL, NULL, 0, NULL}    }
};                else {
























































































}}  return 0;  printf("\n ------ Event located ------ \n");  }     }        itime=kmax;        iloc=i;        corrmax=corrmatrix[i];     if (corrmatrix[i]>corrmax){     #pragma omp critical     corrmatrix[i]=sqrt(stkmax)/((float) nsta);     }              }                   kmax=k;                  stkmax=stk0p*stk0s;                 if (stk0p*stk0s>stkmax){                 }                    stk0s=0. + stk0s;                    stk0p=0. + stk0p;                 else {                 }                    stk0s=stackf_s[j][is] + stk0s;                    stk0p=stackf_p[j][ip] + stk0p;                 if (is<nsamples){              is=ts[j] + k;              ip=tp[j] + k;           for(j=0;j<nsta;j++){           stk0s=0.;           stk0p=0.;       for(k=0;k<nsamples;k++){      }         ts[j] = its[rdist_ind][zdist_ind];         tp[j] = itp[rdist_ind][zdist_ind];         zdist_ind = round((double)zdist / dz);         rdist_ind = round((double)rdist / dx);  /*round the nearest index*/         rdist=sqrt(xdist+ydist); /*this is the distance in meters*/         zdist=sqrt(pow(2,(z[iz]-staz[j])));         ydist=pow(2,(y[iy]-stay[j]));         xdist=pow(2,(x[ix]-stax[j]));       for(j=0;j<nsta;j++){       /*For each grid point compute the distance to each station and extract the associated traveltime*/             stkmax=0.;       int iz = i % 10; /*this is the index of the z coordinate*/       int iy = (i / 10) % 10; /*this is the index of the y coordinate*/       int ix = i / 100;   /*this is the index of the x coordinate*/       printf("\b\b\b\b\b%3ld %%", (100*iter++)/(nxyz-2));    for(i=0;i<nxyz;i++){    #pragma omp parallel for shared(iter) private(ip,is,stkmax,stk0p,stk0s,k,j)     printf(" Location process complete at : %3d %%",0);        omp_set_num_threads(nproc);    int itime = 0;    int iloc = 0;    double corrmax = 0.0;    int kmax = 0;    int dz = z[1] - z[0];   /*traveltime table dz*/    int dx = x[1] - x[0];   /*traveltime table dx*/    iter=0;    double xdist, ydist, rdist, zdist, stk0p, stk0s, stkmax;    int ip, tp, is, ts, rdist_ind, zdist_ind;    long int iter, i, j, k; int stacking(long int nrs, long int nzs, long int nsta, long int nx, long int ny, long int nz, long int nsamples, long int nxyz, int itp[nrs][nzs], int its[nrs][nzs], double stax[nsta], double stay[nsta], double staz[nsta], double x[nx], double y[ny], double z[nz], double stackf_p[nsta][nsamples], double stackf_s[nsta][nsamples], double corrmatrix[nxyz], int nproc){}    return m;    import_array();       return NULL;    if (m==NULL)    m = PyModule_Create(&modlocation);    PyObject *m;PyMODINIT_FUNC PyInit_location(void){};       module_methods       -1,       module_docstring,       "location",       PyModuleDef_HEAD_INIT,static struct PyModuleDef modlocation = {                    stk0p=0. + stk0p;
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