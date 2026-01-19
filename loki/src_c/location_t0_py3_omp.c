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

/* New prototype: sparse output instead of dense corrmatrix */
int stacking(long int nrs, long int nzs, long int nsta, long int nx, long int ny, long int nz,
             long int nsamples, long int nxyz,
             int itp[nrs][nzs], int its[nrs][nzs],
             double stax[nsta], double stay[nsta], double staz[nsta],
             double x[nx], double y[ny], double z[nz],
             double stackf_p[nsta][nsamples], double stackf_s[nsta][nsamples],
             long int *out_indices, double *out_values, long int *out_nnz,
             long int iloc[3], long int *itime, long int nproc);

/* Python wrapper of the C function stacking */
static char module_docstring[] = "Module for computing of the location";
static char stacking_docstring[] = "location through waveform stacking";

/* wrapper */

static PyObject *py_stacking(PyObject *self, PyObject *args){
   PyArrayObject *itp, *its, *stax, *stay, *staz, *x, *y, *z, *stackf_p, *stackf_s;
   long int nrs, nzs, nsta, nx, ny, nz, nsamples, nxyz, nproc;
   long int iloc[3], itime;
   npy_intp dims_sparse[1];

   /* checking the format of the arguments */
   if(!PyArg_ParseTuple(args, "OOOOOOOOOOl",
                        &itp, &its, &stax, &stay, &staz,
                        &x, &y, &z, &stackf_p, &stackf_s, &nproc)){
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
      PyErr_SetString(PyExc_RuntimeError, "tp (itp) is not a contiguous array");
      return NULL; 
   }

   if(!PyArray_Check(its) || !PyArray_ISCONTIGUOUS(its)){
      PyErr_SetString(PyExc_RuntimeError, "ts (its) is not a contiguous array");
      return NULL; 
   }

   /* Dimension checks */

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
      PyErr_SetString(PyExc_RuntimeError, "itp is not a 2D array");
      return NULL; 
   }

   if((PyArray_NDIM(its) != 2)){
      PyErr_SetString(PyExc_RuntimeError, "its is not a 2D array");
      return NULL; 
   }

   /* find the dimension of arrays */
   nsta     = (long int)PyArray_DIM(stackf_p, 0);
   nsamples = (long int)PyArray_DIM(stackf_p, 1);
   nx       = (long int)PyArray_DIM(x, 0);
   ny       = (long int)PyArray_DIM(y, 0);
   nz       = (long int)PyArray_DIM(z, 0);
   nrs      = (long int)PyArray_DIM(itp, 0);
   nzs      = (long int)PyArray_DIM(itp, 1);
   nxyz     = nx * ny * nz;

   /* Allocate maximum possible size for sparse arrays (worst case: all non-zero) */
   dims_sparse[0] = nxyz;

   PyArrayObject *indices_arr = (PyArrayObject *)PyArray_SimpleNew(1, dims_sparse, NPY_LONG);
   PyArrayObject *values_arr  = (PyArrayObject *)PyArray_SimpleNew(1, dims_sparse, NPY_DOUBLE);

   if (!indices_arr || !values_arr) {
      Py_XDECREF(indices_arr);
      Py_XDECREF(values_arr);
      PyErr_SetString(PyExc_RuntimeError, "Failed to allocate sparse output arrays");
      return NULL;
   }

   long int nnz = 0;

   /* call stacking with sparse outputs */
   if (0 != stacking(nrs, nzs, nsta, nx, ny, nz, nsamples, nxyz,
                     (int (*)[nzs])PyArray_DATA(itp), (int (*)[nzs])PyArray_DATA(its),
                     (double *)PyArray_DATA(stax), (double *)PyArray_DATA(stay), (double *)PyArray_DATA(staz),
                     (double *)PyArray_DATA(x), (double *)PyArray_DATA(y), (double *)PyArray_DATA(z),
                     (double (*)[nsamples])PyArray_DATA(stackf_p), (double (*)[nsamples])PyArray_DATA(stackf_s),
                     (long int *)PyArray_DATA(indices_arr), (double *)PyArray_DATA(values_arr), &nnz,
                     iloc, &itime, nproc)) {

      Py_DECREF(indices_arr);
      Py_DECREF(values_arr);
      PyErr_SetString(PyExc_RuntimeError, "Running stacking failed.");
      return NULL;
   }

   /* Trim arrays to nnz (create new arrays sized nnz, copy data) */
   npy_intp dims_trim[1];
   dims_trim[0] = nnz;

   PyArrayObject *indices_trim = (PyArrayObject *)PyArray_SimpleNew(1, dims_trim, NPY_LONG);
   PyArrayObject *values_trim  = (PyArrayObject *)PyArray_SimpleNew(1, dims_trim, NPY_DOUBLE);

   if (!indices_trim || !values_trim) {
      Py_XDECREF(indices_trim);
      Py_XDECREF(values_trim);
      Py_DECREF(indices_arr);
      Py_DECREF(values_arr);
      PyErr_SetString(PyExc_RuntimeError, "Failed to allocate trimmed sparse arrays");
      return NULL;
   }

   memcpy(PyArray_DATA(indices_trim), PyArray_DATA(indices_arr), nnz * sizeof(long int));
   memcpy(PyArray_DATA(values_trim),  PyArray_DATA(values_arr),  nnz * sizeof(double));

   Py_DECREF(indices_arr);
   Py_DECREF(values_arr);

   PyObject *iloctime = Py_BuildValue("(iiii)", (int)iloc[0], (int)iloc[1], (int)iloc[2], (int)itime);

   PyObject *indices_py = Py_BuildValue("O", indices_trim);
   PyObject *values_py  = Py_BuildValue("O", values_trim);

   Py_DECREF(indices_trim);
   Py_DECREF(values_trim);

   /* Return (iloc+time, indices, values, nnz) */
   PyObject *locres = Py_BuildValue("(OOOO)", iloctime, indices_py, values_py, PyLong_FromLong(nnz));

   Py_DECREF(iloctime);
   Py_DECREF(indices_py);
   Py_DECREF(values_py);

   return locres;
}

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
}

/* new section voronoi-inspired grid search (Sambridge, 1999) */

#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>


int stacking(long int nrs, long int nzs, long int nsta, long int nx, long int ny, long int nz,
             long int nsamples, long int nxyz,
             int itp[nrs][nzs], int its[nrs][nzs],
             double stax[nsta], double stay[nsta], double staz[nsta],
             double x[nx], double y[ny], double z[nz],
             double stackf_p[nsta][nsamples], double stackf_s[nsta][nsamples],
             long int *out_indices, double *out_values, long int *out_nnz,
             long int iloc[3], long int *itime, long int nproc)
{
    omp_set_num_threads((int)nproc);

    long int nnz_local = 0;

    /* Step 0: allocate buffers */
    long int *idx_buf = (long int *) malloc((size_t)nxyz * sizeof(long int));
    double  *val_buf = (double  *) malloc((size_t)nxyz * sizeof(double));
    if (!idx_buf || !val_buf) { free(idx_buf); free(val_buf); return -1; }

    /* Step A: generate seeds on coarse grid */
    int nseed = (int)floor(sqrt((double)nxyz))*2;
    if(nseed < 1) nseed = 1;

    int *seed_ix = (int*) malloc((size_t)nseed*sizeof(int));
    int *seed_iy = (int*) malloc((size_t)nseed*sizeof(int));
    int *seed_iz = (int*) malloc((size_t)nseed*sizeof(int));
    double *seed_coherence = (double*) malloc((size_t)nseed*sizeof(double));
    if(!seed_ix || !seed_iy || !seed_iz || !seed_coherence){
        free(seed_ix); free(seed_iy); free(seed_iz); free(seed_coherence);
        free(idx_buf); free(val_buf); return -1;
    }

    int ndiv = (int)round(pow((double)nseed, 1.0/3.0));
    if(ndiv < 1) ndiv = 1;

    long int count = 0;
    for(int ixd=0; ixd<ndiv && count<nseed; ixd++){
        int ix = (int)((long)ixd*(nx-1)/(ndiv-1>0?(ndiv-1):1));
        for(int iyd=0; iyd<ndiv && count<nseed; iyd++){
            int iy = (int)((long)iyd*(ny-1)/(ndiv-1>0?(ndiv-1):1));
            for(int izd=0; izd<ndiv && count<nseed; izd++){
                int iz = (int)((long)izd*(nz-1)/(ndiv-1>0?(ndiv-1):1));
                seed_ix[count] = ix;
                seed_iy[count] = iy;
                seed_iz[count] = iz;
                seed_coherence[count] = 0.0;
                count++;
            }
        }
    }
    nseed = (int)count;

    /* Step B: coarse stacking for seed values */
    #pragma omp parallel for
    for(int s=0; s<nseed; s++){
        int ix = seed_ix[s];
        int iy = seed_iy[s];
        int iz = seed_iz[s];

        long int *tp = (long int*) malloc(nsta*sizeof(long int));
        long int *ts = (long int*) malloc(nsta*sizeof(long int));
        if(!tp || !ts){ if(tp) free(tp); if(ts) free(ts); seed_coherence[s]=-1.0; continue; }

        for(int j=0; j<nsta; j++){
            double xdist=(x[ix]-stax[j])*(x[ix]-stax[j]);
            double ydist=(y[iy]-stay[j])*(y[iy]-stay[j]);
            double zdist=z[iz]-staz[j];
            double rdist=sqrt(xdist+ydist);

            double dx = (nx>1)?(x[1]-x[0]):1.0;
            double dz = (nz>1)?(z[1]-z[0]):1.0;

            /* bilinear interpolation */
            double r_idx = rdist/dx;
            double z_idx = zdist/dz;
            int r0 = (int)floor(r_idx); if(r0<0) r0=0; if(r0>=nrs-1) r0=nrs-2;
            int z0 = (int)floor(z_idx); if(z0<0) z0=0; if(z0>=nzs-1) z0=nzs-2;
            double wr = r_idx - r0;
            double wz = z_idx - z0;

            double tpi = (1.0-wr)*(1.0-wz)*itp[r0][z0] + wr*(1.0-wz)*itp[r0+1][z0]
                       + (1.0-wr)*wz*itp[r0][z0+1] + wr*wz*itp[r0+1][z0+1];
            double tsi = (1.0-wr)*(1.0-wz)*its[r0][z0] + wr*(1.0-wz)*its[r0+1][z0]
                       + (1.0-wr)*wz*its[r0][z0+1] + wr*wz*its[r0+1][z0+1];

            tp[j]=(long int)(tpi+0.5);
            ts[j]=(long int)(tsi+0.5);
        }

        double stkmax=-1.0;
        for(int k=0; k<nsamples; k++){
            double stk0p=0.0, stk0s=0.0;
            for(int j=0; j<nsta; j++){
                int ip=tp[j]+k; int is=ts[j]+k;
                if(ip<0 || ip>=nsamples) continue;
                if(is<0 || is>=nsamples) continue;
                stk0p+=stackf_p[j][ip];
                stk0s+=stackf_s[j][is];
            }
            if(stk0p+stk0s>stkmax) stkmax=stk0p+stk0s;
        }
        seed_coherence[s]=stkmax/(double)nsta;
        free(tp); free(ts);
    }

    /* Step C: pick top seeds */
    int N_top = min(3,nseed);
    int *top_seeds = (int*) malloc(N_top*sizeof(int));
    double *top_vals = (double*) malloc(N_top*sizeof(double));
    for(int t=0; t<N_top; t++){ top_seeds[t]=-1; top_vals[t]=-1.0; }

    for(int s=0; s<nseed; s++){
        double val = seed_coherence[s];
        for(int t=0; t<N_top; t++){
            if(val>top_vals[t]){
                for(int u=N_top-1; u>t; u--){ top_vals[u]=top_vals[u-1]; top_seeds[u]=top_seeds[u-1]; }
                top_vals[t]=val; top_seeds[t]=s; break;
            }
        }
    }

    /* Step D: final refinement with full bilinear interpolation */
    int neigh = 20;
    double corrmax=-1.0;

    for(int t=0; t<N_top; t++){
        int sid = top_seeds[t];
        if(sid<0) continue;

        int ix0=seed_ix[sid], iy0=seed_iy[sid], iz0=seed_iz[sid];

        int ix_min=max(0,ix0-neigh), ix_max=min((int)nx-1,ix0+neigh);
        int iy_min=max(0,iy0-neigh), iy_max=min((int)ny-1,iy0+neigh);
        int iz_min=max(0,iz0-neigh), iz_max=min((int)nz-1,iz0+neigh);

        long int nref = (ix_max-ix_min+1)*(iy_max-iy_min+1)*(iz_max-iz_min+1);

        #pragma omp parallel for
        for(long int ri=0; ri<nref; ri++){
            int ix = ix_min + ri/( (iy_max-iy_min+1)*(iz_max-iz_min+1) );
            int iy = iy_min + (ri/(iz_max-iz_min+1)) % (iy_max-iy_min+1);
            int iz = iz_min + ri % (iz_max-iz_min+1);

            long int idx = (long int)ix*ny*nz + (long int)iy*nz + (long int)iz;

            long int *tp = (long int*) malloc(nsta*sizeof(long int));
            long int *ts = (long int*) malloc(nsta*sizeof(long int));
            if(!tp || !ts){ if(tp) free(tp); if(ts) free(ts); continue; }

            for(int j=0; j<nsta; j++){
                double xdist=(x[ix]-stax[j])*(x[ix]-stax[j]);
                double ydist=(y[iy]-stay[j])*(y[iy]-stay[j]);
                double zdist=z[iz]-staz[j];
                double rdist=sqrt(xdist+ydist);

                double dx = (nx>1)?(x[1]-x[0]):1.0;
                double dz = (nz>1)?(z[1]-z[0]):1.0;

                double r_idx = rdist/dx;
                double z_idx = zdist/dz;
                int r0 = (int)floor(r_idx); if(r0<0) r0=0; if(r0>=nrs-1) r0=nrs-2;
                int z0 = (int)floor(z_idx); if(z0<0) z0=0; if(z0>=nzs-1) z0=nzs-2;
                double wr = r_idx-r0; double wz = z_idx-z0;

                double tpi = (1.0-wr)*(1.0-wz)*itp[r0][z0] + wr*(1.0-wz)*itp[r0+1][z0]
                           + (1.0-wr)*wz*itp[r0][z0+1] + wr*wz*itp[r0+1][z0+1];
                double tsi = (1.0-wr)*(1.0-wz)*its[r0][z0] + wr*(1.0-wz)*its[r0+1][z0]
                           + (1.0-wr)*wz*its[r0][z0+1] + wr*wz*its[r0+1][z0+1];

                tp[j]=(long int)(tpi+0.5);
                ts[j]=(long int)(tsi+0.5);
            }

            double stkmax=-1.0;
            long int kmax=0;
            for(int k=0; k<nsamples; k++){
                double stk0p=0.0, stk0s=0.0;
                for(int j=0; j<nsta; j++){
                    int ip=tp[j]+k; int is=ts[j]+k;
                    if(ip<0||ip>=nsamples) continue;
                    if(is<0||is>=nsamples) continue;
                    stk0p+=stackf_p[j][ip];
                    stk0s+=stackf_s[j][is];
                }
                if(stk0p+stk0s>stkmax){ stkmax=stk0p+stk0s; kmax=k; }
            }

            double val = stkmax/(double)nsta;
            if(val!=0.0){
                #pragma omp critical
                {
                    idx_buf[nnz_local]=idx;
                    val_buf[nnz_local]=val;
                    nnz_local++;
                    if(val>corrmax){ corrmax=val; iloc[0]=idx; iloc[1]=ix; iloc[2]=iy; *itime=kmax; }
                }
            }

            free(tp); free(ts);
        }
    }

    *out_nnz = nnz_local;
    for(long int ii=0; ii<nnz_local; ii++){
        out_indices[ii]=idx_buf[ii];
        out_values[ii]=val_buf[ii];
    }

    free(idx_buf); free(val_buf);
    free(seed_ix); free(seed_iy); free(seed_iz); free(seed_coherence);
    free(top_seeds); free(top_vals);

    return 0;
}

