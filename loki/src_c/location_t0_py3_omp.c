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

    /* ---------------- Progress bookkeeping ---------------- */
    long int work_done = 0;
    int last_pct = -1;

    /* ---------------- Step 0: buffers ---------------- */
    long int *idx_buf = (long int*) malloc((size_t)nxyz*sizeof(long int));
    double  *val_buf = (double *) malloc((size_t)nxyz*sizeof(double));
    if(!idx_buf || !val_buf){ free(idx_buf); free(val_buf); return -1; }

    /* Pre-compute spacing */
    const double dx = (nx>1)?(x[1]-x[0]):1.0;
    const double dz = (nz>1)?(z[1]-z[0]):1.0;
    const double inv_dx = 1.0/dx;
    const double inv_dz = 1.0/dz;

    /* ---------------- Step A: seed generation ---------------- */
    int nseed = (int)floor(sqrt((double)nxyz))*2;
    if(nseed < 1) nseed = 1;

    int *seed_ix = malloc(nseed*sizeof(int));
    int *seed_iy = malloc(nseed*sizeof(int));
    int *seed_iz = malloc(nseed*sizeof(int));
    double *seed_coherence = malloc(nseed*sizeof(double));
    if(!seed_ix||!seed_iy||!seed_iz||!seed_coherence) return -1;

    int ndiv = (int)round(pow((double)nseed,1.0/3.0));
    if(ndiv<1) ndiv=1;

    long int cnt=0;
    for(int ixd=0;ixd<ndiv&&cnt<nseed;ixd++){
        int ix=(int)((long)ixd*(nx-1)/(ndiv-1>0?ndiv-1:1));
        for(int iyd=0;iyd<ndiv&&cnt<nseed;iyd++){
            int iy=(int)((long)iyd*(ny-1)/(ndiv-1>0?ndiv-1:1));
            for(int izd=0;izd<ndiv&&cnt<nseed;izd++){
                int iz=(int)((long)izd*(nz-1)/(ndiv-1>0?ndiv-1:1));
                seed_ix[cnt]=ix; seed_iy[cnt]=iy; seed_iz[cnt]=iz;
                seed_coherence[cnt]=0.0;
                cnt++;
            }
        }
    }
    nseed=(int)cnt;

    /* ---------------- Step B: coarse stacking ---------------- */
    #pragma omp parallel
    {
        long int *tp = malloc(nsta*sizeof(long int));
        long int *ts = malloc(nsta*sizeof(long int));

        #pragma omp for
        for(int s=0;s<nseed;s++){
            int ix=seed_ix[s], iy=seed_iy[s], iz=seed_iz[s];
            double xv=x[ix], yv=y[iy], zv=z[iz];

            for(int j=0;j<nsta;j++){
                double dxs=xv-stax[j], dys=yv-stay[j], dzs=zv-staz[j];
                double r=sqrt(dxs*dxs+dys*dys);
                double ri=r*inv_dx, zi=dzs*inv_dz;
                int r0=(int)floor(ri), z0=(int)floor(zi);
                if(r0<0)r0=0; if(r0>=nrs-1)r0=nrs-2;
                if(z0<0)z0=0; if(z0>=nzs-1)z0=nzs-2;
                double wr=ri-r0, wz=zi-z0;
                double wrc=1.0-wr, wzc=1.0-wz;
                tp[j]=(long int)(wrc*wzc*itp[r0][z0]+wr*wzc*itp[r0+1][z0]+wrc*wz*itp[r0][z0+1]+wr*wz*itp[r0+1][z0+1]+0.5);
                ts[j]=(long int)(wrc*wzc*its[r0][z0]+wr*wzc*its[r0+1][z0]+wrc*wz*its[r0][z0+1]+wr*wz*its[r0+1][z0+1]+0.5);
            }

            double stkmax=-1.0;
            for(int k=0;k<nsamples;k++){
                double sum=0.0;
                for(int j=0;j<nsta;j++){
                    int ip=tp[j]+k, is=ts[j]+k;
                    if(ip>=0&&ip<nsamples&&is>=0&&is<nsamples)
                        sum+=stackf_p[j][ip]+stackf_s[j][is];
                }
                if(sum>stkmax) stkmax=sum;
            }
            seed_coherence[s]=stkmax/(double)nsta;

            /* ---- progress ---- */
            #pragma omp atomic
            work_done++;

            #pragma omp critical
            {
                int pct = (int)(100.0 * work_done / nseed);
                if(pct != last_pct){
                    last_pct = pct;
                    fprintf(stderr,"\rProgress: %d%%",pct);
                    fflush(stderr);
                }
            }
        }
        free(tp); free(ts);
    }

    /* ---------------- Step C: pick top seeds ---------------- */
    int N_top=min(3,nseed);
    int *top_seeds=malloc(N_top*sizeof(int));
    double *top_vals=malloc(N_top*sizeof(double));
    for(int i=0;i<N_top;i++){ top_seeds[i]=-1; top_vals[i]=-1.0; }

    for(int s=0;s<nseed;s++){
        double v=seed_coherence[s];
        for(int t=0;t<N_top;t++){
            if(v>top_vals[t]){
                for(int u=N_top-1;u>t;u--){
                    top_vals[u]=top_vals[u-1];
                    top_seeds[u]=top_seeds[u-1];
                }
                top_vals[t]=v;
                top_seeds[t]=s;
                break;
            }
        }
    }

    /* ---------------- Compute total work for Step D ---------------- */
    int neigh = 20;
    long int total_work = work_done;

    for(int t=0;t<N_top;t++){
        int sid=top_seeds[t];
        if(sid<0) continue;

        int ix0=seed_ix[sid], iy0=seed_iy[sid], iz0=seed_iz[sid];
        int ix_min=max(0,ix0-neigh), ix_max=min((int)nx-1,ix0+neigh);
        int iy_min=max(0,iy0-neigh), iy_max=min((int)ny-1,iy0+neigh);
        int iz_min=max(0,iz0-neigh), iz_max=min((int)nz-1,iz0+neigh);

        total_work +=
            (ix_max-ix_min+1) *
            (iy_max-iy_min+1) *
            (iz_max-iz_min+1);
    }

    /* ---------------- Step D: refinement ---------------- */
    double corrmax=-1.0;
    long int best_idx=0,best_ix=0,best_iy=0,best_time=0;

    for(int t=0;t<N_top;t++){
        int sid=top_seeds[t];
        if(sid<0) continue;

        int ix0=seed_ix[sid], iy0=seed_iy[sid], iz0=seed_iz[sid];
        int ix_min=max(0,ix0-neigh), ix_max=min((int)nx-1,ix0+neigh);
        int iy_min=max(0,iy0-neigh), iy_max=min((int)ny-1,iy0+neigh);
        int iz_min=max(0,iz0-neigh), iz_max=min((int)nz-1,iz0+neigh);
        long int nref=(ix_max-ix_min+1)*(iy_max-iy_min+1)*(iz_max-iz_min+1);

        #pragma omp parallel
        {
            long int *tp=malloc(nsta*sizeof(long int));
            long int *ts=malloc(nsta*sizeof(long int));
            long int *idx_t=malloc(nref*sizeof(long int));
            double *val_t=malloc(nref*sizeof(double));
            long int nnz_t=0;
            double tcorr=-1.0;
            long int tidx=0,tix=0,tiy=0,ttime=0;

            #pragma omp for nowait
            for(long int ri=0;ri<nref;ri++){
                int ix=ix_min+ri/((iy_max-iy_min+1)*(iz_max-iz_min+1));
                int iy=iy_min+(ri/(iz_max-iz_min+1))%(iy_max-iy_min+1);
                int iz=iz_min+ri%(iz_max-iz_min+1);
                long int idx=(long int)ix*ny*nz+(long int)iy*nz+iz;

                double xv=x[ix], yv=y[iy], zv=z[iz];

                for(int j=0;j<nsta;j++){
                    double dxs=xv-stax[j], dys=yv-stay[j], dzs=zv-staz[j];
                    double r=sqrt(dxs*dxs+dys*dys);
                    double ri=r*inv_dx, zi=dzs*inv_dz;
                    int r0=(int)floor(ri), z0=(int)floor(zi);
                    if(r0<0)r0=0; if(r0>=nrs-1)r0=nrs-2;
                    if(z0<0)z0=0; if(z0>=nzs-1)z0=nzs-2;
                    double wr=ri-r0, wz=zi-z0;
                    double wrc=1.0-wr, wzc=1.0-wz;
                    tp[j]=(long int)(wrc*wzc*itp[r0][z0]+wr*wzc*itp[r0+1][z0]+wrc*wz*itp[r0][z0+1]+wr*wz*itp[r0+1][z0+1]+0.5);
                    ts[j]=(long int)(wrc*wzc*its[r0][z0]+wr*wzc*its[r0+1][z0]+wrc*wz*its[r0][z0+1]+wr*wz*its[r0+1][z0+1]+0.5);
                }

                double stkmax=-1.0; long int kmax=0;
                for(int k=0;k<nsamples;k++){
                    double sum=0.0;
                    for(int j=0;j<nsta;j++){
                        int ip=tp[j]+k, is=ts[j]+k;
                        if(ip>=0&&ip<nsamples&&is>=0&&is<nsamples)
                            sum+=stackf_p[j][ip]+stackf_s[j][is];
                    }
                    if(sum>stkmax){ stkmax=sum; kmax=k; }
                }

                double val=stkmax/(double)nsta;
                if(val>0.0){
                    idx_t[nnz_t]=idx;
                    val_t[nnz_t]=val;
                    nnz_t++;
                    if(val>tcorr){ tcorr=val; tidx=idx; tix=ix; tiy=iy; ttime=kmax; }
                }

                #pragma omp atomic
                work_done++;

                #pragma omp critical
                {
                    int pct = (int)(100.0 * work_done / total_work);
                    if(pct != last_pct){
                        last_pct = pct;
                        fprintf(stderr,"\rProgress: %d%%",pct);
                        fflush(stderr);
                    }
                }
            }

            #pragma omp critical
            {
                for(long int i=0;i<nnz_t;i++){
                    idx_buf[nnz_local]=idx_t[i];
                    val_buf[nnz_local]=val_t[i];
                    nnz_local++;
                }
                if(tcorr>corrmax){
                    corrmax=tcorr;
                    best_idx=tidx; best_ix=tix; best_iy=tiy; best_time=ttime;
                }
            }

            free(tp); free(ts); free(idx_t); free(val_t);
        }
    }

    fprintf(stderr,"\n");

    iloc[0]=best_idx;
    iloc[1]=best_ix;
    iloc[2]=best_iy;
    *itime=best_time;

    *out_nnz=nnz_local;
    for(long int i=0;i<nnz_local;i++){
        out_indices[i]=idx_buf[i];
        out_values[i]=val_buf[i];
    }

    free(idx_buf); free(val_buf);
    free(seed_ix); free(seed_iy); free(seed_iz); free(seed_coherence);
    free(top_seeds); free(top_vals);

    return 0;
}
