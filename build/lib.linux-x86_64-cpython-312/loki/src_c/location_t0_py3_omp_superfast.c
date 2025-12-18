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
    nx = (long int)PyArray_DIM(x, 0) - 1;
    ny = (long int)PyArray_DIM(y, 0) - 1;
    nz = (long int)PyArray_DIM(z, 0);


    nrs = (int)PyArray_DIM(itp, 0);
    nzs = (int)PyArray_DIM(itp, 1);

    /* nxyz as declared by user arrays (we'll ensure consistency inside stacking()) */
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


/* new section voronoi-inspired grid search (Sambridge, 1999) */

#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>

/* Helper: deduplicate sorted coordinate array with tolerance.
   Inputs:
     arr_in (len n_in) assumed sorted (monotonic) increasing
   Outputs:
     arr_out (allocated, length returned via *n_out). Caller must free arr_out.
*/
static double *dedup_sorted_coords(const double *arr_in, long int n_in, long int *n_out) {
    if (n_in <= 0) {
        *n_out = 0;
        return NULL;
    }
    double tol = 1e-9;
    double *out = (double*) malloc((size_t)n_in * sizeof(double));
    if (!out) { *n_out = 0; return NULL; }
    long int cnt = 0;
    out[cnt++] = arr_in[0];
    for (long int i = 1; i < n_in; i++) {
        double diff = arr_in[i] - out[cnt-1];
        /* use a tolerance relative to the local value if needed */
        double local_tol = tol * (1.0 + fabs(arr_in[i]));
        if (fabs(diff) > local_tol) {
            out[cnt++] = arr_in[i];
        }
    }
    *n_out = cnt;
    /* optionally shrink */
    out = (double*) realloc(out, (size_t)cnt * sizeof(double));
    return out;
}

int stacking(long int nrs, long int nzs, long int nsta, long int nx_in, long int ny_in, long int nz, long int nsamples, long int nxyz_in,
             int itp[nrs][nzs], int its[nrs][nzs],
             double stax[nsta], double stay[nsta], double staz[nsta],
             double x_in[nx_in], double y_in[ny_in], double z[nz],
             double stackf_p[nsta][nsamples], double stackf_s[nsta][nsamples],
             double corrmatrix[nxyz_in], long int iloc[3], long int *itime, long int nproc)
{
    printf("[stacking] Starting function...\n"); fflush(stdout);

    long int i, j, k, ip, is;
    int ix, iy, iz;
    int rdist_ind, zdist_ind;
    double xdist, ydist, rdist, zdist;
    double stk0p, stk0s, stkmax;
    double corrmax = -1.0;

    /* Deduplicate x and y arrays (in case the Python caller passed repeated coords).
       This often fixes accidental extra boundary points (e.g. 20002 instead of 2001).
       We'll produce local arrays x[] and y[] pointing to deduped memory (or to original if no change). */
    long int nx = nx;
    long int ny = ny;
    double *x = NULL;
    double *y = NULL;
    double *x_alloc = NULL;
    double *y_alloc = NULL;

    /* If input arrays are already unique (nx small), dedup will return same length or slightly smaller */
    x_alloc = dedup_sorted_coords(x_in, nx_in, &nx);
    if (x_alloc == NULL) {
        /* fall back to original pointer but copy into a non-const pointer */
        x = (double*) x_in;
    } else {
        x = x_alloc;
    }

    y_alloc = dedup_sorted_coords(y_in, ny_in, &ny);
    if (y_alloc == NULL) {
        y = (double*) y_in;
    } else {
        y = y_alloc;
    }

    /* recompute nxyz safely from deduped nx, ny, and given nz */
    long int nxyz_expected = (long int) nx * (long int) ny * (long int) nz;
    /* If input corrmatrix passed from Python does not match this new size, we must detect and fail early.
       The Python wrapper created corrmatrix with the original nxyz_in = nx_in*ny_in*nz; so it may be larger.
       To avoid mismatch on return, we'll only use the corrmatrix buffer if its size matches expected.
       Otherwise we will print an informative error and exit (safe behavior). */
    long int nxyz_passed = nxyz_in;

    printf("[stacking] dedup: nx_in=%ld -> nx=%ld, ny_in=%ld -> ny=%ld, nz=%ld\n",
           nx_in, nx, ny_in, ny, nz); fflush(stdout);
    printf("[stacking] computed nxyz_expected=%ld, nxyz_passed=%ld\n", nxyz_expected, nxyz_passed); fflush(stdout);

    if (nxyz_expected != nxyz_passed) {
        /* Attempt to avoid hard failure: if passed buffer is larger, we will use only the first nxyz_expected entries
           and leave the remaining entries untouched. But the Python-side array length will still be larger, and reshaping to
           another shape may fail. So the robust choice is to return an error asking caller to pass coherent x,y sizes.
           However, to help automatic cases where input had exactly repeated coordinates that increased size by a small amount,
           we attempt to continue but will warn strongly. */

        fprintf(stderr,
                "[stacking] ERROR: grid size mismatch after deduplication: expected nxyz=%ld but passed buffer length is %ld.\n"
                "This likely means the Python caller's x/y arrays contain extra points. Please pass consistent x,y arrays.\n",
                nxyz_expected, nxyz_passed);
        /* free allocated temp memory and return failure */
        if (x_alloc) free(x_alloc);
        if (y_alloc) free(y_alloc);
        return -2;
    }

    /* Use deduped nx, ny, nxyz_expected from here onwards */
    nxyz_in = nxyz_expected; /* update local var to represent true used length */

    double dx = (nx > 1) ? (x[1] - x[0]) : 1.0;
    double dz = (nz > 1) ? (z[1] - z[0]) : 1.0;

    printf("[stacking] nx=%ld ny=%ld nz=%ld nsamples=%ld nsta=%ld nproc=%ld dx=%g dz=%g\n",
           nx, ny, nz, nsamples, nsta, nproc, dx, dz); fflush(stdout);

    omp_set_num_threads((int)nproc);

    /* Step A: seeds */
    printf("[stacking] Step A: Seed initialization...\n"); fflush(stdout);
    int nseed = (int)floor(sqrt((double)nxyz_in)) * 2;
    if(nseed < 1) nseed = 1;

    long int *seed_indices = (long int*) malloc((size_t)nseed * sizeof(long int));
    double *seed_coherence = (double*) malloc((size_t)nseed * sizeof(double));
    int *seed_ix = (int*) malloc((size_t)nseed * sizeof(int));
    int *seed_iy = (int*) malloc((size_t)nseed * sizeof(int));
    int *seed_iz = (int*) malloc((size_t)nseed * sizeof(int));

    if(!seed_indices || !seed_coherence || !seed_ix || !seed_iy || !seed_iz){
        fprintf(stderr, "[stacking] Memory allocation failed for seeds\n");
        free(seed_indices); free(seed_coherence); free(seed_ix); free(seed_iy); free(seed_iz);
        if (x_alloc) free(x_alloc);
        if (y_alloc) free(y_alloc);
        return -1;
    }

    int ndiv = (int)round(pow((double)nseed, 1.0/3.0));
    if(ndiv < 1) ndiv = 1;
    long int count = 0;
    for(int ixd = 0; ixd < ndiv && count < nseed; ixd++){
        for(int iyd = 0; iyd < ndiv && count < nseed; iyd++){
            for(int izd = 0; izd < ndiv && count < nseed; izd++){
                /* compute seed positions using deduped nx,ny,nz */
                int s_ix = (int)((long)ixd * (nx - 1) / (ndiv - 1 > 0 ? (ndiv - 1) : 1));
                int s_iy = (int)((long)iyd * (ny - 1) / (ndiv - 1 > 0 ? (ndiv - 1) : 1));
                int s_iz = (int)((long)izd * (nz - 1) / (ndiv - 1 > 0 ? (ndiv - 1) : 1));
                seed_ix[count] = s_ix;
                seed_iy[count] = s_iy;
                seed_iz[count] = s_iz;
                seed_indices[count] = (long int)s_ix * (long int)ny * (long int)nz +
                                      (long int)s_iy * (long int)nz + (long int)s_iz;
                seed_coherence[count] = 0.0;
                count++;
            }
        }
    }
    nseed = (int)count;
    printf("[stacking] Step A done: %d seeds created\n", nseed); fflush(stdout);

    /* Step B: coarse stacking */
    printf("[stacking] Step B: Coarse stacking start...\n"); fflush(stdout);

    #pragma omp parallel for private(ix,iy,iz,j,k,ip,is,xdist,ydist,rdist,zdist,rdist_ind,zdist_ind,stk0p,stk0s,stkmax)
    for(int si = 0; si < nseed; si++){
        if(si % 100 == 0) {
            #pragma omp critical
            printf("[stacking]   Seed %d/%d...\n", si, nseed), fflush(stdout);
        }

        ix = seed_ix[si];
        iy = seed_iy[si];
        iz = seed_iz[si];

        long int *tp = (long int*) malloc((size_t)nsta * sizeof(long int));
        long int *ts = (long int*) malloc((size_t)nsta * sizeof(long int));
        if(!tp || !ts){
            seed_coherence[si] = -1.0;
            if(tp) free(tp);
            if(ts) free(ts);
            continue;
        }

        for(j = 0; j < nsta; j++){
            xdist = (x[ix] - stax[j]) * (x[ix] - stax[j]);
            ydist = (y[iy] - stay[j]) * (y[iy] - stay[j]);
            zdist = z[iz];
            rdist = sqrt(xdist + ydist);

            double r_idx = rdist / dx;
            rdist_ind = (int)floor(r_idx);
            zdist_ind = (int)floor(zdist / dz);

            if(rdist_ind < 0) rdist_ind = 0;
            if(rdist_ind >= nrs - 1) rdist_ind = (int)nrs - 2;
            if(zdist_ind < 0) zdist_ind = 0;
            if(zdist_ind >= nzs - 1) zdist_ind = (int)nzs - 2;

            double w = r_idx - (double)rdist_ind;
            double tpi = (1.0 - w) * (double)itp[rdist_ind][zdist_ind] + w * (double)itp[rdist_ind + 1][zdist_ind];
            double tsi = (1.0 - w) * (double)its[rdist_ind][zdist_ind] + w * (double)its[rdist_ind + 1][zdist_ind];

            tp[j] = (long int)(tpi + 0.5);
            ts[j] = (long int)(tsi + 0.5);
        }

        stkmax = -1.0;
        for(k = 0; k < nsamples; k++){
            stk0p = 0.0; stk0s = 0.0;
            for(j = 0; j < nsta; j++){
                ip = tp[j] + k;
                is = ts[j] + k;
                if(ip < 0 || is < 0) continue;
                if(ip >= nsamples || is >= nsamples) continue;
                stk0p += stackf_p[j][ip];
                stk0s += stackf_s[j][is];
            }
            if(stk0p + stk0s > stkmax) stkmax = stk0p + stk0s;
        }
        seed_coherence[si] = stkmax / ((double)nsta);
        free(tp); free(ts);
    }

    printf("[stacking] Step B complete\n"); fflush(stdout);

    /* Step C: Simplified neighbor marking instead of full Voronoi assignment */
    printf("[stacking] Step C: Selecting top seeds and preparing local refinement...\n");
    fflush(stdout);

    /* Sort seeds by coherence */
    int N_top = 25;
    if (nseed < N_top) N_top = nseed;

    int *top_seeds = (int*) malloc((size_t)N_top * sizeof(int));
    double *top_values = (double*) malloc((size_t)N_top * sizeof(double));
    if (!top_seeds || !top_values) {
        fprintf(stderr, "[stacking] Memory allocation failed for top seeds\n");
        free(top_seeds); free(top_values);
        free(seed_indices); free(seed_coherence); free(seed_ix); free(seed_iy); free(seed_iz);
        if (x_alloc) free(x_alloc);
        if (y_alloc) free(y_alloc);
        return -1;
    }
    for (int t = 0; t < N_top; t++) { top_seeds[t] = -1; top_values[t] = -1.0; }

    /* Find N_top highest coherence seeds */
    for (int s = 0; s < nseed; s++) {
        double val = seed_coherence[s];
        for (int t = 0; t < N_top; t++) {
            if (val > top_values[t]) {
                for (int u = N_top - 1; u > t; u--) {
                    top_values[u] = top_values[u-1];
                    top_seeds[u] = top_seeds[u-1];
                }
                top_values[t] = val;
                top_seeds[t] = s;
                break;
            }
        }
    }

    printf("[stacking] Selected top %d seeds for refinement.\n", N_top); fflush(stdout);

    /* initialize corrmatrix to -1 to mark "uncomputed" cells -- buffer length matches nxyz_expected checked earlier */
    #pragma omp parallel for
    for (long int ii = 0; ii < nxyz_in; ii++) corrmatrix[ii] = -1.0;

    /* Mark neighbor zones (radius 'neigh') around top seeds */
    int neigh = 5;
    long int marked_cells = 0;

    for (int t = 0; t < N_top; t++) {
        int sid = top_seeds[t];
        if (sid < 0) continue;

        int ix0 = seed_ix[sid];
        int iy0 = seed_iy[sid];
        int iz0 = seed_iz[sid];

        long int ix_min = max(0, ix0 - neigh);
        long int ix_max = min((int)nx - 1, ix0 + neigh);
        long int iy_min = max(0, iy0 - neigh);
        long int iy_max = min((int)ny - 1, iy0 + neigh);
        long int iz_min = max(0, iz0 - neigh);
        long int iz_max = min((int)nz - 1, iz0 + neigh);

        for (long int ix2 = ix_min; ix2 <= ix_max; ix2++)
            for (long int iy2 = iy_min; iy2 <= iy_max; iy2++)
                for (long int iz2 = iz_min; iz2 <= iz_max; iz2++) {
                    long int idx = ix2 * (long int)ny * (long int)nz + iy2 * (long int)nz + iz2;
                    corrmatrix[idx] = 0.0;  // just mark as "to be refined"
                    marked_cells++;
                }
    }

    printf("[stacking] Marked %ld cells for refinement (~%.6f%% of grid)\n",
        marked_cells, 100.0 * marked_cells / (double)nxyz_in);
    fflush(stdout);

    /* Cleanup Step C temporary arrays */
    free(top_seeds);
    free(top_values);

    /* Step D: refinement */
    printf("[stacking] Step D: Refinement phase start...\n"); fflush(stdout);

    /* refine only marked cells */
    #pragma omp parallel for private(ix,iy,iz,j,k,ip,is,xdist,ydist,rdist,zdist,rdist_ind,zdist_ind,stk0p,stk0s,stkmax)
    for (long int gi = 0; gi < nxyz_in; gi++) {
        if (corrmatrix[gi] < 0.0) continue;  /* skip unmarked cells */

        ix = (int)(gi / (ny * nz));
        iy = (int)((gi / nz) % ny);
        iz = (int)(gi % nz);

        /* do the same micro stacking as coarse but record best stacking value/time */
        long int *tp = (long int*) malloc((size_t)nsta * sizeof(long int));
        long int *ts = (long int*) malloc((size_t)nsta * sizeof(long int));
        if(!tp || !ts){
            if(tp) free(tp);
            if(ts) free(ts);
            continue;
        }

        for(j = 0; j < nsta; j++){
            xdist = (x[ix] - stax[j]) * (x[ix] - stax[j]);
            ydist = (y[iy] - stay[j]) * (y[iy] - stay[j]);
            zdist = z[iz];
            rdist = sqrt(xdist + ydist);

            double r_idx = rdist / dx;
            rdist_ind = (int)floor(r_idx);
            zdist_ind = (int)floor(zdist / dz);

            if(rdist_ind < 0) rdist_ind = 0;
            if(rdist_ind >= nrs - 1) rdist_ind = (int)nrs - 2;
            if(zdist_ind < 0) zdist_ind = 0;
            if(zdist_ind >= nzs - 1) zdist_ind = (int)nzs - 2;

            double w = r_idx - (double)rdist_ind;
            double tpi = (1.0 - w) * (double)itp[rdist_ind][zdist_ind] + w * (double)itp[rdist_ind + 1][zdist_ind];
            double tsi = (1.0 - w) * (double)its[rdist_ind][zdist_ind] + w * (double)its[rdist_ind + 1][zdist_ind];

            tp[j] = (long int)(tpi + 0.5);
            ts[j] = (long int)(tsi + 0.5);
        }

        /* fine stacking across nsamples */
        stkmax = -1.0;
        long int kmax_local = 0;
        for(k = 0; k < nsamples; k++){
            stk0p = 0.0; stk0s = 0.0;
            for(j = 0; j < nsta; j++){
                ip = tp[j] + k;
                is = ts[j] + k;
                if(ip < 0 || is < 0) continue;
                if(ip >= nsamples || is >= nsamples) continue;
                stk0p += stackf_p[j][ip];
                stk0s += stackf_s[j][is];
            }
            if(stk0p + stk0s > stkmax){
                stkmax = stk0p + stk0s;
                kmax_local = k;
            }
        }

        corrmatrix[gi] = stkmax / ((double)nsta);

        /* update global maximum coherency safely */
        #pragma omp critical
        {
            if(corrmatrix[gi] > corrmax){
                corrmax = corrmatrix[gi];
                iloc[0] = gi;
                iloc[1] = ix;
                iloc[2] = iy;
                *itime = kmax_local;
            }
        }

        free(tp); free(ts);
    }

    printf("[stacking] Finished all steps successfully!\n"); fflush(stdout);

    /* Final diagnostics */
    printf("[stacking] Final dims used: nx=%ld ny=%ld nz=%ld nxyz=%ld (passed buffer size=%ld)\n",
           nx, ny, nz, nxyz_in, nxyz_passed);
    fflush(stdout);

    /* free deduped arrays if allocated */
    if (x_alloc) free(x_alloc);
    if (y_alloc) free(y_alloc);

    /* free seed arrays */
    free(seed_indices); free(seed_coherence); free(seed_ix); free(seed_iy); free(seed_iz);

    return 0;
}
