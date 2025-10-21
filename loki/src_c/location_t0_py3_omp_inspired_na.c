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
  


/* ---- Seed structure for Voronoi refinement ---- */
typedef struct {
    int ix, iy, iz;
    double coh;
} SeedPoint;


/* ---- Progress bar helper ---- */
void print_progress(double progress, const char *stage)
{
    int barWidth = 50;
    int pos = (int)(barWidth * progress);
    printf("\r[%s] [", stage);
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) printf("=");
        else if (i == pos) printf(">");
        else printf(" ");
    }
    printf("] %3.0f%%", progress * 100.0);
    fflush(stdout);
}


/* ---- Stacking function (replacement) ---- */
int stacking(long int nrs, long int nzs, long int nsta, long int nx, long int ny, long int nz, long int nsamples,
             long int nxyz, int itp[nrs][nzs], int its[nrs][nzs],
             double stax[nsta], double stay[nsta], double staz[nsta],
             double x[nx], double y[ny], double z[nz],
             double stackf_p[nsta][nsamples], double stackf_s[nsta][nsamples],
             double corrmatrix[nxyz], long int iloc[3], long int *itime, long int nproc)
{
    long int i, j, k;
    long int tp[nsta], ts[nsta];
    double dx = x[1] - x[0];
    double dz = z[1] - z[0];
    long int best_ix = 0, best_iy = 0, best_iz = 0, best_k = 0;

    omp_set_num_threads(nproc);

    /* initialize corrmatrix */
    for (i = 0; i < nxyz; i++) corrmatrix[i] = 0.0;

    /* ---- COARSE GRID SAMPLING ---- */
    long int stride = nx / 20;
    if (stride < 1) stride = 1;
    double best_coh = -1.0;
    long int total_points = ((nx + stride - 1) / stride) * ((ny + stride - 1) / stride) * ((nz + stride - 1) / stride);
    long int processed = 0;

    #pragma omp parallel for collapse(3) private(j,k,tp,ts) shared(best_coh,best_ix,best_iy,best_iz,best_k,processed)
    for (long int ix = 0; ix < nx; ix += stride) {
        for (long int iy = 0; iy < ny; iy += stride) {
            for (long int iz = 0; iz < nz; iz += stride) {
                double stkmax = -1.0;
                int kmax = 0;

                for (j = 0; j < nsta; j++) {
                    double xdist = pow(x[ix] - stax[j], 2);
                    double ydist = pow(y[iy] - stay[j], 2);
                    double zdist = z[iz];
                    double rdist = sqrt(xdist + ydist);

                    double r_idx = rdist / dx;
                    int rdist_ind = (int)floor(r_idx);
                    int zdist_ind = (int)floor(zdist / dz);
                    int rdist_ind_next = rdist_ind + 1;
                    if (rdist_ind >= nx - 1) { rdist_ind = nx - 2; rdist_ind_next = nx - 1; }
                    if (zdist_ind > nz - 1) zdist_ind = nz - 1;

                    double w = r_idx - rdist_ind;
                    double tpi = (1 - w) * itp[rdist_ind][zdist_ind] + w * itp[rdist_ind_next][zdist_ind];
                    double tsi = (1 - w) * its[rdist_ind][zdist_ind] + w * its[rdist_ind_next][zdist_ind];
                    tp[j] = (long int)(tpi + 0.5);
                    ts[j] = (long int)(tsi + 0.5);
                }

                for (k = 0; k < nsamples; k++) {
                    double stk0p = 0.0, stk0s = 0.0;
                    for (j = 0; j < nsta; j++) {
                        long int ip = tp[j] + k;
                        long int is = ts[j] + k;
                        if (is < nsamples) {
                            stk0p += stackf_p[j][ip];
                            stk0s += stackf_s[j][is];
                        }
                    }
                    if (stk0p + stk0s > stkmax) {
                        stkmax = stk0p + stk0s;
                        kmax = k;
                    }
                }

                double coh = stkmax / (double)nsta;
                corrmatrix[ix * (ny * nz) + iy * nz + iz] = coh;

                #pragma omp critical
                {
                    processed++;
                    if (processed % 1000 == 0 || processed == total_points) {
                        print_progress((double)processed / total_points, "Coarse grid");
                    }
                    if (coh > best_coh) {
                        best_coh = coh;
                        best_ix = ix; best_iy = iy; best_iz = iz;
                        best_k = kmax;
                    }
                }
            }
        }
    }
    print_progress(1.0, "Coarse grid");
    printf("\n");

    /* ---- collect coarse-sampled points (sequential) ---- */
    typedef struct { int ix, iy, iz; double coh; } CoarsePoint;
    long int max_coarse = ((nx + stride - 1) / stride) * ((ny + stride - 1) / stride) * ((nz + stride - 1) / stride);
    CoarsePoint *cp = (CoarsePoint *) malloc(sizeof(CoarsePoint) * (size_t)max_coarse);
    if (!cp) {
        fprintf(stderr, "malloc failed for coarse points\n");
        return -1;
    }
    long int cp_cnt = 0;
    for (long int ix = 0; ix < nx; ix += stride) {
        for (long int iy = 0; iy < ny; iy += stride) {
            for (long int iz = 0; iz < nz; iz += stride) {
                long int idx = ix * (ny * nz) + iy * nz + iz;
                double coh = corrmatrix[idx];
                if (coh != 0.0) {
                    cp[cp_cnt].ix = (int)ix;
                    cp[cp_cnt].iy = (int)iy;
                    cp[cp_cnt].iz = (int)iz;
                    cp[cp_cnt].coh = coh;
                    cp_cnt++;
                }
            }
        }
    }

    /* sort coarse points by coh desc */
    int cmp_coh_desc(const void *a, const void *b) {
        double ca = ((const CoarsePoint *)a)->coh;
        double cb = ((const CoarsePoint *)b)->coh;
        if (ca < cb) return 1;
        if (ca > cb) return -1;
        return 0;
    }
    qsort(cp, (size_t)cp_cnt, sizeof(CoarsePoint), cmp_coh_desc);

    /* pick top seeds from coarse result */
    int top_n_seeds = 10; /* tune: number of starting seeds */
    if (top_n_seeds > cp_cnt) top_n_seeds = (int)cp_cnt;

    /* evaluated mask to avoid recompute */
    unsigned char *evaluated = (unsigned char *) calloc((size_t)nxyz, sizeof(unsigned char));
    if (!evaluated) {
        free(cp);
        fprintf(stderr, "calloc failed for evaluated mask\n");
        return -1;
    }

    /* mark coarse-sampled points evaluated */
    for (long int ii = 0; ii < cp_cnt; ii++) {
        int ix = cp[ii].ix, iy = cp[ii].iy, iz = cp[ii].iz;
        long int idx = (long int)ix * (ny * nz) + iy * nz + iz;
        evaluated[idx] = 1;
    }

    /* build initial seeds array from top coarse points */
    int max_seeds = 200; /* upper cap, tune as needed */
    if (max_seeds < top_n_seeds) max_seeds = top_n_seeds;
    SeedPoint *seeds = (SeedPoint *) malloc(sizeof(SeedPoint) * (size_t)max_seeds);
    if (!seeds) {
        free(cp); free(evaluated);
        fprintf(stderr, "malloc failed for seeds\n");
        return -1;
    }
    int nseeds = top_n_seeds;
    for (int s = 0; s < nseeds; s++) {
        seeds[s].ix = cp[s].ix;
        seeds[s].iy = cp[s].iy;
        seeds[s].iz = cp[s].iz;
        seeds[s].coh = cp[s].coh;
    }

    /* free coarse points buffer (not needed anymore) */
    free(cp);

    /* Evaluation budget to avoid runaway */
    long int eval_count = 0;
    long int max_evals = total_points * 10 + 10000; /* tune: how much extra work allowed */
    if (max_evals > nxyz) max_evals = nxyz;

    /* helper lambda-like inline: evaluate a single grid point (updates corrmatrix and evaluated) */
    /* We'll implement as a small inner function via macros to allow access to outer variables */
#define EVALUATE_POINT_AND_STORE(ix_,iy_,iz_, out_coh, out_k) do {                      \
    long int _ix = (ix_), _iy = (iy_), _iz = (iz_);                                      \
    long int _idx = _ix * (ny * nz) + _iy * nz + _iz;                                    \
    if (!evaluated[_idx]) {                                                              \
        double _stkmax = -1.0; int _kmax = 0;                                            \
        for (j = 0; j < nsta; j++) {                                                     \
            double _xdist = pow(x[_ix] - stax[j], 2);                                   \
            double _ydist = pow(y[_iy] - stay[j], 2);                                   \
            double _zdist = z[_iz];                                                      \
            double _rdist = sqrt(_xdist + _ydist);                                      \
            double _r_idx = _rdist / dx;                                                \
            int _rdist_ind = (int)floor(_r_idx);                                        \
            int _zdist_ind = (int)floor(_zdist / dz);                                   \
            int _rdist_ind_next = _rdist_ind + 1;                                       \
            if (_rdist_ind >= nx - 1) { _rdist_ind = nx - 2; _rdist_ind_next = nx - 1; }\
            if (_zdist_ind > nz - 1) _zdist_ind = nz - 1;                               \
            double _w = _r_idx - _rdist_ind;                                            \
            double _tpi = (1 - _w) * itp[_rdist_ind][_zdist_ind] +                       \
                           _w * itp[_rdist_ind_next][_zdist_ind];                       \
            double _tsi = (1 - _w) * its[_rdist_ind][_zdist_ind] +                       \
                           _w * its[_rdist_ind_next][_zdist_ind];                       \
            tp[j] = (long int)(_tpi + 0.5);                                             \
            ts[j] = (long int)(_tsi + 0.5);                                             \
        }                                                                                \
        for (k = 0; k < nsamples; k++) {                                                 \
            double _stk0p = 0.0, _stk0s = 0.0;                                          \
            for (j = 0; j < nsta; j++) {                                                 \
                long int _ip = tp[j] + k;                                               \
                long int _is = ts[j] + k;                                               \
                if (_is < nsamples) {                                                    \
                    _stk0p += stackf_p[j][_ip];                                         \
                    _stk0s += stackf_s[j][_is];                                         \
                }                                                                        \
            }                                                                            \
            if (_stk0p + _stk0s > _stkmax) {                                             \
                _stkmax = _stk0p + _stk0s;                                               \
                _kmax = k;                                                               \
            }                                                                            \
        }                                                                                \
        double _coh = _stkmax / (double)nsta;                                            \
        corrmatrix[_idx] = _coh;                                                         \
        evaluated[_idx] = 1;                                                             \
        eval_count++;                                                                    \
        out_coh = _coh; out_k = _kmax;                                                   \
    } else {                                                                             \
        out_coh = corrmatrix[_idx]; out_k = 0;                                           \
    }                                                                                    \
} while (0)

    /* ---- Local multiscale hill-climb per seed ---- */
    int iteration = 0;
    const double tol = 1e-8; /* small tolerance to accept tiny improvements */
    int seed_index = 0;

    /* For each seed, perform a multi-scale local search (like original neighbor-guided),
       but limited by overall eval_count and max_evals. We iterate over seeds in priority order. */
    for (seed_index = 0; seed_index < nseeds && eval_count < max_evals; seed_index++) {
        /* start from this seed */
        int sx = seeds[seed_index].ix;
        int sy = seeds[seed_index].iy;
        int sz = seeds[seed_index].iz;
        double scoh = seeds[seed_index].coh;

        /* ensure the seed point itself is evaluated (should be from coarse but double-check) */
        {
            double seed_coh; int seed_k;
            EVALUATE_POINT_AND_STORE(sx, sy, sz, seed_coh, seed_k);
            seeds[seed_index].coh = seed_coh;
            if (seed_coh > best_coh) { best_coh = seed_coh; best_ix = sx; best_iy = sy; best_iz = sz; best_k = seed_k; }
        }

        /* multi-scale steps: start from half the coarse stride down to 1 */
        int step = max(1, stride / 2);
        while (step >= 1 && eval_count < max_evals) {
            int local_improved = 0;
            iteration++;
            /* try neighbors at this step until no improvement */
            int inner_changed = 1;
            while (inner_changed && eval_count < max_evals) {
                inner_changed = 0;
                for (int dix = -1; dix <= 1 && eval_count < max_evals; dix++) {
                    for (int diy = -1; diy <= 1 && eval_count < max_evals; diy++) {
                        for (int diz = -1; diz <= 1 && eval_count < max_evals; diz++) {
                            if (dix == 0 && diy == 0 && diz == 0) continue;
                            int nx_ix = sx + dix * step;
                            int nx_iy = sy + diy * step;
                            int nx_iz = sz + diz * step;
                            if (nx_ix < 0 || nx_ix >= nx || nx_iy < 0 || nx_iy >= ny || nx_iz < 0 || nx_iz >= nz) continue;

                            double coh_val; int kval;
                            EVALUATE_POINT_AND_STORE(nx_ix, nx_iy, nx_iz, coh_val, kval);

                            if (coh_val > seeds[seed_index].coh + tol) {
                                /* move seed */
                                seeds[seed_index].coh = coh_val;
                                seeds[seed_index].ix = nx_ix;
                                seeds[seed_index].iy = nx_iy;
                                seeds[seed_index].iz = nx_iz;
                                sx = nx_ix; sy = nx_iy; sz = nx_iz; /* update center */
                                local_improved = 1;
                                inner_changed = 1;
                                if (coh_val > best_coh) {
                                    best_coh = coh_val;
                                    best_ix = nx_ix; best_iy = nx_iy; best_iz = nx_iz; best_k = kval;
                                }
                            }
                        }
                    }
                }
            } /* while inner_changed */

            /* if we improved locally at this scale, try same step again to fully exploit; else reduce step */
            if (!local_improved) step = step / 2;
            /* loop continues until step becomes 0 -> stop when step < 1 */
        } /* while step >=1 */
    } /* for each seed */

#undef EVALUATE_POINT_AND_STORE

    /* free dynamic buffers */
    free(seeds);
    free(evaluated);

    /* ---- FINAL RESULTS ---- */
    iloc[0] = best_ix;
    iloc[1] = best_iy;
    iloc[2] = best_iz;
    *itime = best_k;

    printf("\n ------ Event located (Adaptive local Voronoi-like refinement) ------ \n");
    return 0;
}
