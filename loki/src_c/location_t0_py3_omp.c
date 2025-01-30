#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#ifndef max
    #define max(a, b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
    #define min(a, b) (((a) < (b)) ? (a) : (b))
#endif

/* Function Prototype */
int stacking(long int nxyz, long int nsta, long int nsamples, void *itp, void *its, 
             double *stalta_p, double *stalta_s, 
             double *corrmatrix, long int *iloc, long int *itime, 
             int nproc, int is_1d_tp, int is_1d_ts);

/* Python Wrapper */
static PyObject *py_stacking(PyObject *self, PyObject *args) {
    PyArrayObject *itp, *its, *stalta_p, *stalta_s, *corrmatrix;
    long int nxyz, nsamples, nsta, nproc;
    long int iloc, itime;
    npy_intp dims[2];
    int is_1d_tp = 0, is_1d_ts = 0;

    /* Parse Python Arguments */
    if (!PyArg_ParseTuple(args, "OOOOi", &itp, &its, &stalta_p, &stalta_s, &nproc)) {
        PyErr_SetString(PyExc_RuntimeError, "Invalid arguments for the C function stacking");
        return NULL;
    }

    /* Ensure Contiguous Arrays */
    stalta_p = (PyArrayObject *)PyArray_FROMANY(stalta_p, NPY_DOUBLE, 2, 2, NPY_ARRAY_IN_ARRAY);
    stalta_s = (PyArrayObject *)PyArray_FROMANY(stalta_s, NPY_DOUBLE, 2, 2, NPY_ARRAY_IN_ARRAY);
    itp = (PyArrayObject *)PyArray_FROMANY(itp, NPY_INT, 1, 2, NPY_ARRAY_IN_ARRAY);
    its = (PyArrayObject *)PyArray_FROMANY(its, NPY_INT, 1, 2, NPY_ARRAY_IN_ARRAY);

    if (!stalta_p || !stalta_s || !itp || !its) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to ensure contiguous arrays");
        return NULL;
    }

    /* Determine Array Dimensionality */
    if (PyArray_NDIM(itp) == 1) is_1d_tp = 1;
    if (PyArray_NDIM(its) == 1) is_1d_ts = 1;

    /* Validate Array Dimensions */
    if (PyArray_NDIM(stalta_p) != 2 || PyArray_NDIM(stalta_s) != 2) {
        PyErr_SetString(PyExc_RuntimeError, "stalta_p or stalta_s is not a 2D array");
        return NULL;
    }

    /* Determine Array Sizes */
    nsta = (long int) PyArray_DIM(stalta_p, 0);
    nsamples = (long int) PyArray_DIM(stalta_p, 1);
    nxyz = (long int) PyArray_DIM(itp, 0);
    dims[0] = nxyz;
    dims[1] = nsamples;

    /* Create Output Array */
    corrmatrix = (PyArrayObject*) PyArray_SimpleNew(2, dims, NPY_DOUBLE);

    /* Call Stacking Function */
    if (stacking(nxyz, nsta, nsamples, PyArray_DATA(itp), PyArray_DATA(its), 
                 PyArray_DATA(stalta_p), PyArray_DATA(stalta_s), PyArray_DATA(corrmatrix), 
                 &iloc, &itime, nproc, is_1d_tp, is_1d_ts) != 0) {
        PyErr_SetString(PyExc_RuntimeError, "running stacking failed.");
        return NULL;
    }

    /* Prepare Return Values */
    PyObject *result = Py_BuildValue("(i,i)", iloc, itime);
    PyObject *corr_matrix = Py_BuildValue("O", corrmatrix);

    /* Cleanup */
    Py_DECREF(stalta_p);
    Py_DECREF(stalta_s);
    Py_DECREF(itp);
    Py_DECREF(its);
    Py_DECREF(corrmatrix);

    return Py_BuildValue("OO", result, corr_matrix);
}

/* Module Initialization */
static PyMethodDef module_methods[] = {
    {"stacking", py_stacking, METH_VARARGS, "location through waveform stacking"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef modlocation_t0 = {
    PyModuleDef_HEAD_INIT,
    "location_t0",
    "Module for computing of the location",
    -1,
    module_methods
};

PyMODINIT_FUNC PyInit_location_t0(void) {
    PyObject *m;
    m = PyModule_Create(&modlocation_t0);
    if (m == NULL) return NULL;
    import_array();
    return m;
}


/* Stacking Function */
int stacking(long int nxyz, long int nsta, long int nsamples, void *itp, void *its, 
             double *stalta_p, double *stalta_s, 
             double *corrmatrix, long int *iloc, long int *itime, 
             int nproc, int is_1d_tp, int is_1d_ts) {

    long int iter, i, j, k, kmax;
    int ip, is;
    double stk0p, stk0s, stkmax, corrmax;

    iter = 0;
    corrmax = -1.0;  // Ensure corrmax is initialized

    omp_set_num_threads(nproc);
    printf(" Location process complete at : %3d %%", 0);

    // Debug: Print input dimensions
    printf("\nDebug: itp dimensions: %ld x %ld\n", nxyz, nsamples);
    printf("Debug: its dimensions: %ld x %ld\n", nxyz, nsamples);
    printf("Debug: stalta_p dimensions: %ld x %ld\n", nsta, nsamples);
    printf("Debug: stalta_s dimensions: %ld x %ld\n", nsta, nsamples);

    #pragma omp parallel for shared(iter, corrmax, iloc, itime) private(ip, is, stkmax, stk0p, stk0s, kmax, k, j)
    for (i = 0; i < nxyz; i++) {
        if (i % 100000 == 0) {  // Print every 10 iterations
            printf("Debug: Processing i = %ld\n", i);
        }

        stkmax = 0.;  // Initialize stkmax
        for (k = 0; k < nsamples; k++) {
            stk0p = 0.;
            stk0s = 0.;
            for (j = 0; j < nsta; j++) {
                // Ensure that 'j' does not exceed nsta (bounds of stalta_p and stalta_s)
                if (j >= nsta) {
                    printf("Error: 'j' index out of bounds: j = %ld, max = %ld\n", j, nsta);
                    continue;
                }
                
                // Debug: print first few indices
                if (i == 0 && k == 0 && j < 5) {
                    printf("Debug: itp[%ld] = %d, its[%ld] = %d\n", i, ((int *)itp)[i], i, ((int *)its)[i]);
                    printf("Debug: stalta_p[%ld, %ld] = %f, stalta_s[%ld, %ld] = %f\n", j, k, stalta_p[j * nsamples + k], j, k, stalta_s[j * nsamples + k]);
                }

                ip = max(0, min(nsamples - 1, ((int *)itp)[i] + k));
                is = max(0, min(nsamples - 1, ((int *)its)[i] + k));

                // Ensure ip and is are within bounds before accessing arrays
                if (ip >= 0 && ip < nsamples && j < nsta) {
                    stk0p += stalta_p[j * nsamples + ip];  // Accessing as 1D array
                }
                if (is >= 0 && is < nsamples && j < nsta) {
                    stk0s += stalta_s[j * nsamples + is];  // Accessing as 1D array
                }
            }
            if (stk0p * stk0s > stkmax) {
                stkmax = stk0p * stk0s;
                kmax = k;
            }
            corrmatrix[i * nsamples + k] = sqrt(stkmax) / ((float)nsta);  // Indexing as 1D array
        }

        #pragma omp critical
        if (corrmatrix[i * nsamples + kmax] > corrmax) {  // Indexing as 1D array
            corrmax = corrmatrix[i * nsamples + kmax];  // Indexing as 1D array
            *iloc = i;
            *itime = kmax;
        }
    }

    printf("\n ------ Event located ------ \n");
    return 0;
}
