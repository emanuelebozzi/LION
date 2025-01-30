#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef max
    #define max(a, b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
    #define min(a, b) (((a) < (b)) ? (a) : (b))
#endif

/* Function Prototype */
int stacking(long int nxyz, long int nsta, long int nsamples, void *itp, void *its, 
             double stalta_p[nsta][nsamples], double stalta_s[nsta][nsamples], 
             double corrmatrix[nxyz][nsamples], long int *iloc, long int *itime, 
             int is_1d_tp, int is_1d_ts);

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
    if (!PyArray_Check(stalta_p) || !PyArray_ISCONTIGUOUS(stalta_p)) {
        PyErr_SetString(PyExc_RuntimeError, "stalta_p is not a contiguous array");
        return NULL;
    }
    if (!PyArray_Check(stalta_s) || !PyArray_ISCONTIGUOUS(stalta_s)) {
        PyErr_SetString(PyExc_RuntimeError, "stalta_s is not a contiguous array");
        return NULL;
    }
    if (!PyArray_Check(itp) || !PyArray_ISCONTIGUOUS(itp)) {
        PyErr_SetString(PyExc_RuntimeError, "tp is not a contiguous array");
        return NULL;
    }
    if (!PyArray_Check(its) || !PyArray_ISCONTIGUOUS(its)) {
        PyErr_SetString(PyExc_RuntimeError, "ts is not a contiguous array");
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
    if (!is_1d_tp && PyArray_NDIM(itp) != 2) {
        PyErr_SetString(PyExc_RuntimeError, "tp is not a 2D array");
        return NULL;
    }
    if (!is_1d_ts && PyArray_NDIM(its) != 2) {
        PyErr_SetString(PyExc_RuntimeError, "ts is not a 2D array");
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
                 &iloc, &itime, is_1d_tp, is_1d_ts) != 0) {
        PyErr_SetString(PyExc_RuntimeError, "running stacking failed.");
        return NULL;
    }

    /* Return Results */
    PyObject *result = Py_BuildValue("(i,i)", iloc, itime);
    PyObject *corr_matrix = Py_BuildValue("O", corrmatrix);
    Py_DECREF(corrmatrix);

    PyObject *final_result = Py_BuildValue("OO", result, corr_matrix);
    Py_DECREF(result);
    Py_DECREF(corr_matrix);

    return final_result;
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
             double stalta_p[nsta][nsamples], double stalta_s[nsta][nsamples], 
             double corrmatrix[nxyz][nsamples], long int *iloc, long int *itime, 
             int is_1d_tp, int is_1d_ts) {

    long int iter, i, j, k, kmax;
    int ip, is;
    double stk0p, stk0s, stkmax, corrmax;

    iter = 0;
    corrmax = -1.0;  // Ensure corrmax is initialized

    // Check for valid dimensions
    if (nxyz <= 0 || nsta <= 0 || nsamples <= 0) {
        return -1;  // Invalid array dimensions
    }

    // Sequential execution, no OpenMP
    for (i = 0; i < nxyz; i++) {
        stkmax = -1e6;  // Initialize to a very low value to find maximum in the next loop

        for (k = 0; k < nsamples; k++) {
            stk0p = 0.;
            stk0s = 0.;

            for (j = 0; j < nsta; j++) {
                ip = max(0, min(nsamples - 1, ((int *)itp)[i] + k));
                is = max(0, min(nsamples - 1, ((int *)its)[i] + k));

                // Ensure ip and is are within the bounds
                if (ip < 0 || ip >= nsamples || is < 0 || is >= nsamples) {
                    printf("Error: Invalid access: ip = %d, is = %d, nsamples = %ld\n", ip, is, nsamples);
                    return -2;  // Return an error for invalid access
                }

                stk0p += stalta_p[j][ip];
                stk0s += stalta_s[j][is];
            }

            // Update the correlation matrix for each k (sample)
            corrmatrix[i][k] = sqrt(stk0p + stk0s) / ((float)nsta);

            // Update stkmax to track the maximum correlation for this i
            if (corrmatrix[i][k] > stkmax) {
                stkmax = corrmatrix[i][k];
                kmax = k;  // Track the index of the max correlation
            }
        }

        // Check if this is the new global max correlation
        if (stkmax > corrmax) {
            corrmax = stkmax;
            *iloc = i;
            *itime = kmax;  // Record the index with max correlation
        }
    }

    return 0;
}
