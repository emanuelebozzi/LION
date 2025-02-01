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

// Define constants for the maximum valid values for tp and ts
#define MAX_VALID_VALUE 100000.0
#define MIN_VALID_VALUE 0.0  // You can change this as needed

// Prototypes
int tt_f2i(double dt, long int nxz, long int nsta, double (*tp)[nsta], double (*ts)[nsta], int (*itp)[nsta], int (*its)[nsta], int nproc);

// Python wrapper of the C function stacking
static char module_docstring[] = "Module for computing traveltime processing";
static char tt_f2i_docstring[] = "Traveltime processing";

// Helper function to validate tp and ts arrays
void validate_array(double *tp, double *ts, long int nxz, long int nsta) {
    for (long int i = 0; i < nxz; i++) {
        for (long int j = 0; j < nsta; j++) {
            // Check and validate tp
            if (isnan(tp[i * nsta + j]) || tp[i * nsta + j] < MIN_VALID_VALUE || tp[i * nsta + j] > MAX_VALID_VALUE) {
                tp[i * nsta + j] = 0.0;  // Default to 0 if invalid
            }
            // Check and validate ts
            if (isnan(ts[i * nsta + j]) || ts[i * nsta + j] < MIN_VALID_VALUE || ts[i * nsta + j] > MAX_VALID_VALUE) {
                ts[i * nsta + j] = 0.0;  // Default to 0 if invalid
            }
        }
    }
}

// Wrapper function to interface with Python
static PyObject *py_tt_f2i(PyObject *self, PyObject *args) {
    double dt;
    PyArrayObject *tp, *ts, *itp, *its;
    long int nxz, nsta, nproc;
    npy_intp dims[2];
    
    // Parse the input arguments
    if (!PyArg_ParseTuple(args, "dOOi", &dt, &tp, &ts, &nproc)) {
        PyErr_SetString(PyExc_RuntimeError, "Invalid arguments for the C function ttprocessing");
        return NULL;
    }

    // Checking that tp and ts are contiguous arrays
    if (!PyArray_Check(tp) || !PyArray_ISCONTIGUOUS(tp)) {
        PyErr_SetString(PyExc_RuntimeError, "tp is not a contiguous array");
        return NULL;
    }
    if (!PyArray_Check(ts) || !PyArray_ISCONTIGUOUS(ts)) {
        PyErr_SetString(PyExc_RuntimeError, "ts is not a contiguous array");
        return NULL;
    }

    // Checking the dimensions of tp and ts arrays
    if (PyArray_NDIM(tp) != 2) {
        PyErr_SetString(PyExc_RuntimeError, "tp is not a 2D array");
        return NULL;
    }
    if (PyArray_NDIM(ts) != 2) {
        PyErr_SetString(PyExc_RuntimeError, "ts is not a 2D array");
        return NULL;
    }

    // Get the dimensions of the tp array
    nxz = dims[0] = (long int)PyArray_DIM(tp, 0);
    nsta = dims[1] = (long int)PyArray_DIM(tp, 1);

    // Create empty arrays for the results (itp and its)
    itp = (PyArrayObject*)PyArray_SimpleNew(2, dims, NPY_INT);
    its = (PyArrayObject*)PyArray_SimpleNew(2, dims, NPY_INT);

    // Convert tp and ts arrays to C-style arrays for processing
    double *tp_data = (double *)PyArray_DATA(tp);
    double *ts_data = (double *)PyArray_DATA(ts);
    
    // Validate the tp and ts arrays before processing
    validate_array(tp_data, ts_data, nxz, nsta);

    // Call the C function tt_f2i
    if (tt_f2i(dt, nxz, nsta, (double (*)[nsta])tp_data, (double (*)[nsta])ts_data, (int (*)[nsta])PyArray_DATA(itp), (int (*)[nsta])PyArray_DATA(its), nproc) != 0) {
        PyErr_SetString(PyExc_RuntimeError, "Running tt_f2i failed.");
        return NULL;
    }

    // Return the results
    PyObject *result = Py_BuildValue("OO", itp, its);
    Py_DECREF(itp);
    Py_DECREF(its);
    return result;
}

// Module methods
static PyMethodDef module_methods[] = {
    {"tt_f2i", py_tt_f2i, METH_VARARGS, tt_f2i_docstring},
    {NULL, NULL, 0, NULL}
};

// Module definition
static struct PyModuleDef modtt_processing = {
    PyModuleDef_HEAD_INIT,
    "tt_processing",
    module_docstring,
    -1,
    module_methods
};

// Module initialization
PyMODINIT_FUNC PyInit_tt_processing(void) {
    PyObject *m;
    m = PyModule_Create(&modtt_processing);
    if (m == NULL)
        return NULL;
    import_array();
    return m;
};

// C function to process the tp and ts arrays
int tt_f2i(double dt, long int nxz, long int nsta, double (*tp)[nsta], double (*ts)[nsta], int (*itp)[nsta], int (*its)[nsta], int nproc) {
    long int i;

    // Set the number of OpenMP threads
    omp_set_num_threads(nproc);

    // Parallel processing loop for calculating itp and its
    #pragma omp parallel for
    for (i = 0; i < nxz; i++) {
        itp[i][0] = (int)lround(tp[i][0] / dt);
        its[i][0] = (int)lround(ts[i][0] / dt);
    }

    return 0;
}
