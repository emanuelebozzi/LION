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

// Python wrapper of the C function
static char module_docstring[] = "Module for computing traveltime processing";
static char tt_f2i_docstring[] = "Traveltime processing";

// Helper function to validate tp and ts arrays
void validate_array(double *tp, double *ts, long int nxz, long int nsta) {
    for (long int i = 0; i < nxz; i++) {
        for (long int j = 0; j < nsta; j++) {
            long int idx = i * nsta + j;

            if (isnan(tp[idx]) || tp[idx] < MIN_VALID_VALUE || tp[idx] > MAX_VALID_VALUE) {
                tp[idx] = 0.0;
            }
            if (isnan(ts[idx]) || ts[idx] < MIN_VALID_VALUE || ts[idx] > MAX_VALID_VALUE) {
                ts[idx] = 0.0;
            }
        }
    }
}

// Wrapper function to interface with Python
static PyObject *py_tt_f2i(PyObject *self, PyObject *args) {
    double dt;
    PyArrayObject *tp, *ts;
    long int nxz, nsta, nproc;
    npy_intp dims[2];

    // Parse the input arguments
    if (!PyArg_ParseTuple(args, "dOOi", &dt, &tp, &ts, &nproc)) {
        PyErr_SetString(PyExc_RuntimeError, "Invalid arguments for the C function tt_f2i");
        return NULL;
    }

    // Ensure inputs are NumPy arrays
    if (!PyArray_Check(tp) || !PyArray_ISCONTIGUOUS(tp) || PyArray_TYPE(tp) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_RuntimeError, "tp must be a contiguous array of type double");
        return NULL;
    }
    if (!PyArray_Check(ts) || !PyArray_ISCONTIGUOUS(ts) || PyArray_TYPE(ts) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_RuntimeError, "ts must be a contiguous array of type double");
        return NULL;
    }

    // Ensure both arrays are 2D
    if (PyArray_NDIM(tp) != 2) {
        PyErr_SetString(PyExc_RuntimeError, "tp must be a 2D array");
        return NULL;
    }
    if (PyArray_NDIM(ts) != 2) {
        PyErr_SetString(PyExc_RuntimeError, "ts must be a 2D array");
        return NULL;
    }

    // Ensure both arrays have the same shape
    if (PyArray_DIM(tp, 0) != PyArray_DIM(ts, 0) || PyArray_DIM(tp, 1) != PyArray_DIM(ts, 1)) {
        PyErr_SetString(PyExc_RuntimeError, "tp and ts must have the same shape");
        return NULL;
    }

    // Get dimensions
    nxz = dims[0] = (long int)PyArray_DIM(tp, 0);
    nsta = dims[1] = (long int)PyArray_DIM(tp, 1);

    // Create output arrays
    PyArrayObject *itp = (PyArrayObject*)PyArray_SimpleNew(2, dims, NPY_INT);
    PyArrayObject *its = (PyArrayObject*)PyArray_SimpleNew(2, dims, NPY_INT);

    if (!itp || !its) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to allocate output arrays");
        return NULL;
    }

    // Get raw data pointers
    double *tp_data = (double *)PyArray_DATA(tp);
    double *ts_data = (double *)PyArray_DATA(ts);
    int *itp_data = (int *)PyArray_DATA(itp);
    int *its_data = (int *)PyArray_DATA(its);

    // Validate arrays before processing
    validate_array(tp_data, ts_data, nxz, nsta);

    // Call the C processing function
    if (tt_f2i(dt, nxz, nsta, (double (*)[nsta])tp_data, (double (*)[nsta])ts_data, (int (*)[nsta])itp_data, (int (*)[nsta])its_data, nproc) != 0) {
        PyErr_SetString(PyExc_RuntimeError, "Error occurred during tt_f2i execution.");
        return NULL;
    }

    // Return the results as a tuple
    PyObject *result = Py_BuildValue("OO", itp, its);

    // Decrease reference counts to avoid memory leaks
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
    PyObject *m = PyModule_Create(&modtt_processing);
    if (m == NULL)
        return NULL;
    import_array();  // Initialize NumPy C API
    return m;
}

// C function to process the tp and ts arrays
int tt_f2i(double dt, long int nxz, long int nsta, double (*tp)[nsta], double (*ts)[nsta], int (*itp)[nsta], int (*its)[nsta], int nproc) {
    omp_set_num_threads(nproc);

    // Parallel loop for processing
    #pragma omp parallel for
    for (long int i = 0; i < nxz; i++) {
        for (long int j = 0; j < nsta; j++) {
            itp[i][j] = (int)lround(tp[i][j] / dt);
            its[i][j] = (int)lround(ts[i][j] / dt);
        }
    }

    return 0;
}
