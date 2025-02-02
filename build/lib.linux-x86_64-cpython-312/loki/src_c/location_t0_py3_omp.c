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

/* Prototypes */
int stacking(long int nxyz, long int nsta, long int nsamples, int *itp, int *its, double *stalta_p, double *stalta_s, double *corrmatrix, int nproc);

/* Python wrapper of the C function stacking */

static char module_docstring[] = "Module for computing the location";
static char stacking_docstring[] = "Location through waveform stacking";

static PyObject *py_stacking(PyObject *self, PyObject *args) {
    PyArrayObject *itp, *its, *stalta_p, *stalta_s, *corrmatrix;
    long int nxyz, nsamples, nsta, nproc;
    npy_intp dims[1];

    if (!PyArg_ParseTuple(args, "O!O!O!O!i",
                          &PyArray_Type, &itp,
                          &PyArray_Type, &its,
                          &PyArray_Type, &stalta_p,
                          &PyArray_Type, &stalta_s,
                          &nproc)) {
        PyErr_SetString(PyExc_RuntimeError, "Invalid arguments for the C function stacking");
        return NULL;
    }

    nsta = (long int)PyArray_DIM(stalta_p, 0);
    nsamples = (long int)PyArray_DIM(stalta_p, 1);
    nxyz = (long int)PyArray_DIM(itp, 0);  // nxyz from itp

    dims[0] = nxyz;

    corrmatrix = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);

    if (corrmatrix == NULL) {
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory for corrmatrix.");
        return NULL;
    }

    if (stacking(nxyz, nsta, nsamples,
                 (int *)PyArray_DATA(itp), (int *)PyArray_DATA(its),
                 (double *)PyArray_DATA(stalta_p), (double *)PyArray_DATA(stalta_s),
                 (double *)PyArray_DATA(corrmatrix), nproc) != 0) {
        PyErr_SetString(PyExc_RuntimeError, "Running stacking failed.");
        Py_DECREF(corrmatrix);
        return NULL;
    }

    Py_INCREF(corrmatrix);
    return PyArray_Return(corrmatrix);
}

static PyMethodDef module_methods[] = {
    {"stacking", py_stacking, METH_VARARGS, stacking_docstring},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef modlocation = {
    PyModuleDef_HEAD_INIT,
    "location_t0",  // Or "location" if you prefer
    module_docstring,
    -1,
    module_methods
};

PyMODINIT_FUNC PyInit_location_t0(void) {  // Or PyInit_location if you prefer
    PyObject *m;
    m = PyModule_Create(&modlocation);
    if (m == NULL)
        return NULL;
    import_array();
    return m;
}

int stacking(long int nxyz, long int nsta, long int nsamples, int *itp, int *its, double *stalta_p, double *stalta_s, double *corrmatrix, int nproc) {
    long int iter, i, j, k;
    int ip, is;
    double stk0p, stk0s, stkmax;
    iter = 0;

    omp_set_num_threads(nproc);

    printf("Location process complete at : %3d %%", 0);
    printf("nxyz: %ld, nsta: %ld, nsamples: %ld\n", nxyz, nsta, nsamples); // Print dimensions

#pragma omp parallel for shared(iter) private(ip, is, stkmax, stk0p, stk0s, k, j) schedule(dynamic) // Or schedule(static)
    for (i = 0; i < nxyz; i++) {
        printf("\b\b\b\b\b%3ld %%", (100 * iter++) / (nxyz - 2));
        stkmax = 0.;
        for (k = 0; k < nsamples; k++) {
            stk0p = 0.;
            stk0s = 0.;
            for (j = 0; j < nsta; j++) {
                ip = itp[i * nsta + j] + k;
                is = its[i * nsta + j] + k;

                printf("i: %ld, j: %ld, k: %ld, ip: %d, is: %d\n", i, j, k, ip, is); // Print indices

                if (ip < nsamples && is < nsamples) {
                    stk0p += stalta_p[j * nsamples + ip];
                    stk0s += stalta_s[j * nsamples + is];
                } else {
                    printf("WARNING: ip or is out of bounds! ip: %d, is: %d, nsamples: %ld\n", ip, is, nsamples); // Print warning
                }
            }
            stkmax = max(stk0p * stk0s, stkmax);
        }
        corrmatrix[i] = sqrt(stkmax) / ((float)nsta);
    }
    printf("\n ------ Event located ------ \n");
    return 0;
}