#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION  // Prevent using deprecated NumPy API
#include <Python.h>  // Include Python header for Python C-API integration
#include <numpy/arrayobject.h>  // Include NumPy header for array handling
#include <math.h>  // For mathematical functions (e.g., sqrt)
#include <stdio.h>  // For standard input/output functions
#include <stdlib.h>  // For standard library functions (e.g., malloc, free)
#include <omp.h>  // For OpenMP (parallel processing)

/* Ensure max and min macros are defined */
#ifndef max
    #define max( a, b ) ( ((a) > (b)) ? (a) : (b) )  // Define max macro if not already defined
#endif
#ifndef min
    #define min( a, b ) ( ((a) < (b)) ? (a) : (b) )  // Define min macro if not already defined
#endif

/* Function prototype for the core stacking function */
int stacking(long int nx, long int nz, long int nxz, long int nxyz, long int nsta, long int nsamples, 
             int itp[nxz], int its[nxz], double stalta_p[nsta][nsamples], 
             double stalta_s[nsta][nsamples], double corrmatrix[nxyz], 
             long int *iloc, long int *itime, int nproc);

/* Python module documentation */
static char module_docstring[] = "Module for computing the location";  // Description for the module
static char stacking_docstring[] = "Location through waveform stacking";  // Description for the stacking function

/* Python wrapper function for the core stacking function */
static PyObject *py_stacking(PyObject *self, PyObject *args) {
    PyArrayObject *itp, *its, *stalta_p, *stalta_s, *corrmatrix;
    long int nx, nz, nxz, nxyz, nsamples, nsta, nproc;
    long int iloc, itime;
    npy_intp dims[1];  // Array to store dimensions of the correlation matrix

    if (!PyArg_ParseTuple(args, "iiOOOOi", &nx, &nz, &itp, &its, &stalta_p, &stalta_s, &nproc)) {
        PyErr_SetString(PyExc_RuntimeError, "Invalid arguments for the C function stacking");
        return NULL;
    }

    /* Ensure arrays are contiguous in memory */
    if (!PyArray_Check(stalta_p) || !PyArray_ISCONTIGUOUS(stalta_p)) {
        PyErr_SetString(PyExc_RuntimeError, "stalta_p is not a contiguous array");
        return NULL;
    }
    if (!PyArray_Check(stalta_s) || !PyArray_ISCONTIGUOUS(stalta_s)) {
        PyErr_SetString(PyExc_RuntimeError, "stalta_s is not a contiguous array");
        return NULL;
    }
    if (!PyArray_Check(itp) || !PyArray_ISCONTIGUOUS(itp)) {
        PyErr_SetString(PyExc_RuntimeError, "itp is not a contiguous array");
        return NULL;
    }
    if (!PyArray_Check(its) || !PyArray_ISCONTIGUOUS(its)) {
        PyErr_SetString(PyExc_RuntimeError, "its is not a contiguous array");
        return NULL;
    }

    /* Get dimensions from input arrays */
    nsta = (long int) PyArray_DIM(stalta_p, 0);  // Number of stations
    nsamples = (long int) PyArray_DIM(stalta_p, 1);  // Number of samples per station
    
    nxz = nx * nz; 
    nxyz = nx * nx * nz;  // Total number of grid points in 3D space

    /* Allocate memory for the correlation matrix */
    dims[0] = nxyz;  // Set dimension of correlation matrix array
    corrmatrix = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_DOUBLE);  // Create 1D NumPy array for the correlation matrix

    /* Call the stacking function */
    if (stacking(nx, nz, nxz, nxyz, nsta, nsamples, 
                (int*) PyArray_DATA(itp), (int*) PyArray_DATA(its), 
                (double (*)[nsamples]) PyArray_DATA(stalta_p), 
                (double (*)[nsamples]) PyArray_DATA(stalta_s), 
                (double*) PyArray_DATA(corrmatrix), &iloc, &itime, nproc) != 0) {
        PyErr_SetString(PyExc_RuntimeError, "Running stacking failed.");
        return NULL;  // Return error if stacking function fails
    }

    /* Prepare output Python objects */
    PyObject *iloctime = Py_BuildValue("(l,l)", iloc, itime);  // Build tuple with location and time
    PyObject *cohermat = Py_BuildValue("O", corrmatrix);  // Build Python object for the correlation matrix
    PyObject *locres = Py_BuildValue("(OO)", iloctime, cohermat);  // Bundle location and matrix into a tuple

    return locres;  // Return the results as a Python tuple
}

/* Module specifications and initialization */
static PyMethodDef module_methods[] = {
    {"stacking", py_stacking, METH_VARARGS, stacking_docstring},  // Register the stacking function
    {NULL, NULL, 0, NULL}  // End of methods list
};

/* Module definition */
static struct PyModuleDef modlocation_t0 = {
    PyModuleDef_HEAD_INIT,
    "location_t0",  // Name of the module
    module_docstring,  // Module description
    -1,  // Size of per-interpreter state of the module, -1 means the module keeps state in global variables
    module_methods  // List of module methods
};

/* Module initialization function */
PyMODINIT_FUNC PyInit_location_t0(void) {
    PyObject *m;
    m = PyModule_Create(&modlocation_t0);  // Create the module
    if (m == NULL)
        return NULL;  // Return NULL if module creation fails
    import_array();  // Initialize NumPy C-API
    return m;  // Return the initialized module
}

/* Stacking Function */
int stacking(long int nx, long int nz, long int nxz, long int nxyz, long int nsta, long int nsamples, 
             int itp[nxz], int its[nxz], double stalta_p[nsta][nsamples], double stalta_s[nsta][nsamples], 
             double corrmatrix[nxyz], long int *iloc, long int *itime, int nproc) {
    
    long int iter = 0, i, j, k, w, kmax;
    int ip, is;
    double stk0p, stk0s, stkmax;
    
    // Set OpenMP threads
    omp_set_num_threads(nproc);
    
    // Initialize tracking variables
    double corrmax = -1.0;
    
    printf("Location process complete at: %3d %%", 0);
    
    #pragma omp parallel for shared(iter, corrmax, iloc, itime) private(ip, is, stkmax, stk0p, stk0s, kmax, k, j)
   
    for (w = 0; w < nxyz; w++) {  // Loop over the 3D location grid points

        // Loop over the 2D travel time lookup tables
        for (i = 0; i < nxz; i++) {
            // Printing progress
            printf("\b\b\b\b\b%3ld %%", (100 * iter++) / (nxyz - 2));

            stkmax = -1.0;  // Reset for each new grid point
            kmax = 0;

            // Iterate over the samples
            for (k = 0; k < nsamples; k++) {
                stk0p = 0.0;
                stk0s = 0.0;

                // Loop over stations and calculate stack contributions
                for (j = 0; j < nsta; j++) {
                    ip = itp[i] + k;  // Just access directly as 1D arrays
                    is = its[i] + k;  // Just access directly as 1D arrays
                   
                    if (ip < nsamples && is < nsamples) {
                        stk0p += stalta_p[j][ip];  // Add positive travel time
                        stk0s += stalta_s[j][is];  // Add negative travel time
                    } else {
                        stk0p = 0.0;
                        stk0s = 0.0;
                    }
                }

                // If the product of stacks is greater than the current max, update it
                if (stk0p * stk0s > stkmax) {
                    stkmax = stk0p * stk0s;
                    kmax = k;  // Store the index corresponding to the max stack value
                }
            }

            // Store the result in the correct position in the correlation matrix
            corrmatrix[i] = sqrt(stkmax) / ((float) nsta);

            #pragma omp critical
            if (corrmatrix[i] > corrmax) {
                corrmax = corrmatrix[i];
                *iloc = i;
                *itime = kmax;
            }
        }
    }

    printf("\n ------ Event located ------ \n");  // Print completion message

    return 0; // Return success
}
