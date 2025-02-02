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

    /* Verify dimensionality of input arrays */
    if (PyArray_NDIM(stalta_p) != 2) {
        PyErr_SetString(PyExc_RuntimeError, "stalta_p is not a 2D array");
        return NULL;
    }
    if (PyArray_NDIM(stalta_s) != 2) {
        PyErr_SetString(PyExc_RuntimeError, "stalta_s is not a 2D array");
        return NULL;
    }
    if (PyArray_NDIM(itp) != 2 || PyArray_DIM(itp, 1) != 1 || PyArray_NDIM(its) != 2 || PyArray_DIM(its, 1) != 1) {
        PyErr_SetString(PyExc_RuntimeError, "itp and its must be 2D arrays with a second dimension of size 1");
        return NULL;
    }

    /* Get dimensions from input arrays */
    nsta = (long int) PyArray_DIM(stalta_p, 0);  // Number of stations
    nsamples = (long int) PyArray_DIM(stalta_p, 1);  // Number of samples per station
    
    nxz = nx*nz; 
    nxyz = nx * nx * nz;  // Total number of grid points in 3D space

    /* Allocate memory for the correlation matrix */
    dims[0] = nxyz;  // Set dimension of correlation matrix array
    corrmatrix = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_DOUBLE);  // Create 1D NumPy array for the correlation matrix

    /* Convert itp and its to 1D arrays (flatten the 2D arrays) */
    int *itp_flat = (int*) PyArray_DATA(itp);  // Extract data as 1D array
    int *its_flat = (int*) PyArray_DATA(its);  // Extract data as 1D array

    /* Call the stacking function */
    if (stacking(nx, nz, nxz, nxyz, nsta, nsamples, itp_flat, its_flat, 
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

    /* Clean up references */
    Py_DECREF(iloctime);  // Decrement reference count of iloctime
    Py_DECREF(cohermat);  // Decrement reference count of cohermat
    Py_DECREF(corrmatrix);  // Decrement reference count of corrmatrix

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

/* Core stacking function with cylindrical symmetry mapping */
int stacking(long int nx, long int nz, long int nxz, long int nxyz, long int nsta, long int nsamples, 
             int itp[nxz], int its[nxz], double stalta_p[nsta][nsamples], 
             double stalta_s[nsta][nsamples], double corrmatrix[nxyz], 
             long int *iloc, long int *itime, int nproc) {

    long int i, j, k, kmax;
    int ip, is;
    double stk0p, stk0s, stkmax, corrmax;
    corrmax = -1.0;  // Initialize correlation maximum to a very low value

    omp_set_num_threads(nproc);  // Set the number of OpenMP threads
    printf("Location process complete at : %3d %%", 0);  // Print the starting point

    /* Parallel loop using OpenMP */
    #pragma omp parallel for shared(corrmax, iloc, itime) private(ip, is, stkmax, stk0p, stk0s, kmax, k, j)
    for (i = 0; i < nxyz; i++) {  // Loop over all grid points in the correlation matrix
        printf("\b\b\b\b\b%3ld %%", (100 * i) / (nxyz - 2));  // Display progress
        stkmax = -1.0;  // Reset the maximum stacking value
        kmax = 0;  // Reset the maximum sample index

        /* Loop over all samples */
        for (k = 0; k < nsamples; k++) {
            stk0p = 0.0;  // Reset the stacking value for P-wave
            stk0s = 0.0;  // Reset the stacking value for S-wave

            /* Loop over all stations */
            for (j = 0; j < nsta; j++) {
                // Map `i` to cylindrical symmetry indices for P and S-wave sample indices
                long int ix = i % nx;  // X index in cylindrical coordinate system
                long int iz = i / nx;  // Z index in cylindrical coordinate system

                ip = itp[ix + iz * nx];  // Get the P-wave sample index based on cylindrical grid
                is = its[ix + iz * nx];  // Get the S-wave sample index based on cylindrical grid

                if (is < nsamples) {  // If sample index is within valid range
                    stk0p += stalta_p[j][ip];  // Add the P-wave value
                    stk0s += stalta_s[j][is];  // Add the S-wave value
                }
            }

            /* Update the maximum stacking value if necessary */
            if (stk0p * stk0s > stkmax) {
                stkmax = stk0p * stk0s;  // Update max stacking value
                kmax = k;  // Update the corresponding sample index
            }
        }

        /* Store the correlation value for the current grid point in the correlation matrix */
        corrmatrix[i] = sqrt(stkmax) / ((float) nsta);

        /* Critical section to update the overall maximum correlation */
        #pragma omp critical
        if (corrmatrix[i] > corrmax) {
            corrmax = corrmatrix[i];  // Update the maximum correlation value
            *iloc = i;  // Update the location index
            *itime = kmax;  // Update the time index
        }
    }

    printf("\n ------ Event located ------ \n");  // Print completion message
    return 0;  // Return success
}
