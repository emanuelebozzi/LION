#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION  // Ensure compatibility with NumPy API version 1.7
#include <Python.h>  // Include Python API
#include <numpy/arrayobject.h>  // Include NumPy array object handling
#include <math.h>  // Include math functions (e.g., sqrt, fmax)
#include <stdio.h>  // Include standard I/O for printing
#include <stdlib.h>  // Standard library for dynamic memory allocation
#include <omp.h>  // OpenMP for parallel processing

// Define max and min macros to avoid conflicts with system-defined ones
#ifndef max
    #define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif
#ifndef min
    #define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

/* Prototypes of functions */

/* The 'stacking' function performs the location calculation through waveform stacking.
   The parameters include grid dimensions, station information, and the number of processors for parallel computation. */
int stacking(long int nxz, long int nxyz, long int nsta, long int nsamples, 
             int *itp, int *its, double *stalta_p, double *stalta_s, 
             double *corrmatrix, int nproc);

/* Python wrapper of the C function 'stacking' */
static char module_docstring[]="Module for computing the location";  // Description of the module
static char stacking_docstring[]="location through waveform stacking";  // Description of the 'stacking' function


/* Wrapper for Python interface to the C function 'stacking' */
static PyObject *py_stacking(PyObject *self, PyObject *args){
    PyArrayObject *itp, *its, *stalta_p, *stalta_s, *corrmatrix;  // Declare NumPy array objects
    long int nxz, nxyz, nsamples, nsta, nproc;  // Declare grid dimensions and other parameters
    npy_intp dims[1];  // Declare dimensions for the output array (correlation matrix)

    // Parse Python arguments: arrays itp, its, stalta_p, stalta_s, and integer nproc
    if(!PyArg_ParseTuple(args, "OOOOi", &itp, &its, &stalta_p, &stalta_s, &nproc)){
        PyErr_SetString(PyExc_RuntimeError, "Invalid arguments for the C function stacking");
        return NULL; 
    }

    // Check if the arrays are contiguous in memory for efficient processing
    if(!PyArray_Check(stalta_p) || !PyArray_ISCONTIGUOUS(stalta_p)){
        PyErr_SetString(PyExc_RuntimeError, "stalta_p is not a contiguous array");
        return NULL; 
    }

    if(!PyArray_Check(stalta_s) || !PyArray_ISCONTIGUOUS(stalta_s)){
        PyErr_SetString(PyExc_RuntimeError, "stalta_s is not a contiguous array");
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

    // Ensure that all arrays are 2D
    if((PyArray_NDIM(stalta_p) != 2)){
        PyErr_SetString(PyExc_RuntimeError, "stalta_p is not a 2D array");
        return NULL; 
    }

    if((PyArray_NDIM(stalta_s) != 2)){
        PyErr_SetString(PyExc_RuntimeError, "stalta_s is not a 2D array");
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

    // Get the dimensions of the station arrays
    nsta = (long int) PyArray_DIM(stalta_p, 0);  // Number of stations
    nsamples = (long int) PyArray_DIM(stalta_p, 1);  // Number of samples
    nxz = (long int) PyArray_DIM(itp, 0);  // Number of travel times in the 2D lookup table in XZ plane

    // Calculate the number of grid points in the 3D space (assuming square grid in Z dimension)
    nxyz = nxz * (long int)sqrt(nxz);  // Grid size (assuming square Z dimension)

    dims[0] = nxyz;  // Set the dimensions of the output correlation matrix
    corrmatrix = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_DOUBLE);  // Create a new empty correlation matrix array

    // Call the stacking function to compute the waveform stacking
    if (0 != stacking(nxz, nxyz, nsta, nsamples, 
                      (int *)PyArray_DATA(itp), 
                      (int *)PyArray_DATA(its), 
                      (double *)PyArray_DATA(stalta_p), 
                      (double *)PyArray_DATA(stalta_s), 
                      (double *)PyArray_DATA(corrmatrix), 
                      nproc)) {
        PyErr_SetString(PyExc_RuntimeError, "running stacking failed.");
        return NULL;
    }

    // Return the resulting correlation matrix as a NumPy array
    Py_INCREF(corrmatrix);
    return PyArray_Return(corrmatrix);
}


/* Module specifications and initialization */
static PyMethodDef module_methods[] = {
    /* {method_name, Cfunction, argument_types, docstring} */
    {"stacking", py_stacking, METH_VARARGS, stacking_docstring},  // Register the stacking function
    {NULL, NULL, 0, NULL}  // End of method list
};

// Define the module
static struct PyModuleDef modlocation_t0 = {
    PyModuleDef_HEAD_INIT,  // Initialize module
    "location_t0",  // Module name (this should match the module you are trying to import in Python)
    module_docstring,  // Module docstring
    -1,  // Size of the module (keeping it negative to indicate no specific size)
    module_methods  // Register the methods of the module
};

// Module initialization function
PyMODINIT_FUNC PyInit_location_t0(void){
    PyObject *m;  // Declare the module object
    m = PyModule_Create(&modlocation_t0);  // Create the module
    if (m == NULL)  // Check for creation errors
        return NULL;
    import_array();  // Initialize NumPy array API
    return m;  // Return the module object
}


// Stacking function for computing location through waveform stacking
int stacking(long int nxz, long int nxyz, long int nsta, long int nsamples, 
             int *itp, int *its, double *stalta_p, double *stalta_s, 
             double *corrmatrix, int nproc) {

    long int iter, i, j, k, idx;  // Loop counters
    int ip, is;  // Indices for P and S wave traveltimes
    double stk0p, stk0s, stkmax;  // Variables to hold temporary stack values
    double stack2D[nxz];  // Array to hold the 2D stack for cylindrical symmetry

    iter = 0;  // Initialize iteration counter
    omp_set_num_threads(nproc);  // Set number of threads for parallel processing

    printf(" Location process complete at : %3d %%", 0);  // Print progress (initially 0%)

    // Step 1: Stack energy for each station in the XZ plane (cylindrical symmetry)
    #pragma omp parallel for private(ip, is, stk0p, stk0s, k, j)  // Parallel loop over all grid points
    for (i = 0; i < nxz; i++) {
        stack2D[i] = 0.0;  // Initialize the stack for the current grid point in XZ plane
        for (k = 0; k < nsamples; k++) {  // Loop over all samples
            stk0p = 0.0;  // Initialize P-wave stack
            stk0s = 0.0;  // Initialize S-wave stack
            for (j = 0; j < nsta; j++) {  // Loop over all stations
                ip = itp[i] + k;  // P-wave traveltime index for the current grid point and sample
                is = its[i] + k;  // S-wave traveltime index for the current grid point and sample
                if (is < nsamples) {  // Ensure S-wave index is within bounds
                    stk0p += stalta_p[j * nsamples + ip];  // Sum STA/LTA values for P-wave
                    stk0s += stalta_s[j * nsamples + is];  // Sum STA/LTA values for S-wave
                }
            }
            stack2D[i] = fmax(stack2D[i], stk0p * stk0s);  // Keep the maximum of P*S for each sample
        }
    }

    // Step 2: Combine 2D stack into 3D grid with interpolation
    #pragma omp parallel for private(stkmax, idx)  // Parallel loop to fill the 3D correlation matrix
    for (i = 0; i < nxyz; i++) {
        printf("\b\b\b\b\b%3ld %%", (100 * iter++) / (nxyz - 2));  // Print progress

        idx = i % nxz;  // Simplified mapping of 3D grid index to 2D cylindrical coordinate index
        stkmax = sqrt(stack2D[idx]);  // Interpolate using square root of the stacked energy

        corrmatrix[i] = stkmax / ((double)nsta);  // Normalize by number of stations
    }

    printf("\n ------ Event located ------ \n");  // Print message indicating location is complete
    return 0;  // Return 0 to indicate success
}
