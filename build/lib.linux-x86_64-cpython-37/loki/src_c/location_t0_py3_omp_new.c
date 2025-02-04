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
// i've added nx, nz, nxz for the cylindrical symmetry of the problem //
// itp and its (traveltime lookup tables are 1D, since nsta = 1 and cyclindric symmetry)

int stacking(long int nx, long int nz, long int nxz, long int nxyz, long int nsta, long int nsamples, 
             int itp[nxz], int its[nxz], double stalta_p[nsta][nsamples], 
             double stalta_s[nsta][nsamples], double corrmatrix[nxyz], 
             long int *iloc, long int *itime, int nproc);

/* Python module documentation */
static char module_docstring[] = "Module for computing the location";  // Description for the module
static char stacking_docstring[] = "Location through waveform stacking";  // Description for the stacking function

/* Python wrapper function for the core stacking function */
// Possible error messages 

static PyObject *py_stacking(PyObject *self, PyObject *args) {
    PyArrayObject *itp, *its, *stalta_p, *stalta_s, *corrmatrix;
    long int nx, nz, nxz, nxyz, nsamples, nsta, nproc;
    long int iloc, itime;
    npy_intp dims[1];  // Array to store dimensions of coherence matrix 

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
    nsta = (long int) PyArray_DIM(stalta_p, 0);  // Number of stations from observed data, NOT traveltime table
    nsamples = (long int) PyArray_DIM(stalta_p, 1);  // Number of samples per station
    
    nxz = nx * nz; //number of points for the 2D traveltime lookup table
    nxyz = nx * nx * nz;  // Total number of grid points in 3D space (even location)

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
    PyObject *iloctime = Py_BuildValue("(i,i)", iloc, itime);  // Build tuple with location and time
    PyObject *cohermat = Py_BuildValue("O", corrmatrix); Py_DECREF(corrmatrix);  // Build Python object for the correlation matrix
    PyObject *locres = Py_BuildValue("(OO)", iloctime, cohermat); Py_DECREF(iloctime); Py_DECREF(cohermat);  // Bundle location and matrix into a tuple

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


// Compute distances, depth on 3D and 2D for the location grid and traveltime lookup tables 

// key question is : are the traveltimes from NonLinLoc ordered the same way???

// Struct for index, distance, and depth (both 3d and 3d)
typedef struct {
    int index;
    double distance;
    double depth;
} GridInfo;

// Function that calculates distance from origin (0,0,0)
double calculate_distance_xyz(int x, int y, int z) {
    return sqrt(x * x + y * y + z * z);
}

// Function that processes the 3D location grid assigning ordered indexes, distances, depth from reference 0,0,0

void process_grid_xyz(int nx, int ny, int nz, GridInfo *grid_matrix) {
    int idx = 0; // For tracking position in the array

    for (int z = 0; z < nz; z++) {
        for (int y = 0; y < ny; y++) {
            for (int x = 0; x < nx; x++) {
                grid_matrix[idx].index = x + y * nx + z * nx * ny;
                grid_matrix[idx].distance = calculate_distance_xyz(x, y, z);
                grid_matrix[idx].depth = z;
                idx++;
            }
        }
    }
}


// the same but 2D (understand if this can be done more efficiently)
double calculate_distance_xz(int x, int z) {
    return sqrt(x * x + z * z);
}

void process_grid_xz(int nx, int nz, GridInfo *grid_matrix) {
    int idx = 0;  // For tracking position in the array

    for (int z = 0; y < nz; z++) {
        for (int x = 0; x < nx; x++) {
            grid_matrix[idx].index = x + z * nx;                  // 2D indexing
            grid_matrix[idx].distance = calculate_distance_xy(x, z);
            grid_matrix[idx].depth = z;                           // Depth as y-coordinate, or remove if not needed
            idx++;
        }
    }
}


// Allocate memory for the grid matrix
GridInfo *grid_loc3d = (GridInfo *)malloc(nxyz * sizeof(GridInfo));
GridInfo *grid_tt2d = (GridInfo *)malloc(nxz * sizeof(GridInfo));

// Process the grid
process_grid_xyz(nx, ny, nz, grid_loc3d);
process_grid_xz(nx, nz, grid_tt2d);


//assign values
index3d = grid_loc3d.index; 
distance3d = grid_loc3d.distance; 
depth3d = grid_loc3d.depth; 

index2d = grid_tt2d.index; 
distance2d = grid_tt2d.distance; 
depth2d = grid_tt2d.depth; 


/* Stacking Function */
int stacking(long int index3d, long int index2d, long int nx, long int nz, long int nxz, long int nxyz, long int nsta, long int nsamples, 
             int itp[nxz], int its[nxz], double distance2d, double distance3d, double depth2d, double depth3d, double stalta_p[nsta][nsamples], double stalta_s[nsta][nsamples], 
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
   

    // OUTHER LOOK OVER 3D LOCATION GRID POINTS 

    for (w = 0; w < nxyz; w++) {  

        printf("\b\b\b\b\b%3ld %%", (100 * iter++) / (nxyz - 2));

        //for each point read from the created grid the current distance and depth in 3D

        curr_distance =  distance3d[w] 
        curr_depth = depth3d[w] 
        
        //find the index of the closest distance-depth point in 2D (cartesian distance from the current 3D distances) 

        //assuming NonLinLoc traveltimes are ordered the same way !!! 

        int closest_tt = 0; for (int i = 1; i < size; i++) if (sqrt(pow(distance2d[i] - curr_distance, 2) + pow(depth2d[i] - curr_depth, 2)) < sqrt(pow(distance[closest_index] - curr_distance, 2) + pow(depth[closest_index] - curr_depth, 2))) closest_index = i;
        

        //initialize the stack energy 

        stkmax = -1.0;  // Reset for each new grid point
        kmax = 0;

        // INNER LOOP FOR THE TIME SAMPLES 
        for (k = 0; k < nsamples; k++) {
            stk0p = 0.0;
            stk0s = 0.0;

            // INNER LOOP FOR ALL THE STATIONS 
            // and calculate stack contributions
            for (j = 0; j < nsta; j++) {
                ip = itp[closest_tt] + k;  // access the traveltime (P)
                is = its[closest_tt] + k;  // access the traveltime (s)
                
                if (ip < nsamples && is < nsamples) {
                    stk0p += stalta_p[j][ip];  // STACK P
                    stk0s += stalta_s[j][is];  // STACK S
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


        // FOR EACH 3D LOCATION GRID POINT STORE TO COHERENCE AND NORMALIZE FOR THE NUMBER OF STATIONS 

        corrmatrix[w] = sqrt(stkmax) / ((float) nsta);

        #pragma omp critical
        if (corrmatrix[w] > corrmax) {
            corrmax = corrmatrix[i];
            *iloc = w;
            *itime = kmax;
        }    
    }

    }

    printf("\n ------ Event located ------ \n");  // Print completion message

    return 0; // Return success
}
