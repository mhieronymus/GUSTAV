# GUSTAV
GpU Spline Transfer And eValuation (GUSTAV) (GUSTAVO -- Optimized)

Loads a B-Spline from an hdf5-file, transfers it to a GPU and evaluates it on 
random points. Compares the execution time to that on a CPU.

## TODO
Test case to see if values are correct

## Prerequisites
 - [CUDA](https://developer.nvidia.com/cuda-zone)
 - [CFITSIO](https://heasarc.gsfc.nasa.gov/fitsio/)
 - [GNU Scientific Library](https://www.gnu.org/software/gsl/)