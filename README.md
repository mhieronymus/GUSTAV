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

## Results
Using 2^23 evaluations on a Xeon E5-2620 and a Tesla M40 we got 
 - loading:                 0.062455 s
 - vanilla on CPU:          14.4715 s
 - simple on CPU:           14.7673 s
 - my version on CPU:       14.834 s
 - host to device copy:     2.85369 s
 - GPU version:             0.63246 s

Given a more realistic amount of 2^27 evaluations, we got:
 - loading:                 0.0601048 s
 - vanilla on CPU:          229.986 s
 - host to device copy:     2.8747 s
 - GPU version:             9.93812 s
  
That is a speedup of 4.2 and 17.9 respectively.