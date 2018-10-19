# GUSTAV
GpU Spline Transfer And eValuation (GUSTAV) (GUSTAVO -- Optimized)

Loads a B-Spline from an hdf5-file, transfers it to a GPU and evaluates it on 
random points. Compares the execution time to that on a CPU.

## TODO
 - Simple minimzer (or scan)
 - Track length optimization
 - Newton Method (Cascade energy optimization)
 - Fit statistics (llh evaluation)

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
  
That is a speedup of 4.2 and 17.9 respectively. The code for evaluating
Splines takes about 30.6 % of total runtime in the reconstruction 
(DRAGON sample and 100 frames).

## Notes
This is a very naive version which does not do anything in a good way.

Debug info tells us, we use 160 registers per thread.
There are 256 KB of registers per SM in a GTX 1080
That is 256*1024/(8*160) = 204,8 many threads that can run on a SM.
Hence we use 192 threads and 40 blocks -> 7680 threads running in parallel.

Quite some registers are spilled.