#include "helper/types.hpp"
#include "helper/load.hpp"
#include "helper/timing.hpp"
#include "cpu_bsplv.hpp"
#include "gpu_bsplv.h"


int main(int argc, char* argv[]) {
    
    value_t * spline_cffs = nullptr; // Coefficients
    value_t * spline_knts = nullptr; // knots
    index_t * spline_dgrs = nullptr; // degrees 
    index_t * n_knts      = nullptr; // Number of knots per dimension
    index_t * n_splines   = nullptr; // Number of splines per dimension
    index_t n_cffs, n_dgrs, n_dims;
    index_t n_evals = 1 << 10;       // Number of times to evaluate the splines

    TIMERSTART(load_from_file)
    // Load a spline file from a file given by the command line
    load_splines(spline_cffs, n_cffs, 
                 spline_knts, n_knts, 
                 spline_dgrs, n_dgrs, 
                 n_dims, *argv);
    TIMERSTOP(load_from_file)

    TIMERSTART(cpu)
    cpu_eval_splines(spline_cffs, n_cffs, 
                     spline_knts, n_knts, 
                     spline_dgrs, n_dgrs, 
                     n_evals, n_dims, 
                     n_splines);
    TIMERSTOP(cpu)

    TIMERSTART(gpu_overall)
    // gpu_eval_splines(spline_cffs, n_cffs, 
    //                  spline_knts, n_knts, 
    //                  spline_dgrs, n_dgrs, 
    //                  n_evals, n_dims);
    TIMERSTOP(gpu_overall)

    return 0;
}