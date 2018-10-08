#ifndef _CPU_BSPLV_H
#define _CPU_BSPLV_H

#include <random>
#include "helper/types.hpp"

void cpu_eval_splines(
    value_t * spline_cffs, index_t n_cffs, 
    value_t * spline_knts, index_t * n_knts, 
    index_t * spline_dgrs, index_t n_dgrs, 
    index_t n_evals, index_t n_dims, 
    index_t * splines_per_dim);

#endif