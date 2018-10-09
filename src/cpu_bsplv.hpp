#ifndef _CPU_BSPLV_H
#define _CPU_BSPLV_H

#include <random>
#include "helper/types.hpp"

// For reference 
#include "helper/photospline/bspline.h"

void cpu_eval_splines(splinetable * table, index_t n_evals);

#endif