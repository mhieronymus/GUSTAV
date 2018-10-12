#ifndef _CPU_BSPLV_H
#define _CPU_BSPLV_H

#include <random>
#include <iostream>


#include "helper/types.hpp"
#include "helper/timing.hpp"
#include "helper/load.hpp"

// For reference 
#include "helper/photospline/bspline.h"

class CPU_BSPLV{
public:
    CPU_BSPLV(char * filename);
    ~CPU_BSPLV();

    void eval_splines(index_t n_evals, value_t * y_array);
    void eval_splines_vanilla(index_t n_evals, value_t * y_array);
    void eval_splines_simple(index_t n_evals, value_t * y_array);
    void load_table(char * filename);
    void free_table();

    splinetable * table;
private:
    void find_left_knot(value_t * knots, index_t n_knts, index_t & left,
        value_t x);
    value_t ndsplineeval_core_2(index_t * centers, index_t maxdegree, 
        value_t * biatx);
    void bsplvb_2(value_t * knots, index_t jhigh, index_t index, value_t x,
        index_t left, value_t * biatx, index_t nknots);
    value_t wiki_bsplines(value_t * knots, value_t x, index_t left, 
        index_t order);

};
#endif