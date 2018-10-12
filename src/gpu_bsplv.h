#ifndef _GPU_BSPLV_H
#define _GPU_BSPLV_H

#include <random>
#include "helper/types.hpp"
#include "helper/timing.hpp"

class GPU_BSPLV{
public:
    GPU_BSPLV(splinetable * table);
    ~GPU_BSPLV();

    // table is only needed for generating points
    void eval_splines(splinetable * table, 
        index_t n_evals, value_t * y_array);

    void print_table();
private:

    Splinetable * Table = nullptr;
    index_t * Order = nullptr;
    value_t * Knots = nullptr, * Coefficients = nullptr;
    long * Nknots = nullptr, * Naxes = nullptr;
    unsigned long * Strides = nullptr;
    index_t ndim, maxdegree;
};
#endif