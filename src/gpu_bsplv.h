#ifndef _GPU_BSPLV_H
#define _GPU_BSPLV_H

#include <random>
#include "helper/types.hpp"
#include "helper/timing.hpp"

class GPU_BSPLV{
public:
    GPU_BSPLV(splinetable * table);
    ~GPU_BSPLV();

    void gpu_eval_splines(Splinetable * Table, index_t ndim, index_t n_evals,
        value_t * y_array);

    void print_table(index_t ndim);
private:
    // Device functions 
    Splinetable * Table = nullptr;
    index_t * Order = nullptr;
    value_t * Knots = nullptr, * Coefficients = nullptr;
    long * Nknots = nullptr, * Naxes = nullptr;
    unsigned long * Strides = nullptr;
};
#endif