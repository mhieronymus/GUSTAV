#include <iostream>
#include "helper/types.hpp"
#include "helper/timing.hpp"
#include "cpu_bsplv.hpp"
#include "gpu_bsplv.h"

const bool DEBUG = false;

int main(int argc, char* argv[]) {
    
    bool verbose = false;            // Print some data about the splines
    index_t n_evals = 1 << 23;       // Number of times to evaluate the splines
    n_evals = 1 << 10;
    value_t * y_array = new value_t[n_evals];
    CPU_BSPLV * cpu = nullptr;
    if(argc > 1) {
        char* filename = argv[1];
        cpu = new CPU_BSPLV(filename);
    } else {
        std::cout << "Please enter a path to a fits file with photosplines ";
        std::cout << "such as\n";
        std::cout << "/gpfs/fs1/home/mhierony/data/splines/LowEnergyMieCascades-ShorterStack-560.abs.fits\n";
        std::cout << "falling back to default\n";
        cpu = new CPU_BSPLV("/gpfs/fs1/home/mhierony/data/splines/LowEnergyMieCascades-ShorterStack-560.abs.fits");
    }

    splinetable * table = cpu->table;
    
    if(verbose)
    {
        std::cout << "Some information about the loaded splines\n";
        std::cout << "Number of dimensions " << table->ndim << "\n";
        std::cout << "Order\n";
        for(index_t i=0; i<table->ndim; ++i) 
        {
            std::cout << table->order[i] << ", ";
        }

        std::cout << "\nNumber of knots\n";
        index_t n_total_knts = 0;
        for(index_t i=0; i<table->ndim; ++i) 
        {
                n_total_knts += table->nknots[i];
                std::cout << table->nknots[i] << ", ";
        }
        std::cout << "\nTotal of " << n_total_knts << " knots\n";

        std::cout << "\nThe knots of each dimension\n";
        for(index_t i=0; i<table->ndim; ++i) 
        {
            std::cout << "Total knots here: " << table->nknots[i] << "\n";
            for(index_t j=0; j<table->nknots[i]; ++j) {
                std::cout << table->knots[i][j] << ", ";
            }
            std::cout << "\n";
        }

        std::cout << "\ncoefficients\n";
        index_t arraysize = 1;
        for(index_t i=0; i<table->ndim-1; ++i)
            arraysize *= table->naxes[i];
        for(index_t i=0; i<arraysize && i<10; ++i)
            std::cout << table->coefficients[i] << ", ";
        std::cout << "\nTotal of " << arraysize << " coefficients\n";

        std::cout << "\nstrides\n";
        for(index_t i=0; i<table->ndim; ++i) 
        {
            std::cout << table->strides[i] << ", ";
        }
        
        std::cout << "\nnaxes\n";
        for(index_t i=0; i<table->ndim; ++i)
            std::cout << table->naxes[i] << ", ";
        std::cout << "\nFinished with those prints\n\n\n" << std::flush;
    }

    start("cpu_vanilla");
    cpu->eval_splines_vanilla(n_evals, y_array);
    stop();
    std::cout << "\nSome results\n";
    for(index_t i=0; i<100 && i<n_evals; i +=10)
        std::cout << y_array[i] << ", ";
    std::cout << "\n\n";

    start("cpu");
    cpu->eval_splines(n_evals, y_array);
    stop();
    std::cout << "\nSome results\n";
    for(index_t i=0; i<100 && i<n_evals; i +=10)
        std::cout << y_array[i] << ", ";
    std::cout << "\n\n";

    start("cpu_simple");
    cpu->eval_splines_simple(n_evals, y_array);
    stop();
    std::cout << "\nSome results\n";
    for(index_t i=0; i<100 && i<n_evals; i +=10)
        std::cout << y_array[i] << ", ";
    std::cout << "\n\n";
    
    start("H2D");
    GPU_BSPLV * gpu = new GPU_BSPLV(table);
    stop();
    if(verbose)
        gpu->print_table();

    start("gpu");
    gpu->eval_splines(table, n_evals, y_array);
    stop();
    // TIMERSTART(gpu_overall)
    // gpu_eval_splines(spline_cffs, n_cffs, 
    //                 spline_knts, n_knts, 
    //                 spline_dgrs, n_dgrs, 
    //                 n_evals, n_dims);
    // TIMERSTOP(gpu_overall)
 
    // Free stuff
    delete cpu;
    delete gpu;
    
    return 0;
}