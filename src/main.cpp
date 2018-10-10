#include <iostream>
#include "helper/types.hpp"
#include "helper/load.hpp"
#include "helper/timing.hpp"
#include "cpu_bsplv.hpp"
#include "gpu_bsplv.h"

const bool DEBUG = false;

int main(int argc, char* argv[]) {
    
    // value_t * spline_cffs = nullptr; // Coefficients
    // value_t * spline_knts = nullptr; // knots
    // index_t * spline_dgrs = nullptr; // degrees 
    // index_t * n_knts      = nullptr; // Number of knots per dimension
    // index_t * n_splines   = nullptr; // Number of splines per dimension
    // index_t n_cffs, n_dgrs, n_dims;
    splinetable * tablestruct_ = new splinetable();
    index_t n_evals = 1 << 10;       // Number of times to evaluate the splines
    n_evals = 4;
    if(argc > 1) {
        char* filename = argv[1];
        std::cout << "loading from " << filename << "\n";
        // TIMERSTART(load_from_file)
        // Load a spline file from a file given by the command line
        // load_splines(spline_cffs, n_cffs, 
        //             spline_knts, n_knts, 
        //             spline_dgrs, n_dgrs, 
        //             n_dims, n_splines, filename);
        if(!load_splines(tablestruct_, filename))
        {
            std::cout << "Failed to load. Abort the mission!\n";
            std::cout << "I repeat: Abort the mission!\n";
            return 0;
        }
        // TIMERSTOP(load_from_file)
    } else {
        std::cout << "Please enter a path to a fits file with photosplines ";
        std::cout << "such as\n";
        std::cout << "/home/mhierony/Documents/Uni/Masterthesis/data/splines/LowEnergyMieCascades-ShorterStack-560.abs.fits\n";
        std::cout << "falling back to default\n";
        // TIMERSTART(load_from_file)
        // Load a spline file from a file given by the command line
        // load_splines(spline_cffs, n_cffs, 
        //             spline_knts, n_knts, 
        //             spline_dgrs, n_dgrs, 
        //             n_dims, n_splines,
        //             "/home/mhierony/Documents/Uni/Masterthesis/data/splines/LowEnergyMieCascades-ShorterStack-560.abs.fits");
        if(!load_splines(tablestruct_, "/home/mhierony/Documents/Uni/Masterthesis/data/splines/LowEnergyMieCascades-ShorterStack-560.abs.fits"))
        {
            std::cout << "Failed to load. Abort the mission!\n";
            std::cout << "I repeat: Abort the mission!\n";
            return 0;
        }
        
        // TIMERSTOP(load_from_file)
    }

    std::cout << "Some information about the loaded splines\n";
    std::cout << "Number of dimensions " << tablestruct_->ndim << "\n";
    std::cout << "Order\n";
    for(index_t i=0; i<tablestruct_->ndim; ++i) {
        std::cout << tablestruct_->order[i] << ", ";
    }

    std::cout << "\nNumber of knots\n";
    index_t n_total_knts = 0;
    for(index_t i=0; i<tablestruct_->ndim; ++i) {
            n_total_knts += tablestruct_->nknots[i];
            std::cout << tablestruct_->nknots[i] << ", ";
    }
    std::cout << "\nTotal of " << n_total_knts << " knots\n";

    std::cout << "\nThe knots of each dimension\n";
    for(index_t i=0; i<tablestruct_->ndim; ++i) 
    {
        std::cout << "Total knots here: " << tablestruct_->nknots[i] << "\n";
        for(index_t j=0; j<tablestruct_->nknots[i]; ++j) {
            std::cout << tablestruct_->knots[i][j] << ", ";
        }
        std::cout << "\n";
    }

    std::cout << "\ncoefficients\n";
    index_t arraysize = 1;
    for(index_t i=0; i<tablestruct_->ndim-1; ++i)
        arraysize *= tablestruct_->naxes[i];
    for(index_t i=0; i<arraysize && i<10; ++i)
        std::cout << tablestruct_->coefficients[i] << ", ";
    std::cout << "\nTotal of " << arraysize << " coefficients\n";

    std::cout << "\nstrides\n";
    for(index_t i=0; i<tablestruct_->ndim; ++i) {
        std::cout << tablestruct_->strides[i] << ", ";
    }
    
    std::cout << "\nnaxes\n";
    for(index_t i=0; i<tablestruct_->ndim; ++i)
        std::cout << tablestruct_->naxes[i] << ", ";
    std::cout << "\nFinished with those prints\n";

    if(!DEBUG) 
    {
        TIMERSTART(cpu)
        cpu_eval_splines(tablestruct_, n_evals);
        TIMERSTOP(cpu)

        TIMERSTART(gpu_overall)
        // gpu_eval_splines(spline_cffs, n_cffs, 
        //                 spline_knts, n_knts, 
        //                 spline_dgrs, n_dgrs, 
        //                 n_evals, n_dims);
        TIMERSTOP(gpu_overall)
    } else {
        cpu_eval_splines(tablestruct_, n_evals);
    }

    
    // Free stuff
    for(index_t i=0; i<tablestruct_->ndim; ++i) 
        free(tablestruct_->knots[i] - tablestruct_->order[i]);
    
    free(tablestruct_->knots);
    free(tablestruct_->order);
    free(tablestruct_->nknots);
    free(tablestruct_->periods);
    free(tablestruct_->coefficients);
    free(tablestruct_->naxes);
    free(tablestruct_->strides);
    free(tablestruct_->extents[0]);
    free(tablestruct_->extents);
    for(index_t i=0; i<tablestruct_->naux; ++i) 
    {
        free(tablestruct_->aux[i][1]);
        free(tablestruct_->aux[i][0]);
        free(tablestruct_->aux[i]);
    }
    free(tablestruct_->aux);
    
    return 0;
}