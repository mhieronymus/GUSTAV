#include "gpu_bsplv.h"

/*
 * Allocate memory on the GPU and move over the splinetable.
 * We ommit stuff such as aux. Hence we don't copy the whole struct but we
 * have to bind the arrays to the struct.
 */
 GPU_BSPLV::GPU_BSPLV(splinetable * table)
{
    cudaMalloc(&Table, sizeof(*Table));    

    // Allocate the arrays
    cudaMalloc(&Order, table->ndim*sizeof(*Order));
    // Get the total number of knots
    index_t total_nknots = 0;
    for(index_t k=0; k<table->ndim; ++k) total_nknots += table->nknots[k];
    cudaMalloc(&Knots, total_nknots*sizeof(*Knots));
    cudaMalloc(&Nknots, table->ndim*sizeof(*Nknots));
    // Extends?
    // Periods?
    // Get the total number of coefficients
    index_t n_coeffs = 0;
    for(index_t c=0; c<table->ndim; ++c) n_coeffs += table->naxes[c];
    cudaMalloc(&Coefficients, n_coeffs*sizeof(*Coefficients));
    cudaMalloc(&Naxes, table->ndim*sizeof(*Naxes));
    cudaMalloc(&Strides, table->ndim*sizeof(*Strides));

    // Copy the data
    cudaMemcpy(Order, 
        table->order, table->ndim*sizeof(table->order[0]), H2D);
    cudaMemcpy(Nknots, 
        table->nknots, n_coeffs*sizeof(table->nknots[0]), H2D);
    total_nknots = 0;
    for(index_t i=0; i<table->ndim; ++i) 
    {
        cudaMemcpy(&(Knots[total_nknots]), 
            table->knots[i], table->nknots[i]*sizeof(table->knots[0]), H2D);
        total_nknots += table->nknots[i];
    }
    cudaMemcpy(Coefficients, 
        table->coefficients, n_coeffs*sizeof(table->coefficients[0]), H2D);
    cudaMemcpy(Naxes, 
        table->naxes, table->ndim*sizeof(table->naxes[0]), H2D);
    cudaMemcpy(Strides, 
        table->strides, table->ndim*sizeof(table->strides[0]), H2D);

    // Bind the data
    cudaMemcpy(&(Table->order), 
        &Order, sizeof(Table->order), H2D);
    cudaMemcpy(&(Table->nknots), 
        &Nknots, sizeof(Table->order), H2D);
    cudaMemcpy(&(Table->knots), 
        &Knots, sizeof(Table->order), H2D);
    cudaMemcpy(&(Table->coefficients), 
        &Coefficients, sizeof(Table->order), H2D);
    cudaMemcpy(&(Table->naxes), 
        &Naxes, sizeof(Table->order), H2D);
    cudaMemcpy(&(Table->strides), 
        &Strides, sizeof(Table->order), H2D);
}


// void GPU_BSPLV::table_to_gpu(
//     splinetable * table)
// {
//     Splinetable * Table = nullptr;
//     cudaMalloc(&Table, sizeof(*Table));
//     index_t * Order = nullptr;
//     value_t * Knots = nullptr, * Coefficients = nullptr;
//     long * Nknots = nullptr, * Naxes = nullptr;
//     unsigned long * Strides = nullptr;

//     // Allocate the arrays
//     cudaMalloc(&Order, table->ndim*sizeof(*Order));
//     // Get the total number of knots
//     index_t total_nknots = 0;
//     for(index_t k=0; k<table->ndim; ++k) total_nknots += table->nknots[k];
//     cudaMalloc(&Knots, total_nknots*sizeof(*Knots));
//     cudaMalloc(&Nknots, table->ndim*sizeof(*Nknots));
//     // Extends?
//     // Periods?
//     // Get the total number of coefficients
//     index_t n_coeffs = 0;
//     for(index_t c=0; c<table->ndim; ++c) n_coeffs += table->naxes[c];
//     cudaMalloc(&Coefficients, n_coeffs*sizeof(*Coefficients));
//     cudaMalloc(&Naxes, table->ndim*sizeof(*Naxes));
//     cudaMalloc(&Strides, table->ndim*sizeof(*Strides));

//     // Copy the data
//     cudaMemcpy(Order, 
//         table->order, table->ndim*sizeof(table->order[0]), H2D);
//     cudaMemcpy(Nknots, 
//         table->nknots, n_coeffs*sizeof(table->nknots[0]), H2D);
//     total_nknots = 0;
//     for(index_t i=0; i<table->ndim; ++i) 
//     {
//         cudaMemcpy(&(Knots[total_nknots]), 
//             table->knots[i], table->nknots[i]*sizeof(table->knots[0]), H2D);
//         total_nknots += table->nknots[i];
//     }
//     cudaMemcpy(Coefficients, 
//         table->coefficients, n_coeffs*sizeof(table->coefficients[0]), H2D);
//     cudaMemcpy(Naxes, 
//         table->naxes, table->ndim*sizeof(table->naxes[0]), H2D);
//     cudaMemcpy(Strides, 
//         table->strides, table->ndim*sizeof(table->strides[0]), H2D);

//     // Bind the data
//     cudaMemcpy(&(Table->order), 
//         &Order, sizeof(Table->order), H2D);
//     cudaMemcpy(&(Table->nknots), 
//         &Nknots, sizeof(Table->order), H2D);
//     cudaMemcpy(&(Table->knots), 
//         &Knots, sizeof(Table->order), H2D);
//     cudaMemcpy(&(Table->coefficients), 
//         &Coefficients, sizeof(Table->order), H2D);
//     cudaMemcpy(&(Table->naxes), 
//         &Naxes, sizeof(Table->order), H2D);
//     cudaMemcpy(&(Table->strides), 
//         &Strides, sizeof(Table->order), H2D);
//     return Table;
//     // // cudaMalloc((void**)&Table, sizeof(Splinetable));
//     // // Table->order = nullptr;
    
//     // // Allocate the arrays
//     // cudaMalloc(&(Table->order), table->ndim*sizeof(table->order[0]));
//     // // Get the total number of knots
//     // index_t total_nknots = 0;
//     // for(index_t k=0; k<table->ndim; ++k) total_nknots += table->nknots[k];
//     // cudaMalloc(&(Table->knots), total_nknots*sizeof(table->knots[0]));
//     // cudaMalloc(&(Table->nknots), table->ndim*sizeof(table->nknots[0]));
//     // // Extends?
//     // // Periods?
//     // // Get the total number of coefficients
//     // index_t n_coeffs = 0;
//     // for(index_t c=0; c<table->ndim; ++c) n_coeffs += table->naxes[c];
//     // cudaMalloc(&(Table->coefficients), n_coeffs*sizeof(table->coefficients[0]));
//     // cudaMalloc(&(Table->naxes), table->ndim*sizeof(table->naxes[0]));
//     // cudaMalloc(&(Table->strides), table->ndim*sizeof(table->strides[0]));

//     // // Copy the arrays
//     // cudaMemcpy(Table->order, 
//     //     table->order, table->ndim*sizeof(table->order[0]), H2D);
//     // total_nknots = 0;
//     // for(index_t i=0; i<table->ndim; ++i) 
//     // {
//     //     cudaMemcpy(&(Table->knots[total_nknots]), 
//     //         table->knots[i], table->nknots[i]*sizeof(table->knots[0]), H2D);
//     //     total_nknots += table->nknots[i];
//     // }

//     // cudaMemcpy(Table->coefficients, 
//     //     table->coefficients, n_coeffs*sizeof(table->coefficients[0]), H2D);
//     // cudaMemcpy(Table->naxes, 
//     //     table->naxes, table->ndim*sizeof(table->naxes[0]), H2D);
//     // cudaMemcpy(Table->strides, 
//     //     table->strides, table->ndim*sizeof(table->strides[0]), H2D);

//     // return Table;
// }

// void GPU_BSPLV::gpu_eval_splines(
//     Splinetable * Table, 
//     index_t ndim,
//     index_t n_evals,
//     value_t * y_array)
// {

// }

GPU_BSPLV::~GPU_BSPLV()
{
    cudaFree(Order);
    cudaFree(Knots);
    cudaFree(Coefficients);
    cudaFree(Naxes);
    cudaFree(Strides);
    cudaFree(Table);
}
