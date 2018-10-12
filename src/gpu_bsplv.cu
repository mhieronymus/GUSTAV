#include "gpu_bsplv.h"
// CUDA libraries have to be imported here
#include <curand.h>
#include <curand_kernel.h>

// Here are device functions. The class is further below.
__global__
void print_table_d(
    index_t ndim, 
    Splinetable * Table)
{
    index_t id = threadIdx.x;
    if(id == 0) 
    {
        printf("Some information\nOrder\n");
        for(index_t i=0; i<ndim; ++i)
            printf("%d, ", Table->order[i]);
        printf("\n");

        printf("\nNumber of knots\n");
        index_t n_total_knts = 0;
        for(index_t i=0; i<ndim; ++i) 
        {
                n_total_knts += Table->nknots[i];
                printf("%d, ", Table->nknots[i]);
        }
        printf("\nTotal of %d knots\n", n_total_knts);

        printf("\nThe knots of each dimension\n");
        index_t offset = 0;
        for(index_t i=0; i<ndim; ++i) 
        {
           printf("Total knots here: %d\n", Table->nknots[i]);
            for(index_t j=0; j<Table->nknots[i]; ++j) {
                printf("%f, ", Table->knots[offset + j]);
            }
            printf("\n");
            offset += Table->nknots[i];
        }

        printf("\ncoefficients\n");
        index_t arraysize = 1;
        for(index_t i=0; i<ndim-1; ++i)
            arraysize *= Table->naxes[i];
        for(index_t i=0; i<arraysize && i<10; ++i)
            printf("%f, ", Table->coefficients[i]);
        printf("\nTotal of %d coefficients\n", arraysize);

        printf("\nstrides\n");
        for(index_t i=0; i<ndim; ++i) 
        {
            printf("%d, ", Table->strides[i]);
        }
        
       printf("\nnaxes\n");
        for(index_t i=0; i<ndim; ++i)
        printf("%d, ", Table->naxes[i]);
        printf("\nFinished with those prints\n");
    }

}

/*
 * Allocate memory on the GPU and move over the splinetable.
 * We ommit stuff such as aux. Hence we don't copy the whole struct but we
 * have to bind the arrays to the struct.
 */
 GPU_BSPLV::GPU_BSPLV(splinetable * table)
{
    cudaMalloc(&Table, sizeof(*Table));                                 CUERR

    // Allocate the arrays
    cudaMalloc(&Order, table->ndim*sizeof(*Order));                     CUERR
    // Get the total number of knots
    index_t total_nknots = 0;
    for(index_t k=0; k<table->ndim; ++k) total_nknots += table->nknots[k];
    cudaMalloc(&Knots, total_nknots*sizeof(*Knots));                    CUERR
    cudaMalloc(&Nknots, table->ndim*sizeof(*Nknots));                   CUERR
    // Extends?
    // Periods?
    // Get the total number of coefficients
    index_t n_coeffs = 0;
    for(index_t c=0; c<table->ndim; ++c) n_coeffs += table->naxes[c];
    cudaMalloc(&Coefficients, n_coeffs*sizeof(*Coefficients));          CUERR
    cudaMalloc(&Naxes, table->ndim*sizeof(*Naxes));                     CUERR
    cudaMalloc(&Strides, table->ndim*sizeof(*Strides));                 CUERR

    // Copy the data
    cudaMemcpy(Order, 
        table->order, table->ndim*sizeof(table->order[0]), H2D);        CUERR
    cudaMemcpy(Nknots, 
        table->nknots, table->ndim*sizeof(table->nknots[0]), H2D);      CUERR
    // TODO: This: Something doesn't work out here
    total_nknots = 0;
    for(index_t i=0; i<table->ndim; ++i) 
    {
        cudaMemcpy(Knots+total_nknots, table->knots[i], 
            table->nknots[i]*sizeof(table->knots[0][0]), H2D);             CUERR
        total_nknots += table->nknots[i];
    }
    cudaMemcpy(Coefficients, table->coefficients, 
        n_coeffs*sizeof(table->coefficients[0]), H2D);                  CUERR
    cudaMemcpy(Naxes, 
        table->naxes, table->ndim*sizeof(table->naxes[0]), H2D);        CUERR
    cudaMemcpy(Strides, 
        table->strides, table->ndim*sizeof(table->strides[0]), H2D);    CUERR

    // Bind the data
    cudaMemcpy(&(Table->order), 
        &Order, sizeof(Table->order), H2D);                             CUERR
    cudaMemcpy(&(Table->nknots), 
        &Nknots, sizeof(Table->nknots), H2D);                           CUERR
    cudaMemcpy(&(Table->knots), 
        &Knots, sizeof(Table->knots), H2D);                             CUERR
    cudaMemcpy(&(Table->coefficients), 
        &Coefficients, sizeof(Table->coefficients), H2D);               CUERR
    cudaMemcpy(&(Table->naxes), 
        &Naxes, sizeof(Table->naxes), H2D);                             CUERR
    cudaMemcpy(&(Table->strides), 
        &Strides, sizeof(Table->strides), H2D);                         CUERR

    ndim = table->ndim;
    maxdegree = 0;
    for(index_t d=0; d<table->ndim; d++) 
        maxdegree = maxdegree > table->order[d] ? maxdegree : table->order[d];
}

__device__
void tablesearchcenters_device(
    Splinetable * Table,
    value_t * par,
    index_t * centers,
    index_t ndim)
{

}

__device__ 
value_t ndsplineeval_core_device(
    index_t * centers,
    index_t maxdegree,
    value_t * biatx,
    index_t ndim) 
{
    return 0;
}

__device__ 
void bsplvb_device(
    value_t * knots,
    index_t jhigh,
    index_t index,
    value_t x,
    index_t left,
    value_t * biatx,
    index_t nknots)
{

}

__global__
void eval_splines_device(
    value_t * Y,
    index_t n_evals,
    value_t * range,
    index_t ndim,
    index_t maxdegree,
    Splinetable * Table,
    value_t * biatx_cache,
    value_t * par_cache,
    index_t * centers_cache)
{
    int id = threadIdx.x + blockDim.x * blockIdx.x;
    curandState state;
    // seed, subsequence to use, offset within the subsequence
    curand_init(42, id, 0, &state);
    for(index_t i=id; id<n_evals; i+=blockDim.x*gridDim.x)
    {
        value_t * biatx = biatx_cache + id*(ndim*(maxdegree+1));
        value_t * par = par_cache + id*ndim;
        index_t * centers = centers_cache + id*ndim;

        index_t offset = 0;
        for(index_t d=0; d<ndim; ++d) 
        {
            par[d] = range[d] * curand_uniform(&state) + Table->knots[offset];
            offset += Table->nknots[d];
        }

        tablesearchcenters_device(Table, par, centers, ndim);
        offset = 0;
        for(index_t d=0; d<ndim; ++d)
        {
            bsplvb_device(
                &(Table->knots[offset]),
                Table->order[d],
                1,
                par[d],
                centers[d],
                &(biatx[d*(maxdegree+1)]),
                Table->nknots[d]);
            offset += Table->nknots[d];
        }
        Y[i] = ndsplineeval_core_device(centers, maxdegree, biatx, ndim);
    }
}

GPU_BSPLV::~GPU_BSPLV()
{
    cudaFree(Order);
    cudaFree(Knots);
    cudaFree(Coefficients);
    cudaFree(Naxes);
    cudaFree(Strides);
    cudaFree(Table);
}

/*
 * A method to check if the table has been properly copied.
 */
void GPU_BSPLV::print_table()
{
    cudaSetDevice(0);                                                   CUERR

    index_t n_threads = 32;
    index_t n_blocks = 1;

    print_table_d<<<n_blocks, n_threads>>>(ndim, Table);                CUERR 
    cudaDeviceSynchronize();
}

void GPU_BSPLV::eval_splines(
    splinetable * table,
    index_t n_evals,
    value_t * y_array)
{
    index_t threads = 1024;
    // There are (ususally) max 64 warps per SM.
    // The GTX 1080 has five SMs per GPC and four GPCs
    index_t blocks = SDIV(n_evals, threads) > 40 ? 40 : SDIV(n_evals, threads);
    value_t * Cache = nullptr;
    index_t * Centers_cache = nullptr;
    index_t par_cache = ndim * blocks * threads;
    index_t biatx_cache = blocks * threads * ndim * (maxdegree+1);
    cudaMalloc(&Cache, (par_cache+biatx_cache) * sizeof(value_t));      CUERR
    cudaMalloc(&Centers_cache, par_cache * sizeof(index_t));            CUERR

    // We add another array that is needed to create random data 
    value_t * Range = nullptr, * range = nullptr, * Y = nullptr;
    cudaMallocHost(&range, ndim*sizeof(value_t));                       CUERR 
    cudaMalloc(&Range, ndim*sizeof(value_t));                           CUERR 
    cudaMalloc(&Y, ndim*sizeof(value_t));                               CUERR
    // Some cache that is needed for par, biatx and centers


    for(index_t d=0; d<table->ndim; d++) 
        range[d] = table->knots[d][table->nknots[d]-1] 
                   - table->knots[d][0];
    cudaMemcpy(Range, range, ndim*sizeof(value_t), H2D);                CUERR 

    
    eval_splines_device<<<blocks, threads>>>(
        Y,                  // The evaluations
        n_evals,            // Length of Y
        Range,              // Range of the knots in each dimension
        ndim,               // Dimensionality
        maxdegree,          // Maximum degree of the splines
        Table,              // The splinetable
        Cache,              // biatx cache
        Cache+biatx_cache,  // par cache
        Centers_cache);                                                 CUERR


    cudaFreeHost(range);                                                CUERR
    cudaFree(Range);                                                    CUERR
    cudaFree(Y);                                                        CUERR
    cudaFree(Cache);                                                    CUERR
    cudaFree(Centers_cache);                                            CUERR
}
