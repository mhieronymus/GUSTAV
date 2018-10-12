#include "gpu_bsplv.h"
// CUDA libraries have to be imported here
#include <curand.h>
#include <curand_kernel.h>

// TODO: Move caches and delta_l and delta_r to shared memory

bool verbose = false; // Print the amount of allocated memory

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

__device__
void tablesearchcenters_device(
    Splinetable * Table,
    value_t * par,
    index_t * centers,
    index_t ndim)
{
    index_t offset = 0;
    for(index_t i=0; i<ndim; ++i)
    {
        // Check if outside of the table
        // Should (!) not happen and needs to be handled by giving an error
        if(par[i] <= Table->knots[offset] ||
           par[i] > Table->knots[offset + Table->nknots[i]-1])
           return;
        if(par[i] < Table->knots[offset + Table->order[i]])
        {
            centers[i] = Table->order[i];
            offset += Table->nknots[i];
            continue;
        } else if(par[i] >= Table->knots[offset + Table->naxes[i]])
        {
            centers[i] = Table->naxes[i]-1;
            offset += Table->nknots[i];
            continue;
        }

        index_t min = Table->order[i];
        index_t max = Table->nknots[i]-2;
        do 
        {
            centers[i] = (max+min)/2;
            if(par[i] < Table->knots[offset + centers[i]])
                max = centers[i]-1;
            else 
                min = centers[i]+1;
        } while(par[i] < Table->knots[offset + centers[i]] ||
            par[i] >= Table->knots[offset + centers[i]+1]);
        
        if(centers[i] == Table->naxes[i]) centers[i]--;
        offset += Table->nknots[i];
    }
}

// TODO: Move basis_tree and decomposed_pos to shared memory
__device__ 
value_t ndsplineeval_core_device(
    Splinetable * Table,
    index_t * centers,
    index_t maxdegree,
    value_t * biatx,
    index_t ndim) 
{
    value_t result = 0;
    value_t basis_tree[7]; // ndim+1
    index_t decomposed_pos[6]; // ndim

    index_t tablepos = 0;
    basis_tree[0] = 1;
    index_t nchunks = 1;
    
    for(index_t i=0; i<ndim; ++i)
    {
        decomposed_pos[i] = 0;
        tablepos += (centers[i]-Table->order[i]) * Table->strides[i];
        basis_tree[i+1] = basis_tree[i] * biatx[i*(maxdegree+1)];
    }
    for(index_t i=0; i<ndim - 1; ++i)
        nchunks *= (Table->order[i] + 1);
    index_t n = 0;
    while(true)
    {
        for(index_t i=0; i<Table->order[ndim-1]+1; ++i)
        {
            result += basis_tree[ndim-1] 
                * biatx[(ndim-1)*(maxdegree+1) + i] 
                * Table->coefficients[tablepos+i];
        }
        if(++n == nchunks) break;

        tablepos += Table->strides[ndim-2];
        decomposed_pos[ndim-2]++;

        // Now to higher dimensions 
        index_t i;
        for(i=ndim-2; decomposed_pos[i]>Table->order[i]; --i)
        {
            decomposed_pos[i-1]++;
            tablepos += (Table->strides[i-1]
                - decomposed_pos[i]*Table->strides[i]);
            decomposed_pos[i] = 0;
        }
        // for(index_t j=i; j < table->ndim-1; ++j)
        for (index_t j = i; j < ndim-1; ++j)
            basis_tree[j+1] = basis_tree[j] 
                * biatx[j*(maxdegree+1) + decomposed_pos[j]];
    }
    return result;
}

/*
 * Evaluate splines but sequentially. There is no parallel execution here.
 * However, in later implementations one could make use of a warp.
 */
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
    index_t j = 0;
    // In case if x is outside of the full support of the spline surface.
    if(left == jhigh) 
    {
        while(left >= 0 && x < knots[left]) left--;
    } else if(left == nknots-jhigh-2) 
    {
        while (left < nknots-1 && x > knots[left+1]) left++;	
    }

    if(index != 2) { 
        biatx[j] = 1;
        // Check if no more columns need to be evaluated.
        if(j >= jhigh) return;
    }
    // For simplicity we assume maxdegree to be less than 6
    value_t delta_r[6];
    value_t delta_l[6];
    do 
    {        
        delta_r[j] = knots[left + j + 1] - x;
        delta_l[j] = x - knots[left-j];
        value_t saved = 0;
        for(index_t i=0; i<=j; ++i) 
        {
            value_t term = biatx[i] / (delta_r[i] + delta_l[j-i]);
            biatx[i] = saved + delta_r[i] * term;
            saved = delta_l[j-i] * term;
        }
        biatx[j+1] = saved;
        j++;
        
    } while(j < jhigh); // shouldn't that be sufficient?

    /* 
	 * If left < (spline order), only the first (left+1)
	 * splines are valid; the remainder are utter nonsense.
     * TODO: Check why
	 */
    index_t i = jhigh-left;
    if (i > 0) {
        for (j = 0; j < left+1; j++)
			biatx[j] = biatx[j+i]; /* Move valid splines over. */
		for ( ; j <= jhigh; j++)
			biatx[j] = 0.0; /* The rest are zero by construction. */
    }
    i = left+jhigh+2-nknots;
    if (i > 0) {
        for (j = jhigh; j > i-1; j--)
			biatx[j] = biatx[j-i];
		for ( ; j >= 0; j--)
			biatx[j] = 0.0;
    }

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
    for(index_t i=id; i<n_evals; i+=blockDim.x*gridDim.x)
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
        Y[i] = ndsplineeval_core_device(Table, centers, 
            maxdegree, biatx, ndim);
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

    if(verbose)
    {
        value_t memory = table->ndim*sizeof(*Order) 
            + total_nknots*sizeof(*Knots)
            + table->ndim*sizeof(*Nknots)
            + n_coeffs*sizeof(*Coefficients)
            + table->ndim*sizeof(*Naxes)
            + table->ndim*sizeof(*Strides);
        memory /= (1024*8);
        std::cout << "\nAllocating " 
            << memory << " kByte on the GPU for the spline\n";
        size_t free_mem, total_mem;
        cudaMemGetInfo(&free_mem, &total_mem);
        std::cout << "Free Memory: " 
            << ((value_t) free_mem) / (1024*1024*8) 
            << " MBytes of " << ((value_t) total_mem) / (1024*1024*8)  << "\n";

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
    // Debug info tells us, we use 160 registers per thread
    // There are 256 KB of registers per SM in a GTX 1080
    // That is 256*1024/(8*160) many threads that can run on a SM
    // That is 204,8 threads
    index_t threads = 192;
    // There are (ususally) max 64 warps per SM.
    // The GTX 1080 has five SMs per GPC and four GPCs
    index_t blocks = SDIV(n_evals, threads) > 40 ? 40 : SDIV(n_evals, threads);
    value_t * Cache = nullptr;
    index_t * Centers_cache = nullptr;
    index_t par_cache = ndim * blocks * threads;
    index_t biatx_cache = blocks * threads * ndim * (maxdegree+1);

    if(verbose)
    {
        value_t memory = (par_cache+biatx_cache) * sizeof(value_t)
            + par_cache * sizeof(index_t);
        memory /= (1024*1024*8);
        std::cout << "\nAllocating " << memory << " MBytes for caching\n";
        value_t memory_2 = (ndim*sizeof(value_t)) / (1024*1024*8);
        std::cout << "\nAllocating " << memory_2 
            << " MBytes for the results\n";
        memory += memory_2 + (ndim*sizeof(value_t))/(1024*1024*8);
        std::cout << "\nTotal allocated for evaluating: " 
            << memory << " MBytes \n" << std::flush;
        std::cout << "Running with " << threads << " threads and " 
            << blocks << " blocks\n";
    }

    cudaMalloc(&Cache, (par_cache+biatx_cache) * sizeof(value_t));      CUERR
    cudaMalloc(&Centers_cache, par_cache * sizeof(index_t));            CUERR

    // We add another array that is needed to create random data 
    value_t * Range = nullptr, * range = nullptr, * Y = nullptr;
    cudaMallocHost(&range, ndim*sizeof(value_t));                       CUERR 
    cudaMalloc(&Range, ndim*sizeof(value_t));                           CUERR 
    cudaMalloc(&Y, n_evals*sizeof(value_t));                            CUERR
    // Some cache that is needed for par, biatx and centers

    for(index_t d=0; d<table->ndim; d++) 
        range[d] = table->knots[d][table->nknots[d]-1] 
                   - table->knots[d][0];
    cudaMemcpy(Range, range, ndim*sizeof(value_t), H2D);                CUERR 
    if(verbose)
    {
        size_t free_mem, total_mem;
        cudaMemGetInfo(&free_mem, &total_mem);
        std::cout << "Free Memory: " 
            << ((value_t) free_mem) / (1024*1024*8) 
            << " MBytes of " << ((value_t) total_mem) / (1024*1024*8)  << "\n";
    }
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
