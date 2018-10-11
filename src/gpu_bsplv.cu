#include "gpu_bsplv.h"

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
                printf("%d, ", Table->knots[offset + j]);
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
        table->nknots, n_coeffs*sizeof(table->nknots[0]), H2D);         CUERR
    total_nknots = 0;
    for(index_t i=0; i<table->ndim; ++i) 
    {
        cudaMemcpy(&(Knots[total_nknots]), table->knots[i], 
            table->nknots[i]*sizeof(table->knots[0]), H2D);             CUERR
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
        &Nknots, sizeof(Table->order), H2D);                            CUERR
    cudaMemcpy(&(Table->knots), 
        &Knots, sizeof(Table->order), H2D);                             CUERR
    cudaMemcpy(&(Table->coefficients), 
        &Coefficients, sizeof(Table->order), H2D);                      CUERR
    cudaMemcpy(&(Table->naxes), 
        &Naxes, sizeof(Table->order), H2D);                             CUERR
    cudaMemcpy(&(Table->strides), 
        &Strides, sizeof(Table->order), H2D);                           CUERR
}

/*
 * A method to check if the table has been properly copied.
 */
void GPU_BSPLV::print_table(
    index_t ndim)
{
    cudaSetDevice(0);                                                   CUERR

    index_t n_threads = 32;
    index_t n_blocks = 1;

    print_table_d<<<n_blocks, n_threads>>>(ndim, Table);                CUERR 
    cudaDeviceSynchronize();
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
