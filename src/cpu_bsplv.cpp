#include "cpu_bsplv.hpp"

CPU_BSPLV::CPU_BSPLV(
    char * filename) 
{
    load_table(filename);
}

CPU_BSPLV::~CPU_BSPLV()
{
    free_table();
}

// Binary search to find the leftmost spline with support.
// Uses left as first guess
void CPU_BSPLV::find_left_knot(
    value_t * knots,
    index_t n_knts,
    index_t & left,
    value_t x) 
{
    if(left == n_knts) {left--; return;}
    if(left == n_knts-1) return;
    if(knots[left] < x && knots[left+1] > x) return;
    index_t high = n_knts;
    do 
    {
        if(knots[left] > x) {
            high = left;
            left = left >> 1;
        } else 
        {
            left += (high - left) >> 1;
        }
    } while(knots[left] > x || knots[left+1] < x);
}

/*
 * Evaluate the result of a full spline basis given a set of knots, a position,
 * an order and the leftmost spline for the position (or -1) that has support.
 */ 
value_t CPU_BSPLV::ndsplineeval_core_2(
    index_t * centers, 
    index_t maxdegree, 
    value_t * biatx) 
{
    value_t result = 0;
    value_t basis_tree[table->ndim+1];
    index_t decomposed_pos[table->ndim];

    index_t tablepos = 0;
    basis_tree[0] = 1;
    index_t nchunks = 1;
    
    for(index_t i=0; i<table->ndim; ++i)
    {
        decomposed_pos[i] = 0;
        tablepos += (centers[i]-table->order[i]) * table->strides[i];
        basis_tree[i+1] = basis_tree[i] * biatx[i*(maxdegree+1)];
    }
    for(index_t i=0; i<table->ndim - 1; ++i)
        nchunks *= (table->order[i] + 1);
    index_t n = 0;
    while(true)
    {
        // We can expect this to be true most of the time, I guess?
        // I can test that against the normal approach later
        for (index_t i = 0; __builtin_expect(i < table->order[table->ndim-1] +
		    1, 1); i++)
        // for(index_t i=0; i<table->order[table->ndim-1]+1; ++i)
        {
            value_t tmp = basis_tree[table->ndim-1] ;
            tmp = biatx[(table->ndim-1)*(maxdegree+1) + i] ;
            tmp = table->coefficients[tablepos+i]; // invalid read
            result += basis_tree[table->ndim-1] 
                * biatx[(table->ndim-1)*(maxdegree+1) + i] 
                * table->coefficients[tablepos+i];
        }
        if (__builtin_expect(++n == nchunks, 0)) break;
        // if(++n == nchunks) break;

        tablepos += table->strides[table->ndim-2];
        decomposed_pos[table->ndim-2]++;

        // Now to higher dimensions 
        index_t i;
        for(i=table->ndim-2; decomposed_pos[i]>table->order[i]; --i)
        {
            decomposed_pos[i-1]++;
            tablepos += (table->strides[i-1]
                - decomposed_pos[i]*table->strides[i]);
            decomposed_pos[i] = 0;
        }
        // for(index_t j=i; j < table->ndim-1; ++j)
        for (index_t j = i; __builtin_expect(j < table->ndim-1, 1); ++j)
            basis_tree[j+1] = basis_tree[j] 
                * biatx[j*(maxdegree+1) + decomposed_pos[j]];
    }
    return result;
}

// See "A Practical Guide to Splines (revisited)" by Carl de Boor ()
// Chapter X, page 112f
/* See "A Practical Guide to Splines (revisited)" by Carl de Boor ()
 * Chapter X, page 112f
 * Input: knots  = knot sequence (non-decreasing)
 *        jhigh  = the number of columns-1 of values that shall be generated
 *                 should agree with the order
 *        index  = order of spline
 *        x      = the point at which the splines shall be evaluated
 *        left   = index with knots[left] <= x <= knots[left+1]
 *        biatx  = help array that stores the evaluations
 */
void CPU_BSPLV::bsplvb_2(
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
   
    value_t delta_r[jhigh];
    value_t delta_l[jhigh];
    do 
    {
        index_t jp1 = j+1;
        
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

value_t CPU_BSPLV::wiki_bsplines(
    value_t * knots, 
    value_t x, 
    index_t left, 
    index_t order)
{
    if(order == 0)
    {
        if(knots[left] <= x && x < knots[left+1]) return 1;
        return 0;
    }
    value_t result = (x-knots[left])  
        * wiki_bsplines(knots, x, left, order-1)
        / (knots[left+order] - knots[left]);
    result += (knots[left+order+1] - x) 
        * wiki_bsplines(knots, x, left+1, order-1)
        / (knots[left+order+1] - knots[left+1]);
    return result;

}


// Create n_evals times a random array of size n_dims and evaluate the spline
void CPU_BSPLV::eval_splines(
    index_t n_evals,
    value_t * y_array) 
{
    std::mt19937 engine(42);
    std::uniform_real_distribution<value_t> normal_dist(0, 1);

    value_t range[table->ndim];
    index_t maxdegree = 0;
    for(index_t d=0; d<table->ndim; d++) 
    {
        range[d] = table->knots[d][table->nknots[d]-1] 
                   - table->knots[d][0];
        maxdegree = maxdegree > table->order[d] ? maxdegree : table->order[d];
    }

    for(index_t i=0; i<n_evals; ++i) 
    {
        value_t biatx[table->ndim*(maxdegree+1)];

        // index_t lefties[table->ndim];
        index_t centers[table->ndim];
        // Get a random parameter within the splines 
        value_t par[table->ndim];

        // That's the way IceCube does that. But why search for tablecenters?
        for(index_t d=0; d<table->ndim; ++d) 
            par[d] = range[d] * normal_dist(engine) + table->knots[d][0];
    
        // We actually search for table centers. Who knows, why
        tablesearchcenters(table, par, centers);
         
        // For every dimension
        for(index_t d=0; d<table->ndim; ++d) 
        {   
            // Get the values
            bsplvb_2(
                table->knots[d], 
                table->order[d], 
                1, 
                par[d], 
                centers[d], 
                &(biatx[d*(maxdegree+1)]),
                table->nknots[d] );
           
        }
        // store y or something
        // y_array[i] = biatx[(table->ndim)*(maxdegree + 1) - 1];
        y_array[i] = ndsplineeval_core_2(centers, maxdegree, biatx); 
    }
}

/*
 * The original IceCube call.
 */ 

void CPU_BSPLV::eval_splines_vanilla(
    index_t n_evals,
    value_t * y_array) 
{
    std::mt19937 engine(42);
    std::uniform_real_distribution<value_t> normal_dist(0, 1);
    value_t range[table->ndim];

    for(index_t d=0; d<table->ndim; d++) 
    {
        range[d] = table->knots[d][table->nknots[d]-1] 
                   - table->knots[d][0];
    }

    for(index_t i=0; i<n_evals; ++i) 
    {
        index_t centers[table->ndim];
        // Get a random parameter within the splines 
        value_t par[table->ndim];
       
        for(index_t d=0; d<table->ndim; d++) 
            par[d] = range[d] * normal_dist(engine) + table->knots[d][0];
        
        // We actually search for table centers. Who knows, why
        tablesearchcenters(table, par, centers);

        y_array[i] = ndsplineeval(table, par, centers, 0);  
    }
}

/*
 * Evaluate the splines as in IceCube but use my multi dimensional overall 
 * evaluation.
 */
void CPU_BSPLV::eval_splines_simple(
    index_t n_evals,
    value_t * y_array) 
{
    std::mt19937 engine(42);
    std::uniform_real_distribution<value_t> normal_dist(0, 1);

    value_t range[table->ndim];
    index_t maxdegree = 0;
    for(index_t d=0; d<table->ndim; d++) 
    {
        range[d] = table->knots[d][table->nknots[d]-1] 
                   - table->knots[d][0];
        maxdegree = maxdegree > table->order[d] ? maxdegree : table->order[d];
    }

    for(index_t i=0; i<n_evals; ++i) 
    {
        value_t biatx[table->ndim*(maxdegree+1)];

        // index_t lefties[table->ndim];
        index_t centers[table->ndim];
        // Get a random parameter within the splines 
        value_t par[table->ndim];
       
    
        for(index_t d=0; d<table->ndim; d++) 
            par[d] = range[d] * normal_dist(engine) + table->knots[d][0];

        tablesearchcenters(table, par, centers);
         
        // For every dimension
        for(index_t d=0; d<table->ndim; d++) 
        {   
            bsplvb_simple(
                table->knots[d],    
                table->nknots[d],   
                par[d],            
                centers[d],              
                table->order[d]+1, 
                &(biatx[d*(maxdegree+1)]) );          
           
        }
        // store y or something
        // y_array[i] = biatx[(table->ndim)*(maxdegree + 1) - 1];
        y_array[i] = ndsplineeval_core_2(centers, maxdegree, biatx);
    }
}

void CPU_BSPLV::load_table(
    char * filename)
{
    if(table != nullptr)
        free_table();
    table = new splinetable();
    std::cout << "loading from " << filename << "\n";
    start("load_from_file");
    // Load a spline file from a file given by the command line
    if(!load_splines(table, filename))
    {
        std::cout << "Failed to load. Abort the mission!\n";
        std::cout << "I repeat: Abort the mission!\n";
    }
    stop();
}

void CPU_BSPLV::free_table()
{
    for(index_t i=0; i<table->ndim; ++i) 
        free(table->knots[i] - table->order[i]);
    
    free(table->knots);
    free(table->order);
    free(table->nknots);
    free(table->periods);
    free(table->coefficients);
    free(table->naxes);
    free(table->strides);
    free(table->extents[0]);
    free(table->extents);
    for(index_t i=0; i<table->naux; ++i) 
    {
        free(table->aux[i][1]);
        free(table->aux[i][0]);
        free(table->aux[i]);
    }
    free(table->aux);
    delete table;
}