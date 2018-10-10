#include "cpu_bsplv.hpp"

#include <iostream>
// Binary search to find the leftmost spline with support.
// Uses left as first guess
void find_left_knot(
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
value_t ndsplineeval_core_2(
    splinetable * table, 
    index_t * centers, 
    index_t maxdegree, 
    value_t * biatx) 
{
    value_t result = 0;
    value_t basis_tree[table->ndim];
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
        // for(index_t i=0; __builtin_expect(i < table->order[table->ndim-1]+1, 
            // 1); ++i)
        for(index_t i=0; i<table->order[table->ndim-1]+1; ++i)
        {
            value_t tmp = basis_tree[table->ndim-1] ;
            tmp = biatx[(table->ndim-1)*(maxdegree+1) + i] ;
            tmp = table->coefficients[tablepos+i]; // invalid read
            result += basis_tree[table->ndim-1] 
                * biatx[(table->ndim-1)*(maxdegree+1) + i] 
                * table->coefficients[tablepos+i];
        }
        // if(__builtin_expect(++n == nchunks, false)) break;
        if(++n == nchunks) break;

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
        // for(index_t j=i; __builtin_expect(j < table->ndim-1, 1); ++j)
        for(index_t j=i; j < table->ndim-1; ++j)
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
void bsplvb_2(
    value_t * knots,
    index_t jhigh,
    index_t index,
    value_t x,
    index_t left,
    value_t * biatx,
    index_t nknots,
    index_t jmax = 20)
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

value_t wiki_bsplines(
    value_t * knots, value_t x, index_t left, index_t order)
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
void cpu_eval_splines(
    splinetable * table,
    index_t n_evals) 
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
        value_t biatx2[table->ndim*(maxdegree+1)];

        // index_t lefties[table->ndim];
        index_t centers[table->ndim];
        // Get a random parameter within the splines 
        value_t par[table->ndim];
       
        // for(index_t d=0; d<table->ndim; d++) 
        // {
        //     value_t exact;
        //     do {
        //         par[d] = range[d] * normal_dist(engine) + table->knots[d][0];
        //         index_t left = SDIV( 
        //         (par[d] - table->knots[d][0])*table->nknots[d], 
        //         range[d]);
            
        //         find_left_knot(
        //             table->knots[d], 
        //             table->nknots[d], 
        //             left, 
        //             par[d]);
        //         centers[d] = left;

        //         exact = (table->knots[d], par[d], centers[d],
        //             table->order[d] );
        //     }while(exact < 0.25);
        //     if(!DEBUG)
        //         std::cout << "\n" << par[d] << "\n" << std::flush;

        // }
        // That's the way IceCube does that. But why search for tablecenters?
        for(index_t d=0; d<table->ndim; d++) 
        {   
            par[d] = range[d] * normal_dist(engine) + table->knots[d][0];

            // index_t left = SDIV( 
            //     (par[d] - table->knots[d][0])*table->nknots[d], 
            //     range[d]);
            
            //     find_left_knot(
            //         table->knots[d], 
            //         table->nknots[d], 
            //         left, 
            //         par[d]);
            //     centers[d] = left;
        }
        // We actually search for table centers. Who knows, why
        tablesearchcenters(table, par, centers);
         
        // For every dimension
        for(index_t d=0; d<table->ndim; d++) 
        {   
            // par[d] = range[d] * normal_dist(engine) + table->knots[d][0];
            // // For every spline s that has support in par...
            // // We do a binary search starting at a knot that might be near.
            // // It is exact if the grid has equal distances everywhere.
            // // std::cout << "Using par[d] = " << par[d] << "\n";
            // // std::cout << table->knots[d][table->nknots[d]-1] << " - " << table->knots[d][0] << " = " << range[d] << "\n";
            // index_t left = SDIV( 
            //     (par[d] - table->knots[d][0])*table->nknots[d], 
            //     range[d]);
            
            // find_left_knot(
            //     table->knots[d], 
            //     table->nknots[d], 
            //     left, 
            //     par[d]);
            // lefties[d] = left;

            // Get the values
            bsplvb_2(
                table->knots[d], 
                table->order[d], 
                1, 
                par[d], 
                centers[d], 
                &(biatx[d*(maxdegree+1)]),
                table->nknots[d] );
       
            // index_t maxdegree = maxorder(table->order, table->ndim)+1;
            // basically just to get the max degree of all dims
            bsplvb_simple(
                table->knots[d],    // double* knots
                table->nknots[d],   // unsigned nknots
                par[d],             // double x
                centers[d],               // int left, centers[n]
                table->order[d]+1,  // int degree
                &(biatx2[d*(maxdegree+1)]) );            // float* biatx
           
        }
        // store y or something
        value_t y = ndsplineeval_core_2(table, centers, maxdegree, biatx);
        value_t y2 = ndsplineeval_core_2(table, centers, maxdegree, biatx2);
        // value_t y = -1;
        // value_t y2 = -1;
        // value_t y3 = -1;

        value_t y1 = ndsplineeval(table, par, centers, 0);

        if(!DEBUG) {

        
            std::cout << "\nbiatx:\n";
            for(index_t d=0; d<table->ndim; ++d)
            {
                std::cout << "Exact ";
                value_t exact = wiki_bsplines(table->knots[d], par[d], centers[d],
                    table->order[d] );
                std::cout << exact << "\n";
                for(index_t j=0; j<maxdegree + 1; ++j) 
                {
                    std::cout << "("  << biatx[d*(maxdegree+1)+j]
                            << ", " << biatx2[d*(maxdegree+1)+j] << "), ";
                }
                std::cout << "\n";
            }
            std::cout << "(my_value, using bsplvb, using ndsplineeval)\n";
            
            std::cout << "(" << y << ", " << y2 << ", " << y1 << ")\n";
        }
        
    }
}

/// From IceCube

// double
// ndsplineeval(const struct splinetable *table, const double *x, const int *centers,
//     int derivatives)
// {
// 	assert(table->ndim>0);
// 	int n;
// 	int maxdegree = maxorder(table->order, table->ndim) + 1; 
// 	float localbasis[table->ndim][maxdegree];
	
// 	for (n = 0; n < table->ndim; n++) {
// 		if (derivatives & (1 << n)) {
// 			bspline_deriv_nonzero(table->knots[n], 
// 			    table->nknots[n], x[n], centers[n],
// 			    table->order[n], localbasis[n]);
// 		} else {
// 			bsplvb_simple(table->knots[n], table->nknots[n],
// 			    x[n], centers[n], table->order[n] + 1,
// 			    localbasis[n]);
// 		}
// 	}

// 	return ndsplineeval_core(table, centers, maxdegree, localbasis);
// }

// static int
// maxorder(int *order, int ndim)
// {
// 	int i, max = 0;
	
// 	for (i = 0; i < ndim; i++)
// 		if (order[i] > max)
// 			max = order[i];
	
// 	return (max);
// }
   
// /*
//  * A brain-dead reimplementation of de Boor's BSPLVB, which generates
//  * the values of the non-zero B-splines at x from the bottom up without
//  * unnecessarily recalculating terms. 
//  * 
//  * NB: for bsplvb_simple(), bspline_nonzero(), and bspline_deriv_nonzero(),
//  * `left' must be the index of the nearest fully-supported knot
//  * span for splines of order n, i.e. n <= left <= nknots-n-2. For bsplvb(),
//  * `left' must be the index of the nonzero 0-th order spline, i.e.
//  * knots[left] <= x < knots[left+1].
//  *
//  * See Chapter X in: 
//  * 
//  * Carl de Boor. A Practical Guide to Splines, volume 27 of Applied
//  *     Mathematical Sciences. Springer-Verlag, 1978.
//  */

// void
// bsplvb_simple(const double *knots, const unsigned nknots,
//     double x, int left, int degree, float *restrict biatx)
// {
// 	assert(degree>0);
// 	int i, j;
// 	double saved, term;
// 	double delta_l[degree], delta_r[degree];
	
// 	biatx[0] = 1.0;
	
// 	/*
// 	 * Handle the (rare) cases where x is outside the full
// 	 * support of the spline surface.
// 	 */
// 	if (left == degree-1)
// 		while (left >= 0 && x < knots[left])
// 			left--;
// 	else if (left == nknots-degree-1)
// 		while (left < nknots-1 && x > knots[left+1])
// 			left++;	
	
// 	/* 
// 	 * NB: if left < degree-1 or left > nknots-degree-1,
// 	 * the following loop will dereference addresses ouside
// 	 * of knots[0:nknots]. While terms involving invalid knot
// 	 * indices will be discarded, it is important that `knots'
// 	 * have (maxdegree-1)*sizeof(double) bytes of padding
// 	 * before and after its valid range to prevent segfaults
// 	 * (see parsefitstable()).
// 	 */
// 	for (j = 0; j < degree-1; j++) {
// 		delta_r[j] = knots[left+j+1] - x;
// 		delta_l[j] = x - knots[left-j];
		
// 		saved = 0.0;
		
// 		for (i = 0; i < j+1; i++) {
// 			term = biatx[i] / (delta_r[i] + delta_l[j-i]);
// 			biatx[i] = saved + delta_r[i]*term;
// 			saved = delta_l[j-i]*term;
// 		}
		
// 		biatx[j+1] = saved;
// 	}
	
// 	/* 
// 	 * If left < (spline order), only the first (left+1)
// 	 * splines are valid; the remainder are utter nonsense.
// 	 */
// 	if ((i = degree-1-left) > 0) {
// 		for (j = 0; j < left+1; j++)
// 			biatx[j] = biatx[j+i]; /* Move valid splines over. */
// 		for ( ; j < degree; j++)
// 			biatx[j] = 0.0; /* The rest are zero by construction. */
// 	} else if ((i = left+degree+1-nknots) > 0) {
// 		for (j = degree-1; j > i-1; j--)
// 			biatx[j] = biatx[j-i];
// 		for ( ; j >= 0; j--)
// 			biatx[j] = 0.0;
// 	}
// }


// /*
//  * The N-Dimensional tensor product basis version of splineeval.
//  * Evaluates the results of a full spline basis given a set of knots,
//  * a position, an order, and a central spline for the position (or -1).
//  * The central spline should be the index of the 0th order basis spline
//  * that is non-zero at the position x.
//  *
//  * x is the vector at which we will evaluate the space
//  */

// static double
// ndsplineeval_core(const struct splinetable *table, const int *centers, int maxdegree,
//     float localbasis[table->ndim][maxdegree])
// {
// 	assert(table->ndim>0);
// 	int i, j, n, tablepos;
// 	float result;
// 	float basis_tree[table->ndim+1];
// 	int nchunks;
// 	int decomposedposition[table->ndim];

// 	tablepos = 0;
// 	for (n = 0; n < table->ndim; n++) {
// 		decomposedposition[n] = 0;
// 		tablepos += (centers[n] - table->order[n])*table->strides[n];
// 	}

// 	basis_tree[0] = 1;
// 	for (n = 0; n < table->ndim; n++)
// 		basis_tree[n+1] = basis_tree[n]*localbasis[n][0];
// 	nchunks = 1;
// 	for (n = 0; n < table->ndim - 1; n++)
// 		nchunks *= (table->order[n] + 1);

// 	result = 0;
// 	n = 0;
// 	while (1) {
// 		for (i = 0; __builtin_expect(i < table->order[table->ndim-1] +
// 		    1, 1); i++) {

// 			result += basis_tree[table->ndim-1]*
// 			    localbasis[table->ndim-1][i]*
// 			    table->coefficients[tablepos + i];
// 		}

// 		if (__builtin_expect(++n == nchunks, 0))
// 			break;

// 		tablepos += table->strides[table->ndim-2];
// 		decomposedposition[table->ndim-2]++;

// 		/* Carry to higher dimensions */
// 		for (i = table->ndim-2;
// 		    decomposedposition[i] > table->order[i]; i--) {
// 			decomposedposition[i-1]++;
// 			tablepos += (table->strides[i-1]
// 			    - decomposedposition[i]*table->strides[i]);
// 			decomposedposition[i] = 0;
// 		}
// 		for (j = i; __builtin_expect(j < table->ndim-1, 1); j++)
// 			basis_tree[j+1] = basis_tree[j]*
// 			    localbasis[j][decomposedposition[j]];
// 	}

// 	return result;
// }


