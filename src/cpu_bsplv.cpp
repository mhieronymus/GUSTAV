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
value_t bsplvb_2(
    value_t * knots,
    index_t jhigh,
    index_t index,
    value_t x,
    index_t left,
    value_t * biatx,
    index_t jmax = 20)
{
    index_t j = 0;
    if(index != 1) { 
        j = 0;
        biatx[j] = 1;
        // Check if no more columns need to be evaluated.
        if(j >= jhigh) return biatx[j];
    }
    value_t * delta_r = new value_t[jmax];
    value_t * delta_l = new value_t[jmax];
    do 
    {
        index_t jp1 = j+1;
        
        delta_r[j] = knots[left + j + 1] - x;
        delta_l[j] = x - knots[left-j];
        value_t saved = 0;
        for(index_t i=0; i<j; ++i) 
        {
            value_t term = biatx[i] / (delta_r[i] + delta_l[jp1-i-1]);
            biatx[i] = saved + delta_r[i] * term;
            saved = delta_l[jp1-i-1] * term;
        }
        biatx[jp1] = saved;
        j = jp1;
        

    } while(j < jhigh);
    delete[] delta_r;
    delete[] delta_l;
    return biatx[j];
}


// Create n_evals times a random array of size n_dims and evaluate the spline
void cpu_eval_splines(
    splinetable * table,
    index_t n_evals) 
{

    std::mt19937 engine(42);
    std::uniform_real_distribution<value_t> normal_dist(0, 1);

    // // Store the range of each dimension given by the difference of its first
    // // and last knot over all splines in that dimension
    // // Also store the lower bounds
    // value_t * range = new value_t[n_dims];
    // value_t * lwr_bnds = new value_t[n_dims];
    // index_t knots_offset = 0;
    // index_t dgrs_offset = 0;
    // index_t * last_knot_idx = new index_t[n_dims];
    // for(index_t d=0; d<n_dims; d++) 
    // {   
    //     index_t curr_last_knot_idx = 0;
    //     for(index_t i=0; i<splines_per_dim[d]; i++) 
    //         curr_last_knot_idx += spline_dgrs[i + dgrs_offset] + 2;

    //     range[d] =  spline_knts[knots_offset + curr_last_knot_idx] 
    //               - spline_knts[knots_offset];
    //     lwr_bnds[d] = spline_knts[knots_offset];
    //     knots_offset += curr_last_knot_idx;
    //     dgrs_offset += splines_per_dim[d];
    //     last_knot_idx[d] = curr_last_knot_idx;
    // }
    value_t * range = new value_t[table->ndim];
    for(index_t d=0; d<table->ndim; d++) 
        range[d] = table->knots[d][table->nknots[d]-1] 
                   - table->knots[d][0];
    for(index_t i=0; i<n_evals; ++i) 
    {
        
        // Get a random parameter within the splines 
        value_t * par = new value_t[table->ndim];
        for(index_t d=0; d<table->ndim; d++) 
            par[d] = range[d] * normal_dist(engine) + table->knots[d][0];
         
        value_t y = 0;
        value_t y2 = 0;
        value_t y3 = 0;
        // For every dimension
        for(index_t d=0; d<table->ndim; d++) 
        {   
            // For every spline s that has support in par...
            // We do a binary search starting at a knot that might be near.
            // It is exact if the grid has equal distances everywhere.
            // std::cout << "Using par[d] = " << par[d] << "\n";
            // std::cout << table->knots[d][table->nknots[d]-1] << " - " << table->knots[d][0] << " = " << range[d] << "\n";
            index_t left = SDIV( 
                (par[d] - table->knots[d][0])*table->nknots[d], 
                range[d]);
            
            find_left_knot(
                table->knots[d], 
                table->nknots[d], 
                left, 
                par[d]);
            
            value_t * biatx = new value_t[table->order[d]+1];
    
            y = bsplvb_2(
                table->knots[d], 
                table->order[d], 
                table->order[d], 
                par[d], 
                left, 
                biatx);
            // std::cout << "\nbiatx\n";
            // for(index_t k=0; k<table->order[d]+1; ++k)
            //     std::cout << biatx[k] << ", ";
            // std::cout << "\n";
            value_t * biatx2 = new value_t[table->order[d]+1]; // max degree?
            // index_t maxdegree = maxorder(table->order, table->ndim)+1;
            // basically just to get the max degree of all dims
            bsplvb_simple(
                table->knots[d],    // double* knots
                table->nknots[d],   // unsigned nknots
                par[d],             // double x
                left,               // int left, centers[n]
                table->order[d]+1,  // int degree
                biatx2);            // float* biatx
            y2 = biatx2[table->order[d]];

            value_t * biatx3 = new value_t[table->order[d]+1];
            value_t * delta_l = new value_t[20];
            value_t * delta_r = new value_t[20];
            bsplvb(
                table->knots[d],
                par[d],
                left,
                table->order[d], 
                table->order[d],
                biatx3,
                delta_l,
                delta_r);
            y3 = biatx3[table->order[d]];
            if(DEBUG)
                std::cout << y << " vs " << y2 << " vs " << y3 << ", ";
            delete[] biatx;
            delete[] biatx2;
            delete[] biatx3;
            delete[] delta_l;
            delete[] delta_r;
        }
        delete[] par;
        // store y or something
        if(DEBUG)
            std::cout << "\n";
    }
    delete[] range;
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


