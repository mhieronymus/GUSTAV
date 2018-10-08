#include "cpu_bsplv.hpp"

// Binary search to find the leftmost spline with support.
// Uses left as first guess
void find_left_knot(
    value_t * knots,
    index_t n_knts,
    index_t & left,
    value_t x) 
{
    if(knots[left] < x && knots[left+1] > x) return;
    do 
    {
        if(knots[left] > x) {
            left = left >> 1;
        } else 
        {
            left = left + (left << 1);
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
    do 
    {
        index_t jp1 = j+1;
        value_t * delta_r = new value_t[jmax];
        delta_r[j] = knots[left + j] - x;
        value_t * delta_l = new value_t[jmax];
        delta_l[j] = x - knots[left+1-j];
        value_t saved = 0;
        for(index_t i=0; i<j; ++i) 
        {
            value_t term = biatx[i] / (delta_r[i] + delta_l[jp1-i]);
            biatx[i] = saved + delta_r[i] * term;
            saved = delta_l[jp1-i] * term;
        }
        biatx[jp1] = saved;
        j = jp1;
    } while(j < jhigh);
    return biatx[j];
}


// Create n_evals times a random array of size n_dims and evaluate the spline
void cpu_eval_splines(
    value_t * spline_cffs, 
    index_t n_cffs, 
    value_t * spline_knts, 
    index_t * n_knts, 
    index_t * spline_dgrs,
    index_t n_dgrs, 
    index_t n_evals,
    index_t n_dims,
    index_t * splines_per_dim) 
{

    std::mt19937 engine(42);
    std::uniform_real_distribution<value_t> normal_dist(0, 1);

    // Store the range of each dimension given by the difference of its first
    // and last knot over all splines in that dimension
    // Also store the lower bounds
    value_t * range = new value_t[n_dims];
    value_t * lwr_bnds = new value_t[n_dims];
    index_t knots_offset = 0;
    index_t dgrs_offset = 0;
    index_t * last_knot_idx = new index_t[n_dims];
    for(index_t d=0; d<n_dims; d++) 
    {   
        index_t curr_last_knot_idx = 0;
        for(index_t i=0; i<splines_per_dim[d]; i++) 
            curr_last_knot_idx += spline_dgrs[i + dgrs_offset] + 2;

        range[d] =  spline_knts[knots_offset + curr_last_knot_idx] 
                  - spline_knts[knots_offset];
        lwr_bnds[d] = spline_knts[knots_offset];
        knots_offset += curr_last_knot_idx;
        dgrs_offset += splines_per_dim[d];
        last_knot_idx[d] = curr_last_knot_idx;
    }

    for(index_t i=0; i<n_evals; ++i) 
    {
        // Get a random parameter within the splines 
        value_t * par = new value_t[n_dims];
        for(index_t d=0; d<n_dims; d++) 
            par[d] = range[d] * normal_dist(engine) + lwr_bnds[d];
         
        value_t y = 0;
        knots_offset = 0;
        dgrs_offset = 0;
        // For every dimension
        for(index_t d=0; d<n_dims; d++) 
        {   
            index_t n_knts = last_knot_idx[d];
            // For every spline s that has support in par...
            // If we assume knots of equal distance for simplicity, we could
            // simply the search for the starting t_i aka left index. 
            // Or we do a binary search starting at a knot that might be near.
            // It is exact if the grid has equal distances everywhere.
            index_t left = SDIV( par[d] - lwr_bnds[d], 
                                 range[d]/last_knot_idx[d] );
            find_left_knot(&(spline_knts[knots_offset]), n_knts, left, par[d]);
            
            value_t * biatx = new value_t[spline_dgrs[dgrs_offset]+1];
            y += bsplvb_2(&(spline_knts[knots_offset]), 
                        spline_dgrs[dgrs_offset]+1, spline_dgrs[dgrs_offset],
                        par[d], left, biatx);

            knots_offset = last_knot_idx[d];
            dgrs_offset += splines_per_dim[d];
        }
        // store y or something
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


