#ifndef _BSPLINE_H
#define _BSPLINE_H

#include "helper/photospline/splinetable.h"
#include "helper/types.hpp"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Compute the value of the ith nth-order basis spline of a set
 * defined by knots at the poindex_t x.
 */

value_t bspline(const value_t *knots, value_t x, index_t i, index_t n);
value_t bspline_deriv(const value_t *knots, value_t x, index_t i, index_t n, unsigned order);

/*
 * A brain-dead reimplementation of de Boor's BSPLVB, which generates
 * the values of the non-zero B-splines at x from the bottom up without
 * unnecessarily recalculating terms. 
 * 
 * NB: for bsplvb_simple(), bspline_nonzero(), and bspline_deriv_nonzero(),
 * `left' must be the index of the nearest fully-supported knot
 * span for splines of order n, i.e. n <= left <= nknots-n-2. For bsplvb(),
 * `left' must be the index of the nonzero 0-th order spline, i.e.
 * knots[left] <= x < knots[left+1].
 *
 * See Chapter X in: 
 * 
 * Carl de Boor. A Practical Guide to Splines, volume 27 of Applied
 *     Mathematical Sciences. Springer-Verlag, 1978.
 */

void bsplvb_simple(const value_t *knots, const unsigned nknots,
    value_t x, index_t left, index_t jhigh, value_t * biatx);
void bsplvb(const value_t *knots, const value_t x, const index_t left, const index_t jlow,
    const index_t jhigh, value_t *biatx,
    value_t *delta_l, value_t *delta_r);
void bspline_nonzero(const value_t *knots, const unsigned nknots,
    const value_t x, index_t left, const index_t n, value_t *values, value_t *derivs);
void bspline_deriv_nonzero(const value_t *knots, const unsigned nknots,
    const value_t x, const index_t left, const index_t n, value_t *biatx);

/*
 * Evaluates the results of a full spline basis given a set of knots,
 * a position, an order, and a central spline for the position (or -1).
 * The central spline should be the index of the 0th order basis spline
 * that is non-zero at the position x.
 */

value_t splineeval(const value_t *knots, const value_t *weights, index_t nknots, value_t x,
    index_t order, index_t center);

/*
 * Spline table based hypersurface evaluation. ndsplineeval() takes a spline
 * coefficient table, a vector at which to evaluate the surface, and a vector
 * indicating the evaluation centers, as for splineeval().
 *
 * tablesearchcenters() provides a method to acquire a centers vector
 * for ndsplineeval() using a binary search. Depending on how the table
 * was produced, a more efficient method may be available.
 */

index_t tablesearchcenters(const struct splinetable *table, const value_t *x, index_t *centers);

value_t ndsplineeval(const struct splinetable *table, const value_t *x, 
    const index_t *centers, index_t derivatives);
// value_t ndsplineeval_linalg(const struct splinetable *table, const value_t *x, 
//     const index_t *centers, index_t derivatives);


// slicing function to get 1-d slice with certain normalizaton properties
void
ndsplineeval_slice_coeffs(const struct splinetable *table, const value_t *x, const index_t *centers, value_t *results,index_t slice_dimension, index_t derivative, index_t area_norm);

/*
* Evaluates the spline surface, optionally differentiated to the given order
* in any dimension.
*/

value_t ndsplineeval_deriv(const struct splinetable *table, const value_t *x, 
    const index_t *centers, const unsigned *derivatives);

/* Evaluate a spline surface and all its derivatives at x */

void ndsplineeval_gradient(const struct splinetable *table, const value_t *x,
    const index_t *centers, value_t *evaluates);

/*
 * Convolve a table with the spline defined on a set of knots along a given 
 * dimension and store the spline expansion of the convolved surface in the
 * table. This will raise the order of the splines in the given dimension by
 * (n_knots - 1), i.e. convolving with an order-0 spline (a box function, 
 * defined on two knots) will raise the order of the spline surface by 1.
 */
index_t splinetable_convolve(struct splinetable *table, const index_t dim,
    const value_t *knots, size_t n_knots);

#ifdef __cplusplus
}
#endif

#endif /* _BSPLINE_H */
