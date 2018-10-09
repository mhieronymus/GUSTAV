#ifndef _TYPES_H
#define _TYPES_H 

#include <stdbool.h> 

typedef double value_t;
typedef int index_t;
// using value_t = double;
// using index_t = int;

// safe division
#define SDIV(x,y)(((x)+(y)-1)/(y))
extern const bool DEBUG;

struct splinetable {
	index_t ndim;
	index_t *order;     // Order the splines (ndim many entries)

	value_t **knots;    // Knots in each dimension
	long *nknots;       // Number of knots in each dimension

	value_t **extents;

	value_t *periods;

	value_t *coefficients;  // All coefficients for each spline
	long *naxes;            // Number of splines for each dimension
	unsigned long *strides; // stride to each dimension

	index_t naux;
	char ***aux;
};

#endif