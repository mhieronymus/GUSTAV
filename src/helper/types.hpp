#ifndef _TYPES_H
#define _TYPES_H 

typedef double value_t;
typedef int index_t;
// using value_t = double;
// using index_t = int;

// safe division
#define SDIV(x,y)(((x)+(y)-1)/(y))

#endif