#ifndef _LOAD_H
#define _LOAD_H
/* The code here is based (or mostly taken) from I3PhotoSplineTable.cxx
 * SetupTable
 * 
 */ 
#include <iostream>
#include "types.hpp"
#include "photospline/geo_type.h"
#include "photospline/splinetable.h"

// int readsplinefitstable(const char *path, struct splinetable *table);
// int splinetable_read_key(const struct splinetable *table, 
//     value_t type, const char *key, void *result);

bool load_splines(splinetable *& tablestruct_, char * filename);

#endif 

// weak scaling linear (problem size per processor is fixed)