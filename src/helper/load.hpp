#ifndef _LOAD_H
#define _LOAD_H
/* The code here is based (or mostly taken) from I3PhotoSplineTable.cxx
 * SetupTable
 * 
 */ 
#include "types.hpp"
#include "photospline/geo_type.h"
#include "photospline/splinetable.h"
// struct splinetable {
// 	index_t ndim;
// 	index_t *order;

// 	value_t **knots;
// 	index_t *nknots;

// 	value_t **extents;

// 	value_t *periods;

// 	value_t *coefficients;
// 	index_t *naxes;
// 	index_t *strides;

// 	index_t naux;
// 	char ***aux;
// };

// int readsplinefitstable(const char *path, struct splinetable *table);
// int splinetable_read_key(const struct splinetable *table, 
//     value_t type, const char *key, void *result);

bool load_splines(
    value_t * spline_cffs,
    index_t n_cffs,
    value_t * spline_knts,
    index_t * n_knts,
    index_t * spline_dgrs,
    index_t n_dgrs,
    index_t n_dims,
    char * filename)
{
    index_t errorvalue_;

	splinetable * tablestruct_ = new splinetable();

    geo_type geometry_;
    int geotype_;

	if (readsplinefitstable(filename, &*tablestruct_) == 0) {
		long geo, geotype, par, err;
		double nGroupTable;

		err = splinetable_read_key(tablestruct_, SPLINETABLE_index_t,
		    "GEOMETRY", &geo);
		if (err)
			geo = CYLINDRICAL;
		geometry_ = (geo_type)geo;

		err = splinetable_read_key(tablestruct_, SPLINETABLE_index_t,
		    "GEOTYPE", &geotype);
		if (err) {
			geotype = 0;
		}
		geotype_ = geotype;

		err = splinetable_read_key(tablestruct_,
		    SPLINETABLE_value_t, "NGROUP", &nGroupTable);
        // See I3PhotonicsServiceCommons.h
        value_t nGroup = 1.35634;
		if (err || nGroupTable == -1) {
			nGroupTable = nGroup;
		}
		value_t nGroupTable_ = nGroupTable;

		err = splinetable_read_key(tablestruct_, SPLINETABLE_index_t,
		    "PARITY", &par);
		if (err) {
			// The first L1 tables were made with positive source
			// angles, but the first level2 tables should have
			// EVEN parity
			if (geotype_ == 0)
				par = ODD;
			else
				par = EVEN;
		}
		parity_type parity_ = (parity_type)par;

        spline_cffs     = tablestruct_->coefficients; 
        n_cffs          = 0;
        for(index_t i=0; i<tablestruct_->ndim; ++i) 
            n_cffs += tablestruct_->order[i]+1;
        // n_cffs          = tablestruct_->ndim * (tablestruct_->order+2);
        n_knts          = (index_t*) tablestruct_->nknots;
        index_t n_total_knts = 0;
        for(index_t i=0; i<tablestruct_->ndim; ++i) 
            n_total_knts += n_knts[i];
        // Flatten the array
        spline_knts = new value_t[n_total_knts];
        index_t offset = 0;
        for(index_t i=0; i<tablestruct_->ndim; ++i) 
        {
            for(index_t j=0; j<n_knts[i]; ++j)
            {
                spline_knts[j + offset] = tablestruct_->knots[i][j];
            }
            offset += n_knts[i];
        }
        // spline_knts     = tablestruct_->knots;
        
        spline_dgrs     = tablestruct_->order;
        n_dgrs          = tablestruct_->ndim;
        n_dims          = tablestruct_->ndim;
		return true;
	} else {
		delete[] tablestruct_; 
		return false;
	}
}

#endif 

// weak scaling linear (problem size per processor is fixed)