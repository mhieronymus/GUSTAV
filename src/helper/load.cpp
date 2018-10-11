#include "load.hpp"

bool load_splines(
    splinetable *& tablestruct_,
    char * filename)
{
    index_t errorvalue_;
    if(DEBUG) {
        std::cout << "Loading from " << filename << "\n";
    }
	// splinetable * tablestruct_ = new splinetable();

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

  
        std::cout << "Loading successful\n";
		return true;
	} else {
        if(DEBUG) {
            std::cout << "Loading failed\n";
        }
		return false;
	}
}