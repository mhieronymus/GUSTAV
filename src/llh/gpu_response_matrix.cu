/* 
 * Every function that is needed to create a response matrix, except for 
 * evaluating the splines, is here.
 * It is basically like GetResponseMatrix()
 */

 // TODO: Implement PhotonicsSource in CUDA

 // TODO:
 // run that in parallel? Only possible by calling that multiple times
 // This version evaluates one llh per block.
 __device__
 value_t scale_light_yield_per_block(
     PhotonicsSource_d & source
     value_t raw_yield,
     value_t energy,
     index_t derivative) 
{
    if(energy < 0) energy = source.energy;
    if(energy < 0.01) energy = 0.01;

    // Cherenkov photons per meter in 300 - 600 nm range
    // assuming beta = 1 and a wavelength dependent index of refraction
    value_t light_factor = 32582.0;

    /* Compute effective track length */
    /* EM Shower: see AMANDA-IR/20020803 Table 2.1 */
    if(source.type == 1) 
    {
        light_factor *= 5.21;
        if(derivative == 0)
        {
            light_factor *= energy;
        } else if(derivative != 1)
        {
            light_factor = 0; // higher derivatives
        } 
    /* Hadronic Shower */
    } else if(source.type == 2)
    {
        light_factor *= 0.860*4.076;
        if(derivative == 0)
        {
            light_factor *= energy;
        } else if(derivative != 1)
        {
            light_factor = 0; // higher derivatives
        } 

    /* Monopole: like a muon, but more */
    } else if(source.type == 9  || source.type == 10 || 
        source.type == 11 || source.type == 12)
    {
        // The n_ice_phase=1.3195 is an I3Constant. This should be given
        // by the host
        value_t beta_n = 1.3195 * source.speed;
        value_t mp_amp = 1.0 - (1.0/powf(beta_n, 2));
        if(derivative==0)
        {
            light_factor *= (1.172 + 0.0324*log(energy));
        } else if(derivative==1)
        {
            light_factor *= 0.0324/energy;
        } else if(derivative==2)
        {
            light_factor *= -0-324/(energy*energy);
        } else 
        {
            light_factor = 0;
        }

    /* Sub-cherenkov monopole */
    } else if(source.type == 13)
    {
        light_factor = 0;

    /* Muon */
    } else 
    {
        if(derivative==0)
        {
            light_factor *= (1.172 + 0.0324*log(energy));
        } else if(derivative==1)
        {
            light_factor *= 0.0324/energy;
        } else if(derivative==2)
        {
            light_factor *= -0-324/(energy*energy);
        } else 
        {
            light_factor = 0;
        }
    }
    return raw_yield * light_factor;
}

// TODO:
// Add xOM, yOM, zOM, max_radius, n_group, c_vacuum, meanPEs, last_source,
// geo_time, emission_point_distance, source_type, source_length,
// emissionpoint_distance, raw_yield, 

// Based on SelectSource in I3PhotoSplineService.cxx
__device__
void select_source(
    value_t & meanPEs,
    value_t * gradient,
    value_t & emission_point_distance,
    value_t & geo_time,
    PhotonicsSource_d const & source,
    bool get_amp)
{
    // check extends 
    value_t ePDist = sqrtf( powf( (xOM_ - source.x), 2) 
                          + powf( (yOM_ - source.y), 2) 
                          + powf( (zOM_ - source.z), 2) );
    if(ePDist > max_radius_)
    {
        meanPEs = 0.0;
        emission_point_distance = ePDist;
        geo_time = emission_point_distance * n_group_/c_vacuum_;
        return;
    }
    // return cahced values for same DOM/source config
    if(meanPEs_ >= 0 && last_source_ == source && gradient == NULL)
    {
        meanPEs = meanPEs_;
        geo_time = geo_time_;
        emission_point_distance = emission_point_distance_;
        return;
    }

    calculate_photonics_input(false, xOM_, yOM_, zOM_, source);

    emission_point_distance = emission_point_distance_;
    source_type_ = source.type;
    source_length_ = source.length;

    value_t tablecoordinates[6], buffer[7];
    fill_table_coordinates(tablecoordinates, false);

    // Check if coordinates are within bounds
    // TODO: Move this function and give amp_table as argument
    bool supported = amp_table_->check_support(tablecoordinates);

    if(geo_type_ == POINTSOURCE)
    {
        // showers and other point-like sources
        // Different c_vacuum! From I3PhotonicsCommons
        geo_time = r_*n_group_/c_vacuum;
    } else if(geo_type_ == INFINITEMUON)
    {
        /**
		 * Calculate "geometric time", which has two options for
		 * infinite muons:
		 *  - time for muon to propagate to emission point + time for
		 *  light to reach OM from there for muons (length >= 0)
		 *  - cascade-like direct time for lightsaber (length = -1)
         */
         if(source.length == -1)
            geo_time = r_*n_group_/c_vacuum;
        else
            geo_time = (emission_point_offset_ + emission_point_distance_)
                / c_vacuum_;
    } else 
    {
        // Log_fatal -> terminate
    }
    if(!supported) 
    {
        raw_yield_ = -1;
    } else if(!get_amp)
    {
        raw_yield = 0;
    } else if(gradient == NULL)
    {
        // TODO: Check how that can be put outside of the current flow.
        if(amp_table->Eval(tablecoordinates, &raw_yield_) == 0)
            raw_yield = exp(raw_yield_);
        else 
            raw_yield_ = -1;
    }
}