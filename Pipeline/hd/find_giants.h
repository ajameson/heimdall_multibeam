/***************************************************************************
 *
 *   Copyright (C) 2012 by Ben Barsdell and Andrew Jameson
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#pragma once

#include <vector>

// TODO: Any way to avoid including this here?
#include <thrust/device_vector.h>

#include <boost/shared_ptr.hpp>

#include "hd/types.h"
#include "hd/error.h"

struct GiantFinder_impl;

struct GiantFinder {
	GiantFinder();
  hd_error exec(const hd_float* d_data, hd_size count, hd_size nbeams,
                hd_float thresh, hd_size merge_dist, hd_size rel_filter_width, 
                hd_size rel_dm_delay,
                thrust::device_vector<hd_float>& d_giant_peaks,
                thrust::device_vector<hd_size>&  d_giant_inds,
                thrust::device_vector<hd_size>&  d_giant_begins,
                thrust::device_vector<hd_size>&  d_giant_ends);
private:
	boost::shared_ptr<GiantFinder_impl> m_impl;
};
