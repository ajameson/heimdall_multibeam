/***************************************************************************
 *
 *   Copyright (C) 2012 by Ben Barsdell and Andrew Jameson
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#pragma once

#include "hd/types.h"
#include "hd/error.h"

#include <boost/shared_ptr.hpp>

struct RemoveBaselinePlan_impl;

struct RemoveBaselinePlan {
	RemoveBaselinePlan();
	hd_error exec(hd_float* d_data,
	              hd_size   beam_stride,
	              hd_size   count,
	              hd_size   smooth_radius,
                hd_size   nbeams);
private:
	boost::shared_ptr<RemoveBaselinePlan_impl> m_impl;
};
