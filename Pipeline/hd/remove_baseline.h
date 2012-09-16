
#pragma once

#include "hd/types.h"
#include "hd/error.h"

#include <boost/shared_ptr.hpp>

struct RemoveBaselinePlan_impl;

struct RemoveBaselinePlan {
	RemoveBaselinePlan();
	hd_error exec(hd_float* d_data,
	              hd_size   count,
	              hd_size   smooth_radius);
private:
	boost::shared_ptr<RemoveBaselinePlan_impl> m_impl;
};
