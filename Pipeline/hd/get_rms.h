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

struct GetRMSPlan_impl;

struct GetRMSPlan {
	GetRMSPlan();
	hd_float exec(hd_float* d_data, hd_size count);

private:
	boost::shared_ptr<GetRMSPlan_impl> m_impl;
};

struct GetRMSPlanMB_impl;

struct GetRMSPlanMB {
  GetRMSPlanMB();
  hd_error exec_multibeam (hd_float* d_data, hd_float * d_rms, hd_size beam_stride, hd_size beam_count, hd_size nbeams);
  hd_float exec_multibeam (hd_float* d_data, hd_size beam_stride, hd_size beam_count, hd_size nbeams);

private:
  boost::shared_ptr<GetRMSPlanMB_impl> m_impl;
};


// Convenience functions for one-off calls
hd_float get_rms(hd_float* d_data, hd_size count);
hd_error normalise(hd_float* d_data, hd_size count);

// For Multibeam access
hd_error get_rms_multibeam(hd_float* d_data, hd_float* d_rms, hd_size beam_stride, hd_size beam_count, hd_size nbeams);
hd_error normalise_multibeam (hd_float* d_data, hd_float block_rms, hd_size count);
hd_error normalise_multibeam(hd_float* d_data, hd_float * d_rms, hd_size beam_stride, hd_size cur_nsamps, hd_size nbeams);
