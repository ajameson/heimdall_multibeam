/***************************************************************************
 *
 *   Copyright (C) 2012 by Ben Barsdell and Andrew Jameson
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#pragma once

#include "hd/error.h"
#include "hd/types.h"

hd_error median_filter3(const hd_float* d_in,
                        hd_size         count,
                        hd_float*       d_out);

hd_error median_filter5(const hd_float* d_in,
                        hd_size         count,
                        hd_float*       d_out);

hd_error median_scrunch3(const hd_float* d_in,
                         hd_size         count,
                         hd_float*       d_out);

hd_error median_scrunch5(const hd_float* d_in,
                         hd_size         count,
                         hd_float*       d_out);

// Note: This can operate 'in-place'
hd_error mean_filter2(const hd_float* d_in,
                      hd_size         count,
                      hd_float*       d_out);

hd_error linear_stretch(const hd_float* d_in,
                        hd_size         in_count,
                        hd_float*       d_out,
                        hd_size         out_count);

// Median-scrunches the corresponding elements from a collection of arrays
// Note: This cannot (currently) handle count not being a multiple of 3
hd_error median_scrunch3_array(const hd_float* d_in,
                               hd_size         array_size,
                               hd_size         count,
                               hd_float*       d_out);

// Median-scrunches the corresponding elements from a collection of arrays [TS]
// Note: This cannot (currently) handle count not being a multiple of 5
hd_error median_scrunch5_array(const hd_float* d_in,
                               hd_size         array_size,
                               hd_size         count,
                               hd_float*       d_out);
// Median-scrunches the corresponding elements from a collection of beams [ST] 
// Note: This cannot (currently) handle count not being a multiple of 5
hd_error median_scrunch5_beam(const hd_float* d_in,
                              hd_size         nbeam,
                              hd_size         count,
                              hd_float*       d_out);


// Mean-scrunches the corresponding elements from a collection of arrays
// Note: This cannot (currently) handle count not being a multiple of 2
hd_error mean_scrunch2_array(const hd_float* d_in,
                             hd_size         array_size,
                             hd_size         count,
                             hd_float*       d_out);
