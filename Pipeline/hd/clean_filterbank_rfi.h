
#pragma once

#include "hd/types.h"
#include "hd/error.h"
#include <dedisp.h>

hd_error clean_filterbank_rfi(dedisp_plan    plan,
                              const hd_byte* h_in,
                              hd_size        nsamps,
                              hd_size        nbits,
                              hd_byte*       h_out,
                              int*           h_killmask,
                              hd_float       dm,
                              hd_float       dt,
                              hd_float       baseline_length,
                              hd_float       rfi_tol,
                              hd_size        rfi_min_beams,
                              hd_size        boxcar_max);
