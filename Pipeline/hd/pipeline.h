/***************************************************************************
 *
 *   Copyright (C) 2012 by Ben Barsdell and Andrew Jameson
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "hd/types.h"
#include "hd/params.h"
#include "hd/error.h"

hd_error hd_create_pipeline(hd_pipeline* pipeline, hd_params params);
hd_error hd_execute(hd_pipeline pipeline,
                    const hd_byte* filterbank, hd_size nsamps, hd_size nbits,
                    hd_size first_idx, hd_size* nsamps_processed);
void     hd_destroy_pipeline(hd_pipeline pipeline);

#ifdef __cplusplus
} // closing brace for extern "C"
#endif
