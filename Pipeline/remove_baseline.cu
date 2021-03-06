/***************************************************************************
 *
 *   Copyright (C) 2012 by Ben Barsdell and Andrew Jameson
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "hd/remove_baseline.h"
#include "hd/median_filter.h"
//#include "hd/write_time_series.h"

#include <thrust/device_vector.h>

#define OPTIMIZED

class RemoveBaselinePlan_impl {
	thrust::device_vector<hd_float> buf1;
	thrust::device_vector<hd_float> buf2;
	thrust::device_vector<hd_float> baseline;
public:
	hd_error exec(hd_float* d_data, hd_size beam_stride, hd_size beam_count,
	              hd_size smooth_radius, hd_size nbeams) {
		
		thrust::device_ptr<hd_float> d_data_begin(d_data);
	
		// This algorithm works by scrunching the data down to a time resolution
		//   representative of the desired smoothing length and then stretching
		//   it back out again. The scrunching is done using the median-of-5
		//   to ensure robustness against outliers (e.g., strong RFI spikes).
    hd_size dm_delay    = beam_stride - beam_count;

    hd_size beam_samples = beam_stride;

    // Note: This parameter allows tuning to match the smoothing length
    //         of the original iterative-clipping algorithm.
    hd_float oversample = 2;
    // Find the desired time resolution
    hd_size  sample_count = (hd_size)(oversample * hd_float(beam_samples)/(2*smooth_radius) + 0.5);
    if( sample_count == 0 ) {
      // Too few samples, no need to baseline
      return HD_NO_ERROR;
    }

    // As we will use median-of-5, round to sample_count times a power of five
    hd_size nscrunches  = (hd_size)(log(beam_samples/sample_count)/log(5.));
    hd_size count_round = pow(5.,nscrunches)*sample_count;

    /*
    std::cerr << "beam_stride=" << beam_stride << " nbeams=" << nbeams
              << " beam_count=" << beam_count << " smooth_radius=" << smooth_radius
              << " sample_count=" << sample_count << " nscrunches=" << nscrunches
              << " count_round=" << count_round << std::endl;
    */

#ifdef OPTIMIZED
    buf1.resize(nbeams * count_round);
    buf2.resize(nbeams * (count_round/5));
#else
    buf1.resize(count_round);
    buf2.resize(count_round/5);
#endif
    hd_float* buf1_ptr = thrust::raw_pointer_cast(&buf1[0]);
    hd_float* buf2_ptr = thrust::raw_pointer_cast(&buf2[0]);

#ifdef OPTIMIZED
    hd_size nbeams_minus = nbeams - 1;
    hd_size nsamps = beam_stride*nbeams - dm_delay;
    baseline.resize (beam_stride*nbeams);
    hd_float* baseline_ptr = thrust::raw_pointer_cast(&baseline[0]);

    // perform resampling of each beam to the rounded size [count_round]
    // first do nbeam-1, which are all the same length, and then the final beam
    linear_stretch_beam (d_data, beam_stride, buf1_ptr, count_round, nbeams_minus);
    linear_stretch (d_data + (nbeams_minus * beam_stride), beam_count, 
                    buf1_ptr + (nbeams_minus * count_round), count_round);

    // then median scrunch untill we reach the sample size, each beam is count_round
    for (hd_size size=count_round; size>sample_count; size/=5) 
    {
      median_scrunch5_beam (buf1_ptr, size, nbeams, buf2_ptr);
      std::swap(buf1_ptr, buf2_ptr);
    }

    // Note: output is now at buf1_ptr
    thrust::device_ptr<hd_float> buf1_begin(buf1_ptr);
    thrust::device_ptr<hd_float> buf2_begin(buf2_ptr);

    linear_stretch_beam (buf1_ptr, sample_count, baseline_ptr, beam_stride, nbeams_minus);
    linear_stretch (buf1_ptr + (nbeams_minus * sample_count), sample_count,
                    baseline_ptr + (nbeams_minus * beam_stride), beam_count);

    // TESTING
    // hd_size to_write = beam_count;
    // hd_size offset = nbeams_minus * beam_stride;
    // write_device_time_series(d_data + offset, to_write, 1.f, "new_pre_baseline.tim");
    // write_device_time_series(baseline_ptr + offset, to_write, 1.f, "new_the_baseline.tim");

    // Now we just subtract it off
    thrust::transform(d_data_begin, d_data_begin + nsamps,
                      baseline.begin(),
                      d_data_begin,
                      thrust::minus<hd_float>());

    // TESTING
    // write_device_time_series(d_data + offset, to_write, 1.f, "new_post_baseline.tim");

#else

    baseline.resize(beam_stride);
    hd_float* baseline_ptr = thrust::raw_pointer_cast(&baseline[0]);

    for (hd_size ibeam=0; ibeam<nbeams; ibeam++)
    {
      hd_size beam_offset = ibeam * beam_stride;
      hd_size count = beam_stride;
      if (ibeam == (nbeams - 1))
      {
        count = beam_count;
      }

      // First we re-sample to the rounded size
      linear_stretch(d_data + beam_offset, beam_count, buf1_ptr, count_round);

      // Then we median scrunch until we reach the sample size
      for( hd_size size=count_round; size>sample_count; size/=5 ) {
        median_scrunch5(buf1_ptr, size, buf2_ptr);
        std::swap(buf1_ptr, buf2_ptr);
      }
      // Note: Output is now at buf1_ptr
      thrust::device_ptr<hd_float> buf1_begin(buf1_ptr);
      thrust::device_ptr<hd_float> buf2_begin(buf2_ptr);
    
      // Then we need to extrapolate the ends
      linear_stretch(buf1_ptr, sample_count, buf2_ptr+1, sample_count*2);
      buf2_begin[0]                = 2*buf2_begin[1] - buf2_begin[2];
      buf2_begin[sample_count*2+1] = (2*buf2_begin[sample_count*2] -
                                      buf2_begin[sample_count*2-1]);
    
      // And finally we stretch back to the original length
      linear_stretch(buf2_ptr, sample_count*2+2, baseline_ptr, count);
    
      // TESTING
      // write_device_time_series(d_data+beam_offset, count, 1.f, "old_pre_baseline.tim");
      // write_device_time_series(baseline_ptr, count, 1.f, "old_the_baseline.tim");
    
      // Now we just subtract it off
      thrust::transform(d_data_begin + beam_offset, d_data_begin+beam_offset+count,
                        baseline.begin(),
                        d_data_begin + beam_offset,
                        thrust::minus<hd_float>());

      // TESTING
      // write_device_time_series(d_data+beam_offset, count, 1.f, "old_post_baseline.tim");
    }
#endif 

    // write_device_time_series(d_data, (beam_stride*nbeams)-dm_delay, 1.f, "post_baselines.tim");

    return HD_NO_ERROR;
	}
};

// Public interface (wrapper for implementation)
RemoveBaselinePlan::RemoveBaselinePlan()
	: m_impl(new RemoveBaselinePlan_impl) {}
hd_error RemoveBaselinePlan::exec(hd_float* d_data, hd_size beam_stride, hd_size count,
                                  hd_size smooth_radius,
                                  hd_size nbeams) {
	return m_impl->exec(d_data, beam_stride, count, smooth_radius, nbeams);
}
