/***************************************************************************
 *
 *   Copyright (C) 2012 by Ben Barsdell and Andrew Jameson
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "hd/get_rms.h"
#include "hd/median_filter.h"
//#include "hd/write_time_series.h"

#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/iterator/constant_iterator.h>

#include <sstream>

template<typename T>
struct absolute_val : public thrust::unary_function<T,T> {
	inline __host__ __device__
	T operator()(T x) const { return abs(x); }
};

class GetRMSPlan_impl {
	thrust::device_vector<hd_float> buf1;
	thrust::device_vector<hd_float> buf2;
public:

  // This algorithm works by taking the absolute values of the data
  //   and then repeatedly scrunching them using median-of-5 in order
  //   to approximate the median absolute deviation. The RMS is then
  //   just 1.4862 times this.

	hd_float exec(hd_float* d_data, hd_size count) {

		thrust::device_ptr<hd_float> d_data_begin(d_data);
		
		buf1.resize(count);
		buf2.resize(count/5);
		hd_float* buf1_ptr = thrust::raw_pointer_cast(&buf1[0]);
		hd_float* buf2_ptr = thrust::raw_pointer_cast(&buf2[0]);

		thrust::transform(d_data_begin, d_data_begin+count,
		                  buf1.begin(),
		                  absolute_val<hd_float>());
		
		for( hd_size size=count; size>1; size/=5 ) {
			median_scrunch5(buf1_ptr, size, buf2_ptr);
			std::swap(buf1_ptr, buf2_ptr);
		}

		// Note: Result is now at buf1_ptr
		thrust::device_ptr<hd_float> buf1_begin(buf1_ptr);
		hd_float med_abs_dev = buf1_begin[0];
		hd_float rms = med_abs_dev * 1.4862;
		
		return rms;
	}
};

class GetRMSPlanMB_impl {
  thrust::device_vector<hd_float> buf1;
  thrust::device_vector<hd_float> buf2;
public:

  hd_error exec_multibeam (hd_float* d_data, hd_float * d_rms, 
                           hd_size beam_stride, hd_size beam_count, 
                           hd_size nbeams) {

    thrust::device_ptr<hd_float> d_data_begin(d_data);

    // beam_count are uncorrupted samples
    hd_size count = nbeams * beam_count;

    buf1.resize(count);
    buf2.resize(count/5);

    hd_float* buf1_ptr = thrust::raw_pointer_cast(&buf1[0]);
    hd_float* buf2_ptr = thrust::raw_pointer_cast(&buf2[0]);

    reblock_abs_beam (d_data, beam_stride, beam_count, buf1_ptr, nbeams);

    for (hd_size size=beam_count; size>1; size/=5)
    {
      median_scrunch5_beam(buf1_ptr, size, nbeams, buf2_ptr);
      std::swap(buf1_ptr, buf2_ptr);
    }

    // Note: Result is now at buf1_ptr
    thrust::device_ptr<hd_float> buf1_begin(buf1_ptr);
    thrust::device_ptr<hd_float> d_rms_begin(d_rms);

    // Convert to RMS 
    thrust::transform (buf1_begin, buf1_begin + nbeams,
                       thrust::make_constant_iterator(1.4862),
                       d_rms_begin,
                       thrust::multiplies<hd_float>());

    return HD_NO_ERROR;
  }

  // compute single RMS across block, discard overlapping regions between beams
  // useful if original timeseries has already been normalised, for wider filters
  hd_float exec_multibeam (hd_float* d_data,
                           hd_size beam_stride, hd_size beam_count,
                           hd_size nbeams) {

    thrust::device_ptr<hd_float> d_data_begin(d_data);
    hd_size count = nbeams * beam_count;

    buf1.resize(count);
    buf2.resize(count/5);

    hd_float* buf1_ptr = thrust::raw_pointer_cast(&buf1[0]);
    hd_float* buf2_ptr = thrust::raw_pointer_cast(&buf2[0]);

    reblock_abs_beam (d_data, beam_stride, beam_count, buf1_ptr, nbeams);

    for (hd_size size=count; size>1; size/=5)
    {
      median_scrunch5(buf1_ptr, size, buf2_ptr);
      std::swap(buf1_ptr, buf2_ptr);
    }

    // Note: Result is now at buf1_ptr
    thrust::device_ptr<hd_float> buf1_begin(buf1_ptr);
    hd_float med_abs_dev = buf1_begin[0];
    hd_float rms = med_abs_dev * 1.4862;

    return rms;
  }

};

// Public interface (wrapper for implementation)
GetRMSPlan::GetRMSPlan()
	: m_impl(new GetRMSPlan_impl) {}
hd_float GetRMSPlan::exec(hd_float* d_data, hd_size count) {
	return m_impl->exec(d_data, count);
}

GetRMSPlanMB::GetRMSPlanMB()
  : m_impl(new GetRMSPlanMB_impl) {}
hd_error GetRMSPlanMB::exec_multibeam(hd_float* d_data, hd_float* d_rms, hd_size beam_stride, hd_size beam_count, hd_size nbeams) {
  return m_impl->exec_multibeam(d_data, d_rms, beam_stride, beam_count, nbeams);
}
hd_float GetRMSPlanMB::exec_multibeam(hd_float* d_data, hd_size beam_stride, hd_size beam_count, hd_size nbeams) {
  return m_impl->exec_multibeam(d_data, beam_stride, beam_count, nbeams);
}

// Convenience functions for one-off calls
hd_float get_rms(hd_float* d_data, hd_size count) {
	return GetRMSPlan().exec(d_data, count);
}

hd_error get_rms_multibeam (hd_float* d_data, hd_float* d_rms, hd_size beam_stride, hd_size beam_count, hd_size nbeams) {
  return GetRMSPlanMB().exec_multibeam(d_data, d_rms, beam_stride, beam_count, nbeams);
}

hd_error normalise(hd_float* d_data, hd_size count)
{
	thrust::device_ptr<hd_float> d_data_begin(d_data);
	thrust::device_ptr<hd_float> d_data_end(d_data + count);
	
	hd_float rms = get_rms(d_data, count);
	thrust::transform(d_data_begin, d_data_end,
	                  thrust::make_constant_iterator(hd_float(1.0)/rms),
	                  d_data_begin,
	                  thrust::multiplies<hd_float>());
	
	return HD_NO_ERROR;
}

struct normalise_beam_kernel
  : public thrust::unary_function<hd_float,hd_float> {
  const hd_float* in;
  const hd_float* rms;
  const hd_size   size;
  normalise_beam_kernel(const hd_float* in_, const hd_float* rms_, hd_size size_)
    : in(in_), rms(rms_), size(size_) {}
  inline __host__ __device__
  hd_float operator()(unsigned int i) const {
    hd_size beam = i / size;
    return in[i] / rms[beam];
  }
};


hd_error normalise_multibeam (hd_float* d_data, hd_float block_rms, hd_size count)
{
  thrust::device_ptr<hd_float> d_data_begin(d_data);
  thrust::device_ptr<hd_float> d_data_end(d_data + count);

  thrust::transform(d_data_begin, d_data_end,
                    thrust::make_constant_iterator(hd_float(1.0)/block_rms),
                    d_data_begin,
                    thrust::multiplies<hd_float>());

  return HD_NO_ERROR;
}

hd_error normalise_multibeam (hd_float* d_data, hd_float * d_rms, hd_size beam_stride, hd_size cur_nsamps, hd_size nbeams)
{
  thrust::device_ptr<hd_float> d_data_begin(d_data);
  using thrust::make_counting_iterator;
  thrust::transform(make_counting_iterator<unsigned int>(0),
                    make_counting_iterator<unsigned int>(cur_nsamps),
                    d_data_begin,
                    normalise_beam_kernel (d_data, d_rms, beam_stride));

  return HD_NO_ERROR;
}

struct reblock_abs_beam_functor
  : public thrust::unary_function<hd_float,hd_float> {
  const hd_float* in;
  hd_size         length;
  hd_float        delta;
  reblock_abs_beam_functor(const hd_float* in_,
                       hd_size in_length,
                       hd_size in_delta)
    : in(in_), length(in_length), delta(in_delta) {}
  inline __host__ __device__
  hd_float operator()(unsigned int o) const {
    hd_size beam = o / length;
    hd_size i = (beam * delta) + o;
    return fabsf(in[i]);
  }
};

// reblock input into output using an input stride and length
hd_error reblock_abs_beam (const hd_float* d_in,
                           hd_size         in_stride,
                           hd_size         in_length,
                           hd_float *      d_out,
                           hd_size         nbeam)
 {
  using thrust::make_counting_iterator;
  hd_size out_count = in_length * nbeam;
  thrust::device_ptr<hd_float> d_out_begin(d_out);
  hd_size delta = in_stride - in_length;

  thrust::transform(make_counting_iterator<unsigned int>(0),
                    make_counting_iterator<unsigned int>(out_count),
                    d_out_begin,
                    reblock_abs_beam_functor(d_in, in_length, delta));
  return HD_NO_ERROR;
}

