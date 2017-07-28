/***************************************************************************
 *
 *   Copyright (C) 2012 by Ben Barsdell and Andrew Jameson
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "hd/find_giants.h"
#include "hd/cached_allocator.cuh"

// TESTING only
#include "hd/stopwatch.h"
#include <iostream>
//#define PRINT_BENCHMARKS

#include <thrust/device_vector.h>
#include <thrust/count.h>
#include <thrust/copy.h>
#include <thrust/adjacent_difference.h>
#include <thrust/scan.h>
#include <thrust/binary_search.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/discard_iterator.h>

#include <thrust/iterator/retag.h>

// Global instance of the custom temporary memory allocator for Thrust
// TODO: We should be calling g_allocator.free_all() somewhere at the end of
//         the application to ensure memory is freed before the underlying
//         device backend (e.g., CUDART) goes out of scope. Not sure exactly
//         where to put it though.
cached_allocator g_allocator;

template<typename T>
struct greater_than_val : public thrust::unary_function<T, bool> {
  T val;
  greater_than_val(T val_) : val(val_) {}
  inline __host__ __device__
  bool operator()(T x) const {
    return x > val;
  }
};


typedef thrust::tuple<hd_float, hd_size> tplfloatint;
struct greater_than_val_mb : public thrust::unary_function<tplfloatint, bool> {
  hd_float  val;
  const hd_size   stride;
  const hd_size   count;
  greater_than_val_mb(hd_float val_, hd_size stride_, hd_size count_) 
    : val(val_), stride(stride_), count(count_) {}
  inline __host__ __device__
  bool operator()(const tplfloatint& x) const {
    hd_size i = thrust::get<1>(x);
    hd_size samp = i % stride;
    hd_float samp_val = thrust::get<0>(x);
    return (samp < count && samp_val > val);
  }
};

template<typename T>
struct maximum_first : public thrust::binary_function<T,T,T> {
  inline __host__ __device__
  T operator()(T a, T b) const {
    return thrust::get<0>(a) >= thrust::get<0>(b) ? a : b;
  }
};

template<typename T>
struct nearby : public thrust::binary_function<T,T,bool> {
  T max_dist;
  nearby(T max_dist_) : max_dist(max_dist_) {}
  inline __host__ __device__
  bool operator()(T a, T b) const {
    return b <= a + max_dist;
  }
};
template<typename T>
struct not_nearby : public thrust::binary_function<T,T,bool> {
  T max_dist;
  not_nearby(T max_dist_) : max_dist(max_dist_) {}
  inline __host__ __device__
  bool operator()(T b, T a) const {
    return b > a + max_dist;
  }
};

template<typename T>
struct plus_one : public thrust::unary_function<T,T> {
  inline __host__ __device__
  T operator()(T x) const { return x+1; }
};

class GiantFinder_impl {
  thrust::device_vector<hd_float> d_giant_data;
  thrust::device_vector<hd_size>  d_giant_data_inds;
  thrust::device_vector<int>      d_giant_data_segments;
  thrust::device_vector<hd_size>  d_giant_data_seg_ids;
public:
  hd_error exec(const hd_float* d_data, hd_size count, hd_size nbeams,
                hd_float thresh, hd_size merge_dist,
                hd_size beam_stride, hd_size beam_count,
                thrust::device_vector<hd_float>& d_giant_peaks,
                thrust::device_vector<hd_size>&  d_giant_inds,
                thrust::device_vector<hd_size>&  d_giant_begins,
                thrust::device_vector<hd_size>&  d_giant_ends) {

    // This algorithm works by extracting all samples in the time series
    //   above thresh (the giant_data), segmenting those samples into
    //   isolated giants (based on merge_dist), and then computing the
    //   details of each giant into the d_giant_* arrays using
    //   reduce_by_key and some scatter operations.
    
    using thrust::copy_if;
    using thrust::make_zip_iterator;
    using thrust::make_tuple;
    using thrust::make_counting_iterator;
  
    typedef thrust::device_ptr<const hd_float> const_float_ptr;
    typedef thrust::device_ptr<hd_float>             float_ptr;
    typedef thrust::device_ptr<hd_size>              size_ptr;
 
    const_float_ptr d_data_begin(d_data);
    const_float_ptr d_data_end(d_data + count);
  
#ifdef PRINT_BENCHMARKS
    Stopwatch timer;
  
    timer.start();
#endif
 
    // Note: The calls to Thrust in this function are retagged to use a
    //         custom temporary memory allocator (cached_allocator.cuh).
    //       This turns out to be critical to performance!
    // Quickly count how much giant data there is so we know the space needed
    hd_size giant_data_count = 
        thrust::count_if(
            make_zip_iterator(make_tuple(thrust::retag<my_tag>(d_data_begin), make_counting_iterator(0u))),
            make_zip_iterator(make_tuple(thrust::retag<my_tag>(d_data_begin), make_counting_iterator(0u)))+count,
            greater_than_val_mb(thresh, beam_stride, beam_count));

    // We can bail early if there are no giants at all
    if( 0 == giant_data_count ) {
      //std::cout << "**** Found ZERO giants" << std::endl;
      return HD_NO_ERROR;
    }
  
#ifdef PRINT_BENCHMARKS
    cudaThreadSynchronize();
    timer.stop();
    std::cout << "count_if time:           " << timer.getTime() << " s" << std::endl;
    timer.reset();
  
    timer.start();
#endif
  
    d_giant_data.resize(giant_data_count);
    d_giant_data_inds.resize(giant_data_count);
  
#ifdef PRINT_BENCHMARKS
    cudaThreadSynchronize();
    timer.stop();
    std::cout << "giant_data resize time:  " << timer.getTime() << " s" << std::endl;
    timer.reset();
  
    // Copy all of the giant data and their locations into one place
  
    timer.start();
#endif
  
    hd_size giant_data_count2 = 
      copy_if(make_zip_iterator(make_tuple(thrust::retag<my_tag>(d_data_begin),
                                           make_counting_iterator(0u))),
              make_zip_iterator(make_tuple(thrust::retag<my_tag>(d_data_begin),
                                           make_counting_iterator(0u)))+count,
              make_zip_iterator(make_tuple(thrust::retag<my_tag>(d_data_begin),
                                           make_counting_iterator(0u))),
              make_zip_iterator(make_tuple(thrust::retag<my_tag>(d_giant_data.begin()),
                                           thrust::retag<my_tag>(d_giant_data_inds.begin()))),
              greater_than_val_mb(thresh, beam_stride, beam_count))
      - make_zip_iterator(make_tuple(thrust::retag<my_tag>(d_giant_data.begin()),
                                     thrust::retag<my_tag>(d_giant_data_inds.begin())));
  
#ifdef PRINT_BENCHMARKS
    cudaThreadSynchronize();
    timer.stop();
    std::cout << "giant_data copy_if time: " << timer.getTime() << " s" << std::endl;
    timer.reset();
  
    timer.start();
#endif
  
    // Create an array of head flags indicating candidate segments
    //thrust::device_vector<int> d_giant_data_segments(giant_data_count);
    d_giant_data_segments.resize(giant_data_count);
    thrust::adjacent_difference(thrust::retag<my_tag>(d_giant_data_inds.begin()),
                                thrust::retag<my_tag>(d_giant_data_inds.end()),
                                thrust::retag<my_tag>(d_giant_data_segments.begin()),
                                not_nearby<hd_size>(merge_dist));
  
    //hd_size giant_count_quick = thrust::count(d_giant_data_segments.begin(),
    //                                          d_giant_data_segments.end(),
    //                                          (int)true);
  
    // The first element is implicitly a segment head
    if( giant_data_count > 0 ) {
      d_giant_data_segments.front() = 0;
      //d_giant_data_segments.front() = 1;
    }
  
    //thrust::device_vector<hd_size> d_giant_data_seg_ids(d_giant_data_segments.size());
    d_giant_data_seg_ids.resize(d_giant_data_segments.size());
    
    thrust::inclusive_scan(thrust::retag<my_tag>(d_giant_data_segments.begin()),
                           thrust::retag<my_tag>(d_giant_data_segments.end()),
                           thrust::retag<my_tag>(d_giant_data_seg_ids.begin()));
  
    // We extract the number of giants from the end of the exclusive scan
    //hd_size giant_count = d_giant_data_seg_ids.back() +
    //  d_giant_data_segments.back() + 1;
    hd_size giant_count = d_giant_data_seg_ids.back() + 1;
    //hd_size giant_count = d_giant_data_seg_ids.back() +
    //  d_giant_data_segments.back();
  
    // Report back the actual number of giants found
    //total_giant_count = giant_count;
    
#ifdef PRINT_BENCHMARKS
    cudaThreadSynchronize();
    timer.stop();
    std::cout << "giant segments time:     " << timer.getTime() << " s" << std::endl;
    timer.reset();
  
    timer.start();
#endif
  
    hd_size new_giants_offset = d_giant_peaks.size();
    // Allocate space for the new giants
    d_giant_peaks.resize(d_giant_peaks.size() + giant_count);
    d_giant_inds.resize(d_giant_inds.size() + giant_count);
    d_giant_begins.resize(d_giant_begins.size() + giant_count);
    d_giant_ends.resize(d_giant_ends.size() + giant_count);
    float_ptr new_giant_peaks_begin(&d_giant_peaks[new_giants_offset]);
    size_ptr  new_giant_inds_begin(&d_giant_inds[new_giants_offset]);
    size_ptr  new_giant_begins_begin(&d_giant_begins[new_giants_offset]);
    size_ptr  new_giant_ends_begin(&d_giant_ends[new_giants_offset]);
  
#ifdef PRINT_BENCHMARKS
    cudaThreadSynchronize();
    timer.stop();
    std::cout << "giants resize time:      " << timer.getTime() << " s" << std::endl;
    timer.reset();
  
    timer.start();
#endif
    
  
    // Now we find the value (snr) and location (time) of each giant's maximum
    hd_size giant_count2 = 
      reduce_by_key(thrust::retag<my_tag>(d_giant_data_inds.begin()), // the keys
                    thrust::retag<my_tag>(d_giant_data_inds.end()),
                    make_zip_iterator(make_tuple(thrust::retag<my_tag>(d_giant_data.begin()),
                                                 thrust::retag<my_tag>(d_giant_data_inds.begin()))),
                    thrust::make_discard_iterator(), // the keys output
                    make_zip_iterator(make_tuple(thrust::retag<my_tag>(new_giant_peaks_begin),
                                                 thrust::retag<my_tag>(new_giant_inds_begin))),
                    nearby<hd_size>(merge_dist),
                    maximum_first<thrust::tuple<hd_float,hd_size> >())
      .second - make_zip_iterator(make_tuple(thrust::retag<my_tag>(new_giant_peaks_begin),
                                             thrust::retag<my_tag>(new_giant_inds_begin)));
  
#ifdef PRINT_BENCHMARKS
    cudaThreadSynchronize();
    timer.stop();
    std::cout << "reduce_by_key time:      " << timer.getTime() << " s" << std::endl;
    timer.reset();
  
    timer.start();
#endif
  
    // Now we make the first segment explicit
    if( giant_count > 0 ) {
      d_giant_data_segments[0] = 1;
    }
  
    // Create arrays of the beginning and end indices of each giant
    thrust::scatter_if(d_giant_data_inds.begin(),
                       d_giant_data_inds.end(),
                       d_giant_data_seg_ids.begin(),
                       d_giant_data_segments.begin(),
                       new_giant_begins_begin);
    thrust::scatter_if(make_transform_iterator(d_giant_data_inds.begin(),
                                                                     plus_one<hd_size>()),
                       make_transform_iterator(d_giant_data_inds.end() - 1,
                                                                     plus_one<hd_size>()),
                       d_giant_data_seg_ids.begin(),
                       d_giant_data_segments.begin() + 1,
                       new_giant_ends_begin);
    
    if( giant_count > 0 ) {
      d_giant_ends.back() = d_giant_data_inds.back() + 1;
    }
  
  
#ifdef PRINT_BENCHMARKS
    cudaThreadSynchronize();
    timer.stop();
    std::cout << "begin/end copy_if time:  " << timer.getTime() << " s" << std::endl;
    timer.reset();
  
    std::cout << "--------------------" << std::endl;
#endif
  
    return HD_NO_ERROR;
  }
  
};

// Public interface (wrapper for implementation)
GiantFinder::GiantFinder()
  : m_impl(new GiantFinder_impl) {}
hd_error GiantFinder::exec(const hd_float* d_data, hd_size count, hd_size nbeams,
                           hd_float thresh, hd_size merge_dist,
                           hd_size beam_stride, hd_size beam_count,
                           thrust::device_vector<hd_float>& d_giant_peaks,
                           thrust::device_vector<hd_size>&  d_giant_inds,
                           thrust::device_vector<hd_size>&  d_giant_begins,
                           thrust::device_vector<hd_size>&  d_giant_ends) {
  return m_impl->exec(d_data, count, nbeams,
                      thresh, merge_dist,
                      beam_stride, beam_count,
                      d_giant_peaks,
                      d_giant_inds,
                      d_giant_begins,
                      d_giant_ends);
}
