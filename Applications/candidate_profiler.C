/***************************************************************************
 *
 *   Copyright (C) 2012 by Ben Barsdell
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <string>
using std::string;
#include <vector>
using std::vector;
#include <fstream>
#include <iomanip>
#include <iterator>
#include <numeric>
#include <cstdlib> // For atoi
#include <cmath>
#include <algorithm>
#include <inttypes.h>

#include "hd/header.h"
#include <dedisp.h>

/*
  #SNR 7.11626 
samp_idx        2872128 
time    183.816 
filter  8 
dm_trial        795 
DM      309.013 
members         40 
begin   2872104 
end     2872200 
nbeams  1 
beam_mask       64 
prim_beam       7 
max_snr         7.11626 
beam    7 



#SNR 7.21726 
samp_idx        6285304 
time    402.259 
filter  7 
dm_trial        797 
DM      310.853 
members         4 
begin   6285300 
end     6285320 
nbeams  1 
beam_mask       64 
prim_beam       7 
max_snr         7.21726 
beam    7 

 */

template<typename InputIterator, typename OutputIterator>
void exclusive_scan(InputIterator begin, InputIterator end,
                    OutputIterator result_begin) {
	*result_begin = 0;
	for( size_t i=1; i<(size_t)(end-begin); ++i ) {
		*(result_begin+i) = *(result_begin+i-1) + *(begin+i-1);
	}
}

struct abs_deviation : public std::unary_function<float,float> {
	float median;
	abs_deviation(float median_) : median(median_) {}
	float operator()(float x) const { return std::fabs(x-median); }
};
struct multiply_by : public std::unary_function<float,float> {
	float val;
	multiply_by(float val_) : val(val_) {}
	float operator()(float x) const { return x * val; }
};
struct add : public std::unary_function<float,float> {
	float val;
	add(float val_) : val(val_) {}
	float operator()(float x) const { return x + val; }
};

int normalise_series_array(float* out, size_t nsamps, size_t nseries) {
	
	vector<float> buf(nsamps);
	for( size_t d=0; d<nseries; ++d ) {
		size_t i = d * nsamps;
		
		// Baseline
		/*
		std::copy(&out[i], &out[i+nsamps], buf.begin());
		std::nth_element(buf.begin(), buf.begin()+buf.size()/2, buf.end());
		float median = buf[buf.size()/2];
		std::transform(&out[i], &out[i+nsamps],
		               &out[i],
		               add(-median));
		*/
		// Note: We allow for a linear ramp in the baseline by finding
		//         the median in the first and second halves separately.
		std::copy(&out[i], &out[i+nsamps], buf.begin());
		std::nth_element(buf.begin(),
		                 buf.begin()+buf.size()/4,
		                 buf.begin()+buf.size()/2);
		float median1 = buf[buf.size()/4];
		std::nth_element(buf.begin()+buf.size()/2,
		                 buf.begin()+3*buf.size()/4,
		                 buf.begin()+buf.size());
		float median2 = buf[3*buf.size()/4];
		double gradient = (median2 - median1) / (nsamps/2);
		float  median0  = median1 - gradient*(nsamps/4);
		for( size_t j=0; j<nsamps; ++j ) {
			float baseline = median0 + j*gradient;
			out[i+j] -= baseline;
		}
		
		// Normalise
		std::transform(&out[i], &out[i+nsamps],
		               buf.begin(),
		               abs_deviation(0.));
		std::nth_element(buf.begin(), buf.begin()+buf.size()/2, buf.end());
		float mad = buf[buf.size()/2];
		float stddev = mad * 1.4826;
		std::transform(&out[i], &out[i+nsamps],
		               &out[i],
		               multiply_by(1./stddev));
	}
	
	return 0;
}

int gen_freq_time_plot(string filename,
                       size_t samp, size_t filter, float dm,
                       size_t out_nsamps, size_t fscrunch,
                       float* out, bool verbose) {
	typedef unsigned int word_type;
	
	std::ifstream in_file(filename.c_str(), std::ios::binary);
	if( !in_file ) {
		cerr << "ERROR: Could not open " << filename << endl;
		return -1;
	}
	
	SigprocHeader header;
	read_header(in_file, header);
	
	// Note: We use a dedisp plan here simply to compute the DM delay
	dedisp_error derror;
	dedisp_plan  plan;
	derror = dedisp_create_plan(&plan,
	                            header.nchans, header.tsamp,
	                            header.fch1, header.foff);
	if( derror != DEDISP_NO_ERROR ) {
		cerr << "dedisp_create_plan failed: "
		     << dedisp_get_error_string(derror) << endl;
		return -1;
	}
	dedisp_set_dm_list(plan, &dm, 1);
	size_t max_delay = dedisp_get_max_delay(plan);
	dedisp_destroy_plan(plan);

  size_t delay = 2;
	size_t tscrunch    = 1<<filter;
	size_t in_nsamps   = out_nsamps * tscrunch;
  if (in_nsamps < max_delay)
  {
    // find smallest power of 2 > max_delay
    while (delay < max_delay)
      delay *= 2;
    in_nsamps = delay;
    tscrunch = in_nsamps / out_nsamps;
  }
	size_t centre_samp = samp + max_delay/2;
	size_t first_samp  = centre_samp - in_nsamps/2;

  if (verbose)
  {	
    cerr << "max delay= " << max_delay<< endl;
    cerr << "delay= " << delay<< endl;
    cerr << "tscrunch = " << tscrunch << endl;
	  cerr << "Input sample range: "
	       << first_samp << " : " << first_samp+in_nsamps << " = " << in_nsamps << " samples" << endl;
    cerr << "DM = " << dm << endl;
  }

  // this is used by the plotter
  cout << "in_nsamps=" << in_nsamps << endl;
	
	size_t chans_per_word = sizeof(word_type)*8/header.nbits;
	size_t mask         = ((unsigned)1<<header.nbits) - 1;//(((unsigned)1<<(header.nbits-1))-1)*2+1;
	
	size_t stride_words = header.nchans/chans_per_word;
	size_t stride_bytes = header.nchans*sizeof(word_type)/chans_per_word;
	in_file.seekg(first_samp*stride_bytes, std::ios::cur);
	std::vector<word_type> packed_data(in_nsamps*stride_words, 0);
	in_file.read((char*)&packed_data[0], in_nsamps*stride_bytes);
	in_file.close();
	
	// Unpack and scrunch
	float peak = 255;
	for( size_t t=0; t<(size_t)in_nsamps; t+=tscrunch ) {
		for( size_t c=0; c<(size_t)header.nchans; c+=fscrunch ) {
			float sum = 0.f;
			for( size_t s=0; s<(size_t)tscrunch; ++s ) {
				for( size_t f=0; f<(size_t)fscrunch; ++f ) {
					size_t w = (c+f) / chans_per_word;
					size_t k = (c+f) % chans_per_word;
					
					size_t nchan_words = header.nchans / chans_per_word;
					word_type x = (packed_data[w + (t+s)*nchan_words]
					               >> (k*header.nbits)) & mask;
					sum += (float)x / mask * peak;
				}
			}
			//out[c/fscrunch + (t/tscrunch)*out_nchans] =

      //cerr << "out[" << (c/fscrunch*out_nsamps + (t/tscrunch)) << "] = " << sum / (tscrunch*fscrunch) << endl;
			out[c/fscrunch*out_nsamps + (t/tscrunch)] =
				sum / (tscrunch*fscrunch);
		}
	}
	
	size_t out_nchans = header.nchans/fscrunch;
	
	return normalise_series_array(out, out_nsamps, out_nchans);
}

int gen_dm_time_plot(string filename,
                     const float* dm_list, size_t dm_count,
                     size_t samp, size_t filter, size_t dm_idx,
                     size_t out_nsamps, size_t out_dm_count,
                     float* out,
                     bool do_filter=true, bool verbose=false) {
	typedef float out_type;
	
	// Scrunching parameters
	float dm_pulse_width = 40;
	float scrunch_tol    = 1.15;
	
	if( dm_idx >= dm_count ) {
		cerr << "Invalid DM index " << dm_idx << endl;
		return -1;
	}
	
	std::ifstream in_file(filename.c_str(), std::ios::binary);
	
	SigprocHeader header;
	read_header(in_file, header);
	
  if (verbose)
  {
	  cerr << "dt = " << header.tsamp << endl;
	  cerr << "f0 = " << header.fch1 << endl;
	  cerr << "df = " << header.foff << endl;
  }
	
	dedisp_error error;
	dedisp_plan  plan;
	error = dedisp_create_plan(&plan,
	                           header.nchans, header.tsamp,
	                           header.fch1, header.foff);
	if( error != DEDISP_NO_ERROR ) {
		cerr << "dedisp_create_plan failed: "
		     << dedisp_get_error_string(error) << endl;
		return -1;
	}
	
	// Set up DM list
	size_t dm_begin = dm_idx >= out_dm_count/2            ? dm_idx - out_dm_count/2 : 0;
	size_t dm_end   = dm_idx <  dm_count - ((out_dm_count-1)/2+1) ? dm_idx + (out_dm_count-1)/2+1 : dm_count;
  if (verbose)
  {
	  cerr << "Target DM      = " << dm_list[dm_idx] << endl;
	  cerr << "DM trial range = " << dm_begin << " : " << dm_end << endl;
  }
	out_dm_count = dm_end - dm_begin;
  if (verbose)
    cerr << "Out DM count = " << out_dm_count << endl;
	error = dedisp_set_dm_list(plan, &dm_list[dm_begin], out_dm_count);
	if( error != DEDISP_NO_ERROR ) {
		cerr << "dedisp_set_dm_list failed: "
		     << dedisp_get_error_string(error) << endl;
		return -1;
	}
	
	error = dedisp_enable_adaptive_dt(plan,
	                                  dm_pulse_width,
	                                  scrunch_tol);
	if( error != DEDISP_NO_ERROR ) {
		cerr << "dedisp_enable_adaptive_dt failed: "
		     << dedisp_get_error_string(error) << endl;
		return -1;
	}
	
	const dedisp_size* scrunch_factors = dedisp_get_dt_factors(plan);
	/*
	for( size_t i=0; i<out_dm_count; ++i ) {
		cerr << scrunch_factors[i] << "\t";
	}
	cerr << endl;
	*/

	size_t width       = 1 << filter;
	size_t nchan_bytes = header.nchans * header.nbits / (8*sizeof(dedisp_byte));
	size_t nsamps;    // nsamps dedispersed
	size_t in_nsamps; // nsamps for dedispersion
	size_t first_samp;

  size_t max_delay = dedisp_get_max_delay(plan);
  size_t delay = 2;

	if( !do_filter ) {
		nsamps     = out_nsamps * width;
		in_nsamps  = out_nsamps * width + dedisp_get_max_delay(plan);
		first_samp = samp >= out_nsamps/2 * width ? samp - out_nsamps/2 * width : 0;
	}
	else {
		nsamps     = out_nsamps + width-1;
		in_nsamps  = out_nsamps + width-1 + dedisp_get_max_delay(plan);
		first_samp = samp > out_nsamps/2 + width/2 ? samp - out_nsamps/2 - width/2 : 0;
	}

	// TODO: Ensure sample index is within upper bound
  if (verbose)
  {
    cerr << "Max dispersion delay = " << dedisp_get_max_delay(plan) << endl;
    cerr << "Input sample range = "
         << first_samp << " : "
         << first_samp+in_nsamps << endl;
    cerr << "(" << in_nsamps << " samples)" << endl;
  }
	
	vector<dedisp_byte> in(in_nsamps * nchan_bytes);
	
	// Read input data
	in_file.seekg(first_samp * nchan_bytes, std::ios::cur);
	in_file.read((char*)&in[0], in.size()*sizeof(dedisp_byte));
	in_file.close();
	
	size_t out_nbits = sizeof(out_type) * 8;
	vector<out_type> dedisped(nsamps * out_dm_count);

	// Dedisperse
	error = dedisp_execute(plan, in_nsamps,
	                       (dedisp_byte*)&in[0], header.nbits,
	                       (dedisp_byte*)&dedisped[0], out_nbits,
	                       0);
	if( error != DEDISP_NO_ERROR ) {
		cerr << "dedisp_execute failed: "
		     << dedisp_get_error_string(error) << endl;
		return -1;
	}
	
	if( !do_filter ) {
		
		// Scrunch in time
		for( size_t d=0; d<out_dm_count; ++d ) {
			for( size_t t=0; t<out_nsamps; ++t ) {
				size_t i = d * out_nsamps + t;
				out[i] = 0.f;
				size_t scrunch = scrunch_factors[d];
				for( size_t j=0; j<width/scrunch; ++j ) {
					size_t k = d*out_nsamps*width + t*(width/scrunch);
					out[i] += dedisped[k + j];
				}
			}
		}
		
	}
	else {
		
		// Un-scrunch in time
		vector<float> unscrunched(dedisped.size());
		for( size_t d=0; d<out_dm_count; ++d ) {
			for( size_t t=0; t<nsamps; ++t ) {
				size_t i = d * nsamps + t/scrunch_factors[d];
				size_t j = d * nsamps + t;
				unscrunched[j] = dedisped[i];
			}
		}
		
		// Filter in time
		vector<float> scanned(nsamps + 1);
		for( size_t d=0; d<out_dm_count; ++d ) {
			size_t i = d * nsamps;
			// Fast boxcar filter via exclusive scan
			// Note: One extra element so that we include the final value
			//exclusive_scan(&dedisped[i], &dedisped[i+nsamps+1],
			//               scanned.begin());
			exclusive_scan(&unscrunched[i], &unscrunched[i+nsamps+1],
			               scanned.begin());
			std::transform(scanned.begin()+width/2 + (width-1)/2+1,
			               scanned.begin()+width/2 + (width-1)/2+1 + out_nsamps,
			               scanned.begin()+width/2 - width/2,
			               &out[d * out_nsamps],
			               std::minus<float>());
		}
		
	}
	
	/*
	// Copy directly
	for( size_t d=0; d<out_dm_count; ++d ) {
	size_t i = d * nsamps;
	std::copy(&dedisped[i + width/2], &dedisped[i + width/2 + out_nsamps],
	&out[d * out_nsamps]);
	}
	*/
	
	return normalise_series_array(out, out_nsamps, out_dm_count);
}

// Helper functors
template<typename T>
struct min_t : public std::binary_function<T,T,T> {
	T operator()(T a, T b) const { return std::min(a,b); }
};
template<typename T>
struct max_t : public std::binary_function<T,T,T> {
	T operator()(T a, T b) const { return std::max(a,b); }
};

void usage(char * binary)
{
  cerr << "Usage: " << binary << " filename.fil dmlist samp filter dm_idx out_nsamps out_dm_count fscrunch [do_filter=0]" << endl;
  cerr << "    sigproc_file        2bit sigproc filterbank file" << endl;
  cerr << "    dmlist              file containing list of DMs to plot" << endl;
  cerr << "    samp                sample index of the candidate" << endl;
  cerr << "    filter              filter (log2) in which candidate was found" << endl;
  cerr << "    dm_idx              line in dmlist which corresponds to nominal DM for candidate" << endl;
  cerr << "    out_nsamps          " << endl;
  cerr << "    out_dm_count        should match the number of lines in dmlist" << endl;
  cerr << "    fscrunch            fscrunching factor for Freq vs Time data" << endl;
  cerr << "    do_filter           not sure - Ben?" << endl;
}

bool IsNan( float f)
{ 
  union { float f; uint32_t x; } u = { f };
  return (u.x << 1) > 0xff000000u;
}

int main(int argc, char* argv[])
{
	if ( argc <= 8 ) 
  {
    cerr << "Error: expected 8 arguments found " << argc << endl;
    usage(argv[0]);
		return 0;
	}
	bool do_filter = false;
  bool verbose = true;
	
	string input_name   = argv[1];
	string dmlist_name  = argv[2];
	size_t samp         = atoi(argv[3]);
	size_t filter       = atoi(argv[4]);
	size_t dm_idx       = atoi(argv[5]);
	size_t out_nsamps   = atoi(argv[6]);//64;
	size_t out_dm_count = atoi(argv[7]);//32;
	size_t fscrunch     = atoi(argv[8]);//16;
	if( argc >= 10 ) {
		do_filter = atoi(argv[9]);
	}
	
  if (verbose)
  {
	  if( do_filter ) {
		  cerr << "Filtering enabled" << endl;
	  }
	  else {
		  cerr << "Filtering disabled" << endl;
	  }
	
  	cerr << "out_nsamps   = " << out_nsamps << endl;
	  cerr << "out_dm_count = " << out_dm_count << endl;
	  cerr << "Reading header info..." << endl;
  }

	std::ifstream in_file(input_name.c_str(), std::ios::binary);
	if( !in_file ) {
		cerr << "ERROR: Could not open " << input_name << endl;
		return -1;
	}
	SigprocHeader header;
	read_header(in_file, header);
	in_file.close();
	
  if (verbose)
	  cerr << "Loading DM list '" << dmlist_name << "'..." << endl;
	vector<float> dm_list;
	std::ifstream dm_file(dmlist_name.c_str());
	if( !dm_file ) {
		cerr << "ERROR: Could not open " << dmlist_name << endl;
		return -1;
	}
	std::copy(std::istream_iterator<float>(dm_file),
	          std::istream_iterator<float>(),
	          std::back_inserter(dm_list));

  if (verbose)
    cerr << "Read " << dm_list.size() << " dms from list" << endl;
	
	vector<float> dm_time_data(out_nsamps * out_dm_count, 0.f);
	
  if (verbose)
	  cerr << "Computing DM-time plot..." << endl;
	int error = gen_dm_time_plot(input_name,
	                             &dm_list[0], dm_list.size(),
	                             samp, filter, dm_idx,
	                             out_nsamps, out_dm_count,
	                             &dm_time_data[0],
	                             do_filter, verbose);
	if( error ) {
		cerr << "gen_dm_time_plot failed" << endl;
		return -1;
	}

  // now generate a 0-dm list
	vector<float> dm0_time_data(out_nsamps, 0.f);
  error = gen_dm_time_plot(input_name,
                           &dm_list[0], dm_list.size(),
                           samp, filter, 0,
                           out_nsamps, 1,
                           &dm0_time_data[0],
                           do_filter, verbose);
  if( error ) {
    cerr << "gen_dm0_time_plot failed" << endl;
    return -1;
  }

	size_t out_nchans = header.nchans / fscrunch;
	vector<float> freq_time_data(out_nsamps * out_nchans, 0.f);
	
  if (verbose)
	  cerr << "Computing freq-time plot..." << endl;
	size_t ft_plot_filter = filter > 2 ? filter - 2 : 0;
	gen_freq_time_plot(input_name,
	                   samp, ft_plot_filter, dm_list[dm_idx],
	                   out_nsamps, fscrunch,
	                   &freq_time_data[0],
                     verbose);
	
  if (verbose)
	  cerr << "Writing freq-time data to freq_time.dat..." << endl;
	//size_t levels = 256;
	std::ofstream freq_time_file("freq_time.dat");
	for( size_t c=0; c<out_nchans; ++c ) {
		for( size_t t=0; t<out_nsamps; ++t ) {
			float  rawval = freq_time_data[c*out_nsamps+t];
      {
        //size_t val;//    = (rawval - minval) / (maxval - minval) * (levels-1);
        /*
        float min_sigma = 0;
        float max_sigma = 8;
        val = std::min(std::max((rawval - min_sigma) /
                                (max_sigma - min_sigma) * (levels-1),
                                 0.f),
                       (float)levels-1);
        */
        //val = rawval >= 3. ? val : 0;
			  //out_file << val << "\t";
        if (IsNan(rawval))
			    rawval = 0;

        // since channels < 150 are always bad for HTRU data
        if ((c * fscrunch) < 150)
          rawval = 0;

			  freq_time_file << rawval  << "\t";
      }
		}
		freq_time_file << "\n";
	}
	freq_time_file.close();
	
  if (verbose)
	  cerr << "Computing result bounds..." << endl;
	float minval = std::accumulate(dm_time_data.begin(),
	                               dm_time_data.end(),
	                               dm_time_data[0],
	                               min_t<float>());
	float maxval = std::accumulate(dm_time_data.begin(),
	                               dm_time_data.end(),
	                               dm_time_data[0],
	                               max_t<float>());
	size_t max_idx = std::max_element(dm_time_data.begin(),
	                                  dm_time_data.end()) - dm_time_data.begin();
  if (verbose)
  {
  	cerr << "Max value        = " << maxval << endl;
	  cerr << "Min value        = " << minval << endl;
  	cerr << "Max data index   = " << max_idx << endl;
	  cerr << "Max sample index = " << max_idx%out_nsamps << endl;
	  cerr << " Dist from given = " << abs(max_idx%out_nsamps - out_nsamps/2) << endl;
  }
	// TODO: Fix this for DM trial bounds
	size_t dm_idx_start = dm_idx > out_dm_count/2 ? dm_idx - out_dm_count/2 : 0;
	size_t max_dm_idx   = dm_idx_start + max_idx/out_nsamps;
  if (verbose)
  {
	  cerr << "Max DM index     = " << max_dm_idx << endl;
	  cerr << " Dist from given = " << abs(max_dm_idx - dm_idx) << endl;
	  cerr << "Extracting SNR vs. DM..." << endl;
  }
	std::ofstream snr_dm_file("snr_dm.dat");
	snr_dm_file << "#dm_trial\tDM\tSNR" << endl;
	for( size_t d=0; d<std::min(out_dm_count,dm_list.size()); ++d ) {
		//size_t t = out_nsamps/2;
		size_t t = max_idx % out_nsamps;
		size_t dd = dm_idx_start + d;
		snr_dm_file << dd << "\t"
		            << std::setprecision(10) << dm_list[dd] << "\t"
		            << dm_time_data[d*out_nsamps + t] << endl;
	}
	snr_dm_file.close();
	
	float width = (1 << filter);
	
  if (verbose)	
	  cerr << "Extracting SNR vs. time..." << endl;
	std::ofstream snr_time_file("snr_time.dat");
	snr_time_file << "#sample\ttime\tSNR\tSNR@DM-1\tSNR@DM+1\tSNR@DM0" << endl;
	for( size_t t=0; t<out_nsamps; ++t ) {
		//size_t d = dm_idx_start + dm_idx;//out_dm_count/2;
		size_t d = max_idx/out_nsamps;
		//size_t tt = samp - out_nsamps/2 + t;
		float tt = (samp + ((float)t - out_nsamps/2)*width) * header.tsamp;
		snr_time_file << tt << "\t"
		              << std::setprecision(10) << tt << "\t"
		              << dm_time_data[d*out_nsamps + t] << "\t"
		              << dm_time_data[(d-1)*out_nsamps + t] << "\t"
		              << dm_time_data[(d+1)*out_nsamps + t] << "\t"
		              << dm0_time_data[t] << endl;
	}
	snr_time_file.close();
	
	size_t max_image_nsamps = 16384; // 4096
	
	if (out_nsamps > max_image_nsamps ) 
  {
    if (verbose)	
    {
		  cerr << "Skipping image output for large data" << endl;
		  cerr << "Done." << endl;
    }
		return 0;
	}

  if (verbose)
	  cerr << "Writing output to 'dm_time.pgm'..." << endl;
	//size_t levels = 16384;
	size_t levels = 256;
	std::ofstream out_file("dm_time.pgm");
	out_file << "P2 "
			 << out_nsamps << " " << out_dm_count << " "
			 << levels-1 << endl;
	for( size_t d=0; d<out_dm_count; ++d ) {
		for( size_t t=0; t<out_nsamps; ++t ) {
			float  rawval = dm_time_data[d*out_nsamps+t];
			size_t val;//    = (rawval - minval) / (maxval - minval) * (levels-1);
			
			float min_sigma = 0;
			float max_sigma = 8;
			val = std::min(std::max((rawval - min_sigma) /
			                        (max_sigma - min_sigma) * (levels-1),
			                         0.f),
			               (float)levels-1);
			
			//val = rawval >= 3. ? val : 0;
			out_file << val << "\t";
		}
		out_file << "\n";
	}
	out_file.close();

  if (verbose) 
	  cerr << "Done" << endl;
}
