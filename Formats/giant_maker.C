/***************************************************************************
 *
 *   Copyright (C) 2012 by Ben Barsdell 
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

/*
  TODO: Sort out quantisation scaling
          How to decide what range to clip to?
        
 */
//#include "hd/header.h"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <fstream>
#include <string>
#include <sstream>
using std::string;
#include <cstdlib>
#include <cmath> // For pow, fmod
#include <vector>
#include <algorithm>


// TODO: Not sure if this works or not
float g_fScale = 2.0f / 0xffffffff;
int g_x1 = 0x67452301;
int g_x2 = 0xefcdab89;
void add_white_noise2(float* _fpDstBuffer,
                      unsigned int _uiBufferSize,
                      float _fLevel) {
	_fLevel *= g_fScale;
	
	while( _uiBufferSize-- ) {
		g_x1 ^= g_x2;
		*_fpDstBuffer++ = g_x2 * _fLevel;
		g_x2 += g_x1;
	}
}

void add_white_noise(float* data, size_t count,
                     float mean=0.0, float rms=1.0) {
	for( size_t i=0; i<count; i+=2 ) {
		// Sample Gaussian dist using Box-Muller transform
		float u1 = drand48();
		float u2 = drand48();
		float R     = sqrt(-2*log(u1));
		float theta = 2*3.141593 * u2;
		float z1 = R * cos(theta);
		float z2 = R * sin(theta);
		
		data[i]   += z1*rms + mean;
		data[i+1] += z2*rms + mean;
		
		//cout << z1*rms + mean << endl;
	}
}

void add_brown_noise(float* data, size_t count,
                     float mean=0.0, float rms=1.0) {
	float z = 0.0;
	for( size_t i=0; i<count; i+=1 ) {
		// Sample Gaussian dist using Box-Muller transform
		float u1    = drand48();
		float u2    = drand48();
		float R     = sqrt(-2*log(u1));
		float theta = 2*3.141593 * u2;
		float z1    = R * cos(theta);
		z += z1*rms;
		
		data[i] += z + mean;
	}
}

// TODO: This is not parameterised very smartly
void add_red_noise(float* data, size_t count,
                   size_t min_k, size_t max_k, float amp) {
	for( size_t k=min_k; k<=max_k; ++k ) {
		float phase = drand48() * 2;
		//float phase = 1.0;
		for( size_t i=0; i<count; i+=1 ) {
			//data[i] += amp*(1+log(k)/log(2)) * sin(((float)i/count+phase) * 3.141593 * k);
			data[i] += amp * sin(((float)i/count+phase) * 3.141593 * k);
		}
	}
}

void add_narrowband_rfi(float* data, size_t count, size_t value) {
	// TODO: Improve this. Try making it intermittent in time
	for( size_t i=0; i<count; ++i ) {
		data[i] += value;
	}
}
/*
struct DelayCalc {
	double f0, df, a, b;
	DelayCalc(double dt, double f0_, double df_)
		: f0(f0_), df(df_), a(4.148808e3/dt), b(1.0/(f0_*f0_)) {}
	double operator()(unsigned v) {
		return a * (1.0 / pow(f0 + v*df, 2.0) - b);
	}
};
*/
// Returns the time delay (EXLUDING THE DM) at a given frequency
float get_delay(unsigned v, double dt, double f0, double df) {
	double k = 4.148808e3;
	return k /* / dt */ * (1.0 / pow(f0 + v*df, 2.0) - 1.0 / (f0*f0));
}
float get_dm_smear(float DM, float f0, float df) {
	double k = 4.148808e3;
	return 2*k * DM * df / (f0*f0*f0);
}

// Adds a top hat signal with the given parameters to a time series
void add_giant(float* data, size_t nsamps, size_t nchans, size_t chan,
               float time, float width, float flux, float DM,
               float dt, float f0, float df)
{
	float start_time = time - width/2;
	float end_time   = time + width/2;
	
	// Delay the times according to the dispersion
	start_time += DM * get_delay(chan, dt, f0, df);
	end_time   += DM * get_delay(chan+1, dt, f0, df);
	
	// The flux is now distributed evenly between start_time and end_time
	
	size_t first_bin = (size_t)(start_time / dt + 0.5);
	size_t last_bin  = (size_t)(end_time / dt + 0.5);
	
	float flux_rate  = flux / (end_time - start_time);
	
	/*
	cout << "width =        " << width << endl;
	cout << "dsprsd width = " << (end_time - start_time) << endl; 
	cout << "start_time =   " << start_time << endl;
	cout << "end_time =     " << end_time << endl;
	cout << "first_bin =    " << first_bin << endl;
	cout << "last_bin =     " << last_bin << endl;
	cout << "flux_rate =    " << flux_rate << endl;
	*/
	
	for( size_t i=first_bin-1; i<=last_bin+1; ++i ) {
		float bin_start_time = (i-0.5) * dt;
		float bin_end_time   = (i+0.5) * dt;
		float overlap_start_time = std::max(start_time, bin_start_time);
		float overlap_end_time   = std::min(end_time,   bin_end_time);
		float overlap_time = std::max(overlap_end_time -
		                              overlap_start_time, 0.f);
		float bin_flux = flux_rate * overlap_time;
		
		// TESTING
		//bin_flux *= sqrt(dt / overlap_time);
		
		// TODO: Is this right or not?
		//bin_flux *= sqrt(overlap_time/dt);
		/*
		cout << "bin_start_time =     " << bin_start_time << endl;
		cout << "bin_end_time =       " << bin_end_time << endl;
		cout << "overlap_start_time = " << overlap_start_time << endl;
		cout << "overlap_end_time =   " << overlap_end_time << endl;
		cout << "overlap_time =       " << overlap_time << endl;
		cout << "bin_flux =           " << bin_flux << endl;
		*/
		//cout << "overlap_time =       " << overlap_time << endl;
		if( i >= 0 && i < nsamps ) {
			data[i] += bin_flux;
		}
	}
}

struct Giant {
	float time;
	float flux;
	float width;
	float DM;
};

int main(int argc, char* argv[])
{
	typedef unsigned int word_type;
	//typedef unsigned char word_type;
	
	float  time    = 15; // s
	//float  time    = 60; // s
	size_t nchans  = 1024;
	size_t nbits   = 2;
	float  dt      = 6.4e-5;
	float  f0      = 1581.804688f;
	float  df      = -0.390625f;
	float  min_val = -7.0;
	float  max_val = +7.0;
	size_t seed    = 1234;
	bool   red_noise = true;
	string out_filename = "giants.fil";
	string giants_filename = "giants.dat";
	//string out_filename = "giants.tim";
	
	// Parse command line arguments
	// ----------------------------
	if( argc == 1 ) {
		cout << "-----------" << endl;
		cout << "Giant Maker" << endl;
		cout << "-----------" << endl;
		cout << "By Ben Barsdell (2012)" << endl;
		cout << "About: Generates filterbank data containing pre-specified giants." << endl;
		cout << "Usage: " << argv[0] << " -i giant_file [options]" << endl;
		cout << "Options:" << endl;
		cout << "-i filename\tInput file containing list of giants" << endl;
		cout << "-o filename\tOutput file for filterbank data" << endl;
		cout << "-n -nbits int\tNumber of bits per sample for output data" << endl;
		cout << "-s -seed int\tSeed value for random number generator" << endl;
		cout << "-t -time float\tDuration of output data in seconds" << endl;
		cout << "-m -range float float\tMin and max with which to scale values" << endl;
		cout << "-r -red int\tAdd red noise [0/1]" << endl;
		cout << "-dt float\tSampling time in seconds" << endl;
		cout << "-f0 float\tFrequency of first channel in MHz" << endl;
		cout << "-df float\tFrequency step between channels in MHz" << endl;
		cout << endl;
		cout << "Input files should list one giant per line as follows:" << endl;
		cout << "Time(secs)\tSNR\twidth(secs)\tDM" << endl;
		cout << "Note that SNR represents the optimal detection SNR." << endl;
		cout << endl;
	}
	int a = 0;
	while( ++a < argc ) {
		if( argv[a] == string("-i") ) { giants_filename = argv[++a]; }
		else if( argv[a] == string("-o") ) { out_filename = argv[++a]; }
		else if( argv[a] == string("-n") || argv[a] == string("-nbits") ) {
			nbits = atoi(argv[++a]);
		}
		else if( argv[a] == string("-s") || argv[a] == string("-seed") ) {
			seed = atoi(argv[++a]);
		}
		else if( argv[a] == string("-t") || argv[a] == string("-time") ) {
			time = atof(argv[++a]);
		}
		else if( argv[a] == string("-dt") ) { dt = atof(argv[++a]); }
		else if( argv[a] == string("-f0") ) { f0 = atof(argv[++a]); }
		else if( argv[a] == string("-df") ) { df = atof(argv[++a]); }
		else if( argv[a] == string("-r") || argv[a] == string("-red") ) {
			red_noise = atoi(argv[++a]);
		}
		else if( argv[a] == string("-m") || argv[a] == string("-range") ) {
			min_val = atof(argv[++a]);
			max_val = atof(argv[++a]);
		}
		else {
			cout << "WARNING: Unknown argument '" << argv[a] << "'" << endl;
		}
	}
	// ----------------------------
	
	size_t nsamps  = (size_t)(time / dt + 0.5);
	size_t stride  = nsamps;
	size_t chans_per_word = sizeof(word_type)*8 / nbits;
	size_t ostride = nchans / chans_per_word;
	
	cout << "nsamps =         " << nsamps << endl;
	cout << "nchans =         " << nchans << endl;
	cout << "stride =         " << stride << endl;
	cout << "nbits =          " << nbits << endl;
	cout << "ostride =        " << ostride << endl;
	cout << "chans_per_word = " << chans_per_word << endl;
	
	std::vector<Giant> giants;
	std::ifstream giants_file(giants_filename.c_str());
	if( !giants_file ) {
		cerr << "ERROR: Could not open '" << giants_filename << "'" << endl;
		return -1;
	}
	string line;
	while( getline(giants_file, line) ) {
		if( line == "" || line[0] == '#' ) continue;
		
		Giant new_giant;
		std::stringstream(line) >> new_giant.time
		                        >> new_giant.flux
		                        >> new_giant.width
		                        >> new_giant.DM;
		if( new_giant.time < 0 || new_giant.time > time ) {
			cout << "WARNING: Giant time out of range" << endl;
		}
		giants.push_back(new_giant);
	}
	giants_file.close();
	
	cout << "Allocating memory..." << endl;
	std::vector<float>     filterbank(stride * nchans, 0.0);
	std::vector<word_type> quantised(ostride * nsamps, 0);
	
	srand48(seed);
	
	cout << "Adding white noise..." << endl;
	for( size_t c=0; c<nchans; ++c ) {
		float* time_series = &filterbank[c*stride];
		// HACK
		//if( c == 0 )
		add_white_noise(time_series, nsamps);
	}
	
	if( red_noise ) {
		cout << "Adding red noise..." << endl;
		std::vector<float> red_noise_vals(nsamps, 0.f);
		add_red_noise(&red_noise_vals[0], nsamps, 0, 8, 0.01);
		
		for( size_t c=0; c<nchans; ++c ) {
			float* time_series = &filterbank[c*stride];
			// TODO: This is not parameterised very smartly
			// HACK
			//if( c == 0 )
			std::transform(time_series, time_series + nsamps,
			               &red_noise_vals[0],
			               time_series,
			               std::plus<float>());
		}
	}
	
	/*
	  Width Avg detection efficiency
	  64us           0.726
	  40us           0.863
	  20us           0.903
	  10us           0.806
	  1us           0.903
	*/
	
	cout << "Adding giants..." << endl;
	for( size_t c=0; c<nchans; ++c ) {
		float* time_series = &filterbank[c*stride];
		
		for( size_t g=0; g<giants.size(); ++g ) {
			
			// Here we adjust the SNR for the effects of finite pulse width,
			//   channel width and sampling time.
			//   We do this because a detection pipeline can't do anything
			//     about them. In contrast, the detection pipeline *can*
			//     choose things like the DM trials and filter widths, so
			//     these are not taken into account here.
			float t_i  = giants[g].width;
			float t_dm = get_dm_smear(giants[g].DM, f0, df);
			
			float observed_width = sqrt(t_i*t_i + t_dm*t_dm + dt*dt);
			//size_t nbins = (size_t)(observed_width / dt + 0.5);
			//float nbins = observed_width / dt;
			// TESTING
			float nbins = observed_width / std::min(dt, t_i);
			
			//cout << giants[g].DM << "\t"
			//     << giants[g].width << "\t"
			//     << t_dm / observed_width * 100 << endl;
			
			add_giant(time_series, nsamps, nchans, c,
			          giants[g].time,
			          giants[g].width,
			          // Note: We normalise the flux such that the given
			          //         values matches what will be detected under
			          //         ideal circumstances.
			          giants[g].flux / sqrt(nchans) * sqrt(nbins),
			          giants[g].DM,
			          dt, f0, df);
		}
	}
	
	/*
	  for( size_t n=0; n<narrowband_rfi_count; ++n ) {
	  add_narrowband_rfi(&filterbank[narrowband_rfi_chan[n]*stride],
	  nsamps, narrowband_rfi_val[n]);
	  }
	*/
	
	double sum    = 0;
	double sum_sq = 0;
	
	cout << "Quantising filterbank data..." << endl;
	for( size_t c=0; c<nchans; ++c ) {
		for( size_t t=0; t<nsamps; ++t ) {
			float  flux = filterbank[c*stride + t];
			flux = std::min(std::max(flux, min_val), max_val);
			double  amp  = (flux - min_val) / (max_val - min_val);
			size_t quant = (size_t)(amp * (1<<nbits));
			//size_t quant = (size_t)(amp * ((1<<nbits)-1) + 0.5);
			quant &= (1<<nbits)-1;
			// TESTING
			sum    += quant;
			double m = 0.5 * ((1<<nbits)-1);
			sum_sq += (quant-m)*(quant-m);
			
			size_t w = c / chans_per_word;
			size_t k = c % chans_per_word;
			//cout << quant << endl;
			quantised[t*ostride + w] |= quant << (k*nbits);
			/*
			if( k > 0 ) {
				cout << quantised[t*ostride + w] << endl;
				cout << "[Paused]" << endl;
				std::cin.get();
			}
			*/
		}
	}
	
	double mean = sum / (nchans*nsamps) / ((1<<nbits)-1);
	double rms  = sqrt(sum_sq / (nchans*nsamps)) * 2 / ((1<<nbits)-1);
	double absrms  = sqrt(sum_sq / (nchans*nsamps)) * 2;
	cout << "Quantised mean   = " << mean << endl;
	cout << "Quantised rms    = " << rms  << endl;
	cout << "Quantised absrms = " << absrms  << endl;
	
	cout << "Writing to file..." << endl;
	std::ofstream out_file(out_filename.c_str());
	
	header_write(out_file, "HEADER_START");
	header_write(out_file, "source_name");
	header_write(out_file, "simulated_source");
	header_write(out_file, "telescope_id", 0);
	header_write(out_file, "machine_id", 0);
	header_write(out_file,
	             0.0, 0.0, 0.0, 0.0);
	header_write(out_file, "data_type", 1); // filterbank
	//header_write(out_file, "data_type", 2); // time series
	header_write(out_file, "refdm", float(0.0));
	header_write(out_file, "fch1", f0);
	header_write(out_file, "foff", df);
	//header_write(out_file, "barycentric", 0);
	header_write(out_file, "nchans", (int)nchans);
	//header_write(out_file, "nchans", (int)1);
	header_write(out_file, "nbits", (int)nbits);
	//header_write(out_file, "nbits", (int)32);
	//header_write(out_file, "tstart", 0.f);
	header_write(out_file, "tsamp", dt);
	//header_write(out_file, "tsamp", float(1));//dt);
	header_write(out_file, "nifs", 1);
	header_write(out_file, "HEADER_END");
	
	out_file.write((char*)&quantised[0], ostride*nsamps*sizeof(word_type));
	//out_file.write((char*)&filterbank[0], nsamps*sizeof(float));
	out_file.close();
	
	cout << "Done." << endl;
}
