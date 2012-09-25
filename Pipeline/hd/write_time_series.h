/***************************************************************************
 *
 *   Copyright (C) 2012 by Ben Barsdell and Andrew Jameson
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include <vector>
#include <string>
using std::string;
#include <fstream>

namespace detail {
// TODO: These were copied from header.hpp. Not sure if this is a good idea.
// Write string value
template<class String, class BinaryStream>
void header_write(BinaryStream& stream, const String& str) {
	std::string s = str;
	int len = s.size();
	// TODO: Apply byte swapping for endian-correctness
	stream.write((char*)&len, sizeof(int));
	// TODO: Apply byte swapping for endian-correctness
	stream.write(s.c_str(), len*sizeof(char));
}

// Write integer value
template<class String, class BinaryStream>
void header_write(BinaryStream& stream, String name, int val) {
	header_write(stream, name);
	// TODO: Apply byte swapping for endian-correctness
	stream.write((char*)&val, sizeof(int));
}

// Write floating-point value
template<class String, class BinaryStream>
void header_write(BinaryStream& stream, String name, double val) {
	header_write(stream, name);
	// TODO: Apply byte swapping for endian-correctness
	stream.write((char*)&val, sizeof(double));
}

// Write coordinates
template<class BinaryStream>
void header_write(BinaryStream& stream,
				  double raj, double dej,
				  double az, double za) {
	header_write(stream, "src_raj",  raj);
	header_write(stream, "src_dej",  dej);
	header_write(stream, "az_start", az);
	header_write(stream, "za_start", za);
}

void write_time_series_header(size_t nbits, float dt, std::ofstream& out_file) {
	// Write the required header information
	header_write(out_file, "HEADER_START");
	//header_write(out_file, "telescope_id", header.telescope_id);
	//header_write(out_file, "machine_id", header.machine_id);
	//header_write(out_file,
	//			 header.src_raj, header.src_dej,
	//			 header.az_start, header.za_start);
	header_write(out_file, "data_type", 2);
	//header_write(out_file, "refdm", header.refdm);
	//header_write(out_file, "fch1", header.f0);
	//header_write(out_file, "barycentric", 6);//header.barycentric);
	header_write(out_file, "nchans", 1);//nbands);
	header_write(out_file, "nbits", (int)nbits);
	//header_write(out_file, "tstart", 0.f);//header.tstart);
	header_write(out_file, "tsamp", dt);
	header_write(out_file, "nifs", 1);//header.nifs);
	header_write(out_file, "HEADER_END");
}
} // namespace detail

// Float data type
void write_host_time_series(const float* data,
							size_t       nsamps,
                            float        dt,
							string       filename)
{
	// Open the output file and write the data
	std::ofstream file(filename.c_str(), std::ios::binary);
	detail::write_time_series_header(32, dt, file);
	size_t size_bytes = nsamps*sizeof(float);
	file.write((char*)data, size_bytes);
	file.close();
}
void write_device_time_series(const float* data,
                              size_t      nsamps,
                              float       dt,
                              string      filename)
{
	std::vector<float> h_data(nsamps);
	cudaError_t error = cudaMemcpy(&h_data[0], data,
	                               nsamps*sizeof(float),
								   cudaMemcpyDeviceToHost);
	if( error != cudaSuccess ) {
		throw std::runtime_error(std::string("write_time_series: cudaMemcpy failed: ") +
								 cudaGetErrorString(error));
	}
	write_host_time_series(&h_data[0], nsamps, dt, filename);
}

// Integer data type
void write_host_time_series(const unsigned int* data,
                            size_t      nsamps,
                            size_t      nbits,
                            float       dt,
                            string      filename)
{
	// Here we convert the data to floats before writing the data
	std::vector<float> float_data(nsamps);
	switch( nbits ) {
	case sizeof(char)*8:
		for( int i=0; i<(int)nsamps; ++i )
			float_data[i] = (float)((unsigned char*)data)[i];
		break;
	case sizeof(short)*8:
		for( int i=0; i<(int)nsamps; ++i )
			float_data[i] = (float)((unsigned short*)data)[i];
		break;
	case sizeof(int)*8:
		for( int i=0; i<(int)nsamps; ++i )
			float_data[i] = (float)((unsigned int*)data)[i];
		break;
	/*
	case sizeof(long long)*8:
		for( int i=0; i<(int)nsamps; ++i )
			float_data[i] = (float)((unsigned long long*)data)[i];
	*/
	default:
		// Unpack to float
		size_t mask = (1 << nbits) - 1;
		size_t spw = sizeof(unsigned int)*8 / nbits; // Samples per word
		for( int i=0; i<(int)nsamps; ++i )
			float_data[i] = (data[i/spw] >> (i % spw * nbits)) & mask;
	}
	write_host_time_series(&float_data[0], nsamps, dt, filename);
}
void write_device_time_series(const unsigned int* data,
                              size_t      nsamps,
                              size_t      nbits,
                              float       dt,
                              string      filename)
{
	size_t nsamps_words = nsamps * nbits/(sizeof(unsigned int)*8);
	std::vector<unsigned int> h_data(nsamps_words);
	cudaError_t error = cudaMemcpy(&h_data[0], data,
	                               nsamps_words*sizeof(unsigned int),
								   cudaMemcpyDeviceToHost);
	if( error != cudaSuccess ) {
		throw std::runtime_error(std::string("write_time_series: cudaMemcpy failed: ") +
								 cudaGetErrorString(error));
	}
	write_host_time_series(&h_data[0], nsamps, nbits, dt, filename);
}
