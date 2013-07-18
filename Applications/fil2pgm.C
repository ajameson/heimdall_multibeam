#include <iostream>
using std::cout;
using std::endl;
#include <fstream>
#include <string>
using std::string;
#include <vector>
#include <cmath>
#include <cstdlib>

// SIGPROC header
#include "hd/header.h"

template<int NBITS, typename T=unsigned int>
struct max_value {
  static const T value = (((unsigned)1<<(NBITS-1))-1)*2+1;
};


typedef unsigned int word_type;
typedef unsigned int size_type;
typedef unsigned char out_type;

int main(int argc, char* argv[])
{
  if( argc <= 2 ) {
    cout << "Usage: " << argv[0] << " infile outfile" << " [nsamps] [skip] [tscrunch] [fscrunch]" << endl;
    return 0;
  }
  string in_filename  = argv[1];
  string out_filename = argv[2];
  size_type tscrunch = 1;
  size_type fscrunch = 1;
  size_type out_nsamps = 4096;
  size_type skip = 0;
  if( argc > 3 ) {
    out_nsamps = atoi(argv[3]);
  }
  if( argc > 4 ) {
    skip = atoi(argv[4]);
  }
  if( argc > 5 ) {
    // Note: Must divide nsamps evenly
    tscrunch = atoi(argv[5]);
  }
  if( argc > 6 ) {
    // Note: Must divide nchans evenly
    fscrunch = atoi(argv[6]);
  }

  if( out_nsamps % tscrunch != 0 ) {
    out_nsamps = out_nsamps/tscrunch*tscrunch;
    cout << "Warning: Adjusting nsamps to multiple of tscrunch: "
   << out_nsamps << endl;
  }

  std::ifstream in_file(in_filename.c_str());

  SigprocHeader header;
  read_header(in_file, header);

  if( header.nchans % fscrunch != 0 ) {
    cout << "Error: fscrunch (" << fscrunch
   << ") does not divide nchans ("
   << header.nchans << ")" << endl;
  }

  size_type chans_per_word = sizeof(word_type)*8/header.nbits;
  size_type mask       = (((unsigned)1<<(header.nbits-1))-1)*2+1;

  //out_nsamps = header.nsamples;

  cout << "Reading file..." << endl;
  size_type stride_words = header.nchans/chans_per_word;
  size_type stride_bytes = header.nchans*sizeof(word_type)/chans_per_word;
  // Skip some samples if requested
  in_file.seekg(skip*stride_bytes, std::ios::cur);
  //std::vector<word_type> packed_data(header.nchans*header.nsamples/chans_per_word);
  //in_file.read((char*)&packed_data[0], header.nchans*header.nsamples/chans_per_word*sizeof(word_type));
  std::vector<word_type> packed_data(out_nsamps*stride_words);
  in_file.read((char*)&packed_data[0], out_nsamps*stride_bytes);
  in_file.close();

  std::vector<out_type> unpacked_data(header.nchans/fscrunch*(out_nsamps/tscrunch));

  //float peak = 255 / sqrt(fscrunch);
  float peak = 255;

  cout << "Unpacking " << header.nbits << "-bit data..." << endl;
  for( size_type t=0; t<(size_type)out_nsamps; t+=tscrunch ) {
    for( size_type c=0; c<(size_type)header.nchans; c+=fscrunch ) {
      float sum = 0.f;
      for( size_type s=0; s<(size_type)tscrunch; ++s ) {
  for( size_type f=0; f<(size_type)fscrunch; ++f ) {
    size_type w = (c+f) / chans_per_word;
    size_type k = (c+f) % chans_per_word;

    word_type x = (packed_data[w + (t+s)*(header.nchans/chans_per_word)]
       >> (k*header.nbits)) & mask;
    sum += (float)x / mask * peak;
    //unpacked_data[c + t*header.nchans] = (float)x / mask * 255;
  }
  //tsum += fsum / fscrunch;
  //unpacked_data[c/fscrunch + t*(header.nchans/fscrunch)] = fsum / fscrunch;
      }
      unpacked_data[c/fscrunch + (t/tscrunch)*(header.nchans/fscrunch)] =
  sum / (tscrunch*fscrunch);
    }
  }

  cout << "Writing output..." << endl;
  std::ofstream out_file(out_filename.c_str());
  out_file << "P5 " << header.nchans/fscrunch << " "
     << out_nsamps/tscrunch << " " << 255 << endl;
  out_file.write((char*)&unpacked_data[0], unpacked_data.size()*sizeof(out_type));
  out_file.close();

  cout << "Done." << endl;
  return 0;
}

