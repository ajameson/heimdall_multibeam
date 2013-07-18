/***************************************************************************
 *
 *   Copyright (C) 2012 by Ben Barsdell and Andrew Jameson
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include "hd/header.h"
#include "hd/SigprocFile.h"

SigprocFile::SigprocFile (const char* filename)
  : m_file_stream(filename, std::ios::binary), DataSource ()
{
  m_error = 0;

  if ( m_file_stream.fail() )
  {
    cerr << "ERROR: Failed to open file '" << filename << "'" << endl;
    m_error = -1;
  }

  SigprocHeader m_header;
  read_header(m_file_stream, m_header);
  if( m_file_stream.fail() )
  {
    cerr << "ERROR: Failed to read from file '" << filename << "'" << endl;
    m_error = -2;
  }

  nchan = m_header.nchans;
  nbit  = m_header.nbits;
  beam  = m_header.ibeam;
  tsamp = m_header.tsamp;  // in seconds
  spectra_rate = 1 / tsamp;

	f0 = m_header.fch1;
	df = m_header.foff;

  utc_start = mjd2utctm (m_header.tstart);
  int buffer_size = 64;
  char buffer[buffer_size];
  strftime (buffer, buffer_size, HD_TIMESTR, localtime (&utc_start));

  stride = (nchan * nbit) / (8 * sizeof(char));

}

SigprocFile::~SigprocFile()
{ 
  m_file_stream.close();
}

size_t SigprocFile::get_data(size_t nsamps, char* data) 
{
  if ( this->get_error() ) {
    return 0;
  }
  size_t nchan_bytes = stride;
  m_file_stream.read((char*)&data[0], nsamps * nchan_bytes);
  size_t bytes_read = m_file_stream.gcount();
  return bytes_read / nchan_bytes;
};

