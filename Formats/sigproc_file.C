
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include "hd/header.h"
#include "hd/sigproc_file.h"

#include "tmutil.h"
#include "dada_def.h"

FileDataSource::FileDataSource (const char* filename)
  : m_file_stream(filename, std::ios::binary)
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

  _nchan = m_header.nchans;
  _nbit  = m_header.nbits;
  _beam  = m_header.ibeam;
  _tsamp = m_header.tsamp;  // in seconds
  _spectra_per_second = 1 / _tsamp;

  _utc_start = mjd2utctm (m_header.tstart);
  int buffer_size = 64;
  char buffer[buffer_size];
  strftime (buffer, buffer_size, DADA_TIMESTR, localtime (&_utc_start));

  _stride = (_nchan * _nbit) / (8 * sizeof(char));

}

FileDataSource::~FileDataSource()
{ 
  m_file_stream.close();
}

size_t FileDataSource::get_data(size_t nsamps, char* data) 
{
  if ( this->error() ) {
    return 0;
  }
  size_t nchan_bytes = _stride;
  m_file_stream.read((char*)&data[0], nsamps*nchan_bytes);
  size_t bytes_read = m_file_stream.gcount();
  return bytes_read / nchan_bytes;
};

