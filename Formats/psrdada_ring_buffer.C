
#include "hd/psrdada_ring_buffer.h"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <errno.h>
#include <string.h>

#include "ascii_header.h"
#include "tmutil.h"

//#define _DEBUG

DadaDataSource::DadaDataSource (key_t dada_id)
{
  dada_error = false;
  printed_first_line = true;

  dada_key = dada_id;

  // create the HDU
  hdu = dada_hdu_create (0);

  // set the shared memory key
  dada_hdu_set_key(hdu, dada_id);

  _nchan = 0;
  _nbit  = 0;
  _npol  = 0;
  _stride = 0;
  _beam = 0;
  _utc_start = 0;
  _tsamp = 0;
  _spectra_per_second = 0;
}

DadaDataSource::~DadaDataSource ()
{
  if (hdu)
  {
		dada_hdu_unlock_write (hdu);
    dada_hdu_disconnect(hdu);
	}
	hdu = 0;
}


bool DadaDataSource::connect()
{
  if (dada_hdu_connect (hdu) < 0)
  {
    dada_error = true;
    return false;
  }

  if (dada_hdu_lock_read (hdu) < 0)
  {
    dada_error = true;
    return false;
  }
	return true;
}

bool DadaDataSource::disconnect()
{
  if (dada_hdu_unlock_read (hdu) < 0)
  {
    dada_error = true;
    return false;
  }

  if (dada_hdu_disconnect (hdu) < 0)
  {
    dada_error = true;
    return false;
  }  
	return true;
}

bool DadaDataSource::read_header()
{
  uint64_t header_size = 0;
  char * header;

  // Wait for the next valid header block
  header = ipcbuf_get_next_read (hdu->header_block, &header_size);
  if (!header)
  {
    cerr << "DadaDataSource::read_header could not get next header" << endl;
    dada_error = true;
    return false;
  }

  //cerr << "==================================================" << endl;
  //cerr << header << endl;
  //cerr << "==================================================" << endl;

  // get the required header params
  if (ascii_header_get (header, "NCHAN", "%d", &_nchan) < 0)
  {
    cerr << "DadaDataSource::read_header could not extract NCHAN from header" << endl;
    dada_error = true;
  }

  if (ascii_header_get (header, "NPOL", "%d", &_npol) < 0)
  {
    cerr << "DadaDataSource::read_header could not extract NPOL from header" << endl;
    dada_error = true; 
  }

  if (ascii_header_get (header, "NBIT", "%d", &_nbit) < 0)
  {
    cerr << "DadaDataSource::read_header could not extract NBIT from header" << endl;
    dada_error = true; 
  }

  if (ascii_header_get (header, "TSAMP", "%f", &_tsamp) < 0)
  {
    cerr << "DadaDataSource::read_header could not extract TSAMP from header" << endl;
    dada_error = true;
  }

  if (ascii_header_get (header, "BEAM", "%d", &_beam) < 0)
  {
    cerr << "DadaDataSource::read_header could not extract BEAM from header, assuming 0" << endl;
		_beam = 0;
  }

  char utc_start_str[64] = "";
  if (ascii_header_get (header, "UTC_START", "%s", utc_start_str) < 0)
  {
    cerr << "DadaDataSource::read_header could not extract UTC_START from header" << endl;
    dada_error = true;
  }
  else
    _utc_start = str2utctime (utc_start_str);

	// TODO fix the fact that we acually have _npol == 2 !!!!!!!!!!!!!!
	_stride = _nchan * (_nbit / 8);
  _spectra_per_second = 1000000 / (double) _tsamp;

#ifdef _DEBUG
  cerr << "DadaDataSource::read_header utc_start=" << _utc_start << endl;
  cerr << "DadaDataSource::read_header tsamp=" << _tsamp << endl;
  cerr << "DadaDataSource::read_header spectra_per_second=" << _spectra_per_second << endl;
#endif

  if (dada_error)
    return false;

  // if we sucessfully parsed the header, mark buffer cleared
  ipcbuf_mark_cleared (hdu->header_block);
  return true;
}

size_t DadaDataSource::get_data(size_t nsamps, char* data)
{
#ifdef _DEBUG
	cerr << "DadaDataSource::get_data: nsamps=" << nsamps << endl;
#endif

  uint64_t bytes_to_read = nsamps * _nchan * (_nbit / 8);
  int64_t  bytes_read = 0;
  
  if (ipcbuf_eod((ipcbuf_t*)hdu->data_block))
    return 0;

#ifdef _DEBUG
	cerr << "DadaDataSource::get_data: ipcio_read for " << bytes_to_read << " bytes" << endl;
#endif

  bytes_read = ipcio_read (hdu->data_block, data, bytes_to_read);

  if (!printed_first_line)
  {
    printed_first_line = true;
    unsigned char * ptr = (unsigned char *) data;
    for (unsigned i=0; i<1024; i++)
    {
      fprintf(stderr, "data[%d]=%d\n",i, ptr[i]);
    }
  }

  if (bytes_read < 0)
  {
    cerr << "DadaDataSource::get_data: ipcio_read error: " <<  strerror(errno) << endl;
    dada_error = true;
    return -1;
  }
  
	// actually return the number of samples read, rounding down
	size_t nsamps_read = bytes_read / (_nchan * (_nbit / 8));

	if (nsamps_read != nsamps)
    if (!ipcbuf_eod((ipcbuf_t*)hdu->data_block))
		  cerr << "DadaDataSource::get_data: returing fewer nsamps than requested!" << endl;

  return nsamps_read;
}
