/***************************************************************************
 *
 *   Copyright (C) 2012 by Andrew Jameson
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <errno.h>
#include <string.h>

#include "hd/PSRDadaRingBuffer.h"

#include "ascii_header.h"
#include "tmutil.h"

//#define _DEBUG

PSRDadaRingBuffer::PSRDadaRingBuffer (key_t dada_id) : DataSource ()
{
  dada_error = false;
  printed_first_line = true;
  resolution = 1;
  curr_block = 0;

  dada_key = dada_id;

  // create the HDU
  hdu = dada_hdu_create (0);

  // set the shared memory key
  dada_hdu_set_key(hdu, dada_id);
}

PSRDadaRingBuffer::~PSRDadaRingBuffer ()
{
  if (hdu)
  {
		dada_hdu_unlock_write (hdu);
    dada_hdu_disconnect(hdu);
	}
	hdu = 0;
}


bool PSRDadaRingBuffer::connect()
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
  curr_block = 0;
  bytes_read = 0;
	return true;
}

bool PSRDadaRingBuffer::disconnect()
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

bool PSRDadaRingBuffer::read_header()
{
  uint64_t header_size = 0;
  char * header;

  // Wait for the next valid header block
  header = ipcbuf_get_next_read (hdu->header_block, &header_size);
  if (!header)
  {
    cerr << "PSRDadaRingBuffer::read_header could not get next header" << endl;
    dada_error = true;
    return false;
  }

#ifdef DEBUG
  cerr << "==================================================" << endl;
  cerr << header << endl;
  cerr << "==================================================" << endl;
#endif

  // get the required header params
  if (ascii_header_get (header, "NCHAN", "%d", &nchan) < 0)
  {
    cerr << "PSRDadaRingBuffer::read_header could not extract NCHAN from header" << endl;
    dada_error = true;
  }

	float bandwidth = 0;
  if (ascii_header_get (header, "BANDWIDTH", "%f", &bandwidth) < 0)
  {
    if (ascii_header_get (header, "BW", "%f", &bandwidth) < 0)
    {
      cerr << "WARNING: PSRDadaRingBuffer::read_header could not extract BANDWIDTH from header" << endl;
      dada_error = true;
    }
  }

	float cfreq = 0;
  if (ascii_header_get (header, "CFREQ", "%f", &cfreq) < 0)
  {
    if (ascii_header_get (header, "FREQ", "%f", &cfreq) < 0)
    {
      cerr << "WARNING: PSRDadaRingBuffer::read_header could not extract CFREQ from header" << endl;
      dada_error = true;
    }
  }

	if (!dada_error)
	{
		float start_freq = cfreq - (bandwidth / 2);
		float chan_width = bandwidth / nchan;
		float first_chan_cfreq = start_freq + (chan_width / 2);

		// set the cfreq of first chan and the delta freq between chans
		f0 = first_chan_cfreq;
		df = chan_width;
	}

  if (ascii_header_get (header, "NPOL", "%d", &npol) < 0)
  {
    cerr << "PSRDadaRingBuffer::read_header could not extract NPOL from header" << endl;
    dada_error = true; 
  }

  if (ascii_header_get (header, "NBIT", "%d", &nbit) < 0)
  {
    cerr << "PSRDadaRingBuffer::read_header could not extract NBIT from header" << endl;
    dada_error = true; 
  }
#ifdef _DEBUG
  cerr << "PSRDadaRingBuffer::read_header nbit=" << nbit << endl;
#endif

  if (ascii_header_get (header, "TSAMP", "%f", &tsamp) < 0)
  {
    cerr << "PSRDadaRingBuffer::read_header could not extract TSAMP from header" << endl;
    dada_error = true;
  }

  if (ascii_header_get (header, "NBEAM", "%d", &nbeams) < 0)
  {
    cerr << "PSRDadaRingBuffer::read_header could not extract NBEAM from header, assuming 1" << endl;
    nbeams = 1;
  }

  if (ascii_header_get (header, "BEAM", "%d", &beam) < 0)
  {
    if (nbeams == 1)
      cerr << "PSRDadaRingBuffer::read_header could not extract BEAM from header, assuming 0" << endl;
    beam = 0;
  }

  if (ascii_header_get (header, "RESOLUTION", "%u", &resolution) < 0)
  {
    resolution = nchan * npol * (nbit/8);
    cerr << "PSRDadaRingBuffer::read_header could not extract RESOLUTION from header, using " << resolution << endl;
  }

  if (ascii_header_get (header, "ORDER", "%s", order) < 0)
  {
    cerr << "PSRDadaRingBuffer::read_header could not extract ORDER from header" << endl;
    dada_error = true;
  }

  if ((nbeams == 1) && (strcmp(order, "TF") != 0))
  {
    cerr << "PSRDadaRingBuffer::read_header 1 beam but data order not TF" << endl;
    dada_error = true;
  }

  if ((nbeams > 1) && (strcmp(order, "STF") != 0))
  {
    cerr << "PSRDadaRingBuffer::read_header nbeam=" << nbeams << " but data order not STF" << endl;
    dada_error = true;
  }

  nsamps_block = (size_t) resolution / (nbeams * nchan * npol * (nbit/8));
  //cerr << "PSRDadaRingBuffer::read_header resolution=" << resolution << " nsamps_block=" << nsamps_block << endl;

  char utc_start_str[64] = "";
  if (ascii_header_get (header, "UTC_START", "%s", utc_start_str) < 0)
  {
    cerr << "PSRDadaRingBuffer::read_header could not extract UTC_START from header" << endl;
    dada_error = true;
  }
  else
    utc_start = str2utctime (utc_start_str);

	stride = nchan * (nbit / 8);
  spectra_rate = 1000000 / (double) tsamp;

  // convert tsamp from usecs (DADA DEFAULT) to seconds
  tsamp /= 1000000;

  // get the DADA buffer size
  buf_sz = ipcbuf_get_bufsz ( (ipcbuf_t *) hdu->data_block);

  // get the samples per DADA buffer
  samples_per_buf = buf_sz / (npol * nchan * (nbit/8));

  if (buf_sz != resolution && nbeams > 1)
  {
    cerr << "PSRDadaRingBuffer::read_header buf_sz != resolution and multibeam data" << endl;
    dada_error = true;
  }

#ifdef _DEBUG
  cerr << "PSRDadaRingBuffer::read_header utc_start_str=" << utc_start_str << endl;
  cerr << "PSRDadaRingBuffer::read_header utc_start=" << utc_start << endl;
  cerr << "PSRDadaRingBuffer::read_header tsamp=" << tsamp << endl;
  cerr << "PSRDadaRingBuffer::read_header spectra_rate =" << spectra_rate << endl;
  cerr << "PSRDadaRingBuffer::read_header resolution=" << resolution << endl;
  cerr << "PSRDadaRingBuffer::read_header nbeams=" << nbeams << endl;
  cerr << "PSRDadaRingBuffer::read_header nsamps_block=" << nsamps_block << endl;
#endif

  if (dada_error)
    return false;

  // if we sucessfully parsed the header, mark buffer cleared
  ipcbuf_mark_cleared (hdu->header_block);
  return true;
}

size_t PSRDadaRingBuffer::get_data(size_t nsamps, char* data)
{
#ifdef _DEBUG
	cerr << "PSRDadaRingBuffer::get_data: nsamps=" << nsamps << " nbeams=" << nbeams << endl;
#endif

  if (nbeams > 1)
  {
    return get_data_block (nsamps, data);
  }

  uint64_t bytes_to_read = nsamps * nchan * (nbit / 8);
  int64_t  bytes_read = 0;
  
  if (ipcbuf_eod((ipcbuf_t*)hdu->data_block))
    return 0;

#ifdef _DEBUG
	cerr << "PSRDadaRingBuffer::get_data: ipcio_read for " << bytes_to_read << " bytes" << endl;
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
    cerr << "PSRDadaRingBuffer::get_data: ipcio_read error: " <<  strerror(errno) << endl;
    dada_error = true;
    return -1;
  }
  
	// actually return the number of samples read, rounding down
	size_t nsamps_read = bytes_read / (nchan * (nbit / 8));

	if (nsamps_read != nsamps)
    if (!ipcbuf_eod((ipcbuf_t*)hdu->data_block))
		  cerr << "PSRDadaRingBuffer::get_data: returing fewer nsamps than requested!" << endl;

  return nsamps_read;
}

// read a single blocks worth of data, nsamps MUST be equal to block size 
size_t PSRDadaRingBuffer::get_data_block (size_t nsamps, char* data)
{
  size_t bytes_to_read = nsamps * nchan * (nbit/8);
  uint64_t block_id, bytes_in_block;
  unsigned ibeam;

#ifdef _DEBUG
    cerr << "PSRDadaRingBuffer::get_data_block: nsamps_requested=" << nsamps << " bytes_to_read=" << bytes_to_read << endl;
#endif

  // open block if necessary
  if (!curr_block)
  {
    curr_block = ipcio_open_block_read (hdu->data_block, &bytes_in_block, &block_id);
    if (!curr_block)
    {
      if (ipcbuf_eod((ipcbuf_t*)hdu->data_block))
      {
        cerr << "PSRDadaRingBuffer::get_data_block: EOD" << endl;
        return 0;
      }
      else
      {
        cerr << "PSRDadaRingBuffer::get_data: ipcio_open_block_read failed" << endl;
        return -1;
      }
    }
  }

#ifdef _DEBUG
    cerr << "PSRDadaRingBuffer::get_data_block: bytes_in_block=" << bytes_in_block << endl;
#endif

  if (bytes_in_block != bytes_to_read)
  {
    // cerr << "PSRDadaRingBuffer::get_data: bytes avaiable in data block less than 1 full block" << endl;
  }

  memcpy ((void *) data, curr_block, bytes_in_block);

  // close the data block
  ipcio_close_block_read (hdu->data_block, bytes_in_block);
  curr_block = 0;

  uint64_t nsamps_read = bytes_in_block / (nchan * (nbit/8));

  if (bytes_in_block < bytes_to_read)
    if (!ipcbuf_eod((ipcbuf_t*)hdu->data_block))
      cerr << "PSRDadaRingBuffer::get_data: end of data" << endl;

#ifdef _DEBUG
    cerr << "PSRDadaRingBuffer::get_data_block: nsamps_read=" << nsamps_read<< endl;
#endif

  return nsamps_read;
}

size_t PSRDadaRingBuffer::open_data_block (char ** buf_ptr)
{
  uint64_t block_id, bytes_in_block;

  if (!curr_block)
  {
    curr_block = ipcio_open_block_read (hdu->data_block, &bytes_in_block, &block_id);
    if (!curr_block)
    {
      if (ipcbuf_eod((ipcbuf_t*)hdu->data_block))
      {
        cerr << "PSRDadaRingBuffer::get_data_block: EOD" << endl;
        return 0;
      }
      else
      {
        cerr << "PSRDadaRingBuffer::get_data: ipcio_open_block_read failed" << endl;
        return -1;
      }
    }
  }
  else
  {
    cerr << "PSRDadaRingBuffer::get_data_block block already open" << endl;
    return -1;
  }

  *buf_ptr = curr_block;

  if (bytes_in_block != buf_sz)
  {
    cerr << "PSRDadaRingBuffer::open_data_block: bytes avaiable in data block [" << bytes_in_block << "] less than 1 full block [" << buf_sz << "]" << endl;
  }

  size_t nsamps_read = bytes_in_block / (nchan * (nbit / 8));
  return nsamps_read;
}

void PSRDadaRingBuffer::close_data_block (size_t nsamps_read)
{
  uint64_t bytes_in_block = nsamps_read * (nchan * (nbit / 8));

  // close the data block
  ipcio_close_block_read (hdu->data_block, bytes_in_block);
  curr_block = 0;

  if (bytes_in_block != buf_sz)
    if (!ipcbuf_eod((ipcbuf_t*)hdu->data_block))
      cerr << "PSRDadaRingBuffer::close_data_block end of data" << endl;

#ifdef _DEBUG
    cerr << "PSRDadaRingBuffer::close_data_block nsamps_read=" << nsamps_read<< endl;
#endif

}
