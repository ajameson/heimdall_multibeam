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
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>

#include "hd/parse_command_line.h"
#include "hd/default_params.h"
#include "hd/pipeline.h"
#include "hd/error.h"

// input formats supported
#include "hd/DataSource.h"
#include "hd/SigprocFile.h"
#ifdef HAVE_PSRDADA
#include "hd/PSRDadaRingBuffer.h"
#endif

#include "hd/stopwatch.h"

int main(int argc, char* argv[]) 
{
  hd_params params;
  hd_set_default_params(&params);
  int ok = hd_parse_command_line(argc, argv, &params);
  size_t nsamps_gulp = params.nsamps_gulp;

  if (ok < 0)
    return 1;
  
  DataSource* data_source = 0;

#ifdef HAVE_PSRDADA
  if( params.dada_id != 0 ) {

    if (params.verbosity)
      cerr << "Createing PSRDADA client" << endl;

    PSRDadaRingBuffer * d = new PSRDadaRingBuffer(params.dada_id);

    // Read from psrdada ring buffer
    if( !d || d->get_error() ) {
      cerr << "ERROR: Failed to initialise connection to psrdada" << endl;
      return -1;
    }

    if (params.verbosity)
      cerr << "Connecting to ring buffer" << endl;
    // connect to PSRDADA ring buffer
    if (! d->connect())
    {
       cerr << "ERROR: Failed to connection to psrdada ring buffer" << endl;
      return -1;
    }

    if (params.verbosity)
      cerr << "Waiting for next header / data" << endl;

    // wait for and then read next PSRDADA header/observation
    if (! d->read_header())
    {
       cerr << "ERROR: Failed to connection to psrdada ring buffer" << endl;
      return -1;
    }

    data_source = (DataSource *) d;
    if (!params.override_beam)
      params.beam = d->get_beam() - 1;
  }
  else 
  {
#endif
    // Read from filterbank file
    data_source = new SigprocFile(params.sigproc_file);
    if( !data_source || data_source->get_error() ) {
      cerr << "ERROR: Failed to open data file" << endl;
      return -1;
    }
#ifdef HAVE_PSRDADA
  }
#endif

  if (!params.override_beam)
    if (data_source->get_beam() > 0)
      params.beam = data_source->get_beam() - 1;
    else
      params.beam = 0;

  params.f0 = data_source->get_f0();
  params.df = data_source->get_df();
  params.dt = data_source->get_tsamp();

  if ( params.verbosity > 0)
    cout << "processing beam " << (params.beam+1)  << endl;

  float tsamp = data_source->get_tsamp() / 1000000;
  size_t stride = data_source->get_stride();
  size_t nbits  = data_source->get_nbit();

  params.nchans = data_source->get_nchan();
  params.utc_start = data_source->get_utc_start();
  params.spectra_per_second = data_source->get_spectra_rate();

  if ( params.verbosity >= 2)
    cout << "allocating filterbank data vector for " << nsamps_gulp << " samples with size " << (nsamps_gulp * stride) << " bytes" << endl;
  std::vector<hd_byte> filterbank(nsamps_gulp * stride);
  
  bool stop_requested = false;
  
  // Create the pipeline object
  // --------------------------
  hd_pipeline pipeline;
  hd_error error;
  error = hd_create_pipeline(&pipeline, params);
  if( error != HD_NO_ERROR ) {
    cerr << "ERROR: Pipeline creation failed" << endl;
    cerr << "       " << hd_get_error_string(error) << endl;
    return -1;
  }
  // --------------------------
  
  if( params.verbosity >= 1 ) {
    cout << "Beginning data processing, requesting " << nsamps_gulp << " samples" << endl;
  }


  // start a timer for the whole pipeline
  //Stopwatch pipeline_timer;

  size_t total_nsamps = 0;
  size_t nsamps_read = data_source->get_data (nsamps_gulp, (char*)&filterbank[0]);
  size_t overlap = 0;
  while( nsamps_read && !stop_requested ) {
    
    if ( params.verbosity >= 1 ) {
      cout << "Executing pipeline on new gulp of " << nsamps_read
           << " samples..." << endl;
    }
    //pipeline_timer.start();
      
    hd_size nsamps_processed;
    error = hd_execute(pipeline, &filterbank[0], nsamps_read+overlap, nbits,
                       total_nsamps, &nsamps_processed);
    if (error == HD_NO_ERROR)
    {
      if (params.verbosity >= 1)
        cout << "Processed " << nsamps_processed << " samples." << endl;
    }
    else if (error == HD_TOO_MANY_EVENTS) 
    {
      if (params.verbosity >= 1)
        cerr << "WARNING: hd_execute produces too many events, some data skipped" << endl;
    }
    else 
    {
      cerr << "ERROR: Pipeline execution failed" << endl;
      cerr << "       " << hd_get_error_string(error) << endl;
      hd_destroy_pipeline(pipeline);
      return -1;
    }

    //pipeline_timer.stop();
    //cout << "pipeline time: " << pipeline_timer.getTime() << " of " << (nsamps_read+overlap) * tsamp << endl;
    //pipeline_timer.reset();

    total_nsamps += nsamps_processed;
    // Now we must 'rewind' to do samples that couldn't be processed
    // Note: This assumes nsamps_gulp > 2*overlap
    std::copy(&filterbank[nsamps_processed * stride],
              &filterbank[(nsamps_read+overlap) * stride],
              &filterbank[0]);
    overlap += nsamps_read - nsamps_processed;
    nsamps_read = data_source->get_data(nsamps_gulp - overlap,
                                        (char*)&filterbank[overlap*stride]);

    // at the end of data, never execute the pipeline
    if (nsamps_read < nsamps_gulp - overlap)
      stop_requested = 1;

  }
    
  if( params.verbosity >= 1 ) {
    cout << "Successfully processed a total of " << total_nsamps
         << " samples." << endl;
  }
    
  if( params.verbosity >= 1 ) {
    cout << "Shutting down..." << endl;
  }
  
  hd_destroy_pipeline(pipeline);
  
  if( params.verbosity >= 1 ) {
    cout << "All done." << endl;
  }
}
