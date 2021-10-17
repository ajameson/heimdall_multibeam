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
#include <stdio.h>
#include <inttypes.h>

#include "hd/parse_command_line.h"
#include "hd/default_params.h"
#include "hd/pipeline.h"
#include "hd/error.h"

// input formats supported
#include "hd/DataSource.h"
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

  int i;

  if (ok < 0)
    return 1;
  
#ifdef HAVE_PSRDADA
  PSRDadaRingBuffer * data_source  = new PSRDadaRingBuffer(params.dada_id);

  if ( params.dada_id != 0 ) 
  {
    if (params.verbosity)
      cerr << "Creating PSRDADA client" << endl;

    // Read from psrdada ring buffer
    if( !data_source || data_source->get_error() ) {
      cerr << "ERROR: Failed to initialise connection to psrdada" << endl;
      return -1;
    }

    if (params.verbosity)
      cerr << "Connecting to ring buffer" << endl;
    // connect to PSRDADA ring buffer
    if (! data_source->connect())
    {
       cerr << "ERROR: Failed to connection to psrdada ring buffer" << endl;
      return -1;
    }

    if (params.verbosity)
      cerr << "Waiting for next header / data" << endl;

    // wait for and then read next PSRDADA header/observation
    if (! data_source->read_header())
    {
       cerr << "ERROR: Failed to connection to psrdada ring buffer" << endl;
      return -1;
    }

    if (!params.override_beam)
      params.beam = data_source->get_beam() - 1;
    params.nbeams = data_source->get_nbeams();
    nsamps_gulp = params.nbeams * data_source->get_nsamps_block();
  }
  else 
#endif
  {
    cerr << "ERROR: only PSRDADA data source supported" << endl;
    return -1;
  }

  if (params.nbeams == 1) 
  {
    cerr << "ERROR: only multiple beam processing supported" << endl;
    return -1;
  }

  if (!params.override_beam)
    if (data_source->get_beam() > 0)
      params.beam = data_source->get_beam() - 1;
    else
      params.beam = 0;

  if (  params.verbosity > 0)
    cout << "Have nbeams in heimdall " << params.nbeams << endl;

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

  // Preallocate memory in the pipeline
  // --------------------------
  error = hd_preallocate(pipeline, nsamps_gulp, params.nbeams);
  if( error != HD_NO_ERROR ) {
    cerr << "ERROR: Pipeline pre-allocation failed" << endl;
    return -1;
  }
  
  size_t max_overlap = hd_get_max_overlap (pipeline, params.nbeams);
  size_t filter_process_size = (nsamps_gulp + max_overlap) * stride;
  if ( params.verbosity >= 1 )
    cout << "allocating filterbank_process(" << filter_process_size << ")" << endl;
  std::vector<hd_byte> filterbank_process(filter_process_size);

  if( params.verbosity >= 1 ) {
    cout << "Beginning data processing, requesting " << nsamps_gulp << " samples" << endl;
  }

  // start a timer for the whole pipeline
  //Stopwatch pipeline_timer;

  // acquire the first block of data from the ring buffer
  // read the first block o
  size_t total_nsamps = 0;
  char * filterbank = NULL;
  size_t nsamps_read = data_source->open_data_block (&filterbank);
  
  // for the first iteration only, can use the raw filterbank
  hd_byte * fb_curr = (hd_byte *) filterbank;
  hd_byte * fb_next = &filterbank_process[0];

  data_source->close_data_block (nsamps_read);
  size_t overlap = 0;

  while( nsamps_read && !stop_requested )
  {
    // if processing data with more than one beam
    hd_size nsamps_processed;
    hd_size nsamps_to_process;
    nsamps_to_process = nsamps_gulp + (overlap * params.nbeams);

    if ( params.verbosity >= 1 ) 
    {
      cout << "Executing pipeline on new gulp of " << nsamps_read
           << " samples, nsamps_to_process=" << nsamps_to_process << endl;
    }     
    if (params.verbosity >= 2)
      fprintf (stderr, "filterbank=%p, nsamps_to_process=%d, nbits=%d "
               "total_nsamps=%d params.nbeams=%d\n", (void *) fb_curr,
               nsamps_to_process, nbits, total_nsamps, params.nbeams);

    // overlap is overlap per beam
    error = hd_execute(pipeline, &fb_curr[0], 
                       nsamps_to_process, nbits, total_nsamps, 
                       params.nbeams, &nsamps_processed);
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

    if (params.verbosity >= 1)
      cout << "Main: nsamps_processed=" << nsamps_processed << endl;

    // determine the per_beam overlap
    overlap = (nsamps_to_process-nsamps_processed)/params.nbeams;
    // find total nsamps per beam
    total_nsamps += nsamps_processed/params.nbeams;

    // check that the containers are large enough to contain the nsamps_gulp + overlap

    // copy overlaps in filterbank_process
    if (params.verbosity >= 1)
      cout << "Rewinding and reading new block of data, overlap=" << overlap << endl;

    // copy the overlapping data from the current buffer to the next buffer
    //   stride is the width of 1 sample from 1 beam (i.e. nchan * nbit)
    //   
    // for each beam (STF ordered)
    //   copy from = nchan * (ibeam+1) * 

    // number of samples per beam in the input (fb_curr)
    const size_t nsamps_to_process_per_beam = nsamps_to_process / params.nbeams;

    // number of samples per beam processed in last iteration
    const size_t nsamps_processed_per_beam  = nsamps_to_process_per_beam - overlap;

    // number of samples per beam in the output (fb_next)
    const size_t nsamps_to_process_next_per_beam = overlap + (nsamps_gulp / params.nbeams);

    // strides for each beam 
    const size_t curr_beam_stride = nsamps_to_process_per_beam * stride;
    const size_t next_beam_stride = nsamps_to_process_next_per_beam * stride;

    // copy from the curr to next filterbank pointers
    for (i=0;i<params.nbeams;i++)
    {
      const size_t curr_from = i * nsamps_to_process_per_beam + nsamps_processed_per_beam;
      const size_t curr_to   = i * nsamps_to_process_per_beam + nsamps_to_process_per_beam;
      const size_t next_from = i * nsamps_to_process_next_per_beam;
#ifdef _DEBUG
      cerr << "[" << i << "] std::copy(" << curr_from << ", " << curr_to << ", " << next_from << ")" << endl;
#endif
      // note: on the first iteration, fb_curr points to the raw ring buffer,
      //       on subsequent iterations, fb_curr and fb_next point to the same buffer
      std::copy (fb_curr + (stride * curr_from),
                 fb_curr + (stride * curr_to),
                 fb_next + (stride * next_from));
    }
   
    // obtain a pointer to the next block of data
    nsamps_read = data_source->open_data_block (&filterbank);

    // check that a full block could be read
    if (nsamps_read < nsamps_gulp)
    {
      stop_requested = 1;
    }
    // pad out the next FB buffer with new data
    else
    {
      // copy new data into filterbank_process
      for (i=0;i<params.nbeams;i++) 
      {
        const size_t fb_from = (i * nsamps_gulp) / params.nbeams;
        const size_t fb_to   = ((i+1) * nsamps_gulp) / params.nbeams;
        const size_t next_from = overlap + (i * (nsamps_gulp/params.nbeams+overlap));

#ifdef _DEBUG
        cerr << "[" << i << "] std::copy(" << fb_from << ", " << fb_to << ", " << next_from << endl;
#endif
        std::copy (filterbank + (stride * fb_from),
                   filterbank + (stride * fb_to),
                   fb_next + (stride * next_from));
      }
    }

    data_source->close_data_block (nsamps_read);

    // after the first iteration fb_curr now points to fb_next
    fb_curr = fb_next;

    if (params.verbosity >= 1)
      cerr << "end of loop nsamps_read=" << nsamps_read << " nsamps_gulp=" << nsamps_gulp << endl;

    // at the end of data, never execute the pipeline
    if (nsamps_read < nsamps_gulp)
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
