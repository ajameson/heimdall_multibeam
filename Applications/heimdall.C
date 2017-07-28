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

  int i;

  if (ok < 0)
    return 1;
  
  DataSource* data_source = 0;

#ifdef HAVE_PSRDADA
  if ( params.dada_id != 0 ) 
  {

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

    params.nbeams = d->get_nbeams();
    nsamps_gulp = params.nbeams * d->get_nsamps_block();
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

  if ( params.verbosity >= 2 )
    cout << "allocating filterbank data vector for " << nsamps_gulp 
         << " samples with size 2 * " << (nsamps_gulp * stride) << " bytes" << endl;
  std::vector<hd_byte> filterbank(nsamps_gulp * stride);

  if ( params.verbosity >= 2)
    cout << "allocating twice-as-large filterbank data vector..." << endl;
  std::vector<hd_byte> filterbank_process1(2 * nsamps_gulp * stride);
  std::vector<hd_byte> filterbank_process2(2 * nsamps_gulp * stride);

  hd_byte * fb_curr = &filterbank_process1[0];
  hd_byte * fb_next = &filterbank_process2[0];
  hd_byte * fb_temp = 0;
  
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

  // copy to filterbank_process if needed
  if (params.nbeams > 1) 
  {
    std::copy(&filterbank[0],&filterbank[stride*nsamps_gulp], fb_curr);
    //cerr << "copying " << (&filterbank[stride*nsamps_gulp] - &filterbank[0]) << " bytes to fb_curr" << endl;
  }
  size_t overlap = 0;
  unsigned first = 2;

  while( nsamps_read && !stop_requested )
  {
    if (params.nbeams == 1) 
    {
      if ( params.verbosity >= 1 )
      {
        cout << "Executing pipeline on new gulp of " << nsamps_read
             << " samples, total_samples=" << nsamps_read+overlap << endl;
      }
      //pipeline_timer.start();

      hd_size nsamps_processed;
      error = hd_execute(pipeline, &filterbank[0], nsamps_read+overlap, nbits,
                         total_nsamps, params.nbeams, &nsamps_processed);
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

    // if processing data with more than one beam
    else if (params.nbeams > 1) 
    {

      hd_size nsamps_processed;
      hd_size nsamps_to_process;
      nsamps_to_process = nsamps_gulp + (overlap*params.nbeams);

      if ( params.verbosity >= 1 ) 
      {
        cout << "Executing pipeline on new gulp of " << nsamps_read
             << " samples, nsamps_to_process=" << nsamps_to_process << endl;
      }     
      if (params.verbosity >= 2)
        fprintf (stderr, "filterbank=%p, nsamps_to_process=%d, nbits=%d "
                 "total_nsamps=%d params.nbeams=%d\n", &filterbank[0], 
                 nsamps_to_process, nbits, total_nsamps, params.nbeams);
      // overlap is overlap per beam
      error = hd_execute(pipeline, fb_curr, 
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

      // number of samples per veam in the output (fb_next)
      const size_t nsamps_to_process_next_per_beam = overlap + (nsamps_gulp / params.nbeams);

      // ensure containers are large enough for this
      if (filterbank_process1.size() < nsamps_to_process_next_per_beam * stride * params.nbeams)
      {
        bool fb_same = (fb_curr == &filterbank_process1[0]);

        filterbank_process1.resize (nsamps_to_process_next_per_beam * stride * params.nbeams);
        filterbank_process2.resize (nsamps_to_process_next_per_beam * stride * params.nbeams);
        if (fb_same)
        {
          fb_curr = &filterbank_process1[0];
          fb_next = &filterbank_process2[0];
        }
        else
        {
          fb_curr = &filterbank_process2[0];
          fb_next = &filterbank_process1[0];
        }
      }
     
      // strides for each beam 
      const size_t curr_beam_stride = nsamps_to_process_per_beam * stride;
      const size_t next_beam_stride = nsamps_to_process_next_per_beam * stride;

      for (i=0;i<params.nbeams;i++) 
      {
        const size_t curr_from = i * nsamps_to_process_per_beam + nsamps_processed_per_beam;
        const size_t curr_to   = i * nsamps_to_process_per_beam + nsamps_to_process_per_beam;
        const size_t next_from = i * nsamps_to_process_next_per_beam;

        const size_t to_copy = curr_to - curr_from;

#ifdef _DEBUG
        cerr << "std::copy [" << i << "] ("
             << curr_from << ", "
             << curr_to << ", " 
             << next_from << ")" << endl;
#endif
        std::copy (fb_curr + (stride * curr_from),
                   fb_curr + (stride * curr_to),
                   fb_next + (stride * next_from));
      }
     
      // read new gulp of data
      nsamps_read = data_source->get_data(nsamps_gulp,(char*)&filterbank[0]);

      // check that a full block could be read
      if (nsamps_read < nsamps_gulp)
      {
        stop_requested = 1;
      }
      // pad out the next FB buffer with new data
      else
      {
        //cerr << "copying new data from FB into FB_NEXT" << endl;
        // copy new data into filterbank_process
        for (i=0;i<params.nbeams;i++) 
        {
#ifdef _DEBUG
          cerr << "std::copy [" << i << "] ("
               << ((i)*nsamps_gulp/params.nbeams) << ", "
               << ((i+1)*nsamps_gulp/params.nbeams) << ", "
               << overlap+i*(nsamps_gulp/params.nbeams+overlap) << ")" << endl;
#endif

          std::copy(&filterbank[stride*((i)*nsamps_gulp/params.nbeams)],
                    &filterbank[stride*((i+1)*nsamps_gulp/params.nbeams)],
                    fb_next + stride*(overlap+i*(nsamps_gulp/params.nbeams+overlap)));
        }
      }

      fb_temp = fb_curr;
      fb_curr = fb_next;
      fb_next = fb_temp;

      //first--;
      if (params.verbosity >= 1)
        cerr << "end of loop nsamps_read=" << nsamps_read << " nsamps_gulp=" << nsamps_gulp << endl;

      // at the end of data, never execute the pipeline
      if (nsamps_read < nsamps_gulp)
        stop_requested = 1;
    }
  }
 
  // final iteration for nsamps which is not a multiple of gulp size - overlap
  if (stop_requested && params.nbeams==1) 
  {
    if (params.verbosity >= 1)
      cout << "Final sub gulp: nsamps_read=" << nsamps_read << " nsamps_gulp=" << nsamps_gulp << " overlap=" << overlap << endl;
    hd_size nsamps_processed;
    hd_size nsamps_to_process = nsamps_read + (overlap * 2 - params.boxcar_max);
    if (nsamps_to_process > nsamps_gulp)
      nsamps_to_process = nsamps_gulp;
    error = hd_execute(pipeline, &filterbank[0], nsamps_to_process, nbits, 
                       total_nsamps, params.nbeams, &nsamps_processed);
    if (params.verbosity >= 1)
      cout << "Final sub gulp: nsamps_processed=" << nsamps_processed << endl;

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
    }
    total_nsamps += nsamps_processed;
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
