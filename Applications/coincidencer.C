/***************************************************************************
 *
 *   Copyright (C) 2012 by Ben Barsdell and Andrew Jameson
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>

using namespace std;

#include <unistd.h>
#include <inttypes.h>
#include <stdlib.h>
#include <libgen.h>
#include <math.h>
#include <signal.h>
#include <string.h>

#include "hd/Candidates.h"
#include "hd/ServerSocket.h"
#include "hd/SocketException.h"

int quit_threads = 0;

void usage(void)
{
  fprintf(stdout, "coincidencer [options] [beam_candidate_files]\n");
  fprintf(stdout, "  -h          print this help text\n");
  fprintf(stdout, "  -a address  interface for candidate events\n");
  fprintf(stdout, "  -n nbeams   number of beams to expect data from\n");
  fprintf(stdout, "  -p port     port for candidate events\n");
  fprintf(stdout, "  -v          verbose output\n");
}

void signal_handler(int signalValue);

int main(int argc, char* argv[])
{
  int arg = 0;

  int port = 0;

  ServerSocket * server = NULL;

  // list of chunks of observations
  std::vector<CandidateChunk *> chunks;

  std::string curr_first_sample_utc;

  int curr_chunk = 0;

  bool persist = false;

  const char * address = "any";

  unsigned int total_beams = 1;

  unsigned int max_chunks_to_wait = 4;

  unsigned int verbose = 0;

  while ((arg = getopt (argc, argv, "a:hn:p:v")) != -1)
  {
    switch (arg)
    {

      case 'a':
        address = strdup (optarg);
        break;

      case 'h':
        usage();
        return 0;   
        break;

      case 'n': 
        total_beams = atoi(optarg);
        break;

      case 'p':
        port = atoi(optarg);
        persist = true;
        break;

      case 'v':
         verbose ++;
         break;

      default:
        usage();
        return 1;
    }
  }

  signal(SIGINT, signal_handler);

  bool continue_processing = true;
  int nfiles = (argc - optind);

  if (port > 0)
    server = new ServerSocket( (char *) address, port );

  while ( continue_processing )
  {

    int beam_number = 0;
    unsigned i = 0;

    if (nfiles > 0)
    {
      chunks.resize(1);

      // create a canidate chunk from a bunch o cand files
      chunks[0] = new CandidateChunk (argc, optind, argv);

      // compute coincidence information from loaded files
      chunks[0]->compute_coincidence();

      // write the output
      chunks[0]->write_coincident_candidates();
    }
    else if ( port > 0)
    {
      bool waiting_for_beams = true;
      int connections_waiting = 0;
      unsigned int curr_beams = 0;

      while ( waiting_for_beams )
      {
        // create the conversational socket
        ServerSocket new_sock;
        connections_waiting = 0;

        // wait for some activity on the socket
        while (!quit_threads && !connections_waiting)
        {
          connections_waiting = server->select (1.0);
        }

        if ( quit_threads )
        {
          waiting_for_beams = false;
          continue_processing = false;
          cerr << "quit_threads now true" << endl;
          continue;
        }

        // wait for a client connection
        server->accept ( new_sock );

        try
        {
          std::string socket_data;
          std::ostringstream oss;
          try
          {
            while (true)
            {
              new_sock >> socket_data;
              oss << socket_data;
            }
          }
          catch ( SocketException& e)
          {
            if (verbose > 1)
            {
              cerr << "main: end of data read:" << endl;
              cerr << oss.str();
              cerr << "==================" << endl;
            }
  
          }

          std::istringstream iss(oss.str());

          string tmp_str;
          string utc_start;
          string first_sample_utc;
          uint64_t first_sample_index;
          int beam_number;
          uint64_t num_events;

          iss >> utc_start >> first_sample_utc >> first_sample_index >> beam_number >> num_events >> ws;

          if (verbose)
          {
            cerr << "main: UTC_START=" << utc_start << " SAMPLE_UTC=" << first_sample_utc 
                 << " SAMPLE_IDX=" << first_sample_index 
                 << " BEAM=" << beam_number << " NUM_EVENTS=" << num_events << endl;
          }

          // check first_sample_utc to see if it matches an existing chunk
          if ( curr_first_sample_utc.empty() || curr_first_sample_utc.compare(first_sample_utc) != 0 )
          {
            curr_chunk = -1;
            time_t youngest = 0;
            for (unsigned ichunk=0; ichunk < chunks.size(); ichunk++)
            {
              // get the relative age between this new beam/chunk and the existing ones
              time_t relative_age = chunks[ichunk]->get_relative_age (first_sample_utc);
              if (relative_age < youngest)
                youngest = relative_age;
              if (relative_age == 0)
                curr_chunk = ichunk;
            }

            // if the new beam's data preceeds the youngest chunk, discard it
            if (youngest < 0)
              curr_chunk = -2;
            
            // if we found a matching chunk
            if (curr_chunk >= 0)
            {
              if (verbose > 1)
                cerr << "main: found existing chunk for " << first_sample_utc << endl;
            }
            // if not, create a new chunk           
            if (curr_chunk == -1)
            {
              if (verbose > 1)
                cerr << "main: creating new chunk for " << first_sample_utc << endl;
              chunks.push_back (new CandidateChunk());
              curr_chunk = chunks.size() - 1; 
            }
            else
              cerr << "main: new beam " << first_sample_utc << " arrived too late" << endl;
          }
      
          // record this for a little efficiency
          if (curr_chunk >= 0)
            curr_first_sample_utc = first_sample_utc;

          // now parse the rest of this beam's data into the chunk
          if (curr_chunk >= 0)
          {
            if (verbose > 1)
              cerr << "main: chunks[" << curr_chunk <<"]->addBeam(" << beam_number << ")" << endl;
            chunks[curr_chunk]->addBeam(utc_start, first_sample_utc, first_sample_index, beam_number, num_events, iss);
            
            // if we have reached the specified number of beams for this 
            // chunk, or if more than 3 chunks are stored, dump a chunk
            if (chunks[curr_chunk]->get_n_beams() == total_beams || chunks.size() > max_chunks_to_wait)
              waiting_for_beams = false;
          }
          iss.str("");
          oss.str("");
        }
        catch ( SocketException& e ) 
        {
          cerr << "caught socket exception!" << endl;
        }
      }

      if (!quit_threads)
      {
        // if we have too many chunks, dump the first one
        if (chunks.size() > max_chunks_to_wait)
          curr_chunk = 0;

        // perform coincidence calculations
        if (verbose > 1)
          cerr << "main: chunks["<<curr_chunk<<"]->compute_coincidence()" << endl;
        chunks[curr_chunk]->compute_coincidence();

        // write the output
        if (verbose > 1)
          cerr << "main: chunks["<<curr_chunk<<"]->write_coincident_candidates()" << endl;
        chunks[curr_chunk]->write_coincident_candidates();

        // discard the curr element on the vector
        delete chunks[curr_chunk];
        chunks.erase(chunks.begin() + curr_chunk);

        curr_chunk = 0;
        curr_first_sample_utc.clear();
      }
    }
    else
    {
      cerr << "ERROR!" << endl;
      persist = false;
    }

    if (! persist )
    {
      continue_processing = false;
    }
  }
  return 0;
}

void signal_handler(int signalValue) 
{
  cerr << "receiged SIGNIT" << endl;
  quit_threads = 1;
}

