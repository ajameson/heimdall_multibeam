/***************************************************************************
 *
 *   Copyright (C) 2012 by Ben Barsdell and Andrew Jameson
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include <vector>
#include <iostream>
#include <sstream>
#include <string>

#include <inttypes.h>
#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE /* glibc2 needs this for strptime  */
#endif
#include <time.h>

class Candidate 
{
  public:

    Candidate (); 
    Candidate (const char * line, unsigned _beam_number);
    ~Candidate ();

    void header();

    bool is_coincident(const Candidate * c);

    friend std::ostream& operator<<(std::ostream& os, const Candidate * c);

    float         snr;
    int64_t       sample_idx;
    float         sample_time;
    unsigned int  filter;
    unsigned int  dm_trial;
    float         dm;
    unsigned int  members;
    int64_t       begin;
    int64_t       end;
    unsigned      nbeams;
    unsigned      beam_mask;
    unsigned int  primary_beam;
    float         max_snr;
    unsigned int  beam;
};

class CandidateChunk 
{
  public:

    CandidateChunk();
    CandidateChunk(int argc, int optind, char ** argv);
    ~CandidateChunk();

    int addBeam(std::string _utc_start, std::string _first_sample_utc, uint64_t _first_sample, unsigned int beam, uint64_t num_events, std::istringstream& ss);


    void resize (unsigned _n_beams);

    unsigned int get_n_beams() const;

    void compute_coincidence();

    void write_coincident_candidates();

    // returns true if utc matches this chunk's utc_start
    bool matches (std::string utc);

    // return the relative age (in seconds) for the specified utc
    time_t get_relative_age (std::string utc);

  private:

    std::vector<std::vector<Candidate *> > cands;

    std::vector<unsigned int> beam_numbers;

    unsigned int n_beams;

    uint64_t first_sample;

    std::string first_sample_utc;

    std::string   utc_start;

    int verbose;

    time_t str2utctime (const char* str);
    time_t str2utctm (struct tm* time, const char* str);

};

