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
    Candidate (const char * line);
    ~Candidate ();

    void header();

    bool is_coincident(const Candidate * c, const unsigned sep_time, const unsigned sep_filter);

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
    unsigned int  primary_beam;
    float         max_snr;
    unsigned int  beam;
};

class CandidateChunk 
{
  public:

    CandidateChunk();
    CandidateChunk(unsigned nbeams);
    CandidateChunk(unsigned nbeams, int argc, int optind, char ** argv);
    ~CandidateChunk();

    void setFirstSampleUTC (const char * str);

    int addBeam(std::string _utc_start, std::string _first_sample_utc, uint64_t _first_sample, unsigned int beam, uint64_t num_events, std::istringstream& ss);

    void resize (unsigned nbeams);

    unsigned int get_nbeam_size() const;
    unsigned int get_nbeam() const;

    void compute_coincidence (unsigned sep_time, unsigned sep_filter);

    void write_coincident_candidates();

    // returns true if utc matches this chunk's utc_start
    bool matches (std::string utc);

    // return the relative age (in seconds) for the specified utc
    time_t get_relative_age (std::string utc);

    std::string get_utc () const { return first_sample_utc; };

  private:

    std::vector<std::vector<Candidate *> > cands;

    // number of beams added to the chunk
    unsigned int nbeam;

    // total number of beams for the chunk
    unsigned int nbeam_size;

    uint64_t first_sample;

    std::string first_sample_utc;

    std::string   utc_start;

    int verbose;

    time_t str2utctime (const char* str);

    time_t str2utctm (struct tm* time, const char* str);

    bool configured;
};

