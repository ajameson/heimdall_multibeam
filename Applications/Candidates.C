/***************************************************************************
 *
 *   Copyright (C) 2012 by Ben Barsdell and Andrew Jameson
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "hd/Candidates.h"

#include "tmutil.h"

#include <fstream>

#include <math.h>
#include <stdlib.h>
#include <libgen.h>
#include <sys/stat.h>
#include <string.h>

using namespace std;

Candidate::Candidate ()
{
  snr = sample_idx = sample_time = filter = dm_trial = dm = members = 0;
  begin = end = nbeams = primary_beam = max_snr = 0;
}

Candidate::Candidate (const char * line)
{
  istringstream iss (line, istringstream::in);
  iss >> snr;
  iss >> sample_idx;
  iss >> sample_time;
  iss >> filter;
  iss >> dm_trial;
  iss >> dm;
  iss >> members;
  iss >> begin;
  iss >> end;
  iss >> beam;

  nbeams = 1;

  primary_beam = beam;
  max_snr = snr;
  
  iss >> ws;

  if (!iss.eof())
    cout << "Candiate::Candidate too many params on input line [" << line << "]" << endl;
}

Candidate::~Candidate ()
{
  //cerr << "Candiate::~Candidate beam=" << beam << " sample_idx=" << sample_idx << endl;
}

void Candidate::header()
{
  cout << "SNR\tsamp_idx\ttime\tfilter\tdm_trial\tDM\tmembers\tbegin\t";
  cout << "end\tnbeams\ttprim_beam\tmax_snr\tbeam" << endl;
}

bool Candidate::is_coincident(const Candidate * c, const unsigned sep_time, const unsigned sep_filter)
{
  //const int64_t sep_time = 3;
  //const uint64_t sep_filter = 4;
  const uint64_t sep_dm = 9999;
  const float    sep_snr = 0.1;
  const int64_t tol = sep_time * powf(2,max(c->filter,filter));

  // change temporal coincidence on bens suggestion 6/8/2012
  return ( (abs(c->sample_idx - sample_idx) <= tol) &&
           (abs(c->dm_trial - dm_trial) <= sep_dm) &&
           (abs(c->filter - filter) <= sep_filter) &&
           ((fabsf(c->snr - snr) / (c->snr + snr)) <= sep_snr));
}

std::ostream& operator<<(std::ostream& os, const Candidate * c) {
  os << c->snr << "\t" << c->sample_idx << "\t" << c->sample_time << "\t"
     << c->filter << "\t" << c->dm_trial << "\t" << c->dm  << "\t"
     << c->members << "\t" << c->begin << "\t" << c->end << "\t"
     << c->nbeams << "\t" << "\t" << c->primary_beam << "\t"
     << c->max_snr <<"\t" << c->beam;
  return os;
}

CandidateChunk::CandidateChunk () {
  first_sample = 0;
  resize(0);
  verbose = 0;
  configured = false;
  nbeam = 0;
  nbeam_size = 0;
}

CandidateChunk::CandidateChunk (unsigned nbeams)
{
  verbose = 0;
  nbeam = 0;
  nbeam_size = 0;
  // resize internal storage
  resize (nbeams);
  configured = false;
}


CandidateChunk::CandidateChunk(unsigned nbeams, int argc, int optind, char ** argv)
{
  nbeam = nbeams;
  nbeam_size = 0;

  // resize internal storage
  resize (nbeams);
  verbose = 0;

  int nfiles = argc - optind;

  char line[1024];
  string beam;

  // parse the date/time stamp from the first file
  if (nfiles > 0)
  {
    char * first_file = strdup (argv[optind]);
    char * pch;
    pch = strtok (first_file, "_");
    if (pch != NULL)
      setFirstSampleUTC(pch);
  }
    
  for (int i=0; i<nfiles; i++)
  {
    if (verbose)
      cerr << "CandidateChunk::CandidateChunk opening file " << argv[optind+i] << endl;


    std::ifstream ifs (argv[optind+i], std::ios::in);
    while (ifs.good())
    {
      ifs.getline(line, 1024, '\n');
      if (!ifs.eof())
      {
        if (verbose)
          cerr << "Adding candidate from " << line << endl;
        Candidate * cand = new Candidate(line);
        cands[(cand->beam-1)].push_back (cand);
      }
    }
  }
  configured = true;
}

CandidateChunk::~CandidateChunk ()
{
  if (verbose)
    cerr << "CandidateChunk::~CandidateChunk" << endl;
  for (unsigned i=0; i<nbeam_size; i++)
  {
    for (unsigned j=0; j<cands[i].size(); j++)
      delete cands[i][j];
    cands[i].erase(cands[i].begin(), cands[i].end());
  }
  cands.erase(cands.begin(), cands.end());
}

void CandidateChunk::setFirstSampleUTC (const char * str)
{
  first_sample_utc.assign (str);
}

// add beam
int CandidateChunk::addBeam (string _utc_start, string _first_sample_utc,
                             uint64_t _first_sample, unsigned int nbeams,
                             uint64_t num_events, std::istringstream& ss)
{
  if (!configured)
  {
    first_sample = _first_sample;
    first_sample_utc = _first_sample_utc;
    utc_start = _utc_start;
    configured = true;
  }
  else
  {
    if (first_sample != _first_sample)
      cerr << "CandidateChunk::addBeam sample mismatch [" << first_sample << " != " << _first_sample << "]" << endl;
    if (utc_start != _utc_start)
      cerr << "CandidateChunk::addBeam utc_start mismatch" << endl;
  }
 
  // increment the number of beams integrated into this chunk
  nbeam += nbeams;

  if (verbose > 1)
    cerr << "CandidateChunk::addBeam adding " << nbeams << " containing " << num_events << " to chunk" << endl;

  char cand_line[1024];
  for (unsigned ievent=0; ievent < num_events; ievent++)
  {
    // read a line from the ss
    ss.getline (cand_line, 1024, '\n');
    if ((!ss.eof()) && strlen(cand_line) > 10)
    {
      // build a candiate
      Candidate * cand = new Candidate(cand_line);
      // determine beam of this candidate
      int ibeam = cand->beam - 1;
      if (ibeam >= 0 && ibeam < nbeam_size)
      {
        if (verbose > 1)
          cerr << "cands[" << ibeam <<"].push_back (" << cand << ") [" << ievent << "]" << endl;
        // append the candidate to the right beam
        cands[ibeam].push_back (cand);
      }
    }
  }
}

void CandidateChunk::resize (unsigned nbeams)
{
  nbeam_size = nbeams;
  cands.resize(nbeam_size);
}

unsigned int CandidateChunk::get_nbeam_size() const
{
  return nbeam_size;  
}

unsigned int CandidateChunk::get_nbeam() const
{
  return nbeam;
}

void CandidateChunk::compute_coincidence (unsigned sep_time, unsigned sep_filter)
{
  float max_snr_j;
  float snr_l;
  unsigned i, j, k, l;

  unsigned int members_tol = 3;

  if (verbose)
    cerr << "CandidateChunk::compute_coincidence nbeams=" << nbeam_size
         << " cands.size()=" << cands.size() << endl;

  // compute coincidence information for beam [i]
  for (i=0; i < cands.size(); i++)
  {
    if (verbose)
      cerr << "compute_coincidence assessing " << cands[i].size() << " candidates from beam " << i << " of " << nbeam << endl;
    for (j=0; j<cands[i].size(); j++)
    {
      max_snr_j = cands[i][j]->snr;

      // narrow, 0 DM stuff has very few members
      //if (cands[i][j]->members < members_tol && cands[i][j]->dm_trial > 0)
      //  continue;

      // compare canidiate [i][j] vs every other candidate
      for (k=0; k<cands.size(); k++)
      {
        // exclude self comparison of beams
        if (i != k)
        {
          for (l=0; l<cands[k].size(); l++)
          {
            snr_l = cands[k][l]->snr;

            //if (cands[k][l]->members < members_tol && cands[i][j]->dm_trial > 0)
            //  continue;

            if ( cands[i][j]->is_coincident (cands[k][l], sep_time, sep_filter ) )
            {
              cands[i][j]->nbeams ++;

              if (snr_l > max_snr_j)
              {
                cands[i][j]->primary_beam = cands[k][l]->beam;
                cands[i][j]->max_snr = snr_l;
                max_snr_j = snr_l;
              }
              break;
            }
          }
        }
      }
    }
  }
}

void CandidateChunk::write_coincident_candidates()
{
  std::string * filename;
  if ( utc_start == "")
    filename = new std::string(first_sample_utc + "_all.cand");
  else
  {
    struct stat dir_stat;
    if ((stat(utc_start.c_str(), &dir_stat) == 0) && (((dir_stat.st_mode) & S_IFMT) == S_IFDIR))
      filename = new std::string(utc_start + "/" + first_sample_utc + "_all.cand");
    else
    {
      cerr << "directory [" << utc_start << "] did not exist, not writing candidate file" << endl;
      return;
    }
  }
  std::ofstream ofs (filename->c_str(), std::ios::out);

  if (verbose)
    cerr << "CandidateChunk::write_coincident_candidates: output_file=" << filename->c_str() << endl;

  for (unsigned i=0; i<cands.size(); i++)
  {
    for (unsigned j=0; j<cands[i].size(); j++)
    {
      ofs << cands[i][j] << endl;
    }
  }
  ofs.close();
}

bool CandidateChunk::matches (std::string utc)
{
  return (first_sample_utc.compare(utc) == 0);
}

time_t CandidateChunk::get_relative_age (std::string utc)
{
  time_t self_age = str2utctime (first_sample_utc.c_str());
  time_t utc_age = str2utctime (utc.c_str());
  return (utc_age - self_age);
}

time_t CandidateChunk::str2utctime (const char* str)
{
  struct tm time;
  return str2utctm (&time, str);
}

time_t CandidateChunk::str2utctm (struct tm* time, const char* str)
{
  
  /* append the GMT+0 timeszone information */
  char * str_utc = (char *) malloc(sizeof(char) * (strlen(str) + 4 + 1));
  sprintf(str_utc, "%s UTC",str);

  const char * format = "%Y-%m-%d-%H:%M:%S %Z";
  
  strptime(str_utc, format, time);

  free(str_utc);
  return timegm(time);
}
