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
  begin = end = nbeams = beam_mask = primary_beam = max_snr = 0;
}

Candidate::Candidate (const char * line, unsigned _beam_number)
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

  nbeams = 1;
  beam_mask = 1 << (_beam_number-1);
  primary_beam = _beam_number;
  beam = _beam_number;
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
  cout << "end\tnbeams\tbeam_mask\tprim_beam\tmax_snr\tbeam" << endl;
}

bool Candidate::is_coincident(const Candidate * c)
{
  const int64_t sep_time = 3;
  const uint64_t sep_filter = 4;
  const uint64_t sep_dm = 9999;
  const float    sep_snr = 0.30;
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
     << c->nbeams << "\t" << c->beam_mask << "\t" << c->primary_beam << "\t"
     << c->max_snr <<"\t" << c->beam;
  return os;
}

CandidateChunk::CandidateChunk () {
  first_sample = 0;
  resize(0);
  verbose = 0;
}

CandidateChunk::CandidateChunk(int argc, int optind, char ** argv)
{
  // resize internal storage
  resize (argc - optind);

  char line[1024];
  string beam;
  int beam_number = 0;

  for (unsigned int i=0; i < n_beams; i++)
  {
    if (verbose)
      cerr << "CandidateChunk::CandidateChunk opening file " << argv[optind+i] << endl;

    // determine beam number from filename
    stringstream ss(basename(argv[optind+i]));
    getline(ss, first_sample_utc, '_');  // candidates
    getline(ss, beam, '.');    // beam number
    beam_number = atoi(beam.c_str());
    beam_numbers[i] = beam_number;

    if (verbose)
      cerr << "CandidateChunk::CandidateChunk parsed beam number as " << beam << endl;

    std::ifstream ifs (argv[optind+i], std::ios::in);
    while (ifs.good())
    {
      ifs.getline(line, 1024, '\n');
      if (!ifs.eof())
      {
        cands[i].push_back( new Candidate(line, beam_number));
      }
    }
  }
}

CandidateChunk::~CandidateChunk ()
{
  if (verbose)
    cerr << "CandidateChunk::~CandidateChunk" << endl;
  for (unsigned i=0; i<n_beams; i++)
  {
    for (unsigned j=0; j<cands[i].size(); j++)
      delete cands[i][j];
    cands[i].erase(cands[i].begin(), cands[i].end());
  }
  cands.erase(cands.begin(), cands.end());
}

// add beam
int CandidateChunk::addBeam (string _utc_start, string _first_sample_utc, 
                             uint64_t _first_sample, unsigned int beam, 
                             uint64_t num_events, std::istringstream& ss)
{
  unsigned int ibeam = n_beams;

  if (n_beams == 0)
  {
    first_sample = _first_sample;
    first_sample_utc = _first_sample_utc;
    utc_start = _utc_start;
  }
  else
  {
    if (first_sample != _first_sample)
      cerr << "CandidateChunk::addBeam sample mismatch [" << first_sample << " != " << _first_sample << "]" << endl;
    if (utc_start != _utc_start)
      cerr << "CandidateChunk::addBeam utc_start mismatch" << endl;
  }

  // resize storage for this new beam
  resize(n_beams + 1);
  beam_numbers[ibeam] = beam;

  if (verbose > 1)
    cerr << "CandidateChunk::addBeam resized to " << n_beams << " beams with beam " << beam << endl;

  char cand_line[1024];
  cands[ibeam].resize(num_events);
  for (unsigned ievent=0; ievent < num_events; ievent++)
  {
    ss.getline(cand_line, 1024, '\n');
    cands[ibeam][ievent] = new Candidate(cand_line, beam);
    //cerr << "cands[" << ibeam << "][" << ievent << "]=" << cands[ibeam][ievent] << endl;
    //Candidate c(cand_line, beam);
    //cands[ibeam].push_back(c);
  }
}

void CandidateChunk::resize (unsigned _n_beams)
{
  n_beams = _n_beams;
  cands.resize(_n_beams);
  beam_numbers.resize(_n_beams);
}

unsigned int CandidateChunk::get_n_beams() const
{ 
  return n_beams;
}

void CandidateChunk::compute_coincidence()
{
  float max_snr_j;
  float snr_l;
  unsigned i, j, k, l;

  unsigned int members_tol = 3;
  unsigned int rfi_mask = 1 << 16;
  unsigned int beam_thresh = 2;

  // compute coincidence information
  for (i=0; i < n_beams; i++)
  {
    for (j=0; j<cands[i].size(); j++)
    {
      max_snr_j = cands[i][j]->snr;

      if (cands[i][j]->members < members_tol)
        continue;

      for (k=0; k < n_beams; k++)
      {
        if (i != k)
        {
          for (l=0; l<cands[k].size(); l++)
          {
            snr_l = cands[k][l]->snr;

            if (cands[k][l]->members < members_tol)
              continue;

            if ( cands[i][j]->is_coincident (cands[k][l]) )
            {
              cands[i][j]->nbeams ++;
              cands[i][j]->beam_mask |= 1 <<  k;

              if (cands[i][j]->nbeams >= beam_thresh + 1)
                cands[i][j]->beam_mask |= rfi_mask;

              if (snr_l > max_snr_j)
              {
                cands[i][j]->primary_beam = beam_numbers[k];
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

  for (unsigned i=0; i< n_beams; i++)
    for (unsigned j=0; j<cands[i].size(); j++)
      ofs << cands[i][j] << endl;
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
