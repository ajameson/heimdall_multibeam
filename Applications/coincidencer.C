
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

#include "hd/server_socket.h"
#include "hd/socket_exception.h"

int verbose = 0;

void usage(void)
{
  fprintf(stdout, "coincidencer [options] [beam_candidate_files]\n");
  fprintf(stdout, "  -h          print this help text\n");
  fprintf(stdout, "  -p port     listen for candidate events on port\n");
  fprintf(stdout, "  -v          verbose output\n");
}

class Candidate {

  public:

    Candidate () 
    {
      snr = sample_idx = sample_time = filter = dm_trial = dm = members = begin = end = nbeams = beam_mask = primary_beam = max_snr = 0;
    }

    Candidate (const char * line, unsigned _beam_number)
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

      if (!iss.eof()) {
        cout << "Candiate::Candidate: too many params on input line [" << line << "]" << endl;
      }
    }

    ~Candidate () { ; }

    friend ostream& operator<<(ostream& os, const Candidate& c);

    void header()
    {
      cout << "SNR\tsamp_idx\ttime\tfilter\tdm_trial\tDM\tmembers\tbegin\tend\tnbeams\tbeam_mask\tprim_beam\tmax_snr\tbeam" << endl;
    }

    bool is_coincident(const Candidate& c)
    {
      const int64_t sep_time = 3;
      const uint64_t sep_filter = 4;
      const uint64_t sep_dm = 9999;
      const float    sep_snr = 0.30;
      const int64_t tol = sep_time * powf(2,max(c.filter,filter));

      // change temporal coincidence on bens suggestion 6/8/2012
      return ( (abs(c.sample_idx - sample_idx) <= tol) &&
               (abs(c.dm_trial - dm_trial) <= sep_dm) &&
               (abs(c.filter - filter) <= sep_filter) &&
               ((fabsf(c.snr-snr) / (c.snr + snr)) <= sep_snr));
    }

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

ostream& operator<<(ostream& os, const Candidate& c) {
    os << c.snr << "\t" << c.sample_idx << "\t" << c.sample_time << "\t"
       << c.filter << "\t" << c.dm_trial << "\t" << c.dm  << "\t"
       << c.members << "\t" << c.begin << "\t" << c.end << "\t"
       << c.nbeams << "\t" << c.beam_mask << "\t" << c.primary_beam << "\t"
       << c.max_snr <<"\t" << c.beam;
    return os;
}

class CandidateChunk {

    public:

      CandidateChunk() {
        first_sample = 0;
        resize(0);
      }

      // create candidate chunk from a group of files
      CandidateChunk(int argc, int optind, char ** argv)
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
              Candidate c(line, beam_number);
              if (verbose > 1)
                cerr << c;
              cands[i].push_back(c);
            }
          }
        }
      }

      // add 
      int addBeam(string _utc_start, string _first_sample_utc, uint64_t _first_sample, unsigned int beam, uint64_t num_events, ServerSocket socket)
      {
        unsigned int this_beam = n_beams;

        if (n_beams == 0)
        {
          first_sample = _first_sample;
          first_sample_utc = _first_sample_utc;
          utc_start = _utc_start;
        }
        else
        {
          if (first_sample != _first_sample)
            cerr << "CandidateChunk::addBeam sample mismatch" << endl;
          if (utc_start != _utc_start)
            cerr << "CandidateChunk::addBeam utc_start mismatch" << endl;
        }

        resize(n_beams + 1);
        beam_numbers[this_beam] = beam;
        std::string cand;

        cerr << "CandidateChunk::addBeam adding cands to beam " << this_beam << endl;

        for (unsigned i=0; i<num_events; i++)
        {
          socket >> cand;
          Candidate c(cand.c_str(), beam);
          cerr << "CandidateChunk: " << c << endl;
          cands[this_beam].push_back(c);
        }

        cerr << "CandidateChunk::addBeam n_beams=" << n_beams << endl;
      }

      void resize (unsigned _n_beams)
      {
        n_beams = _n_beams;
        cands.resize(_n_beams);
        beam_numbers.resize(_n_beams);
      }

      unsigned int get_n_beams() const
      { 
        return n_beams;
      }

      void compute_coincidence()
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
            max_snr_j = cands[i][j].snr;

            if (cands[i][j].members < members_tol)
              continue;

            for (k=0; k < n_beams; k++)
            {
              if (i != k)
              {
                for (l=0; l<cands[k].size(); l++)
                {
                  snr_l = cands[k][l].snr;

                  if (cands[k][l].members < members_tol)
                    continue;

                  if ( cands[i][j].is_coincident (cands[k][l]) )
                  {
                    cands[i][j].nbeams ++;
                    cands[i][j].beam_mask |= 1 <<  k;

                    if (cands[i][j].nbeams >= beam_thresh + 1)
                      cands[i][j].beam_mask |= rfi_mask;

                    if (snr_l > max_snr_j)
                    {
                      cands[i][j].primary_beam = beam_numbers[k];
                      cands[i][j].max_snr = snr_l;
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

      void write_coincident_candidates()
      {
        std::string * filename;
        if ( utc_start == "")
          filename = new std::string(first_sample_utc + "_all.cand");
        else
          filename = new std::string(utc_start + "/" + first_sample_utc + "_all.cand");
        std::ofstream ofs (filename->c_str(), std::ios::out);

        if (verbose)
        {
          cerr << "CandidateChunk::write_coincident_candidates: output_file=" << filename->c_str() << endl;
          ofs << "SNR\tsamp_idx\ttime\tfilter\tdm_trial\tDM\tmembers\tbegin\tend\tnbeams\tbeam_mask\tprim_beam\tmax_snr\tbeam" << endl;
        }
        for (unsigned i=0; i< n_beams; i++)
          for (unsigned j=0; j<cands[i].size(); j++)
            ofs << cands[i][j] << endl;
        ofs.close();
      }

    private:

      std::vector<std::vector<Candidate> > cands;

      std::vector<unsigned int> beam_numbers;

      unsigned int n_beams;

      uint64_t first_sample;

      string   first_sample_utc;

      string   utc_start;

};

int main(int argc, char* argv[])
{
  int arg = 0;

  int port = 0;

  ServerSocket * server = NULL;

  // list of chunks of observations
  std::vector<CandidateChunk *> chunks;

  std::string curr_utc_start;

  std::string curr_first_sample_utc;

  int curr = 0;

  bool persist = false;

  while ((arg = getopt (argc, argv, "hp:v")) != -1)
  {
    switch (arg)
    {
      case 'h':
        usage();
        return 0;   
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

  bool continue_processing = true;
  int nfiles = (argc - optind);

  if (port > 0)
    server = new ServerSocket(port);

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
      cerr << "doing over socket!" << endl;
      
      bool waiting_for_beams = true;
      unsigned int curr_beams = 0;

      chunks.resize(1);
      chunks[0] = new CandidateChunk();

      while ( waiting_for_beams )
      {
        //create the conversational socket
        ServerSocket new_sock;

        // wait for a client connection
        server->accept ( new_sock );

        try
        {
          std::string metadata;

          new_sock >> metadata;
          cerr << "main: received '" << metadata << "'" << endl;

          stringstream ss(metadata, stringstream::in);

          string tmp_str;
          string utc_start;
          string first_sample_utc;
          uint64_t first_sample_index;
          int beam_number;
          uint64_t num_events;

          ss >> utc_start >> first_sample_utc >> first_sample_index >> beam_number >> num_events >> ws;
          if (!ss.eof())
            cerr << "main: ss parse error" << endl;

          if (verbose)
          {
            cerr << "main: utc_start: " << utc_start << endl;
            cerr << "main: first_sample_utc: " << first_sample_utc << endl;
            cerr << "main: first_sample_index: " << first_sample_index << endl;
            cerr << "main: beam_number: " << beam_number << endl;
            cerr << "main: num_events: " << num_events << endl;
          }

          if (!curr_utc_start.empty() && curr_utc_start.compare(utc_start) != 0)
          {
            cerr << "main: UTC_START mismatch!" << endl;
            curr++;
            chunks.resize(curr+1);
            waiting_for_beams = false;
          }
          curr_utc_start = utc_start;

          if (!curr_first_sample_utc.empty() && curr_first_sample_utc.compare(first_sample_utc) != 0)
          {
            cerr << "main: first_sample_utc mismatch!" << endl;
            curr++;
            chunks.resize(curr+1);
            waiting_for_beams = false;
          }
          curr_first_sample_utc = first_sample_utc;

          // now parts the rest of this beam's data into the chunk
          cerr << "main: addBeam to chunk " << curr << endl;
          chunks[curr]->addBeam(utc_start, first_sample_utc, first_sample_index, beam_number, num_events, new_sock);
          cerr << "main: beam added" << endl;

          curr_beams++;
          if (curr_beams == 13)
          {
            waiting_for_beams = false;
          }
        }
        catch ( SocketException& e ) 
        {
          cerr << "caught socket exception!" << endl;
        }
      }

      // perform coincidence calculations
      cerr << "main: chunks[curr]->compute_coincidence()" << endl;
      chunks[curr]->compute_coincidence();

      // write the output
      cerr << "main: chunks[curr]->>write_coincident_candidates()" << endl;
      chunks[curr]->write_coincident_candidates();

      // discard the first element on the vector
      chunks.erase(chunks.begin());
      curr--;

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
