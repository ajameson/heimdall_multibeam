/***************************************************************************
 *
 *   Copyright (C) 2012 by Ben Barsdell
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <string>
using std::string;
#include <vector>
using std::vector;
#include <fstream>
#include <iomanip>
#include <iterator>
#include <numeric>
#include <cstdlib> // For atoi
#include <cmath>
#include <algorithm>
#include <inttypes.h>

#include "hd/header.h"
#include "hd/error.h"
#include "hd/default_params.h"
#include "hd/params.h"

#include <dedisp.h>

void usage(char * binary, hd_params params)
{
  cerr << "Usage: " << binary << " [options]" << endl;
  cerr << "  generate the list of DM trials" << endl;
  cerr << "    -f0 f               centre frequency of highest channel [default " << params.f0 << "] MHz" << endl;
  cerr << "    -df f               width of frequecy channel [default " << params.df << "] MHz" << endl;
  cerr << "    -nchan n            number of frequency channels [default " << params.nchans << "]" << endl;
  cerr << "    -dt s               sampling time [default " << params.dt << "] s" << endl;
  cerr << "    -dm min max         DM range [default " << params.dm_min << " " << params.dm_max << "]" << endl;
  cerr << "    -dm_pulse_width w   [default " << params.dm_pulse_width << "] " << endl;
  cerr << "    -dm_tol t           tolerance between DM trials [default " << params.dm_tol << "]" << endl;
}

int main(int argc, char* argv[])
{
  bool verbose = false;
  dedisp_error derror;
  hd_params   params;
  hd_params   default_params;
  hd_set_default_params (&default_params);
  hd_set_default_params (&params);

  size_t i=0;
  while( ++i < (size_t)argc ) 
  {
    if( argv[i] == string("-h") ) {
      usage (argv[0], default_params);
      return -1;
    }

    if (argv[i] == string("-v")) {
      params.verbosity = std::max (params.verbosity, 1);
    }
    else if (argv[i] == string("-f0")) {
      params.f0 = atof(argv[++i]);
    }
    else if (argv[i] == string("-df")) {
      params.df = atof(argv[++i]);
    }
    else if (argv[i] == string("-nchan")) {
      params.nchans = atoi(argv[++i]);
    }
    else if (argv[i] == string("-dt")) {
      params.dt = atof(argv[++i]);
    }
    else if( argv[i] == string("-dm")) {
      params.dm_min = atof(argv[++i]);
      params.dm_max = atof(argv[++i]);
    }
    else if( argv[i] == string("-dm_pulse_width") ) {
      params.dm_pulse_width = atof(argv[++i]);
    }
    else if( argv[i] == string("-dm_tol") ) {
      params.dm_tol = atof(argv[++i]);
    }
    else {
      cerr << "WARNING: Unknown parameter '" << argv[i] << "'" << endl;
    }
  }

  dedisp_size dm_count;
  const float * dm_list = dedisp_generate_dm_list_guru (params.dm_min,
                                    params.dm_max,
                                    params.dt,
                                    params.dm_pulse_width,
                                    params.f0,
                                    params.df,
                                    params.nchans,
                                    params.dm_tol,
                                    &dm_count);
  if( derror != DEDISP_NO_ERROR ) {
    throw_dedisp_error(derror);
    return 1;
  }

  for( hd_size i=0; i<dm_count; ++i ) 
  {
    cout << dm_list[i] << endl;
  }
}
