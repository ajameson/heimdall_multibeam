#!/usr/bin/env python
#
# Filename: trans_cand_server.py
#
#  * return plots as binary data over socket 

import Dada
import sys
import time
import math
import os
import tempfile

DL = 1
#TSAMP = 0.00032768
verbose = False

def log(lvl, dlvl, message):
  message = message.replace("`","'")
  if (lvl <= dlvl):
    time = Dada.getCurrentDadaTimeUS()
    if args.logfile == None:
      Dada.logMsg(lvl, dlvl, message)
    else:
      fptr = open(args.logfile, "a")
      if (lvl == -1):
        fptr.write("[" + time + "] WARN " + message + "\n")
      elif (lvl == -2):
        fptr.write("[" + time + "] ERR  " + message + "\n")
      else:
        fptr.write("[" + time + "] " + message + "\n")
      fptr.close()

###########################################################################

def plotCandDspsr (fil_file, sample, filter, dm, beam, snr):
  
  # determine filterbank sampling time in micro seconds
  cmd = "header " + fil_file + " | grep 'Sample time' | awk '{print $NF}'"
  p = os.popen(cmd)
  response = p.readline().strip()
  p.close()
  tsamp = float(response) / 1e6;

  # determine the number of channesl
  cmd = "header " + fil_file + " | grep 'Number of channels' | awk '{print $NF}'"
  p = os.popen(cmd)
  response = p.readline().strip()
  p.close()
  nchan_fb = int(response)

  cand_time = (tsamp * sample)
  cmd = "dmsmear -f 835 -b 31.25 -n 320 -d " + str(dm) + " -q 2>&1 "
  log(2, DL, "plotCandDspsr: " + str(cmd))
  p = os.popen(cmd)
  cand_band_smear = p.readline().strip()
  p.close()
  log(2, DL, "plotCandDspsr: cand_band_smear='" + cand_band_smear + "'")

  cand_filter_time = (2 ** filter) * tsamp

  cand_smearing = float(cand_band_smear) + float(cand_filter_time)

  cand_start_time = cand_time - (2.0 * cand_smearing)
  cand_tot_time   = 5.0 * cand_smearing
  if cand_start_time < 0:
    cand_start_time = 0
  #cand_start_time = cand_time
  #cand_tot_time   = cand_smearing

  log(2, DL, "plotCandDspsr: cand_time=" + str(cand_time) + " cand_start_time=" + str(cand_start_time))

  bin_width = tsamp * ( 2 ** filter)
  nbin = int (cand_tot_time / bin_width)

  if nbin < 64:
    nbin = 64

  cmd = "dspsr -q " + fil_file + " -S " + str(cand_start_time) + \
        " -b " + str(nbin) + \
        " -T " + str(cand_tot_time) + \
        " -c " + str(cand_tot_time) + \
        " -D " + str(dm) + \
        " -U 0.01" + \
        " -cepoch=start" + \
        " |& grep unloading | awk '{print $NF}'"

  # create a temporary working directory
  workdir = tempfile.mkdtemp()

  os.chdir(workdir)

  log(2, DL, "plotCandDspsr: " + cmd)
  p = os.popen(cmd)
  response = p.readline().strip()
  log(3, DL, "plotCandDspsr: response=" + response)
  p.close()

  archive = response + ".ar"
  count = 10 
  while ((not os.path.exists(archive)) and count > 0):
    log(1, DL, "plotCandDspsr: archive file [" + archive + "] did not exist")
    time.sleep(1)
    count = count - 1

  binary_data = []

  width = "{0:.2f}".format(round(cand_filter_time*1000,2))
  title = beam + "  DM " + str(dm) + "  Width " + width + "ms"

  # determine the ideal number of channels for the plot
  nchan = int(round(math.pow(float(snr)/4.0,2)))

  log(2, DL, "ideal nchan_plot=" + str(nchan))

  # upper limit on channelisation
  if nchan > nchan_fb:
    nchan = nchan_fb

  # ensure some minimal amount of channelisation
  while (nchan < 64 or nchan_fb % nchan != 0) and nchan_fb > nchan:
    nchan += 1

  log(2, DL, "nchan_fb=" + str(nchan_fb) + " nchan=" + str(nchan))

  cmd = "psrplot -c flux:below:l='' -c above:l='" + title + "' -c x:unit=ms -p freq+ ./" + archive + " -j 'F " + str(nchan) + "' -x -D -/PNG";
  log(2, DL, "plotCandDspsr: " + cmd)
  p = os.popen(cmd)
  binary_data = p.read()
  p.close()
  
  if os.path.exists(archive):
    os.remove(archive)
  os.chdir ("/")
  os.rmdir(workdir)

  return binary_data

############################################################################### 
#
# main
#
if __name__ == "__main__":

  import argparse

  parser = argparse.ArgumentParser(description="Plot a transient event to stdout")
  parser.add_argument("fil_file", help="filterbank file to process")
  parser.add_argument("sample", help="sample number of event start", type=int)
  parser.add_argument("dm", help="dm of event", type=float)
  parser.add_argument("filter", help="filter of event", type=int)
  parser.add_argument("snr", help="snr of event", type=float)
  parser.add_argument('-beam', help="name of beam", type=str)
  parser.add_argument('-logfile', help="file to log output to", type=str, default=None)
  parser.add_argument('-verbose', action="store_true")
  args = parser.parse_args()

  fil_file = args.fil_file
  sample = args.sample
  dm = args.dm
  filter = args.filter
  verbose = args.verbose
  beam = args.beam
  snr = args.snr

  if not os.path.exists(fil_file):
    sys.exit('ERROR: Filterbank file %s was not found!' % fil_file)

  try:
    fptr = open(fil_file, 'r')
  except IOError:
    log(1, DL, "main: could not open filterbank file for reading")
    fil_file = ""
  else:
    fptr.close()

  if len(fil_file) == 0:
    log(1, DL, "main: could not find filterbank file")
  else:
    binary_data = []

    log(2, DL, "main: plotCandDspsr()")
    binary_data = plotCandDspsr(fil_file, sample, filter, dm, beam, snr)
    binary_len = len(binary_data)
    log(3, DL, "main: sending binary data len="+str(binary_len))

  sys.stdout.write(binary_data)

  # exit
  sys.exit(0)
