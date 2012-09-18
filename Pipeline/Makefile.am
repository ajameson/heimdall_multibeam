
include $(top_srcdir)/cuda.mk

INCLUDES = -I$(top_srcdir)/src

lib_LTLIBRARIES = libhdpipeline.la

libhdpipeline_la_SOURCES = default_params.C error.C parse_command_line.C clean_filterbank_rfi.cu get_rms.cu matched_filter.cu remove_baseline.cu find_giants.cu label_candidate_clusters.cu merge_candidates.cu pipeline.cu measure_bandpass.cu median_filter.cu matched_filter.cu 

include_HEADERS = 

# To find the dedisp lib from bbarsdell
AM_CPPFLAGS = -I../ -I../Network @DEDISP_CFLAGS@ @THRUST_CFLAGS@ @CUDA_CFLAGS@

include $(top_srcdir)/config/Makefile.targets

