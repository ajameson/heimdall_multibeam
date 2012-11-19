
INCLUDES = -I$(top_srcdir)/src

lib_LTLIBRARIES = libhdformats.la

libhdformats_la_SOURCES = SigprocFile.C 

include_HEADERS = 

bin_PROGRAMS = 

if HAVE_PSRDADA

libhdformats_la_SOURCES += PSRDadaRingBuffer.C

AM_CPPFLAGS = @PSRDADA_CFLAGS@

endif

include $(top_srcdir)/config/Makefile.targets

