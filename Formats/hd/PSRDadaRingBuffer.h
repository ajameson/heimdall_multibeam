/***************************************************************************
 *
 *   Copyright (C) 2012 by Andrew Jameson
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "dada_hdu.h"

#include "hd/DataSource.h"

class PSRDadaRingBuffer : public DataSource {

  public:

     PSRDadaRingBuffer (key_t dada_id);
    ~PSRDadaRingBuffer ();

    bool connect ();
    bool disconnect ();
    bool read_header ();

    bool   get_error()  const { return dada_error; }
    size_t get_data (size_t nsamps, char* data);
    size_t get_data_block (size_t nsamps, char* data);

    size_t open_data_block (char ** buf_ptr);
    void   close_data_block (uint64_t bytes_in_block);

  private:

    key_t  dada_key;

    dada_hdu_t* hdu;

    bool dada_error;

    bool printed_first_line;

    unsigned resolution;

    char order[4];

    // size of the data block buffer
    uint64_t buf_sz;

    // number of samples per buffer
    uint64_t samples_per_buf;

    // number of bytes read in the currently open buffer
    uint64_t bytes_read;

    // pointer to the current buffer
    char * curr_block;

};

