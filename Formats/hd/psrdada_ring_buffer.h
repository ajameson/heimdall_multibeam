#include "dada_hdu.h"

#include "hd/data_source.h"

class DadaDataSource : public IDataSource {

	public:

 		DadaDataSource (key_t dada_id);
  	~DadaDataSource ();

		bool connect ();
		bool disconnect ();
	  bool read_header ();

	  size_t nchans() const { return _nchan; }
  	size_t nbits()  const { return _nbit; }
	  size_t stride() const { return _stride; }
  	size_t beam()   const { return _beam; }
  	time_t utc_start()  const { return _utc_start; }
  	double spectra_per_second()  const { return _spectra_per_second; }
  	float  tsamp()  const { return _tsamp; }
  	bool   error()  const { return dada_error; }
	  size_t get_data(size_t nsamps, char* data);

	private:

		key_t	dada_key;

  	dada_hdu_t* hdu;

	  bool dada_error;
    bool printed_first_line;

};

