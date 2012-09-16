
#include <fstream>

#include "hd/data_source.h"

class FileDataSource : public IDataSource 
{
  public:

    FileDataSource (const char* filename);
    ~FileDataSource ();

    size_t nchans() const { return _nchan; }
    size_t nbits()  const { return _nbit; }
    size_t stride() const { return _stride; }
    size_t beam()   const { return _beam; } // Note: This isn't meaningful
    time_t utc_start()  const { return _utc_start; }
    double spectra_per_second()  const { return _spectra_per_second; }
    float tsamp()  const { return _tsamp; }
    bool   error()  const { return m_error != 0; }
    size_t get_data(size_t nsamps, char* data);

  private:

    std::ifstream m_file_stream;
    int           m_error;

};
