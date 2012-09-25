
#include <fstream>

#include "hd/DataSource.h"

class SigprocFile: public DataSource 
{
  public:

    SigprocFile (const char* filename);
    ~SigprocFile ();

    bool   get_error() const { return m_error != 0; }
    size_t get_data (size_t nsamps, char* data);

  private:

    std::ifstream m_file_stream;
    int           m_error;

};
