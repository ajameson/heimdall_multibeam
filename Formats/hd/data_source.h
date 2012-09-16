
#ifndef __IDataSource_hpp
#define __IDataSource_hpp

class IDataSource {

	public:
		virtual ~IDataSource() {}
		virtual size_t nchans() const = 0;
		virtual size_t nbits()  const = 0;
		virtual size_t stride() const = 0;
		virtual size_t beam()   const = 0;
    virtual time_t utc_start() const = 0;
		virtual double spectra_per_second() const = 0;
		virtual float tsamp () const = 0;
		virtual bool   error()  const = 0;
		virtual size_t get_data(size_t nsamps, char* data) = 0;

  protected:
    size_t _nchan;
    size_t _nbit;
    size_t _npol;
    size_t _stride;
    size_t _beam;
    time_t _utc_start;
    float _tsamp;
    double _spectra_per_second;

};

#endif
