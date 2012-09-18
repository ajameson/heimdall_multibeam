
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

		static time_t mjd2utctm (double mjd)
		{
			const int seconds_in_day = 86400;
			int days = (int) mjd;
			double fdays = mjd - (double) days;
			double seconds = fdays * (double) seconds_in_day;
			int secs = (int) seconds;
			double fracsec = seconds - (double) secs;
			if (fracsec - 1 < 0.0000001)
				secs++;

			int julian_day = days + 2400001;

			int n_four = 4  * (julian_day+((6*((4*julian_day-17918)/146097))/4+1)/2-37);
			int n_dten = 10 * (((n_four-237)%1461)/4) + 5;

			struct tm gregdate;
			gregdate.tm_year = n_four/1461 - 4712 - 1900; // extra -1900 for C struct tm
			gregdate.tm_mon  = (n_dten/306+2)%12;         // struct tm mon 0->11
			gregdate.tm_mday = (n_dten%306)/10 + 1;

			gregdate.tm_hour = secs / 3600;
			secs -= 3600 * gregdate.tm_hour;


			gregdate.tm_min = secs / 60;
			secs -= 60 * (gregdate.tm_min);

			gregdate.tm_sec = secs;

			gregdate.tm_isdst = -1;
			time_t date = mktime (&gregdate);

			return date;
		}

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
