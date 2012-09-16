// Definition of the ClientSocket class

#ifndef ClientSocket_class
#define ClientSocket_class

#include "hd/socket.h"


class ClientSocket : private Socket
{
 public:

  ClientSocket ( std::string host, int port );
  virtual ~ClientSocket(){};

  const ClientSocket& operator << ( const std::string& ) const;
  const ClientSocket& operator << ( const char * ) const;
  const ClientSocket& operator << ( const size_t n ) const;
  const ClientSocket& operator << ( const float f ) const;
  const ClientSocket& operator >> ( std::string& ) const;

};


#endif
