// Definition of the ServerSocket class

#ifndef ServerSocket_class
#define ServerSocket_class

#include "hd/socket.h"


class ServerSocket : private Socket
{
 public:

  ServerSocket ( char * address, int port );
  ServerSocket (){};
  virtual ~ServerSocket();

  const ServerSocket& operator << ( const std::string& ) const;
  const ServerSocket& operator >> ( std::string& ) const;
  const ServerSocket& operator >> ( std::stringstream& ) const;

  void accept ( ServerSocket& );
  int select (float sleep_seconds );


};


#endif
