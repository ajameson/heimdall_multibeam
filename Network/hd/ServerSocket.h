/***************************************************************************
 *
 *   Copyright (C) 2012 by Ben Barsdell and Andrew Jameson
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#ifndef hd_ServerSocket
#define hd_ServerSocket

#include "hd/Socket.h"

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
