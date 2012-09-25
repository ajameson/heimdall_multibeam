/***************************************************************************
 *
 *   Copyright (C) 2012 by Ben Barsdell and Andrew Jameson
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include "hd/ServerSocket.h"
#include "hd/SocketException.h"


ServerSocket::ServerSocket ( char * address, int port )
{
  if ( ! Socket::create() )
    throw SocketException ( "Could not create server socket." );

  if ( ! Socket::bind ( address, port ) )
    throw SocketException ( "Could not bind to port." );

  if ( ! Socket::listen() )
    throw SocketException ( "Could not listen to socket." );
}

ServerSocket::~ServerSocket()
{
}

const ServerSocket& ServerSocket::operator << ( const std::string& s ) const
{
  if ( ! Socket::send ( s ) )
    throw SocketException ( "Could not write to socket." );

  return *this;
}

const ServerSocket& ServerSocket::operator >> ( std::string& s ) const
{
  if ( ! Socket::recv ( s ) )
    throw SocketException ( "Could not read from socket." );

  return *this;
}

const ServerSocket& ServerSocket::operator >> ( std::stringstream& ss ) const
{
  std::string s;
  if ( ! Socket::recv ( s ) )
    throw SocketException ( "Could not read from socket." );
  return *this;
}


void ServerSocket::accept ( ServerSocket& sock )
{
  if ( ! Socket::accept ( sock ) )
    throw SocketException ( "Could not accept socket." );
}

int ServerSocket::select (float sleep_seconds )
{
  return Socket::select_timeout (sleep_seconds );
}
