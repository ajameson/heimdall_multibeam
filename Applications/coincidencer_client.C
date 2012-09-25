/***************************************************************************
 *
 *   Copyright (C) 2012 by Ben Barsdell and Andrew Jameson
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

#include <iostream>
#include <string>

#include "hd/ClientSocket.h"
#include "hd/SocketException.h"

using namespace std;

int main ( int argc, char * argv[] )
{
  try
    {
      ClientSocket client_socket ( "192.168.0.10", 51041 );

      std::string reply;

      try
      {
        client_socket << "2012-08-20-03:05:25 2012-08-20-03:15:25 1000 12 2\r\n";
        std::cerr << "2012-08-20-03:05:25 2012-08-20-03:15:25 1000 12 2" << endl;
        sleep(1);

        client_socket << "6.75968 5500778 352.05  3       2       0.360012        6       5500778 5500780\r\n";
        std::cerr << "6.75968 5500778 352.05  3       2       0.360012        6       5500778 5500780" << endl;
        sleep(1);

        client_socket << "7.39276 5550148 355.209 4       1       0.180007        10      5550144 5550156\r\n";
        std::cerr << "7.39276 5550148 355.209 4       1       0.180007        10      5550144 5550156" << endl;
        sleep(1);
/*
        client_socket << "2012-08-20-03:05:25" << " ";
        std::cerr << "2012-08-20-03:05:25" << " ";
        sleep(1);

        client_socket << "2012-08-20-03:15:25" << " ";
        std::cerr << "2012-08-20-03:15:25" << " ";
        sleep(1);

        client_socket << "1000" << " ";
        std::cerr << "1000" << " ";
        sleep(1);

        client_socket << "12" << " ";
        std::cerr << "12" << " ";
        sleep(1);

        client_socket << "2" << "\r\n";
        std::cerr << "2" << "\r\n";
        sleep(1);

        client_socket << "6.75968 5500778 352.05  3       2       0.360012        6       5500778 5500780\r\n";
        std::cerr << "6.75968 5500778 352.05  3       2       0.360012        6       5500778 5500780\r\n";
        sleep(1);

        client_socket << "7.39276 5550148 355.209 4       1       0.180007        10      5550144 5550156\r\n";
        std::cerr << "7.39276 5550148 355.209 4       1       0.180007        10      5550144 5550156\r\n";
        sleep(1);
        */
      }
      catch ( SocketException& ) {}

      std::cout << "We received this response from the server:\n\"" << reply << "\"\n";;

    }
  catch ( SocketException& e )
    {
      std::cout << "Exception was caught: " << e.description() << "\n";
    }

  return 0;
}
