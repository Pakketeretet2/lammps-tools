#ifndef COMMUNICATION_H
#define COMMUNICATION_H

/*!
  \file communication.h
  @brief Contains some abstractions for IPC between python and the C++ lib

  \ingroup cpp_lib
*/


extern "C" {
	int get_pipe( const char *pname );

	int close_pipe( const char *pname );
}



#endif // COMMUNICATION_H
