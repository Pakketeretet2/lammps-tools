#include "communication.h"

#include <cstdlib>
#include <iostream>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


int get_pipe( const char *pname )
{
	// Construct a nice name for the pipe:
	int status = mkfifo( pname, 0666 );
	if( status < 0 ){
		unlink(pname);
		int status = mkfifo( pname, 0666 );
		if( status < 0 ){
			std::cerr << "Error " << status
			          << " opening pipe named " << pname
			          << "!\n";
			return status;
		}  
	}
	
	// From here on out there is a FIFO named pname.
	return status;
}



int close_pipe( const char *pname )
{
	int status = unlink(pname);
	return status;
}
