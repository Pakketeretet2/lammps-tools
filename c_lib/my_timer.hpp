#ifndef MY_TIMER_HPP
#define MY_TIMER_HPP


#include <sys/time.h>
#include <iostream>


/*!
  @brief A simple timer class based on sys/time

  This class is a fairly accurate timer. Useful for timing parts of programs
  like loops and stuff.

  \ingroup misc
*/
class my_timer {
public:
	/// Default constructor, no output.
	my_timer() : out( nullptr )
	{ init_tic_toc(); }

	/// Constructor that takes an std::ostream to which stuff is
	/// occasionally printed.
	my_timer(std::ostream &out_stream ) : out( &out_stream )
	{ init_tic_toc(); }

	/// Empty destructor
	~my_timer(){}
	
	/// Sets the "tic"-time to current time.
	void tic()
	{ gettimeofday(&t_tic, nullptr); }


	/*!
	  @brief Computes difference between the "tic"-time and current time.

	  @param msg   A message to print to out in addition
	               to the elapsed time (optional)
	  @returns     The difference between tic-time and
	               current time in milliseconds.
	*/
	double toc( const char *msg = nullptr )
	{
		gettimeofday(&t_toc, nullptr);
		double diff_msec = (t_toc.tv_usec - t_tic.tv_usec)*1e-3 +
			(t_toc.tv_sec  - t_tic.tv_sec)*1000.0;
		double diff_sec  = diff_msec*1e-3;
		if( out ){
			if( msg ){
				*out << msg << ": ";
			}
			*out << diff_msec << " ms elapsed ("
			     << diff_sec << " s).\n";
		}
		return diff_msec;
	}

	/// Enables output stream and sets it to o.
	void enable_output( std::ostream* o )
	{ out = o; }

	/*!
	  @brief Disables the output stream.

	  @warning After calling this, out is lost!
	*/
	void disable_output()
	{ out = nullptr; }
	

private:
	std::ostream *out;
	timeval t_tic, t_toc;

	void init_tic_toc()
	{
		gettimeofday(&t_tic, nullptr);
		gettimeofday(&t_toc, nullptr);
	}
};

#endif // MY_TIMER_HPP
