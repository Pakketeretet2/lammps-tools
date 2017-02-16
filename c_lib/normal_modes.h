#ifndef NORMAL_MODES_H
#define NORMAL_MODES_H

/*!
  \file normal_modes.h
  @brief Contains some functions used for normal mode analysis.

  \ingroup cpp_lib
*/

#include "types.h"

#include <list>
#include <vector>
#include <array>


extern "C" {

	void normal_mode_analysis( const char *pname, void *eigenvalues,
	                           void *eigenvectors, py_int N );

}

void normal_mode_analysis( class dump_reader &r, py_int N, void *eigenvalues,
                           void *eigenvectors );


void get_normal_modes( const std::list< std::vector<std::array< double, 3> > > &block_data,
                       py_int N, void *eigenvalues, void *eigenvectors,
                       std::vector<double> *X_avg_ptr = nullptr );



#endif // NORMAL_MODES_H
