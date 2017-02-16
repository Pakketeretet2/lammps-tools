#ifndef RMSD_FLUCTUATIONS_H
#define RMSD_FLUCTUATIONS_H

#include "types.h"
#include "dump_reader.h"

#include <list>



py_float rmsd_fluctuations_impl( dump_reader &r, std::vector<std::array<double,4> > &rmsd_flucts,
                                 const std::list<long int> *ids, std::list<double> &rmsd_time );

py_float get_msd( const block_data &b1, const block_data &b0, py_float *msd );




#endif // RMSD_FLUCTUATIONS_H
