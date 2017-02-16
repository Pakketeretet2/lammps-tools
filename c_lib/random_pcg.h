/* ----------------------------------------------------------------------
   PCG random number generator

   Paper:
   - http://www.pcg-random.org/pdf/toms-oneill-pcg-family-v1.02.pdf
   - https://web.archive.org/web/20150817092619/http://www.pcg-random.org/pdf/toms-oneill-pcg-family-v1.02.pdf
   According to website submitted to ACM Transactions on Mathematical
   Software and as of August 2015 under review.

   Code contributed by Stefan Paquay and Maarten Sebregts
   (Eindhoven University of Technology).
------------------------------------------------------------------------- */

#ifndef RANPCG_H
#define RANPCG_H

#include "stdint.h"

class RanPCG {
 public:
  RanPCG(int seed, int stream=0);
  ~RanPCG();

  void reset(int seed, int stream=0);

  double uniform() { return uniform12() - 1.0; }
  double uniform_lo_hi(double lo, double hi) { return uniform()*(hi-lo)+lo; }
  double uniform12();
  uint32_t uniformi();

  double gaussian();

private:
  // for gaussian ziggurat algorithm
  double gaussian_zigg_tail();
  void generate_ziggurat();
  void set_layers(double);
  double my_inv_gauss(double);
  double my_erfc(double);

  // pcg
  uint64_t state;
  uint64_t inc;

  // ziggurat
  int nlayers;
  double *layer_x, *layer_y; // Layer width and height
  double A; // Layer area
};

#endif
