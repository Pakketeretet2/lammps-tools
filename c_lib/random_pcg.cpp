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

#include "math.h"

#include "random_pcg.h"

static const double MY_PI  = 3.14159265358979323846; // pi
static const double MY_SQRT2 = 1.41421356237309504880; // sqrt(2)
static const double MY_PIS = 1.77245385090551602729; // sqrt(pi)
static const double MY_PI2 = 1.57079632679489661923; // pi/2
static const double MY_2PI = 6.28318530717958647692; // 2pi

template <typename real_type> inline
double sign(real_type x) { return (x < 0) ? -1.0 : 1.0; }

/* ---------------------------------------------------------------------- */

RanPCG::RanPCG(int seed, int stream)
{
  reset(seed, stream);
  nlayers = 0;
}

/* ---------------------------------------------------------------------- */

RanPCG::~RanPCG()
{
  if (nlayers) {
    if (layer_x) delete [] layer_x;
    if (layer_y) delete [] layer_y;
  }
}

/* ---------------------------------------------------------------------- */

void RanPCG::reset(int seed, int stream)
{
  state = static_cast<uint64_t>(seed);
  inc = static_cast<uint64_t>(stream);
  // make sure the increment is odd:
  inc = 2*inc + 1;
  // more or less overcome the fact that the seed will likely be a small number:
  uniform();
  state += seed;
  uniform();
}

/* ----------------------------------------------------------------------
   32-bits uniform RN
------------------------------------------------------------------------- */

inline uint32_t rotate32(uint32_t v, uint32_t r)
{
  // rotate v by r bits
  // https://en.wikipedia.org/wiki/Bitwise_operation#Circular_rotates_in_C_and_C.2B.2B
  return (v >> r) | (v << ((-r) & 31));
}

uint32_t RanPCG::uniformi()
{
  uint64_t oldstate = state;
  // advance state, constant taken from implementation by Melissa O'Neill:
  // https://github.com/imneme/pcg-c/blob/master/include/pcg_variants.h#L253
  // (PCG_DEFAULT_MULTIPLIER_64)
  state = oldstate * 6364136223846793005ULL + inc;
  // calculate output, section 6.3, PCG-XSH-RR
  return rotate32((oldstate ^ (oldstate >> 18)) >> 27, oldstate >> 59);
}

/* ----------------------------------------------------------------------
   uniform RN [1,2)
------------------------------------------------------------------------- */

double RanPCG::uniform12()
{
  // (ab)use IEEE 754 binary64 format
  // Note: this is not portable to systems not using the IEEE754 binary64 format!
  union { uint64_t u; double d; };
  u = static_cast<uint64_t>(uniformi());
  u = (u << 20) | 0x3FF0000000000000ULL;
  return d;
}

/* ----------------------------------------------------------------------
   gaussian RN, mean = 0, var = 1
------------------------------------------------------------------------- */

double RanPCG::gaussian()
{
  double x, y;
  int n;

  if( !nlayers ) generate_ziggurat();

  // For details about the Ziggurat algoritm, see Marsaglia & Tang's paper
  // from 2000 or the Wikipedia article.
  n = static_cast<int>( uniform_lo_hi(0, static_cast<double>(nlayers)) );

  if (n == 0) {
    // Tail of distribution, need something special
    const double bound_x = A / layer_y[0];
    x = uniform_lo_hi( -bound_x, bound_x );

    if (fabs(x) < layer_x[0])
      return x;
    else
      return gaussian_zigg_tail();
  } else if (n == nlayers-1) {
    // Top layer, just generate random x and y and do check on y:
    const double bound_x = layer_x[n-1];
    x = uniform_lo_hi(-bound_x, bound_x);
    y = uniform_lo_hi(layer_y[n-1], layer_y[n]);

    if (y < exp(-0.5*x*x))
      return x;
    else
      return gaussian();
  } else {
    // Generate x and y, and do normal test.
    double bound_x = layer_x[n];
    x = uniform_lo_hi(-bound_x, bound_x);
    if (fabs(x) < layer_x[n+1])
      return x;
    y = uniform_lo_hi(layer_y[n], layer_y[n+1]);
    if (y < exp(-0.5*x*x))
      return x;
    else
      return gaussian();
  }
}

/* ----------------------------------------------------------------------
   Algorithm for the rare case a point from the tail is needed.
------------------------------------------------------------------------- */

double RanPCG::gaussian_zigg_tail()
{
  // See Marsaglia paper for details.
  double s  = sign( -0.5 + uniform() );
  double x0 = layer_x[0];
  double x  = -log(uniform()) / x0;
  double y  = -log(uniform());
  if (2*y > x*x){
    return s*(x0 + x);
  } else {
    return gaussian_zigg_tail();
  }
}

/* ----------------------------------------------------------------------
  Sets up the ziggurat layers for use later.
------------------------------------------------------------------------- */

void RanPCG::generate_ziggurat()
{
  // *** Generation of the ziggurat, as explained by Marsaglia & Tang: ***
  nlayers = 256;
  layer_x = new double[nlayers];
  layer_y = new double[nlayers];

  // From paper of Marsaglia:
  double r = 3.6541528853610088;
  A        = 0.00492867323399;
  set_layers(r);

  layer_x[nlayers-1] = 0;
  layer_y[nlayers-1] = 1.0;
}


/* ----------------------------------------------------------------------
    Sets values for all layers based on given (guessed) x[0].
------------------------------------------------------------------------- */
void RanPCG::set_layers(double xn)
{
  int N = nlayers;
  layer_x[0] = xn;
  layer_y[0] = exp(-xn*xn*0.5);
  // Recursive relations up to layer 0:
  for (int i = 1; i < N; ++i) {
    layer_y[i] = A/layer_x[i-1] + layer_y[i-1];
    const double tmp = -2.0*log( layer_y[i] );
    layer_x[i] = sqrt( tmp );
    double At = layer_x[i-1]*(layer_y[i] - layer_y[i-1]);
  }
}


/* ----------------------------------------------------------------------
   Computes int( exp(-t^2/2) from t = x to t = inf
------------------------------------------------------------------------- */
double RanPCG::my_erfc(double x)
{
  const double erfc_conv_fac = MY_PIS / MY_SQRT2;
  return erfc_conv_fac * erfc( x / MY_SQRT2 );
}


/* ----------------------------------------------------------------------
   Inverse of exp(-x^2/2);
------------------------------------------------------------------------- */
double RanPCG::my_inv_gauss(double x)
{
  return MY_SQRT2 * sqrt( -2.0*log(x) );
}
