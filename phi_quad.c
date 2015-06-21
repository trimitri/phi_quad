/**
 * This script does some simple Markov Chain Monte Carlo simulations on a Phi^4
 * model. The metropolis Markov Chain algorithm is used.
 */
#include <stdlib.h>       // calloc()
#include <math.h>         // sqrt()
#include "dSFMT/dSFMT.h"  // random numbers
#include "phi_quad.h"     // my header

// length of a stack array; won't work on heap arrays!
#define LENGTH(stack_array) (sizeof(stack_array) / sizeof(stack_array[0]))

typedef unsigned long int ulong;
typedef unsigned long int ullong;
typedef unsigned int uint;
typedef unsigned short int ushort;

static const ulong kSeed = 123456;         // it's seed value
static const uint kDim[3] = {16, 16, 16};  // lattice dimensions
static const ullong kSweeps = (ullong)1e4;
static const ulong kThermalisationSweeps = 1000;
// scale the normal distribution used for Metropolis proposals. The value is
// chosen (by hand) so that acceptance rate is about 0.5
static const double kProposalScale = 1.25;
static const double kLambda = 1.;
static const double k2Kappa = 2.;

static dsfmt_t rng;                   // random number generator
static double *lat;                   // main system
static double *A_1, *A_2;             // observables
static ulong accepted_proposals = 0;  // # of acc. metropolis proposals
static ulong volume;                  // number of lattice sites
static ulong *steps;                  // {L_1, L_1*L_2, L_1*L_2*L_3*...}

int main() {
  SeedRNG(kSeed);
  steps = (ulong *)malloc(LENGTH(kDim) * sizeof(ulong));
  PopulateStepSizes(&steps);

  // calculate lattice site count
  volume = 1;
  for (uint i = 0; i < LENGTH(kDim); i++) {
    volume *= kDim[i];
  }

  // "calloc" automatically initializes to zero
  lat = (double *)calloc(volume, sizeof(double));
  A_1 = (double *)calloc(kSweeps, sizeof(double));
  A_2 = (double *)calloc(kSweeps, sizeof(double));

  // thermalisation (1000 is far more than necessary, 100 would be enough)
  for (uint i = 0; i < kThermalisationSweeps; i++) {
    SweepSequential();
  }

  // main simulation
  for (uint i = 0; i < kSweeps; i++) {
    SweepSequential();
    Observe(i);
  }

  printf("acceptance rate after %.1e sweeps: %f\n", (double)kSweeps,
         (double)accepted_proposals / 2 / kSweeps);

  free(steps);
  free(lat);
  free(A_1);
  free(A_2);
  return 0;
}

// As there are only two lattice sites, we hand-unroll the loop that one would
// expect here.
inline void SweepSequential() {
  for (ulong i = 0; i < volume; i++) {
    Propagate(i);
  }
}

inline void Propagate(ulong site) {
  double u = lat[site];

  // make a proposal on how to alter lattice site state
  double delta = GaussRandomNumber() * kProposalScale;  // change in phi
  double n = u + delta;                                 // new value

  // decide if we want to accept it

  // gauss term
  double gauss_diff = delta * (k2Kappa * SumNeighbors(site) - n - u);

  // quartic term (not dependend on neigbors)
  double tmp = u * u - 1;
  double quartic_diff = tmp * tmp;
  tmp = n * n - 1;
  quartic_diff -= tmp * tmp;
  quartic_diff *= kLambda;

  double energy_diff = gauss_diff + quartic_diff;

  // update site with probablility min(1, exp(energy_diff));
  if (RandDouble() < exp(energy_diff)) {
    lat[site] = n;
    accepted_proposals++;
  }
}

inline void Observe(ulong index) {
  // TODO
}

// As we just plain sum the elements (no Kahan summation or anything), care
// needs to be taken on large systems.
inline double Average(const double *array, const ulong length) {
  double sum = 0.;
  for (ulong i = 0; i < length; i++) {
    sum += array[i];
  }
  return sum / length;
}

double BinnedStatisticalError(const double *array, const ulong length,
                              const ulong bin_size) {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wbad-function-cast"
  ulong n_bins = (ulong)ceil((double)length / bin_size);
#pragma GCC diagnostic pop
  double *local_averages = (double *)calloc(n_bins, sizeof(double));
  for (ulong i = 0; i < length; i++) {
    local_averages[i / bin_size] += array[i];
  }
  for (ulong i = 0; i < n_bins; i++) {
    local_averages[i] /= bin_size;
  }
  double err = Variance(local_averages, n_bins) / n_bins;
  free(local_averages);
  return sqrt(err);
}

double Variance(const double *array, const ulong length) {
  double average = Average(array, length);
  double sum_of_squares = 0.;
  for (ulong i = 0; i < length; i++) {
    sum_of_squares += array[i] * array[i];
  }
  sum_of_squares /= length;
  return sum_of_squares - average * average;
}

void SeedRNG(uint seed) { dsfmt_init_gen_rand(&rng, seed); }

// we use the Marsaglia polar method, taken from
// https://en.wikipedia.org/wiki/Marsaglia_polar_method
double GaussRandomNumber() {
  static char hasSpare = 0;
  static double spare;

  if (hasSpare == 1) {
    hasSpare = 0;
    return spare;
  }

  hasSpare = 1;
  static double u, v, s;
  do {
    u = dsfmt_genrand_close1_open2(&rng) * 2 - 3;
    v = dsfmt_genrand_close1_open2(&rng) * 2 - 3;
    s = u * u + v * v;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
  } while (s >= 1 || s == 0);
#pragma GCC diagnostic pop

  s = sqrt(-2.0 * log(s) / s);
  spare = v * s;
  return u * s;
}

// uses dSFMT from included library
inline double RandDouble() { return dsfmt_genrand_close_open(&rng); }

// Does no checks at all.
void PersistArray(double *array, long long unsigned int element_count,
                  const char *file_name) {
  printf("persisting %i doubles\n", (int)element_count);
  FILE *file = fopen(file_name, "wb");
  fwrite(array, sizeof(double), sizeof(double) * element_count, file);
  fclose(file);
}

inline void PopulateStepSizes(ulong **step_sizes) {
  ulong tmp_step = 1;
  for (uint i = 0; i < LENGTH(kDim); i++) {
    tmp_step *= kDim[i];
    *step_sizes[i] = tmp_step;
  }
}

inline double SumNeighbors(ulong site) {
  double sum = 0;
  sum += lat[site - site % steps[0] + (site + steps[0] - 1) % steps[0]];
  sum += lat[site - site % steps[0] + (site + 1) % steps[0]];
  sum += lat[site - site % steps[1] + (site + steps[1] - steps[0]) % steps[1]];
  sum += lat[site - site % steps[1] + (site + steps[0]) % steps[1]];
  sum += lat[site - site % steps[2] + (site + steps[2] - steps[1]) % steps[2]];
  sum += lat[site - site % steps[2] + (site + steps[1]) % steps[2]];
  return sum;
}
