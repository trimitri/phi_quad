/**
 * This script does some simple Markov Chain Monte Carlo simulations on a 1x2
 * Phi^4 model. The metropolis Markov Chain algorithm is used.
 */
#include <stdlib.h>       // calloc()
#include <math.h>         // sqrt()
#include "dSFMT/dSFMT.h"  // random numbers
#include "phi_quad.h"     // my header

typedef unsigned long int ulong;
typedef unsigned int uint;

static double *lat;  // main system, has only two lattice points
static double *A_1, *A_2, *A_3, *A_4, *A_5;  // observables
static dsfmt_t rng;                          // the random number generator
static const ulong kSeed = 123456;           // it's seed value
static const unsigned long long kSweeps = (unsigned long long)1e9;

// scale the normal distribution used for Metropolis proposals. The value is
// chosen (by hand) so that acceptance rate is 0.5
static double kProposalScale = 1.246;
static ulong accepted_proposals = 0;

int main(/*int arg_count, char **arg_values*/) {
  /*
  if (arg_count == 2) {
    kProposalScale = strtod(arg_values[1], (char**)NULL);
  }
  */

  // "calloc" automatically initializes to zero
  lat = (double *)calloc(2, sizeof(double));
  /*
  A_1 = (double *)calloc(kSweeps, sizeof(double));
  A_2 = (double *)calloc(kSweeps, sizeof(double));
  A_3 = (double *)calloc(kSweeps, sizeof(double));
  A_4 = (double *)calloc(kSweeps, sizeof(double));
  A_5 = (double *)calloc(kSweeps, sizeof(double));
  */

  SeedRNG(kSeed);

  // thermalisation (1000 is far more than necessary, 100 would be enough)
  for (uint i = 0; i < 1000; i++) {
    SweepSequential();
  }

  // main simulation
  for (uint i = 0; i < kSweeps; i++) {
    SweepSequential();
//    Observe(i);
  }

  printf("acceptance rate after %.1e sweeps: %f\n", (double)kSweeps,
         (double)accepted_proposals / 2 / kSweeps);
  /*
  puts("                      A_1     A_2     A_3     A_4     A_5");
  printf("measured values:   %.5f %.5f %.5f %.5f %.5f\n", Average(A_1, kSweeps),
         Average(A_2, kSweeps), Average(A_3, kSweeps), Average(A_4, kSweeps),
         Average(A_4, kSweeps));
         */

  /*
    double *bin_errors = (double *)malloc(5 * 1000 * sizeof(double));
  for (uint i = 0; i < 100; i++) {
    printf("%i at %i\n", 10 * (i + 1), 5 * i);
    bin_errors[5 * i + 0] = BinnedStatisticalError(A_1, kSweeps, 10 * (i + 1));
    bin_errors[5 * i + 1] = BinnedStatisticalError(A_2, kSweeps, 10 * (i + 1));
    bin_errors[5 * i + 2] = BinnedStatisticalError(A_3, kSweeps, 10 * (i + 1));
    bin_errors[5 * i + 3] = BinnedStatisticalError(A_4, kSweeps, 10 * (i + 1));
    bin_errors[5 * i + 4] = BinnedStatisticalError(A_5, kSweeps, 10 * (i + 1));
  }
  PersistArray(bin_errors, 5000, "binned_errors.bin");
  free(bin_errors);
  */

//  SolveNumerically();
/*
  double stat_error[3];
  double tau_int[3];

  double raw_error[3];
  raw_error[0] = Variance(A_1, kSweeps);
  raw_error[1] = Variance(A_2, kSweeps);
  raw_error[2] = Variance(A_4, kSweeps);

  ulong bin_sizes[3] = {10, 100, 1000};
  for (uint i = 0; i < 3; i++) {
    stat_error[0] = BinnedStatisticalError(A_1, kSweeps, bin_sizes[i]);
    stat_error[1] = BinnedStatisticalError(A_2, kSweeps, bin_sizes[i]);
    stat_error[2] = BinnedStatisticalError(A_4, kSweeps, bin_sizes[i]);
    for (uint j = 0; j < 3; j++) {
      tau_int[j] = .5 * kSweeps * pow(stat_error[j] / raw_error[j], 2);
    }
    printf("bin size: %lu\n", bin_sizes[i]);
    printf("stat. error:       %.5f %.5f %.5f %.5f %.5f\n", stat_error[0],
           stat_error[1], stat_error[1], stat_error[2], stat_error[2]);
    printf("tau_int estimate:  %.5f %.5f %.5f %.5f %.5f\n", tau_int[0],
           tau_int[1], tau_int[1], tau_int[2], tau_int[2]);
  }
  */
  free(lat);
  free(A_1);
  free(A_2);
  free(A_3);
  free(A_4);
  free(A_5);
  return 0;
}

// As there are only two lattice sites, we hand-unroll the loop that one would
// expect here.
inline void SweepSequential() {
  Propagate(0);
  Propagate(1);
}

inline void Propagate(ulong site) {
  double u, t;  // phi values "us" (this site) and "them" (other site)
  if (site == 0) {
    u = lat[0];
    t = lat[1];
  } else {
    u = lat[1];
    t = lat[0];
  }

  // make a proposal on how to alter lattice site state
  double delta = GaussRandomNumber() * kProposalScale;  // change in phi
  double n = u + delta;                                 // new value

  // decide if we want to accept it

  // the gauss term simplifies a lot:
  double gauss_diff = delta * (t - n - u);  // if we had gauss model

  // the quartic term not so much:
  double tmp = u * u - 1;
  double quartic_diff = tmp * tmp;
  tmp = n * n - 1;
  quartic_diff -= tmp * tmp;

  double energy_diff = gauss_diff + quartic_diff;

  // update site with probablility min(1, exp(energy_diff));
  if (RandDouble() < exp(energy_diff)) {
    lat[site] = n;
    accepted_proposals++;
  }
}

inline void Observe(ulong index) {
  A_1[index] = lat[0] * lat[1];
  A_2[index] = lat[0] * lat[0];
  A_3[index] = lat[1] * lat[1];
  A_4[index] = A_2[index] * A_2[index];
  A_5[index] = A_3[index] * A_3[index];
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

// Omit A3 and A5, as the problem is symmetric!
// For the same reason, only integrating half the space is also possible, but
// less readable, so it is not done here.
void SolveNumerically() {
  double expect_norm = 0, expect_int_A1 = 0, expect_int_A2 = 0,
         expect_int_A4 = 0;
  //   int sample_bins=10000;
  // loop through all reasonable values for phi_1 and phi_2
  for (int i = 0; i < 5; i++) {
    expect_norm = 0, expect_int_A1 = 0, expect_int_A2 = 0, expect_int_A4 = 0;
    for (double phi_1 = -5 + i; phi_1 <= 5 - i;
         phi_1 += 10. / 60) {  //(100-i*10)) {
      for (double phi_2 = -5 - i; phi_2 <= 5 - i;
           phi_2 += 10. / 60) {  //(100-i*10)) {
        double A1 = phi_1 * phi_2;
        double A2 = phi_1 * phi_1;
        double A3 = phi_2 * phi_2;
        double A4 = A2 * A2;
        double boltzmann_factor =
            exp(A1 - A2 - A3 - (A2 - 1) * (A2 - 1) - (A3 - 1) * (A3 - 1));

        expect_norm += boltzmann_factor;
        expect_int_A1 += boltzmann_factor * A1;
        expect_int_A2 += boltzmann_factor * A2;
        expect_int_A4 += boltzmann_factor * A4;
        // 	printf("test %g \n",phi_2);
      }
    }
    expect_int_A1 /= expect_norm;
    expect_int_A2 /= expect_norm;
    expect_int_A4 /= expect_norm;
    //     printf("Number of bins %i X %i     A1		A2	    A3
    //     A4   	      A5 \n",(100-i*10),(100-i*10) );
    printf(
        "range of phi [%i , %i]   	  A1		A2	    A3		"
        "  A4   	      A5 \n",
        -5 + i, 5 - i);

    printf("	numerical result:  %.11f %.11f %.11f %.11f %.11f\n",
           expect_int_A1, expect_int_A2, expect_int_A2, expect_int_A4,
           expect_int_A4);
  }
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
