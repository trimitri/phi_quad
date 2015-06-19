#pragma once

// Returns a number in [-inf, inf]. Multiple calls will obey gaussian normal
// distribution.
double GaussRandomNumber(void);

// returns a pseudo-random number in [0,1)
double RandDouble(void);

// calculates the arithmetic mean
double Average(const double *array, const unsigned long length);

// calculates the variance (square of standard deviation)
double Variance(const double *array, const unsigned long length);

// the variance of the bin means, normalized by bin size
double BinnedStatisticalError(const double *array,
                              const unsigned long int length,
                              const unsigned long int bin_size);

// Sweep the lattice typewriter style.
void SweepSequential(void);

// Do a Markov chain step on the given lattice point.
void Propagate(unsigned long int lattice_site);

// solve the two-site system numerically
void SolveNumerically(void);

// Measure observables and store to the given place in observable array.
void Observe(unsigned long int observation_index);

// Seeds the pseudo-random number generator.
void SeedRNG(unsigned int seed);

// write double array to disk as matlab-compatible file
void PersistArray(double *array, long long unsigned int element_count,
                  const char *file_name);
