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

// the stat. error of the mean
double BinnedStatisticalError(const double *array,
                              const unsigned long int length,
                              const unsigned long int bin_size);

// Do some basic analysis, printing average, stat. error and tau_int estimate
void BinningAnalysis(const double *data, const unsigned long int length,
                     const unsigned long int bin_size);

// Sweep the lattice typewriter style.
void SweepSequential(void);

// Do a Markov chain step on the given lattice point.
void Propagate(unsigned long int lattice_site);

// Measure observables and store to the given place in observable array.
void Observe(unsigned long int observation_index);

// Seeds the pseudo-random number generator.
void SeedRNG(unsigned int seed);

// write double array to disk as matlab-compatible file
void PersistArray(double *array, long long unsigned int element_count,
                  const char *file_name);

// Sum sites directly neigboring the given site, using periodic boundary
// conditions.
double SumNeighbors(unsigned long int site);

// Sum sites directly neigboring the given site, only taking the neighbors in
// positive direction of each dimension into account. This is to avoid
// double-counting of neigbour pairs.
double SumPositiveNeighbors(unsigned long int site);
