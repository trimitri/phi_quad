Guohua Maier (538652), Franz Gutsch 534348

# CP III Uebung 8: Metropolis on 3D lattice
In this exercise we should calculate the expectation value of given observables
using the Metropolis algorithm.
Those values should be calculated for $λ=1.1$ and different $κ$.
Then we should compare those values with a series expansion code given on the
course web page.
We should find out, for which kappa, then sum is still a good approximation.

In comparison to the last exercise, here we are in a 3D lattice, so we have
$16^3$ lattice points instead of two.
This time we have 6	neighbours insted of 1.
The expectationvalue calculated are also different ones, see exercise sheet.

## Our Results
Below are printed some results for $10^6$ measuring sweeps after $10^5$
thermalisation sweeps on a $16^3$ lattice. The parameters were $λ=1.1$ and
$κ=0.1$.

    acceptance rate after 1.0e+06 sweeps: 0.551963

    A_1:analysing at bin size 10...
    value, error and tau_int est. : 0.086571 0.000023 1169403
    A_1:analysing at bin size 100...
    value, error and tau_int est. : 0.086571 0.000024 1317348
    A_1:analysing at bin size 1000...
    value, error and tau_int est. : 0.086571 0.000025 1435412
    A_1:analysing at bin size 5000...
    value, error and tau_int est. : 0.086571 0.000027 1599951
    A_1:analysing at bin size 10000...
    value, error and tau_int est. : 0.086571 0.000026 1494324

    A_2:analysing at bin size 10...
    value, error and tau_int est. : 0.769726 0.001830 1411599
    A_2:analysing at bin size 100...
    value, error and tau_int est. : 0.769726 0.002006 1696694
    A_2:analysing at bin size 1000...
    value, error and tau_int est. : 0.769726 0.002042 1758566
    A_2:analysing at bin size 5000...
    value, error and tau_int est. : 0.769726 0.002149 1947515
    A_2:analysing at bin size 10000...
    value, error and tau_int est. : 0.769726 0.002024 1727588

    real    4m46.333s
    user    4m46.300s

Although the binning error estimate converges as expected, the approximated
integrated autocorrelation times $τ_\mathrm{int}$ seem to be off.

## Series Expansion vs. Monte-Carlo Simulation
Here we compare the expectationvalues calculated using metropolis with the
expectationvalue using the sum for different values of kappa.

*$κ=0.01$*

    Simulation:
    acceptance rate after 1.0e+04 sweeps: 0.699786
    measured values: A2= 0.564316
    Series approximation:
    last sum-term 3.72781e-26 
    A2 = 0.566723 

*$κ=0.05$*

    Simulation:
    acceptance rate after 1.0e+04 sweeps: 0.692052
    measured values: A2= 0.779794
    Series approximation:
    last sum-term 3.55512e-12 
    A2 = 0.767858 

*$κ=0.1$*

    Simulation:
    acceptance rate after 1.0e+04 sweeps: 0.672449
    measured values:  A2= 1.3062
    Series approximation:
    last sum-term 3.72781e-06 
    A2 = 1.31429 

*$κ=0.15$*

    Simulation:
    acceptance rate after 1.0e+04 sweeps: 0.642251
    measured values:    A2= 3.75251
    Series approximation:
    last sum-term 0.0123959 
    A2 = 3.57713 

*$κ=0.2$*

    Simulation:
    acceptance rate after 1.0e+04 sweeps: 0.525654
    measured values:      A2= 1448.81
    Series approximation:
    last sum-term 3.9089 
    A2 = 40.6817

From the results above, one can see that the sum is useful up to for the value
of $κ=0.15$.
For this value, the error is about 5%. One also notices the last sum-term is
still very small ($10^{-2}$).
For $κ=0.2$ the results of the sum are useless (the last sum-term ist
only one magnitude smaller than the result).





