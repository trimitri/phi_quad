Guohua Maier (538652), Franz Gutsch 534348

# CP III Uebung 8: Monte-Carlo Simulation with Metropolis algorithm for 3 periodic 3d lattice

## Problem Description
In this exercise we should calculate the expectation value usinig the Metropolis
algorithm. Those values should calculated for lambda=1.1 and different kappa.
Then we should compare those values with the values calculated using the same value 
of kappa using the code giving by Dr. Hasenbusch. We should find out, for which kappa, then
sum still makes sense.



## Solution
In comparison to the last exercise, here we are in a 3d  lattice so we have instead of 
lattice points 16^3 lattice points. This time we have 6	neighbours insted of 1. The expectationvalue
calculated are also different ones, see exercise sheet.








## comparison of the values
Here we compare the expectationvalues calculated using metropolis with the expectationvalue using 
the sum for different values of kappa.

guohua@pool11:~/cp3/vl8/UE08>gcc_release.elf 
give kappa 
0.01
acceptance rate after 1.0e+04 sweeps: 0.699786
measured values: A2= 0.564316
guohua@pool11:~/cp3/vl8/UE08>a.out 
give kappa 
0.01
last sum-term 3.72781e-26 
A2 = 0.566723 


guohua@pool11:~/cp3/vl8/UE08>gcc_release.elf 
give kappa 
0.05
acceptance rate after 1.0e+04 sweeps: 0.692052
measured values: A2= 0.779794
guohua@pool11:~/cp3/vl8/UE08>a.out 
give kappa 
0.05
last sum-term 3.55512e-12 
A2 = 0.767858 

guohua@pool11:~/cp3/vl8/UE08>gcc_release.elf 
give kappa 
0.1
acceptance rate after 1.0e+04 sweeps: 0.672449
measured values:  A2= 1.3062
guohua@pool11:~/cp3/vl8/UE08>a.out 
give kappa 
0.1
last sum-term 3.72781e-06 
A2 = 1.31429 

guohua@pool11:~/cp3/vl8/UE08>gcc_release.elf 
give kappa 
0.15
acceptance rate after 1.0e+04 sweeps: 0.642251
measured values:    A2= 3.75251
guohua@pool11:~/cp3/vl8/UE08>a.out 
give kappa 
0.15
last sum-term 0.0123959 
A2 = 3.57713 



guohua@pool11:~/cp3/vl8/UE08>gcc_release.elf 
give kappa 
0.2
acceptance rate after 1.0e+04 sweeps: 0.525654
measured values:      A2= 1448.81
guohua@pool11:~/cp3/vl8/UE08>a.out 
give kappa 
0.2
last sum-term 3.9089 
A2 = 40.6817

From the results above, one can see that the sum is useful up to for the value of kappa=0.15.
For this value, the error is about 5%. One also notice the last sum-term is still very small( 10^-2)
For kappa=0.2 the results of the sum are useless(see also the last sum-term ist only 1 magnetude smaller than the result)





