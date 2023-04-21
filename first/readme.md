# prediction at the fixed point

## Files

run1.m  first.m  runinv1.m  invfirst.m  firstA.m  invfirstP2.m

## Compute trial function of tatal cases

We decompose the time domain into training and test intervals, and develope a trial function for the total number
of cases using NN. The parameters involved in NN are estimated to minimize the error function over the training domain.
Thereafter, we calculate the derivatives of that trial function analytically.

** run **  run1.m to get the output p_T_test.txt 

p_T_test.txt : contains estimated values of the NN parameters involved in the trial function 


##  Outcomes of Inverse Problem

Solving the inverse problem for the differential equation, second order differential equation discussed in the model section, using the trial function 
and its derivatives we obtain trial functions for the time dependent model parameters.

** run **  runinv1.m to get the output  p_T_inv.txt

p_T_inv.txt : contains estimated values of the NN parameters invloved in the model parameters and error functions.

** run ** firstA.m to generate Figure 1 of the article

** run ** invfirstP2.m to generate Figures 2, 3, 4 of the article
