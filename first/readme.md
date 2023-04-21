# prediction at the fixed point



## Files

**run1.m**  : run the function first.m 

**first.m**  : Estimating NN parameters of the trial function associated with total cases

**runinv1.m**  : run the function invfirst.m

**invfirst.m**  : Estamating NN parameters of the trial functions associated with model parameters

**firstA.m**  : generating figure

**invfirstP2.m** : generating figures



## Compute trial function of tatal cases

We decompose the time domain into training and test intervals, and develope a trial function for the total number
of cases using NN. The parameters involved in NN are estimated to minimize the error function over the training domain.
Thereafter, we calculate the derivatives of that trial function analytically.

**run** run1.m to get the output p_T_test.txt 

p_T_two.txt : contains estimated values of the NN parameters involved in the trial function 



##  Outcomes of Inverse Problem

Solving the inverse problem for the differential equation, second order differential equation discussed in the model section, using the trial function 
and its derivatives we obtain trial functions for the time dependent model parameters.

**run**  runinv1.m to get the output  p_T_inv.txt

p_inv_two.txt : contains estimated values of the NN parameters invloved in the model parameters and error functions.



## Figures

**run** firstA.m to generate Figure 1 of the article

**run** invfirstP2.m to generate Figures 2, 3, 4 of the article
