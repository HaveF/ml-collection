Matlab implementation for nonlinear robust fuzzy PCA algorithm which was first introduced
in publication:

P. Luukka, 'A New Nonlinear Fuzzy Robust PCA Algorithm and Similarity
classifier in classification of medical data sets', 
International Journal of Fuzzy Systems, Vol. 13, No. 3, September 2011
pp. 153-162.

Function call:

[scores,loads]=nfrpca(x,LV,method)

INPUTS:
1)x: your data matrix samples in rows and variables in columns.
2)LV: How many variables to use from data, if not specified all variables are used.
3)method:  Optional, will be selected automatically. By writing 
 'nfrpca1' current method is selected. In current version only this one is implemented and will be selected
 automatically.

OUTPUTS:

scores:    The principal component scores; that is, the representation of X in the principal component space. 
	   Rows of SCORE correspond to observations, columns to components.
loads:     principal component coefficients also known as loadings.

One can use this function in a following way:

1) load your data in your matlab

For example in this case write load exampledata3.txt
which should give you following data:

exampledata3 =

    0.4600    0.3400    0.1400    0.0300
    0.5000    0.3400    0.1500    0.0200
    0.4400    0.2900    0.1400    0.0200
    0.7600    0.3000    0.6600    0.2100
    0.4900    0.2500    0.4500    0.1700
    0.7300    0.2900    0.6300    0.1800

  
there we have an artificial data where in four columns we
have data measurements from six samples. 

2) run the program by writing i.e.

[scores,loads]=nfrpca(exampledata3,4,'nfrpca1')

You will get scores and loads i.e. (depend on parameters you chooce)

scores =

   -0.5171   -0.2720   -0.0730   -0.0308
   -0.5499   -0.2841   -0.0733   -0.0067
   -0.4866   -0.2394   -0.0602   -0.0019
   -1.0710    0.0095    0.0141   -0.0023
   -0.7283   -0.0001   -0.0002   -0.0602
   -1.0228    0.0031   -0.0029    0.0099


loads =

   -0.7113   -0.4603    0.2561    0.4654
   -0.2898   -0.4887   -0.5115   -0.6447
   -0.6134    0.7343   -0.2857   -0.0541
   -0.1837    0.1011    0.7689   -0.6040


Notice that in lines 40 to 43 you have parameters which you can vary to
alter the results. Defaults now being:
a=10;     % parameter value in monotonic function y=tanh(a*x)
ALFA0=1; % learning coefficient (0,1]
ETA=.1;  % soft threshold, a small positive number
EXPO=1.5;  % fuzziness variable
 

