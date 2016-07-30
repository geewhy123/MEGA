/****************************************
 *file:GA.hh
 * header containing objective and constraint
 * function prototypes
 * 
 * Gary Yan, Peter Zhang
 * June 13,2012
 ****************************************/


#ifndef GA_H
#define GA_H



double f1(vec x, void *P);

double f2(vec x, void *P);

double g1(vec x, void *P);

double g2(vec x, void *P);

double h1(vec x, void *P);

double h2(vec x, void *P);


#endif
