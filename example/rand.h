#ifndef RAND_H
#define RAND_H

#include <gsl/gsl_rng.h>

extern int nrep;
extern gsl_rng *rg;
extern gsl_rng **rg_replica;

void init_rng(long unsigned int seed);

#endif
