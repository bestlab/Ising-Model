#include "rand.h"
#include <gsl/gsl_rng.h>

void init_rng(long unsigned int seed){
  
	const gsl_rng_type *T;
	T = gsl_rng_default;
	rg = gsl_rng_alloc(T);
	srand((unsigned) seed);
	gsl_rng_env_setup();
	gsl_rng_set(rg, seed);
}
