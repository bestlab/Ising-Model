#include "rand.h"
#include <gsl/gsl_rng.h>

void init_rng(long unsigned int seed){
  
	const gsl_rng_type *T;
	T = gsl_rng_default;
	rg = gsl_rng_alloc(T);
	srand((unsigned) seed);
	gsl_rng_env_setup();
	gsl_rng_set(rg, seed);

	rg_replica=new gsl_rng*[nrep];
	for (int r = 0; r< nrep; r++) {
		rg_replica[r] = gsl_rng_alloc(T);
		gsl_rng_set(rg_replica[r], seed + r*10000); // perhaps there is a better way!
	}
}
