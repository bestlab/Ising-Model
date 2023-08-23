#ifndef RUN_DCA_H
#define RUN_DCA_H

#include <armadillo>

extern gsl_rng *rg;

void init_default(model &mymodel, arma::vec &msa_freq);

void fit(model &mymodel, arma::vec &msa_freq, arma::vec &msa_corr, int nrep);

void read_inputs(std::string in_name);

void read_params(model &mymodel, std::string in_name);

void write_params(model &mymodel, std::string out_name);

void write_restart(model &mymodel, std::string out_name);

void write_seqs(model &mymodel, std::string out_name);

arma::mat get_norm(model &mymodel);

#endif
