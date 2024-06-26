#ifndef RUN_ISING_H
#define RUN_ISING_H

#include <armadillo>

extern gsl_rng *rg;

void init_default(model &mymodel, arma::mat &msa_freq);

void fit(model &mymodel, arma::mat &msa_freq, arma::mat &msa_corr);

void read_inputs(std::string in_name);

void read_params(model &mymodel, std::string in_name);

void write_params(model &mymodel, std::string out_name);

void write_restart(model &mymodel, std::string out_name);

void write_seqs(model &mymodel, std::string out_name);

arma::mat get_norm(model &mymodel);

#endif
