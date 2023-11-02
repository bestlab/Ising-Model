/* mc.cpp Created by Layne Frechette on Feb. 17, 2021
 *
 * Implements Monte Carlo moves for sampling under given parameters.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <armadillo>
#include <sstream>

//User-defined headers
#include "mc.h"
#include "msa.h"
#include "model.h"
#include "rand.h"

// number of times to run sweep for equilibrium sequences
int nequil=1000;

/***************/
/*** Methods ***/
/***************/

void run_mc_traj(model &model, int n){
  //dump_freq = mc_steps/ number sequences
  //sequence number % dump_freq gives which sequence to save
  int dump_freq = n/model.num_seqs;
  int nblock=10;
  model.avg_ene=0;

  model.mom1.zeros();
  model.mom2.zeros();
  model.mom1_err=0;
  model.mom2_err=0;

  std::vector<arma::vec> mom1_temp(nblock);
  std::vector<arma::vec> mom2_temp(nblock);
  for(int i=0; i<nblock; i++){
    mom1_temp[i] = arma::zeros(model.N);
    mom2_temp[i] = arma::zeros(model.N*(model.N-1)/2);
  }

  int block_cnt=0;
  std::vector<int> seq(model.N,0.0);

  //Initialize a random sequence
  if(model.mc_init=="random"){
    for(int i=0; i<model.N; i++){
        seq[i] = (int)gsl_rng_uniform_int(rg, 2);
    }
  }

  //Equilibrate ***NOTE: IF MC_INIT NOT "random", WILL NOT EQUILIBRATE!!!!!***
  if(model.mc_init=="random"){
    for(int t=0;  t<nequil; t++){
        do_sweep(model, seq, 1.0, rg);
    }
  }

  //Production
  int dump_cnt=0;
  std::vector<int> my_seq = seq;

  // Perform n monte carlo steps
  for(int t=0; t<n; t++){

    // Writes monte carlo sequences every dump_freq iterations
    if(t%dump_freq==0){
    std::stringstream seqstr;
    for(int j=0; j<my_seq.size(); j++){
        seqstr << number_to_letter(my_seq[j]);
    }
    model.seqs[dump_cnt] = seqstr.str();
    dump_cnt++;
    }

    // calculate new sequence
    do_sweep(model, my_seq, 1.0, rg);
    model.avg_ene += model.get_energy(my_seq);

	      seq = my_seq;
    //Updates moments over sequence
    for(int i=0; i<model.N; i++){
        if (seq[i]==1){
            model.mom1(i)++;
            mom1_temp[block_cnt](i)++;
        }
    }
    for(int i=0; i<(model.N-1); i++){
        for(int j=(i+1); j<model.N; j++){
            int index = (model.N-1)*i-i*(i+1)/2+j-1;
            if (seq[i]==1&&seq[j]==1){
                model.mom2(index)++;
                mom2_temp[block_cnt](index)++;
            }
        }
    }

    if((t+1)%(n/nblock)==0){
        block_cnt++;
    }
  }

  // Compute average of summed moments
  model.mom1 /= n;
  model.mom2 /= n;
  model.avg_ene /= n;

  // Compute error of moments
  for(int i=0; i<nblock; i++){
    mom1_temp[i] /= (n/nblock);
    mom2_temp[i] /= (n/nblock);

    model.mom1_err += arma::accu(arma::square(mom1_temp[i]-model.mom1))/model.mom1.n_elem;
    model.mom2_err += arma::accu(arma::square(mom2_temp[i]-model.mom2))/model.mom2.n_elem;
  }
  model.mom1_err /= (nblock-1);
  model.mom2_err /= (nblock-1);
  model.mom1_err = sqrt(model.mom1_err);
  model.mom2_err = sqrt(model.mom2_err);
}

void do_sweep(model &model, std::vector<int> &seq, int T, gsl_rng *rng){
  // loops over entire list of contacts
  for(int i=0; i<model.N; i++){
    unsigned long int s = gsl_rng_uniform_int(rng, (unsigned long int)model.N);
    // initialize contact to be present
    unsigned long int r = 1;
    // if contact present in sequence change contact to not present
    if (seq[s]==1) r = 0;
    //calculate change in energy as a result of contact switch
    double dE = model.get_delta_energy(seq, (int)s, (int)r);
    double prob = std::min(1.0, exp(-dE/T));
    double xsi = gsl_rng_uniform(rng);
    // acceptance probability = prob
    // if random number xsi is within acceptance probability then accept
    if(prob>xsi) seq[s] = (int)r;
  }
}

double get_boltzmann(model &model, std::vector<int> &seq, int m, int r){

  double dE = model.get_delta_energy(seq, m, r);
  return exp(-dE);
}
