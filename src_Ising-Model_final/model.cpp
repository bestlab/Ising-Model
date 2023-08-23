/* model.cpp Created by Layne Frechette on Feb. 17, 2021
 *
 * Implements DCA model class.
 */

#include "model.h"
#include <stdlib.h>
#include <vector>


model::model(){
  lambda = 0.01;
}

model::model(int N1, double lam, bool symon, std::string mc_init1){
  N = N1;
  lambda = lam;
  avg_ene = 0;
  num_seqs = 100;
  seqs.resize(num_seqs);
  h.zeros(N);
  J.zeros(N*(N-1)/2);
  mom1.zeros(N);
  mom2.zeros(N*(N-1)/2);
  mom1_err=0;
  mom2_err=0;
  mc_init = mc_init1;
  symmetrize_on=symon;
}

model::~model(){
}

double model::get_energy(std::vector<int> &seq){

  double energy = 0;

  // loop through all contacts
  for(int i=0; i<N; i++){
    if (seq[i] == 1){
      energy += h(i,0);
    }
  }

  // loop through all unique contact pairs
  for(int i=0; i<(N-1); i++){
    for(int j=(i+1); j<N; j++){
      int index = (N-1)*i-i*(i+1)/2+j-1;
      if(seq[i]==1&&seq[j]==1){
       energy += J(index);
      }
    }
  }
  return energy;
}

double model::get_delta_energy(std::vector<int> &seq, int m, int r){

  double dE = 0.;

  // dE = E - E0 = h - 0 = h
  if (r == 1){
    dE += h(m);
  // dE = E - E0 = 0 - h = -h
  } else {
    dE -= h(m);
  }

  // loop through all contacts before site m
  // where switched contact resides
  for (int i=0; i<m; i++) {
      int index = (N-1)*i-i*(i+1)/2+m-1;
      // dE = E - E0 = (1,1) - (0,1) = J - 0 = J
      if(r == 1&&seq[i]==1){
	      dE += J(index);
      // dE = E - E0 = (0,1) - (1,1) = 0 - J = -J
      } else if (r==0&&seq[i]==1) {
	      dE -= J(index);
      }
  }

  // loop through all contacts after site m
  // where switched contact resides
  for(int j=(m+1); j<N; j++){
      int index = (N-1)*m-m*(m+1)/2+j-1;
      // dE = E - E0 = (1,1) - (0,1) = J - 0 = J
      if(r == 1&&seq[j]==1){
	      dE += J(index);
      // dE = E - E0 = (0,1) - (1,1) = 0 - J = -J
      } else if (r==0&&seq[j]==1){
	      dE -= J(index);
      }
  }

  return dE;
}

double model::get_Z(){

  //Uses code from https://www.geeksforgeeks.org/print-all-sequences-of-given-length/

  double Z=0;
  std::vector<int> seq(N);

  while(true){
    double ene = get_energy(seq);
    Z += exp(-ene);
    for(int i=0; i<N; i++){
      printf("%d ", seq[i]);
    }
    printf("\n");
    int p=N-1;
    while(seq[p]==1) p--;
    if(p<0) break;
    seq[p]=seq[p]+1;
    for(int i=p+1; i<N; i++) seq[i]=0;
  }

  return Z;
}
