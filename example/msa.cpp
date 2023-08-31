#include "msa.h"
#include "model.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <numeric>
#include <string>
#include <iostream>
#include <fstream>
#include <armadillo>
#include <mcheck.h>

std::string number_to_letter(int a){

  std::string s;
  if(a==0) s='H';
  else if (a==1) s='P';
  else{
    printf("Error: Too many states for an Ising Model. Exiting.\n");
    exit(-1);
  }
  return s;
}

void read_Nq(std::string msafile, int *N){

  int index = msafile.find(".");
  std::string extension = msafile.substr(index+1);

  if(extension=="dat"){
    std::ifstream stream(msafile);
    std::string firstSeq;
    std::getline(stream,firstSeq);
    *N = firstSeq.length();
    std::cout << "Sequence length: " << *N << std::endl;
  }
  else{
    std::cout << "Error: file extension not recognized. Exiting." << std::endl;
    exit(-1);
  }
}

std::vector<int> read_msa(std::string msafile, int N){

  int index = msafile.find(".");
  std::string extension = msafile.substr(index+1);
  int nseq = 0;

  // Tests if file extension type is correct
  if(extension=="dat"){
    std::string line;
    std::ifstream myfile(msafile);
    while(std::getline(myfile,line)) nseq++;
    std::cout << "No. of sequences: " << nseq << std::endl;

    std::vector<int> msa_seqs(nseq*N);

    int linecount=0;
    std::ifstream file(msafile);
    std::string str;
    while(std::getline(file, str)){
      for(int i=0; i<N; i++){
        msa_seqs[linecount*N+i] = int(str.at(i)-'0');
      }
      linecount++;
    }
    return msa_seqs;
  }
  else{
    printf("Error: unsupported file extension. Exiting.\n");
    exit(-1);
  }
}

std::vector<double> get_weights(std::vector<int> &msa_seqs, int nseq, double delta){

  std::cout << "Constructing sequence weights..." << std::endl;
  std::vector<double> weights(nseq, 0.0);
  int N = msa_seqs.size()/nseq;
  double Meff=0;
  for(int m1=0; m1<nseq; m1++){
    weights[m1] += 1.0;
    if(delta>0.0){
    std::vector<int> seq1(msa_seqs.begin()+m1*N, msa_seqs.begin()+(m1+1)*N);
    for(int m2=0; m2<nseq; m2++){
      if(m1!=m2){
        std::vector<int> seq2(msa_seqs.begin()+m2*N, msa_seqs.begin()+(m2+1)*N);
        double dist = get_hamming_dist(seq1, seq2);
        if(dist<delta*N) weights[m1] += 1.0;
      }
    }
    }
    weights[m1] = 1.0/weights[m1];
  }

  for(int m1=0; m1<nseq; m1++) Meff += weights[m1];
  std::cout << "Computed weights. Meff=" << Meff << std::endl;
  return weights;
}

void lowmem_fill_freq(std::string msafile, arma::mat &msa_freq, arma::mat &msa_corr, int N){

  int linecount=0;
  std::ifstream file(msafile);
  std::string str;
  while(std::getline(file, str)){
    if(linecount%10000==0) std::cout << linecount << std::endl;
    for(int i=0; i<N; i++){
      int q1 = int(str.at(i)-'0');
      if(q1 == 1){
        msa_freq(i,0) += 1.0;
      }
    }
    int l=0;
    for(int i=0; i<N-1; i++){
      for(int j=i+1; j<N; j++){
        int q1 = int(str.at(i)-'0');
        int q2 = int(str.at(j)-'0');
        if (q1==1&&q2==1){
          msa_corr(l,0) += 1.0;
        }
        l++;
      }
    }
    linecount++;
  }
  msa_freq = msa_freq/linecount;
  msa_corr = msa_corr/linecount;
  std::cout << "No. seqs: " << linecount << std::endl;
  std::cout << "Should be N: " << arma::accu(msa_freq) << std::endl;
}

void fill_freq(std::vector<int> &msa_seqs, std::vector<double> &weights,
		arma::mat &msa_freq, arma::mat &msa_corr,int nseq, bool symon){

  //double Meff = std::accumulate(weights.begin(), weights.end(), 0);
  std::cout << nseq << std::endl;
  double Meff=0;
  for(int i=0; i<nseq; i++) Meff+=weights[i];
  std::cout << "Checking Meff: " << Meff << std::endl;
  int N = msa_seqs.size()/nseq;

  //Set freqs and corrs to zero
  msa_freq.zeros();
  msa_corr.zeros();

  //Compute frequencies
  for(int m=0; m<nseq; m++){
    if(m%10000==0) std::cout << m << std::endl;
    for(int i=0; i<N; i++){
      msa_freq(msa_seqs[m*N+i],0) += weights[m];
    }
    for(int i=0; i<(N-1); i++){
      for(int j=(i+1); j<N; j++){
        int index = (N-1)*i-i*(i+1)/2+j-1;
        msa_corr(index,0) += weights[m];
      }
    }
  }

  //Normalize
  msa_freq /= Meff;
  msa_corr /= Meff;
  std::cout << "Computed frequencies." << std::endl;
}

double get_hamming_dist(std::vector<int> v1, std::vector<int> v2){

  if(v1.size()!=v2.size()){
    std::cout << "ERROR: vectors need the same size!" << std::endl;
    exit(-1);
  }
  double dist=0;
  for(int i=0; i<v1.size(); i++){
    if(v1[i]!=v2[i]) dist+=1;
  }
  return dist;
}
