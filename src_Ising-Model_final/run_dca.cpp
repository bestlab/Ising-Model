/* run_dca.cpp Created by Layne Frechette on Feb. 17, 2021
 *
 * Run DCA to fit MSA to Potts model.
 */

//Headers
#include "mc.h"
#include "rand.h"
#include "model.h"
#include "msa.h"
#include "run_dca.h"
#include "omp.h"

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <filesystem>
//#include <experimental/filesystem>
//namespace fs = std::experimental::filesystem; //use this for C++14
namespace fs = std::filesystem;
using std::istringstream;

//Global variables

std::string msa_name;
std::string freq_dir="none";
std::string out_name;
std::string input_name;
std::string conf_name;
std::string mc_init;

std::string scratch_dir;
std::string scratch_dir_orig;
std::string output_dir;
std::string folder_name;
std::string folder_name_orig;

gsl_rng *rg;
gsl_rng **rg_replica;

int max_iter;
int mc_steps;
int nrep;
double eps0_h;
double eps0_J;
double eps_inc;
double eps_dec;
double eps_min_h;
double eps_max_h;
double eps_min_J;
double eps_max_J_N;
double lambda;
double gamma_mom;
double cutoff_freq;
double cutoff_rmsd=2e-3;
bool adaptive_stepsize_on;
bool adaptive_sampling_on;
bool symmetrize_on;
bool verbose;
double delta=0.2; //reweighting threshold
int do_lowmem_fill=0;
std::string reg_type="L2";
std::string conv_type="rmsd";
long unsigned int myseed;
int record_freq=10;
int num_restarts=0;
bool is_restart=0;

/*** Main ***/
int main(int argc, char *argv[]){

  //Read command line arguments
  if(argc == 2){
    //Try reading in parameters from conf file
    conf_name = argv[1];
    if(fs::exists(conf_name)){
      read_inputs(conf_name);
    }
    else{
      printf("ERROR: conf file does not exist! Exiting.\n");
      exit(-1);
    }
  }
  else{
    printf("ERROR: conf file required (1 argument).\n");
    exit(-1);
  }

  if(num_restarts!=0) is_restart=1;

  //Set up random number generator
  init_rng(myseed);

  //Memory allocation
  int N;

  std::cout << msa_name << std::endl;
  read_Nq(msa_name, &N);
  model model(N, lambda, symmetrize_on, mc_init);
  std::cout << model.mc_init << std::endl;
  std::cout << "conv type: " << conv_type << std::endl;
  arma::vec msa_freq(model.N,arma::fill::zeros);
  arma::vec msa_corr(model.N*(model.N-1)/2,arma::fill::zeros);

  //Define data output directory
  folder_name_orig = folder_name;
  scratch_dir_orig = scratch_dir;
  output_dir = "data/";
  folder_name += "/lambda=" + std::to_string(lambda);
  folder_name += "_mc_steps_init=" + std::to_string(mc_steps);
  if(adaptive_sampling_on) folder_name += "_as";
  if(adaptive_stepsize_on) folder_name += "_al";
  if(gamma_mom>0.0) folder_name += "_gamma=" + std::to_string(gamma_mom);
  if(reg_type!="L2") folder_name += "_reg_type=" + reg_type;
  if(nrep!=1) folder_name += "_nrep=" + std::to_string(nrep);
  if(conv_type=="max"){
    folder_name += "_cutoff_freq=" + std::to_string(cutoff_freq);
  }
  else if(conv_type=="rmsd"){
    folder_name += "_cutoff_rmsd=" + std::to_string(cutoff_rmsd);
  }
  if(input_name!="none" && !(is_restart)){
    int index1 = input_name.find_last_of("/");
    std::string thefile = input_name.substr(index1+1);
    int index2 = thefile.find(".");
    std::string thename = thefile.substr(0,index2);
    std::cout << "check: " << thename << std::endl;

    folder_name += "_input=" + thename;
  }
  if(is_restart){
    folder_name += "_restart" + std::to_string(num_restarts);
  }

  output_dir += folder_name + "/";
  scratch_dir += folder_name + "/";

  fs::create_directories(output_dir);
  fs::create_directories(scratch_dir);

  //Get frequencies and pair correlations from MSA
  if(freq_dir!="none"){
    //Try reading in 1- and 2-point frequencies
    //from specifically named files in freq_dir
    std::cout << "Reading in MSA frequencies..." << std::endl;
    msa_freq.load(freq_dir + "/stat_align_1p.txt", arma::arma_ascii);
    msa_corr.load(freq_dir + "/stat_align_2p.txt", arma::arma_ascii);
  }
  else{
    if(do_lowmem_fill){
      std::cout << "Doing low-memory frequency fill..." << std::endl;
      lowmem_fill_freq(msa_name, msa_freq, msa_corr, N);
    }
    else{
      std::vector<int> msa_seqs = read_msa(msa_name, N);
      int nseq = msa_seqs.size()/N;
      std::cout << nseq << " sequences" << std::endl;
      std::vector<double> weights = get_weights(msa_seqs, nseq, delta);
      fill_freq(msa_seqs, weights, msa_freq, msa_corr, nseq, symmetrize_on);
    }
    msa_freq.save(scratch_dir + "stat_align_1p.txt", arma::arma_ascii);
    msa_corr.save(scratch_dir + "stat_align_2p.txt", arma::arma_ascii);
  }

  //Initialize parameters
  if(input_name=="none") init_default(model, msa_freq);
  else if(!(is_restart)) read_params(model,input_name);
  else read_params(model,input_name + "/params_curr.txt");

  //Run DCA
  num_restarts += 1;
  fit(model,msa_freq,msa_corr,nrep);

  //Print parameters to file
  write_params(model,output_dir+out_name);

  return 1;
}

/***********
 * METHODS *
 ***********/

void init_default(model &mymodel, arma::vec &msa_freq){
  //Initialize h's to be consistent with single-site frequencies
  for(int i=0; i<mymodel.N; i++){
    mymodel.h(i) = -1*log(msa_freq(i)+0.0001); //Small regularization to avoid divergence
  }
}

void fit(model &mymodel, arma::vec &msa_freq, arma::vec &msa_corr, int nr){

  bool converged = false;
  int niter=0;

  //Initialize adaptive learning rates
  arma::vec alpha_h(mymodel.N,arma::fill::ones);
  arma::vec alpha_J(mymodel.N*(mymodel.N-1)/2,arma::fill::ones);
  alpha_h = eps0_h*alpha_h;
  alpha_J = eps0_J*alpha_J;

  //Derivatives
  arma::vec dh(mymodel.N,arma::fill::zeros);
  arma::vec dJ(mymodel.N*(mymodel.N-1)/2,arma::fill::zeros);

  //Stored derivatives of previous
  arma::vec dh_prev(mymodel.N,arma::fill::ones);
  arma::vec dJ_prev(mymodel.N*(mymodel.N-1)/2,arma::fill::ones);

  //Changes in parameters (for momentum)
  arma::vec change_h(mymodel.N,arma::fill::zeros);
  arma::vec change_J(mymodel.N*(mymodel.N-1)/2,arma::fill::zeros);

  //Average energy
  arma::vec avg_energies(max_iter, arma::fill::zeros);

  //Cross entropy
  arma::vec entropy(max_iter, arma::fill::zeros);

  /* Check whether this run is a restart. If so,
   * read in values of arma objects.
   */
  if(is_restart){
    alpha_h.load(input_name + "/alpha_h_curr.txt", arma::arma_ascii);
    alpha_J.load(input_name + "/alpha_J_curr.txt", arma::arma_ascii);
    dh_prev.load(input_name + "/dh_prev_curr.txt", arma::arma_ascii);
    dJ_prev.load(input_name + "/dJ_prev_curr.txt", arma::arma_ascii);
    change_h.load(input_name + "/change_h_curr.txt", arma::arma_ascii);
    change_J.load(input_name + "/change_J_curr.txt", arma::arma_ascii);
  }

  //Fit
  std::cout << "Fitting model." << std::endl;
  while(!converged){

    std::cout << std::endl << "iteration " << niter << std::endl;
    //Print current parameters to file
    write_params(mymodel,scratch_dir + "params_curr.txt");
    mymodel.mom1.save(scratch_dir + "stat_MC_1p_curr.txt", arma::arma_ascii);
    mymodel.mom2.save(scratch_dir + "stat_MC_2p_curr.txt", arma::arma_ascii);
    alpha_h.save(scratch_dir + "alpha_h_curr.txt", arma::arma_ascii);
    alpha_J.save(scratch_dir + "alpha_J_curr.txt", arma::arma_ascii);
    dh_prev.save(scratch_dir + "dh_prev_curr.txt", arma::arma_ascii);
    dJ_prev.save(scratch_dir + "dJ_prev_curr.txt", arma::arma_ascii);
    change_h.save(scratch_dir + "change_h_curr.txt", arma::arma_ascii);
    change_J.save(scratch_dir + "change_J_curr.txt", arma::arma_ascii);
    write_restart(mymodel,output_dir + "restart.conf");

    run_mc_traj(mymodel,mc_steps,nr);

    avg_energies[niter] = mymodel.avg_ene;
    std::cout << "Avg. energy: " << avg_energies[niter] << std::endl;
    write_seqs(mymodel, scratch_dir + "seqs_" + std::to_string(niter) + ".txt");

    /*** Compute derivatives ***/
    dh = msa_freq-mymodel.mom1;
    dJ = msa_corr-mymodel.mom2;

    // regularization types
    if(reg_type=="L2"){
      dh += 2*mymodel.lambda*mymodel.h;
      dJ += 2*mymodel.lambda*mymodel.J;
    }
    else if(reg_type=="L1"){
      dh += mymodel.lambda*sign(mymodel.h);
      dJ += mymodel.lambda*sign(mymodel.J);
    }

    change_h = gamma_mom*change_h + alpha_h%dh; //% = element-wise matrix multiplication
    change_J = gamma_mom*change_J + alpha_J%dJ;

    //Dump derivatives to file
    if(verbose && niter%record_freq==0){
	  dh.save(scratch_dir + "dh_" + std::to_string(niter) + ".txt", arma::arma_ascii);
	  dJ.save(scratch_dir + "dJ_" + std::to_string(niter) + ".txt", arma::arma_ascii);
    }

    //Compute RMS difference between MSA and MC
    double diff1 = sqrt(arma::accu(arma::square(dh))/mymodel.mom1.n_elem);
    double diff2 = sqrt(arma::accu(arma::square(dJ))/mymodel.mom2.n_elem);
    double max1 = arma::abs(dh).max();
    double max2 = arma::abs(dJ).max();
    std::cout << "RMS of dh, dJ: " << diff1 << " " << diff2 << std::endl;
    std::cout << "Max of dh, dJ: " << max1 << " " << max2 << std::endl;
    if(mymodel.N==2){
      std::cout << "h:" << std::endl;
      std::cout << mymodel.h << std::endl;
      std::cout << "J:" << std::endl;
      std::cout << mymodel.J << std::endl;
      std::cout << "dh:" << std::endl;
      std::cout << dh << std::endl;
      std::cout << "dJ:" << std::endl;
      std::cout << dJ << std::endl;
      std::cout << "model 1p freq:" << std::endl;
      std::cout << mymodel.mom1 << std::endl;
      std::cout << "model 2p freq:" << std::endl;
      std::cout << mymodel.mom2 << std::endl;
    }

    if(niter>max_iter ||
      (conv_type=="max" && max1<cutoff_freq && max2<cutoff_freq) ||
      (conv_type=="rmsd" && (diff1+diff2)<cutoff_rmsd)
      ){
      std::cout << "Converged!" << std::endl;
      std::cout << "RMS of dh, dJ, tot: " << diff1 << " " << diff2 << " " << diff1+diff2 << std::endl;
      std::cout << "Max of |dh|, |dJ|: " << max1 << " " << max2 << std::endl;
      converged=true;
    }

    if(adaptive_sampling_on){
      //Increase amount of sampling if rms derivatives are within
      //standard error of MC sampling.
      if(diff1<mymodel.mom1_err && diff2<mymodel.mom2_err){
        std::cout << "Std. errors: " << mymodel.mom1_err << " " << mymodel.mom2_err << std::endl;
        std::cout << "Model difference is within sampling error. Increasing sampling..." << std::endl;
        mc_steps *= 2;
        std::cout << "MC steps per iteration: " << mc_steps << std::endl;
      }
    }
    //Dump current statistics and parameters to file
    if(verbose && niter%record_freq==0){
	  mymodel.mom1.save(scratch_dir + "stat_MC_1p_" + std::to_string(niter) + ".txt", arma::arma_ascii);
	  mymodel.mom2.save(scratch_dir + "stat_MC_2p_" + std::to_string(niter) + ".txt",arma::arma_ascii);
	  mymodel.h.save(scratch_dir + "h_" + std::to_string(niter) + ".txt", arma::arma_ascii);
	  mymodel.J.save(scratch_dir + "J_" + std::to_string(niter) + ".txt", arma::arma_ascii);
    }

    //Update parameters
    mymodel.h -= change_h;
    mymodel.J -= change_J;

    //Update learning rates **MAKE THIS INTO A SEPARATE FUNCTION**
    if(adaptive_stepsize_on){
      for(int i=0; i<mymodel.N; i++){
          double prod = dh_prev(i)*dh(i);
          if(prod>0 && alpha_h(i)<eps_max_h) alpha_h(i) = eps_inc*alpha_h(i);
          else if(prod<0 && alpha_h(i)>eps_min_h) alpha_h(i) = eps_dec*alpha_h(i);
          else {};
          if(alpha_h(i)>eps_max_h) alpha_h(i)=eps_max_h;
          if(alpha_h(i)<eps_min_h) alpha_h(i)=eps_min_h;
      }
      for(int i=0; i<mymodel.N-1; i++){
        for(int j=i+1; j<mymodel.N; j++){
          int index = (mymodel.N-1)*i-i*(i+1)/2+j-1;
              double prod = dJ_prev(index)*dJ(index);
              if(prod>0 && alpha_J(index)<eps_max_J_N/mymodel.N) alpha_J(index) = eps_inc*alpha_J(index);
              else if(prod<0 && alpha_J(index)>eps_min_J) alpha_J(index) = eps_dec*alpha_J(index);
              else {};
              if(alpha_J(index)>eps_max_J_N/mymodel.N) alpha_J(index)=eps_max_J_N/mymodel.N;
              if(alpha_J(index)<eps_min_J) alpha_J(index)=eps_min_J;
        }
      }
    }

    //Update previous step gradient
    dh_prev = dh;
    dJ_prev = dJ;

    niter++;
  }

  mymodel.mom1.save(output_dir + "stat_MC_1p.txt", arma::arma_ascii);
  mymodel.mom2.save(output_dir + "stat_MC_2p.txt",arma::arma_ascii);

  //Dump energies to file
  avg_energies.save(output_dir + "avg_ene.txt", arma::arma_ascii);

}

void read_inputs(std::string name){

  // default values...
  nrep = 1;
  verbose = 0;

  //Adapted from https://stackoverflow.com/questions/7868936/read-file-line-by-line-using-ifstream-in-c
  std::ifstream file(name);
  if(file.is_open()){
    std::string line;
    std::cout << "INPUT PARAMETERS:" << std::endl;
    while(std::getline(file,line)){
      int index = line.find(":");
      std::string value = line.substr(index+1);
      std::string key = line.substr(0,index);
      std::cout << key << " " << value << std::endl;
      if(key.compare("msa_name")==0) msa_name = value;
      if(key.compare("freq_dir")==0) freq_dir = value;
      if(key.compare("out_name")==0) out_name = value;
      if(key.compare("input_name")==0) input_name = value;
      if(key.compare("folder_name")==0) folder_name = value;
      if(key.compare("scratch_dir")==0) scratch_dir = value;
      if(key.compare("lambda")==0) lambda = std::stof(value);
      if(key.compare("myseed")==0) myseed = (long unsigned)std::stoi(value);
      if(key.compare("max_iter")==0) max_iter = std::stoi(value);
      if(key.compare("mc_steps")==0) mc_steps = std::stoi(value);
      if(key.compare("eps0_h")==0) eps0_h = std::stof(value);
      if(key.compare("eps0_J")==0) eps0_J = std::stof(value);
      if(key.compare("adaptive_stepsize_on")==0) adaptive_stepsize_on = std::stoi(value);
      if(key.compare("adaptive_sampling_on")==0) adaptive_sampling_on = std::stoi(value);
      if(key.compare("eps_inc")==0) eps_inc = std::stof(value);
      if(key.compare("eps_dec")==0) eps_dec = std::stof(value);
      if(key.compare("eps_min_h")==0) eps_min_h = std::stof(value);
      if(key.compare("eps_max_h")==0) eps_max_h = std::stof(value);
      if(key.compare("eps_min_J")==0) eps_min_J = std::stof(value);
      if(key.compare("eps_max_J_N")==0) eps_max_J_N = std::stof(value);
      if(key.compare("gamma")==0) gamma_mom = std::stof(value);
      if(key.compare("cutoff_freq")==0) cutoff_freq = std::stof(value);
      if(key.compare("cutoff_rmsd")==0) cutoff_rmsd = std::stof(value);
      if(key.compare("delta")==0) delta = std::stof(value);
      if(key.compare("do_lowmem_fill")==0) do_lowmem_fill = std::stoi(value);
      if(key.compare("mc_init")==0) mc_init = value;
      if(key.compare("symmetrize_on")==0) symmetrize_on = std::stoi(value);
      if(key.compare("nrep")==0) nrep = std::stoi(value);
      if(key.compare("verbose")==0) verbose = std::stoi(value);
      if(key.compare("record_freq")==0) record_freq = std::stoi(value);
      if(key.compare("reg_type")==0) reg_type = value;
      if(key.compare("conv_type")==0) conv_type = value;
      if(key.compare("num_restarts")==0) num_restarts = std::stoi(value);
    }
    std::cout << std::endl;
    file.close();
  }
}

void write_restart(model &mymodel, std::string name){

  std::ofstream ofile;
  ofile.open(name);
  ofile << "msa_name:" << msa_name << std::endl;
  ofile << "freq_dir:" << freq_dir << std::endl;
  ofile << "out_name:" << out_name << std::endl; //modify this
  ofile << "input_name:" << scratch_dir << std::endl; //for restart file, input_name is a directory
  ofile << "folder_name:" << folder_name_orig << std::endl;
  ofile << "scratch_dir:" << scratch_dir_orig << std::endl;
  ofile << "lambda:" << mymodel.lambda << std::endl;
  ofile << "myseed:" << myseed << std::endl;
  ofile << "max_iter:" << max_iter << std::endl;
  ofile << "mc_steps:" << mc_steps << std::endl;
  ofile << "eps0_h:" << eps0_h << std::endl;
  ofile << "eps0_J:" << eps0_J << std::endl;
  ofile << "adaptive_stepsize_on:" << adaptive_stepsize_on << std::endl;
  ofile << "adaptive_sampling_on:" << adaptive_sampling_on << std::endl;
  ofile << "eps_inc:" << eps_inc << std::endl;
  ofile << "eps_dec:" << eps_dec << std::endl;
  ofile << "eps_min_h:" << eps_min_h << std::endl;
  ofile << "eps_max_h:" << eps_max_h << std::endl;
  ofile << "eps_min_J:" << eps_min_J << std::endl;
  ofile << "eps_max_J_N:" << eps_max_J_N << std::endl;
  ofile << "gamma:" << gamma_mom << std::endl;
  ofile << "cutoff_freq:" << cutoff_freq << std::endl;
  ofile << "cutoff_rmsd:" << cutoff_rmsd << std::endl;
  ofile << "delta:" << delta << std::endl;
  ofile << "do_lowmem_fill:" << do_lowmem_fill << std::endl;
  ofile << "mc_init:" << mc_init << std::endl;
  ofile << "symmetrize_on:" << symmetrize_on << std::endl;
  ofile << "nrep:" << nrep << std::endl;
  ofile << "verbose:" << verbose << std::endl;
  ofile << "record_freq:" << record_freq << std::endl;
  ofile << "reg_type:" << reg_type << std::endl;
  ofile << "conv_type:" << conv_type << std::endl;
  ofile << "num_restarts:" << num_restarts << std::endl;
  ofile.close();
}

void read_params(model &mymodel, std::string in_name){

  std::ifstream file(in_name);
  if(file.is_open()){
    std::string line;
    while(std::getline(file,line)){
      std::string::size_type pos;
      pos = line.find(' ',0);
      std::string data = line.substr(pos+1);
      std::string id = line.substr(0,pos);
      istringstream my_stream(data);
      if(id=="J"){
        int i;
        int j;
        my_stream >> i;
        my_stream >> j;
        int index = (mymodel.N-1)*i-i*(i+1)/2+j-1;
        my_stream >> mymodel.J(index);
      }
      else if(id=="h"){
        int i;
        my_stream >> i;
        my_stream >> mymodel.h(i);
      }
      else{
        std::cout << "Error: unrecognized specifier. Exiting." << std::endl;
        exit(-1);
      }
    }
  }
}

void write_params(model &mymodel, std::string name){

  std::ofstream ofile;
  ofile.open(name);
  for(int i=0; i<mymodel.N-1; i++){
    for(int j=i+1; j<mymodel.N; j++){
      int index = (mymodel.N-1)*i-i*(i+1)/2+j-1;
      ofile << "J " << i << " " << j << " " << mymodel.J(index) << std::endl;
    }
  }

  for(int i=0; i<mymodel.N; i++){
    ofile << "h " << i << " " << mymodel.h(i) << std::endl;
  }
  ofile.close();
}

void write_seqs(model &mymodel, std::string name){

  std::ofstream ofile;
  ofile.open(name);
  for(int i=0; i<mymodel.num_seqs; i++){
    ofile << mymodel.seqs[i] << std::endl;
  }
  ofile.close();
}
