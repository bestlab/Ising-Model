#ifndef MSA_H
#define MSA_H

#include "model.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <armadillo>

/** Convert amino acid character to integer index.
 * @param a, Amino acid character.
 * @param q, Number of states. */
int letter_to_number(char a, int q);

std::string number_to_letter(int a, int q);

void read_Nq(std::string msafile, int *N, int *q);

/** Read in an MSA and convert to ints.
 * @param msafile, Name of the MSA file.
 * @param *N, Pointer to sequence length.
 * @param *q, Pointer to number of states. */
std::vector<int> read_msa(std::string msafile, int N, int q);

std::vector<double> get_weights(std::vector<int> &msa_seqs, int nseq, double delta);

void lowmem_fill_freq(std::string msafile, arma::mat &msa_freq, arma::cube &msa_corr, int N);

/** Compute single-site frequencies and pair correlations in MSA.
 * @param msa_seqs, Vector of MSA sequences.
 * @param msa_freq, Arma matrix of frequencies.
 * @param msa_corr, Arma matrix of correlations.
 * @param N, Sequence length.
 * @param q, Number of states.
 * @param nseq, Number of sequences in MSA.
 * @param delta, Reweighting threshold. */
void fill_freq(std::vector<int> &msa_seqs, std::vector<double> &weights,
		arma::mat &msa_freq, arma::cube &msa_corr, int nseq, bool symon);

/** Get Hamming distance between two sequences.
 * @param v1, First sequence.
 * @param v2, Second sequence. */
double get_hamming_dist(std::vector<int> v1, std::vector<int> v2);

#endif
