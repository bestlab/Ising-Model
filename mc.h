/* mc.h Created by Layne Frechette on Feb. 17, 2021
 *
 * Implements Monte Carlo moves for sampling under given parameters.
 */

#ifndef MC_H
#define MC_H

//Headers
#include <stdlib.h>
#include <armadillo>
#include <gsl/gsl_rng.h>
#include "model.h"

//Variables
extern long unsigned int myseed;
extern int print_freq;

//Functions

/** Run an MC Trajectory.
 * @param model, The model object.
 * @param nstp, Number of steps in trajectory.*/
void run_mc_traj(model &model, int nstp, int nrep);

/** Perform a sweep (model.N attempted MC moves)
 * @param model, The model object.
 * @param seq, Protein sequence. */
//void do_sweep(model &model, std::vector<int> &seq);
void do_sweep(model &model, std::vector<int> &seq, int T, gsl_rng *rng);

/** Get the Boltzmann factor associated with a single-site mutation.
 * @param model, The model object.
 * @param seq, Protein sequence.
 * @param m, Index of site.
 * @param r, Index of proposed new residue. */
double get_boltzmann(model &model, std::vector<int> &seq, int m, int r);

#endif
