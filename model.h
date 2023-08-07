/* model.h Created by Layne Frechette on Feb. 17, 2021
 *
 * Implements DCA model class.
 */

#ifndef MODEL_H
#define MODEL_H

//Headers
#include <armadillo>

//Class for model
class model {
  public:

    int N; //Length of sequences
    double lambda; //l2 regularization constant
    double avg_ene; //Average energy estimated by most recent MC run
    bool symmetrize_on;
    int num_seqs;
    std::vector<std::string> seqs; //Sample of sequences obtained by most recent MC run
    arma::mat h; //Single-site propensities
    arma::cube J; //Coupling parameters
    arma::mat mom1; //First moment
    arma::cube mom2; //Second moment
    double mom1_err; //Std error in first moment via block averaging
    double mom2_err; //Std error in second moment via block averaging
    std::string mc_init;

    //Constructors
    model();
    model(int N1, double lam, bool symon, std::string mc_init1);
    //Destructor
    ~model();

    //Other Methods

    /** Get the energy of a sequence.
     * @param seq, Protein sequence. */
    double get_energy(std::vector<int> &seq);

    double get_delta_energy(std::vector<int> &seq, int m, int r);

    double get_energy_single(std::vector<int> &seq, int well);

    /** Get the energy of a selected site within a sequence.
     * @param seq, Protein sequence.
     * @param m, Index of site. */
    double get_single_site_energy(std::vector<int> &seq, int m);

    /** Compute partition function. */
    /* WARNING: ONLY USE FOR TOY MODELS! Computing Z is wicked expensive. */
    double get_Z();

    /** Compute energies of each sequence in an MSA.
     *
    void calc_energy_dists(const char *name);
    */

    /** Convert to zero-sum gauge. */
    void convert_to_zero_sum();

    /** Print current parameters to file. */
    void dump_params(int iter);
};

#endif
