/* Declaration for nsgaii params */

# ifndef NSGAII_PARAMS_HPP
# define NSGAII_PARAMS_HPP

# include <vector>

struct params
{
    int popsize;
    int ngen;
    int nobj;
    int ncon;
    int nreal;
    double min_realvar;
    double max_realvar;
    double pcross_real;
    double pmut_real;
    double eta_c;
    double eta_m;
    int nbin;
    std::vector<int> nbit;
    std::vector<double> min_binvar;
    std::vector<double> max_binvar;
    double pcross_bin;
    double pmut_bin;
};


# endif // NSGAII_PARAMS_HPP