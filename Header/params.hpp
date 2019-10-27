/* Declaration for nsgaii params */

# ifndef NSGAII_HEADER_PARAMS_HPP
# define NSGAII_HEADER_PARAMS_HPP

# include <vector>

struct params
{
    int popsize;
    int ngen;
    int nobj;
    int ncon;
    int nreal;
    std::vector<double> min_realvar;
    std::vector<double> max_realvar;
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

    std::string IPAddress;
    int port; 
}; // struct params

struct individual{
    int rank; // 支配等级
    /*  constr_violation 该个体违反约束的情况
        < 0 说明该个体超出了限制条件
        = 0 说明该个体没有超出限制条件
    */
    double constr_violation;
    std::vector<double> xreal;
    std::vector<std::vector<int> > gene;
    std::vector<double> xbin;
    std::vector<double> obj;
    std::vector<double> constr;
    double crowd_dist;
    double power_dist;
    int cnt;
}; // struct individual

# endif // NSGAII_HEADER_PARAMS_HPP