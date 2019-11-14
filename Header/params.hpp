/* Declaration for nsgaii params */

# ifndef NSGAII_HEADER_PARAMS_HPP
# define NSGAII_HEADER_PARAMS_HPP

# include <vector>
# include <string>

struct params
{
    // 种群的大小
    int popsize;
    // 迭代的次数
    int ngen;
    // 目标函数的个数
    int nobj;
    // 约束的数目
    int ncon;
    // 实数编码变量的个数
    int nreal;
    // 实数变量的最小值
    std::vector<double> min_realvar;
    // 实数变量的最大值
    std::vector<double> max_realvar;
    // 实数交叉概率
    double pcross_real;
    // 实数变异概率
    double pmut_real;
    // 实数编码的交叉的参数
    double eta_c;
    // 实数编码的多项式变异的参数
    double eta_m;
    // 二进制编码变量的个数
    int nbin;
    // 二进制交叉的次数
    std::vector<int> nbits;
    // 二进制变量的最小值
    std::vector<double> min_binvar;
    // 二进制变量的最大值
    std::vector<double> max_binvar;
    // 二进制交叉概率
    double pcross_bin;
    // 二进制变异概率
    double pmut_bin;

    std::string IPAddress;
    int port; 
}; // Struct params

# endif // NSGAII_HEADER_PARAMS_HPP