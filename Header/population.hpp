/* Definition of population class */

# ifndef NSGAII_HEADER_POPULATION_HPP
# define NSGAII_HEADER_POPULATION_HPP

# include <string>
# include <vector>
# include "./params.hpp"

// 个体类
class Individual{
public:
    Individual(params p);

    int rank; // 支配等级
    // constr_violation 该个体违反约束的情况
    // < 0 说明该个体超出了限制条件
    // = 0 说明该个体没有超出限制条件
    double constr_violation;
    // 一维数组，染色体实数列表
    std::vector<double> xreal;
    // 二维数组，二进制染色体
    std::vector<std::vector<int> > gene;
    // 一维数组，染色体解码值列表
    std::vector<double> xbin;
    // 一维数组，目标函数值
    std::vector<double> obj;
    // 约束值
    std::vector<double> constr;
    // 个体拥挤距离
    double crowd_dist;
    // 种群中每个个体的距离
    double power_dist;
    int cnt;
private:
    // 对gene进行译码得到xreal
    void _decode(std::vector<int> &nbit, std::vector<double> &min_binvar, std::vector<double> &max_binvar);
}; // Class Individual

// 种群类， 由一组个体类组成
class Population{
public:
    Population(params p);
    // 种群进化一次
    void Evolution();
    // 将种群信息写入文件
    void ReportPop(std::string path);

    // 种群中的个体
    std::vector<Individual*> ind;
    // 当前代数
    int currentGen;
    // 变异概率
    double pcross_real;
    double pcross_bin;
    // 变量的最小值
    std::vector<double> *min_realvar;
    std::vector<double> *min_binvar;
    // 变量的最大值
    std::vector<double> *max_realvar;
    std::vector<double> *max_binvar;

    std::vector<int> *nbits;
    double eta_c; 
private:
    // 计算种群中个体违反约束的情况
    void _computeViolation();
    // 计算支配等级及拥挤度
    void _assign_rank_and_crowding_distance();
    // 种群选择操作
    void _selection(std::vector<Individual*> &childInd);
    // 种群交叉
    void _crossover(Individual *p1, Individual *p2, Individual *c1, Individual *c2, int flag);
}; // Struct population

// 公共方法
// 确定两个个体间的支配关系
int Check_dominance(Individual *a, Individual *b);
// 选出两个个体间处于支配地位的个体
Individual* FetchInd(Individual *a, Individual *b);

# endif // NSGAII_HEADER_POPULATION_HPP