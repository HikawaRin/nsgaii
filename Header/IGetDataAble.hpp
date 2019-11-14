// Definition of interface IGetDataAble

# ifndef NSGAII_HEADER_IGETDATAABLE_HPP
# define NSGAII_HEADER_IGETDATAABLE_HPP

# include <vector>

# include "./population.hpp"

struct DataForm{
    // 该个体所属的种群
    int popIndex;
    // 种群中个体所属的编号
    int indIndex;
    // 实数值
    std::vector<double> xreal;
    // 解码值
    std::vector<double> xbin;
    // 二进制染色体
    std::vector<std::vector<int> > gene;
    // 目标函数值
    std::vector<double> obj;
    // 约束值
    std::vector<double> constr;
}; // Struct DataForms

// Interface IGetDataAble
// 通过此接口可获得DataForm中的目标函数值
// 具有此接口的类需显式实现GetData函数
class IGetDataAble{
public:
    // 获取对应输入的目标函数值
    // 接受一个DataForm数组作为参数
    virtual void GetData(std::vector<DataForm> &dataForms)=0;
}; // Interface IGetDataAble

# endif // NSGAII_HEADER_IGETDATAABLE_HPP
