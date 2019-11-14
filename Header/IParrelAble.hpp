// Definition of Interface IParrelAble

# ifndef NSGAII_HEADER_IPARRELABLE_HPP
# define NSGAII_HEADER_IPARRELABLE_HPP

# include <vector>
# include "./population.hpp"

// Interface IParrelAble
// 该接口提供并行计算的功能
class IParrelAble
{
public:
    // 并行计算多个种群
    // 返回值表示计算的状态
    // 0： 并行计算正常且全部完成
    int ParrelCompute(std::vector<Population>& populations);
}; // Interface IParrelAble

# endif // NSGAII_HEADER_IPARRELABLE_HPP