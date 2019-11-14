/* Declaration of NSGAII class */

# ifndef NSGAII_HEADER_NSGAII_HPP
# define NSGAII_HEADER_NSGAII_HPP

# include "./global.hpp"
# include "./params.hpp"
# include "./population.hpp"
# include "./IGetDataAble.hpp"

class NSGAII{
    public:
        NSGAII(IGetDataAble *dataBase);
        // 计算所有种群中每个个体的目标函数
        void ComputeObj();
        // 种群进化
        void Evolution();

        params Params;
    private:
        int nbinmut;
        int nrealmut;
        int nbincross;
        int nrealcross;

        // 算法使用的种群,不同的种群间相互隔离
        std::vector<Population*> _populations;
        IGetDataAble *_dataBase;
        // Return True if params legal
        // 检查params是否符合要求
        bool _checkparams(params _params);
}; // Class NSGAII

# endif // NSGAII_HEADER_NSGAII_HPP