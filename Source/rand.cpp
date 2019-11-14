// Declaration of random function

# ifndef NSGAII_SOURCE_RAND_CPP
# define NSGAII_SOURCE_RAND_CPP

# include <random>
// 将用于获得随机数引擎的种子
std::random_device rd;
// 以 rd() 播种的标准 mersenne_twister_engine
std::mt19937 gen(rd());

#endif // NSGAII_SOURCE_RAND_CPP