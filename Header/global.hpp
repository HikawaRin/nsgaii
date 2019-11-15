/* This file contains global variable*/

# ifndef NSGAII_HEADER_GLOBAL_HPP
# define NSGAII_HEADER_GLOBAL_HPP
// 文件路径

// 需手动创建输出文件夹，代码无法自动创建文件夹
// 设置参数文件路径
# define SettingPath    "../nsgaii/Input/zdt5.in"
// 存放初始种群数据
# define InitialPopPath "./Out/initial_pop.csv"
// 存放最后一代的种群数据
# define FinalPopPath   "./Out/final_pop.csv"
// 存储最好的一代种群数据
# define BestPopPath    "./Out/best_pop.csv"
// 存储所有的种群数据
# define AllPopPath     "./Out/all_pop.csv"
// 存储该次算法使用的参数 
# define ParamsPath     "./Out/params.csv"
// 存储LOG信息
# define LogPath        "./Out/Log.txt"

// 常数定义

// 使用的计算单元个数
# define CoreNum 1
# define INF 1.0e14
# define EPS 1.0e-14
# define E   2.71828182845905
# define PI  3.14159265358979

# endif // NSGAII_HEADER_GLOBAL_HPP