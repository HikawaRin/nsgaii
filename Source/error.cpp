/* Declaration for handle error*/
# ifndef NSGAII_SOURCE_ERROR_CPP
# define NSGAII_SOURCE_ERROR_CPP

# include <iostream>
// Interupt program 
// 用于在程序中出现异常时输出异常信息
// 根据实际需求更改函数中的内容
void RaiseError(const char* s){
    std::cout << s << std::endl;

    exit(1);
}

# endif // NSGAII_SOURCE_ERROR_CPP