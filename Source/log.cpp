// Declaration of log function

# ifndef NSGAII_SOURCE_LOG_CPP
# define NSGAII_SOURCE_LOG_CPP

# include <iostream>
# include <fstream>
# include <string>

# include "../Header/global.hpp"
# include "./error.cpp"

using namespace std;
bool isExist = false;
void Log(string s){
    ofstream file;
    if (!isExist){
        file.open(LogPath, ios::out | ios::trunc);
        isExist = true;
    }else{
        file.open(LogPath, ios::app);
    }
    if (file.is_open()){
        file << s << endl;
    }else{
        RaiseError("Log file open fail, please check log path");
    }

    file.close();
}

#endif // NSGAII_SOURCE_LOG_CPP