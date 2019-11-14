/* Declaration of function for parse file*/

# ifndef NSGAII_SOURCE_PARSEFILE_CPP
# define NSGAII_SOURCE_PARSEFILE_CPP

# include <fstream>
# include <string>
# include <regex>

# include "../Header/params.hpp"
# include "./error.cpp"

using namespace std;
params LoadJson(string path);
params Loadin(string path);

// Return params filled with data in file
// Support .json .in
params LoadSetting(string path){
    size_t found = path.find_last_of(".");
    string _filetype = path.substr(found + 1);
    
    if (_filetype == "json"){
        // 配置文件为json格式
        return LoadJson(path);
    }else if(_filetype == "in"){
        // 配置文件为in格式
        return Loadin(path);
    }else{
        string s = "Unsupport Setting File Type: ";
        s += _filetype;
        RaiseError(s.c_str());
    }
    params p;
    return p;
}

// Return params filled with json data
params LoadJson(string path){
    params _params;

    return _params;
}

// Return params filled with data write in .in
params Loadin(string path){
    ifstream file;
    file.open(path);
    if ((file.rdstate() & std::ifstream::failbit) != 0){
        RaiseError("Open Setting File Failure, Check Setting File");
    }
    
    params _params;
    char buffer[256];
    int numLines = 0;
    while (!file.eof()){
        file.getline(buffer, 100);
        switch (numLines)
        {
        case 0:
            _params.popsize = atoi(buffer);
            break;
        case 1:
            _params.ngen = atoi(buffer);
            break;
        case 2:
            _params.nobj = atoi(buffer);
            break;
        case 3:
            _params.ncon = atoi(buffer);
            break;
        case 4:
            _params.nreal = atoi(buffer);
            break;
        case 5:
            if (_params.nreal != 0){
                for (int i = 0; i < _params.nreal; i++){
                    if (i != 0){
                        file.getline(buffer, 100);
                    }
                    string s = buffer;
                    size_t found = s.find(" ");
                    double min, max;
                    min = stod(s.substr(0, found));
                    s = s.substr(found + 1);
                    found = s.find(" ");
                    max = stod(s.substr(0, found));
                    if (min > max){
                        RaiseError("Wrong limits entered for the min and max bounds of real variable, hence exiting");
                    }
                    _params.min_realvar.push_back(min);
                    _params.max_realvar.push_back(max);
                } // for (int i = 0; i < _params.nreal; i++)
                break;
            }else{
                numLines++;
            } // if (_params.nreal != 0)
        case 6:
            if (_params.nreal == 0){
                _params.pcross_real = 0.0;
                numLines++;
            }else{
                _params.pcross_real = strtod(buffer, NULL);
                break;
            }
        case 7:
            if (_params.nreal == 0){
                _params.pmut_real = 0.0;
                numLines++;
            }else{
                _params.pmut_real = strtod(buffer, NULL);
                break;
            }
        case 8:
            if (_params.nreal == 0){
                _params.eta_c = 0.0;
                numLines++;
            }else{
                _params.eta_c = strtod(buffer, NULL);
                break;
            }
        case 9:
            if (_params.nreal == 0){
                _params.eta_m = 0.0;
                numLines++;
            }else{
                _params.eta_m = strtod(buffer, NULL);
                break;
            }
        case 10:
            _params.nbin = atoi(buffer);
            break;
        case 11:
            if (_params.nbin == 0){
                numLines++;
            }else{
                for (int i = 0; i < _params.nbin; i++){
                    if (i != 0){
                        file.getline(buffer, 100);
                    }
                    string s = buffer;
                    size_t found = s.find(" ");
                    int bit;
                    double min, max;
                    bit = stoi(s.substr(0, found));
                    if (bit < 1){
                        RaiseError("Wrong number of bits for binary variable entered, hence exiting");
                    }
                    s = s.substr(found + 1);
                    found = s.find(" ");
                    min = stod(s.substr(0, found));
                    s = s.substr(found + 1);
                    found = s.find(" ");
                    max = stod(s.substr(0, found));
                    if (min > max){
                        RaiseError("Wrong limits entered for the min and max bounds of binary variable entered, hence exiting");
                    }
                    _params.nbits.push_back(bit);
                    _params.min_binvar.push_back(min);
                    _params.max_binvar.push_back(max);
                } // for (int i = 0; i < _params.nbin; i++)
                break;
            } // if (_params.nbin == 0)
        case 12:
            if (_params.nbin == 0){
                _params.pcross_bin = 0.0;
                numLines++;
            }else{
                _params.pcross_bin = strtod(buffer, NULL);
                break;
            }
        case 13:
            if (_params.nbin == 0){
                _params.pmut_bin = 0.0;
                numLines++;
            }else{
                _params.pmut_bin = strtod(buffer, NULL);
                break;
            }
        case 14:
            _params.IPAddress = buffer;
            break;
        case 15:
            _params.port = atoi(buffer);
            break;
        default:
            break;
        } // switch (numLines)
        numLines++;
        if (numLines > 15) break;
    } // while (!file.eof())

    file.close();
    return _params;
}

# endif // NSGAII_SOURCE_PARSEFILE_CPP