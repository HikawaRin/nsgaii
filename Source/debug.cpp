// Declaration for debug file

# ifndef NSGAII_SOURCE_DEBUG_CPP
# define NSGAII_SOURCE_DEBUG_CPP

# include <iostream>

# include "../Header/params.hpp"

using namespace std;

template<typename T> void ToString(vector<T> v){
    for (int i = 0; i < v.size(); i++){
        cout << v[i];
        if (i != v.size() - 1){
            cout << ", ";
        }
    }
    cout << endl;
}

void ToString(params _params){
    cout << "popsize: " << _params.popsize << endl;
    cout << "ngen: " << _params.ngen << endl;
    cout << "nobj: " << _params.nobj << endl;
    cout << "ncon: " << _params.ncon << endl;
    cout << "nreal: " << _params.nreal << endl;
    cout << "min_realvar: ";
    ToString<double>(_params.min_realvar);
    cout << "max_realvar: ";
    ToString<double>(_params.max_realvar);
    cout << "pcross_real: " << _params.pcross_real << endl;
    cout << "pmut_real: " << _params.pmut_real << endl;
    cout << "eta_c: " << _params.eta_c << endl;
    cout << "eta_m: " << _params.eta_m << endl;
    cout << "nbin:" << _params.nbin << endl;
    cout << "nbits: ";
    ToString<int>(_params.nbits);
    cout << "min_binvar: ";
    ToString<double>(_params.min_binvar);
    cout << "max_binvar: ";
    ToString<double>(_params.max_binvar);
    cout << "pcross_bin: " << _params.pcross_bin << endl;
    cout << "pmut_bin: " << _params.pmut_bin << endl;

    cout << "ipaddress: " << _params.IPAddress << endl;
    cout << "port: " << _params.port << endl; 
}

# endif // NSGAII_SOURCE_DEBUG_CPP