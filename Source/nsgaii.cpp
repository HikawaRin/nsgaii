/* Declaration of NSGAII */

# include <iostream>
# include <cstdlib>

# include "../Header/nsgaii.hpp"
# include "../Header/params.hpp"
# include "./parsefile.cpp"
# include "./error.cpp"

# include "./debug.cpp"

using namespace std;

NSGAII::NSGAII(){
    params p = LoadSetting(SettingPath);
    if (_checkparams(p)){
        NSGAII::_params = p;
    }
    ToString(NSGAII::_params);
} // NSGAII::NSGAII()

// Return True if params legal 
bool NSGAII::_checkparams(params _params){
    if (_params.popsize < 4 || (_params.popsize % 4) != 0){
        RaiseError("Wrong population size entered, hence exiting");
    }
    if (_params.ngen < 1){
        RaiseError("Wrong nuber of generations entered, hence exiting");
    }
    if (_params.nobj < 1){
        RaiseError("Wrong number of objectives entered, hence exiting");
    }
    if (_params.ncon < 0){
        RaiseError("Wrong number of constraints enetered, hence exiting");
    }
    if (_params.nreal < 0){
        RaiseError("Wrong number of variables entered, hence exiting");
    }
    if (_params.pcross_real < 0.0 || _params.pcross_real > 1.0){
        RaiseError("Entered value of probability of crossover of real variables is out of bounds, hence exiting");
    }
    if (_params.pmut_real < 0.0 || _params.pmut_real > 1.0){
        RaiseError("Entered value of probability of mutation of real variables is out of bounds, hence exiting");
    }
    if (_params.nreal != 0 && _params.eta_c <= 0){
        RaiseError("Wrong value of distribution index for crossover entered, hence exiting");
    }
    if (_params.nreal != 0 && _params.eta_m <= 0){
        RaiseError("Wrong value of distribution index for mutation entered, hence exiting");
    }
    if (_params.nbin < 0){
        RaiseError("Wrong number of binary variables entered, hence exiting");
    }
    if (_params.nbin != 0 && (_params.pcross_bin < 0.0 || _params.pcross_bin > 1.0)){
        RaiseError("Entered value of probability of crossover of binary variables is out of bounds, hence exiting");
    }
    if (_params.nbin != 0 && (_params.nreal == 0 && _params.nbin == 0)){
        RaiseError("Entered value of probability  of mutation of binary variables is out of bounds, hence exiting");
    }
    if (_params.nreal == 0 && _params.nbin == 0){
        RaiseError("Number of real as well as binary variables, both are zero, hence exiting");
    }

    return true;
} // bool NSGAII::_checkparams(params _params)