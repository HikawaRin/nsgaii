/* Declaration of NSGAII */

# ifndef NSGAII_SOURCE_NSGAII_CPP
# define NSGAII_SOURCE_NSGAII_CPP

# include <iostream>
# include <cstdlib>

# include "../Header/nsgaii.hpp"
# include "./population.cpp"
# include "./parsefile.cpp"
# include "./error.cpp"

# include "./debug.cpp"

using namespace std;

void NSGAII::ComputeObj(std::vector<Individual*> &inds){
    // 初始化dataForms
    std::vector<DataForm> dataForms;
    for (int i = 0; i < inds.size(); i++){
        DataForm d;
        d.xreal = inds[i]->xreal;
        d.xbin = inds[i]->xbin;
        d.gene = inds[i]->gene;

        dataForms.push_back(d);
    }
    // 计算目标值
    NSGAII:_dataBase->GetData(dataForms);
    // 将目标值写入种群
    for (int i = 0; i < dataForms.size(); i++){
        inds[i]->obj = dataForms[i].obj;
        inds[i]->constr = dataForms[i].constr;
    }
} // void NSGAII::ComputeObj()

void NSGAII::Evolution(){
    std::vector<Individual*> inds;
    for (int i = 0; i < this->_populations.size(); i++){
        for (int j = 0; j < this->_populations[i]->ind.size(); j++){
            inds.push_back(this->_populations[i]->ind[j]);
        }
    }
    // 计算Obj及Constr_violation
    NSGAII::ComputeObj(inds);
    inds.clear();
    if (CoreNum == 1){
        // 单线程计算
        std::vector<Individual*> childInd;
        for (int i = 0; i < this->_populations[0]->ind.size(); i++){
            Individual *child = new Individual(*(this->_populations[0]->ind[i]));
            childInd.push_back(child);
        }
        for (int i = 0; i < this->Params.ngen; i++){
            cout << "Gen: " << i << endl;
            NSGAII::_populations[0]->Evolution(childInd);
            
            for (auto i : childInd){
                inds.push_back(i);
            }
            this->ComputeObj(inds);
            inds.clear();
            
            this->_populations[0]->RefreshPop(childInd);
            this->_populations[0]->ReportPop(AllPopPath);
        }
        this->_populations[0]->RecordBest(BestPopPath);
    }else{
        // 创建CoreNum个线程，每个线程分配一个种群进行隔离进化
    }
    
} // void NSGAII::Evolution(Population &p)

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

NSGAII::NSGAII(IGetDataAble *dataBase){
    NSGAII::_dataBase = dataBase;

    params p = LoadSetting(SettingPath);
    if (_checkparams(p)){
        NSGAII::Params = p;
    }
    ToString(NSGAII::Params);

    for (int i = 0; i < CoreNum; i++){
        Population *pop = new Population(p);
        NSGAII::_populations.push_back(pop);
    }
    // 初始化结果文件
    ofstream file;
    string s = "";
    s+= "objectives,constraints,real_var,bin_var,constr_violation,rank,crowding_distance";
    file.open(InitialPopPath, ios::out | ios::trunc);
    if (file.is_open()){
        file << s << endl;
    }else{
        RaiseError("InitialPop file open fail, please check log path");
    }
    file.close();
    file.open(FinalPopPath, ios::out | ios::trunc);
    if (file.is_open()){
        file << s << endl;
    }else{
        RaiseError("FinalPop file open fail, please check log path");
    }
    file.close();
    file.open(BestPopPath, ios::out | ios::trunc);
    if (file.is_open()){
        file << s << endl;
    }else{
        RaiseError("BestPop file open fail, please check log path");
    }
    file.close();
    file.open(AllPopPath, ios::out | ios::trunc);
    if (file.is_open()){
        file << s << endl;
    }else{
        RaiseError("AllPop file open fail, please check log path");
    }
    file.close();
} // NSGAII::NSGAII()

# endif // NSGAII_SOURCE_NSGAII_CPP