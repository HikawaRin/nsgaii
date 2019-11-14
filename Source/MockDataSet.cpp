// Declaration of MockDataSet
// this file contain a Mock data set for test
# ifndef NSGAII_SOURCE_MOCKDATASET_CPP
# define NSGAII_SOURCE_MOCKDATASET_CPP

# include <cmath>

# include "../Header/IGetDataAble.hpp"

class MockDataSet : public IGetDataAble{
public:
    void GetData(std::vector<DataForm> &dataForms);
};

void MockDataSet::GetData(std::vector<DataForm> &dataForms){
    double g;
    for (int i = 0; i < dataForms.size(); i++){
        g = 1.0 + dataForms[i].xreal[1];
        dataForms[i].obj.push_back(dataForms[i].xreal[0]);
        dataForms[i].obj.push_back(g*exp(-dataForms[i].obj[0]/g));
        dataForms[i].constr.push_back(dataForms[i].obj[1]/(0.858*exp(-0.541*dataForms[i].obj[0]))-1.0);
        dataForms[i].constr.push_back(dataForms[i].obj[1]/(0.728*exp(-0.295*dataForms[i].obj[0]))-1.0);
    }
}

# endif // NSGAII_SOURCE_MOCKDATASET_CPP