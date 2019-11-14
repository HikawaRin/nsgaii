// Declaration of MockDataSet
// this file contain a Mock data set for test
# ifndef NSGAII_SOURCE_MOCKDATASET_CPP
# define NSGAII_SOURCE_MOCKDATASET_CPP

# include "../Header/IGetDataAble.hpp"

class MockDataSet : public IGetDataAble{
public:
    void GetData(std::vector<DataForm> &dataForms);
};

void MockDataSet::GetData(std::vector<DataForm> &dataForms){

}

# endif // NSGAII_SOURCE_MOCKDATASET_CPP