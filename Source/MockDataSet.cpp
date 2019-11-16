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
    // ctp1
    double g;
    for (int i = 0; i < dataForms.size(); i++){
        g = 1.0 + dataForms[i].xreal[1];
        dataForms[i].obj.push_back(dataForms[i].xreal[0]);
        dataForms[i].obj.push_back(g*exp(-dataForms[i].obj[0]/g));
        dataForms[i].constr.push_back(dataForms[i].obj[1]/(0.858*exp(-0.541*dataForms[i].obj[0]))-1.0);
        dataForms[i].constr.push_back(dataForms[i].obj[1]/(0.728*exp(-0.295*dataForms[i].obj[0]))-1.0);
    }
    
    // zdt5
    // int i, j;
    // int u[11];
    // int v[11];
    // double f1, f2, g, h;
    // for (auto f : dataForms){
    //     for (i=0; i<11; i++)
    //     {
    //         u[i] = 0;
    //     }
    //     for (j=0; j<30; j++)
    //     {
    //         if (f.gene[0][j] == 1)
    //         {
    //             u[0]++;
    //         }
    //     }
    //     for (i=1; i<11; i++)
    //     {
    //         for (j=0; j<4; j++)
    //         {
    //             if (f.gene[i][j] == 1)
    //             {
    //                 u[i]++;
    //             }
    //         }
    //     }
    //     f1 = 1.0 + u[0];
    //     for (i=1; i<11; i++)
    //     {
    //         if (u[i] < 5)
    //         {
    //             v[i] = 2 + u[i];
    //         }
    //         else
    //         {
    //             v[i] = 1;
    //         }
    //     }
    //     g = 0;
    //     for (i=1; i<11; i++)
    //     {
    //         g += v[i];
    //     }
    //     h = 1.0/f1;
    //     f2 = g*h;
    //     f.obj.push_back(f1);
    //     f.obj.push_back(f2);
    // }
}

# endif // NSGAII_SOURCE_MOCKDATASET_CPP