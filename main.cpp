# include <iostream>

# include "Source/nsgaii.cpp"
# include "Source/parsefile.cpp"
# include "Source/MockDataSet.cpp"
# include "Source/debug.cpp"
# include "Source/log.cpp"

using namespace std;

int main(){
    Log("Program Start");
    MockDataSet mock;
    NSGAII nsga2(&mock);
    Log("Load complete");
    nsga2.Evolution();

    return 0;
}
