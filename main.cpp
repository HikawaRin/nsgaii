# include <iostream>

# include "Header/nsgaii.hpp"
# include "Source/parsefile.cpp"
# include "Source/debug.cpp"

using namespace std;

int main(){
    params p = Loadin(SettingPath);
    ToString(p);
    
    return 0;
}