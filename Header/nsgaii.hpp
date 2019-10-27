/* Declaration of NSGAII class */

# ifndef NSGAII_HEADER_NSGAII_HPP
# define NSGAII_HEADER_NSGAII_HPP

# include "global.hpp"
# include "population.hpp"

class NSGAII{
    public:
        NSGAII();
    private:
        int nbinmut;
        int nrealmut;
        int nbincross;
        int nrealcross;
        
        params _params;

        bool _checkparams(params _params);
};

# endif // NSGAII_HEADER_NSGAII_HPP