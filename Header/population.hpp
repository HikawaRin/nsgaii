/* Definition of population class */

# ifndef NSGAII_HEADER_POPULATION_HPP
# define NSGAII_HEADER_POPULATION_HPP

# include "params.hpp"

class Population{
public:
    Population(double seed);
private:
    std::vector<individual> ind;
}; // struct population

# endif // NSGAII_HEADER_POPULATION_HPP