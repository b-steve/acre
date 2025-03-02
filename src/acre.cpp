#include <fstream>
#include <vector>
#include <iostream>
#include <TMB.hpp>
#include <numeric>
#include <cmath>
#include <limits>
#include "init.h"
#include "acreTMB.h"



template<class Type>
Type objective_function<Type>::operator() ()
{
   return acreTMB(this);

   return 0;
}

  
