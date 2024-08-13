#include <TMB.hpp>
#include <vector>
#include <iostream>
#include <numeric>
#include <cmath>
#include "init.h"
#include "acreTMB.h"



template<class Type>
Type objective_function<Type>::operator() ()
{
   return acreTMB(this);

   return 0;
}

  
