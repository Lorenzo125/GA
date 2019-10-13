#ifndef Config_h
#define Config_h

#include "ParametersDomain.h"

class Config {
public:
  int PopulationSize;
  double MutationRate;
  double KeepFraction;
  ParametersDomain ParDomain;
  
  int NumberOfParameters() const {
    return ParDomain.NumberOfParameters();
  }
};

#endif
