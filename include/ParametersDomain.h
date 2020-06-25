#ifndef Parameter_h
#define Parameter_h

#include <string>
#include <vector>

class Parameter {

public:
   Parameter()
   {
      m_name    = "";
      m_min_val = 0.;
      m_max_val = 0.;
   };

   std::string m_name;
   double      m_min_val;
   double      m_max_val;
};

class ParametersDomain {
public:
   ParametersDomain();

   ParametersDomain(size_t);

   void setParameterDomain(size_t, std::string, double, double);

   size_t setNumberOfParameters() const;

   size_t getNumberOfParameters() const;

   double getParamaterMin(size_t i);

   double getParamaterMax(size_t i);

   Parameter &operator[](size_t i) { return m_par[i]; };

private:
   std::vector<Parameter> m_par;
};

#endif
