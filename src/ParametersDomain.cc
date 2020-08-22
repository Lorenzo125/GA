#include "ParametersDomain.h"

ParametersDomain::ParametersDomain() : m_par(0){};

ParametersDomain::ParametersDomain(size_t s) : m_par(s){};

void ParametersDomain::setParameterDomain(size_t i, std::string vname, double vmin, double vmax)
{
   m_par[i].m_name    = vname;
   m_par[i].m_min_val = vmin;
   m_par[i].m_max_val = vmax;
};

double ParametersDomain::getParameterMin(size_t i)
{
   return m_par[i].m_min_val;
};

double ParametersDomain::getParameterMax(size_t i)
{
   return m_par[i].m_max_val;
};

size_t ParametersDomain::getNumberOfParameters() const
{
   return m_par.size();
};
