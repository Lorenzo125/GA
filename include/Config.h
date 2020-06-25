#ifndef Config_h
#define Config_h

#include "ParametersDomain.h"
#include <cmath>

class Config {

public:
   void setPopulationSize(int t) { m_population_size = t; };

   void setMutationRate(double t) { m_mutation_rate = t; };

   void setKeepFraction(double t) { m_keep_fraction = t; };

   void setParameterDomain(const ParametersDomain &t) { m_par_domain = t; };

   int getPopulationSize() const { return m_population_size; };

   double getMutationRate() { return m_mutation_rate; };

   double getKeepFraction() { return m_keep_fraction; };

   ParametersDomain getParamaterDomain() { return m_par_domain; };

   int getKeep() { return m_keep; };

   double getProb(int t) { return m_prob[t]; };

   int getNumberOfParameters() const { return m_par_domain.getNumberOfParameters(); }

   void setConfigCrossover()
   {
      int pool = m_keep_fraction * m_population_size;
      m_keep   = floor(pool); /// arrotonda per difetto
      double k = 0.0;
      // calcolo le probabilità sequenziali
      k = m_keep * (m_keep + 1) / 2;
      m_prob.push_back(m_keep / k);
      for (int i = 2; i <= m_keep; i++) {
         m_prob.push_back((m_keep - i + 1) / k + m_prob[i - 2]);
      };
   };

private:
   int                 m_keep;
   int                 m_population_size;
   std::vector<double> m_prob;
   double              m_mutation_rate;
   double              m_keep_fraction;
   ParametersDomain    m_par_domain;
};

#endif
