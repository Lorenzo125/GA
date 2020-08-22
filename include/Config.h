#ifndef Config_h
#define Config_h

#include "ParametersDomain.h"
#include <cmath>

class Config {
   /// Config default constructor
public:
   /// Set the number of chromosomes in the population
   void setPopulationSize(int t) { m_population_size = t; };

   /// Set the mutation rate
   void setMutationRate(double t) { m_mutation_rate = t; };

   /// Set the selection rate
   void setKeepFraction(double t) { m_keep_fraction = t; };

   /// Set the parameters' search domain
   void setParameterDomain(const ParametersDomain &t) { m_par_domain = t; };

   /// Return the number of chromosomes in the population <br> <br>
   /// Referenced by Population::mutation()
   int getPopulationSize() const { return m_population_size; };

   /// Return the mutation rate <br> <br>
   /// Referenced by Population::mutation()
   double getMutationRate() { return m_mutation_rate; };

   /// Return the selection rate <br> <br>
   /// Referenced by Population::crossover()
   double getKeepFraction() { return m_keep_fraction; };

   /// Returnthe domain of parameters <br> <br>
   /// Referenced by Population::init(), Population::mutation()
   ParametersDomain getParameterDomain() { return m_par_domain; };

   /// Return m_keep <br> <br>
   /// Referenced by Population::crossover()
   int getKeep() { return m_keep; };

   /// Return m_prob <br> <br>
   /// Referenced by Population::crossover()
   double getProb(int t) { return m_prob[t]; };

   /// Return the number of parameters involved <br> <br>
   /// Referenced by Population::init(), Population::crossover(), Population::mutation()
   int getNumberOfParameters() const { return m_par_domain.getNumberOfParameters(); }

   /// Set the values for m_keep and m_prob
   void setConfigCrossover()
   {
      int pool = m_keep_fraction * m_population_size;
      m_keep   = floor(pool);
      double k = 0.0;
      k        = m_keep * (m_keep + 1) / 2;
      m_prob.push_back(m_keep / k);
      for (int i = 2; i <= m_keep; i++) {
         m_prob.push_back((m_keep - i + 1) / k + m_prob[i - 2]);
      };
   };

private:
   int                 m_keep;            /**< Number of chromosomes that survives to selection */
   int                 m_population_size; /**< Number of chromosomes in the population */
   std::vector<double> m_prob;            /**< Vector of probabilities used in Roulette Wheel selection */
   double              m_mutation_rate;   /**< Mutation rate */
   double              m_keep_fraction;   /**< Selection rate */
   ParametersDomain    m_par_domain;      /**< search domain */
};

#endif
