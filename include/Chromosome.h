#ifndef Chromosome_h
#define Chromosome_h

#include "TF1.h"
#include <vector>
#include <iostream>
#include <iomanip>

class Chromosome {
public:
   Chromosome(size_t);

   Chromosome(const std::vector<double> &);

   Chromosome(const Chromosome &);

   void setGene(size_t, double);

   std::vector<double> &Genes() { return m_genes; };

   size_t size() const;

   void setCost(double);

   void updateModel(TF1 *);

   double getGene(size_t i);

   double getIndicator();

   void setIndicatorUp();

   void setIndicatorDown();

   double getCost();

   double &operator[](size_t i) { return m_genes[i]; };

   friend std::ostream &operator<<(std::ostream &os, Chromosome &rhs)
   {
      for (size_t i = 0; i < rhs.size(); ++i) {
         os << std::setw(10) << std::right << rhs[i];
         if (i < rhs.size() - 1) os << ", ";
      }
      os << " --> " << rhs.m_cost << "\n";
      return os;
   };

   friend bool operator<(const Chromosome &l, const Chromosome &r) { return l.m_cost < r.m_cost; }

private:
   std::vector<double> m_genes;
   double              m_cost;
   int                 m_indicator = 0; // ==0, non è ancora stato calcolato il costo, ==1 sì
};

#endif
