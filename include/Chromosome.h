#ifndef Chromosome_h
#define Chromosome_h

#include "TF1.h"
#include <vector>
#include <iostream>
#include <iomanip>

class Chromosome {
public:
   /// Chromosome default constructor for an _n_ - dimensional search space (- _n_ fitting parameters)
   Chromosome(size_t n);

   /// Chromosome constructor with an assigned vector of genes
   Chromosome(const std::vector<double> &n_gene);

   /// Set the chromosome's _t_ - gene <br> <br>
   /// Referenced by Population::init(), Population::crossover(), Population::muation()
   void setGene(size_t, double);

   /// Return m_genes
   std::vector<double> &Genes() { return m_genes; };

   /// Return the number of genes <br> <br>
   /// Referenced by Population::init(), Population::crossover(), Population::mutation()
   size_t size() const;

   /// Set the cost of the chromosomes <br> <br>
   /// Referenced by Evaluate::computeCostFit()
   void setCost(double);

   /// Set the parameters of the model equal to genes of the chromosome <br> <br>
   /// Referenced by Evaluate::computeCostFit()
   void setModel(TF1 *);

   /// Return the chromosome's _t_ - gene
   double getGene(size_t i);

   /// Return the value of m_indicator <br> <br>

   /// Referenced by Evaluate::computeCostFit()
   double getIndicator();

   /// Set the value of the indicator equal to 1; <br> <br>
   /// Referenced by Evaluate::computeCostFit()
   void setIndicatorUp();

   /// Set the value of the indicator equal to 0
   /// Referenced by Evaluate::computeCostFit()
   void setIndicatorDown();

   /// Return the cost of the chromosome
   double getCost();

   /// Return the chromosome's _t_ - gene
   double &operator[](size_t t) { return m_genes[t]; };

   /// Show on video the characteristics of the population
   friend std::ostream &operator<<(std::ostream &os, Chromosome &rhs)
   {
      for (size_t i = 0; i < rhs.size(); ++i) {
         os << std::setw(10) << std::right << rhs[i];
         if (i < rhs.size() - 1) os << ", ";
      }
      os << " --> " << rhs.getCost() << "\n";
      return os;
   };

   /// Compare two chromosomes' cost
   friend bool operator<(const Chromosome &l, const Chromosome &r) { return l.m_cost < r.m_cost; }

private:
   std::vector<double> m_genes;         /**<Vector storing the genes of the chromosome*/
   double              m_cost;          /**<Cost of the chromosome*/
   int                 m_indicator = 0; /**<When the indicator equals 0, the cost needs to be revaluated */
};

#endif
