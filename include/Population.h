#ifndef Population_h
#define Population_h

#include "Chromosome.h"
#include "Config.h"
#include "TRandom.h"
#include <iostream>
#include <numeric>

#include "TH1F.h"

class Population {
public:
   /// Population default constructor, starting assigning the configuration Config
   Population(const Config &);

   /// Population initialization with random starting genes
   void init();

   /// Return the number of chromosomes in the population <br> <br>
   /// Referenced by Evaluate::computeCostFit()
   size_t size() const;

   /// Sort chromosomes according to their cost
   void sort();

   /// Start the mutation process
   void mutation();

   /// Start the crossover process
   void crossover();

   /// Return the configuration Config
   Config getConfig();

   /// Show on display the model with the best set of parameters found and the normalized data in the background
   void draw(TH1F *t_data, TF1 *t_model);

   /// Return the _t_ - chromosome in the population
   Chromosome &getChromosome(size_t t) { return m_chrom[t]; };

   /// Return the _t_ - chromosome in the population
   Chromosome &operator[](size_t t) { return m_chrom[t]; };

   /// Show on video the sequence of genes for the best chromosome
   friend std::ostream &operator<<(std::ostream &os, Population &rhs)
   {
      os << "--- Population ---\n";
      /*for (size_t i = 0; i < rhs.Size(); ++i) {
        os << rhs[i];
      };*/
      os << rhs[0];

      os << "------------------\n";
      return os;
   };

private:
   std::vector<Chromosome> m_chrom; /**<Vector of chromosomes that compose the population*/
   Config                  m_conf;  /**<Configuration used in the population*/
};

#endif
