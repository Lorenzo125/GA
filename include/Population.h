#ifndef Population_h
#define Population_h

#include "Chromosome.h"
#include "Config.h"
#include "TRandom.h"
#include "Mating.h"
#include "Mutation.h"
#include <iostream>
#include <numeric>

#include "TH1F.h"

class Population {
public:
  Population(const Config&);

  Population(const Population&);

  void Init();

  size_t Size() const;

  void Sort();

  void PairAndMate();

  void PairAndMate_Beta();

  void Evolve();

  Config Configuration();

  Chromosome& AccessChromosome(size_t i) {
    return m_chrom[i];
  };

  Chromosome& operator[](size_t i) {
    return m_chrom[i];
  };

  friend std::ostream& operator<<(std::ostream& os, Population& rhs) {
    os << "--- Population ---\n";
    /*for (size_t i = 0; i < rhs.Size(); ++i) { //stampa tutta la popolazione
      os << rhs[i];
    };*/
    os << rhs[0]; //stampa migliore

    os << "------------------\n";
    return os;
  };

private:
  std::vector<Chromosome> m_chrom;
  Config m_conf;
};

#endif
