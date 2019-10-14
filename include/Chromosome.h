#ifndef Chromosome_h
#define Chromosome_h

#include "TF1.h"
#include <vector>
#include <iostream>
#include <iomanip>

class Chromosome {
public:
  Chromosome(size_t);

  Chromosome(const std::vector<double>&);

  Chromosome(const Chromosome&);

  void Init(double, double);

  void SetGene(size_t, double);

  std::vector<double>& Genes() {
    return m_genes;
  };

  size_t Size() const;

  void Cost(double);

  /// Set the parameters of the model according
  /// to the genes of the chromosome
  void UpdateModel(TF1*);

  double& operator[](size_t i) {
    return m_genes[i];
  };

  friend std::ostream& operator<<(std::ostream& os, Chromosome& rhs) {
    for (size_t i = 0; i < rhs.Size(); ++i) {
      os << std::setw(10) << std::right << rhs[i];
      if (i < rhs.Size()-1) os << ", ";
    }
    os << " --> " << rhs.m_cost << "\n";
    return os;
  };

  friend bool operator<(const Chromosome& l, const Chromosome& r) {
    return l.m_cost < r.m_cost;
  }

private:
  std::vector<double> m_genes;
  double m_cost;
};

#endif
