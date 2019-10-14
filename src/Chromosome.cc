#include "Chromosome.h"
#include  "TRandom.h"

Chromosome::Chromosome(size_t ngenes) :
m_genes(ngenes, 0), m_cost(0.) {};

Chromosome::Chromosome(const std::vector<double>& v) :
m_genes(v), m_cost(0.) {};

Chromosome::Chromosome(const Chromosome& c) :
m_genes(c.m_genes), m_cost(c.m_cost) {};

void Chromosome::Init(double vmin, double vmax) {
  for (size_t i = 0; i < m_genes.size(); ++i)
  m_genes[i] = gRandom->Uniform(vmin, vmax);
};

void Chromosome::SetGene(size_t i, double v) {
  m_genes[i] = v;
};

size_t Chromosome::Size() const {
  return m_genes.size();
};

void Chromosome::Cost(double v) {
  m_cost = v;
};

void Chromosome::UpdateModel(TF1* f) { 
  for (size_t i = 0; i < m_genes.size(); ++i) {
    f->SetParameter(i, m_genes[i]);
  }
};