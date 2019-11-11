#include "Chromosome.h"
#include  "TRandom.h"

Chromosome::Chromosome(size_t ngenes) :
m_genes(ngenes, 0), m_cost(0.) {};

Chromosome::Chromosome(const std::vector<double>& v) :
m_genes(v), m_cost(0.) {};

Chromosome::Chromosome(const Chromosome& c) :
m_genes(c.m_genes), m_cost(c.m_cost) {};

void Chromosome::SetGene(size_t i, double v) {
  m_genes[i] = v;
};

size_t Chromosome::Size() const {
  return m_genes.size();
};

double Chromosome::ViewGene(size_t i) {
  return m_genes[i];
};

double Chromosome::ViewIndicator() {
  return m_indicator;
};

void Chromosome::UpIndicator() {
  m_indicator = 1;
};

void Chromosome::DownIndicator() {
  m_indicator = 0;
};

double Chromosome::ViewCost(){
  return m_cost;
};

void Chromosome::Cost(double v) {
  m_cost = v;
};

void Chromosome::UpdateModel(TF2* f) {
  for (size_t i = 0; i < m_genes.size(); ++i) {
    f->SetParameter(i, m_genes[i]);
  }
};
