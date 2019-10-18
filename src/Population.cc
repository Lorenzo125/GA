#include "Population.h"
#include "Mutation.h"
#include "Hybrid.h"

#include "TH1F.h"
#include <random>

Population::Population(const Config& conf) :
m_chrom(conf.PopulationSize, conf.NumberOfParameters()), m_conf(conf) {};

Population::Population(const Population& p) :
m_chrom(p.m_chrom), m_conf(p.m_conf) {};

void Population::Init() {
  for (size_t i = 0; i < m_chrom.size(); ++i) {
    size_t csize = m_chrom[i].Size();
    for (size_t j = 0; j < csize; ++j) {
      m_chrom[i][j] = gRandom->Uniform(m_conf.ParDomain[j].min_val, m_conf.ParDomain[j].max_val);
    }
  }
};

size_t Population::Size() const {
  return m_chrom.size();
};

void Population::Sort() {
  std::sort(m_chrom.begin(), m_chrom.end());
};

void Population::PairAndMate() {
  // number of chromosomes to keep
  int keep = floor(m_conf.KeepFraction*m_conf.PopulationSize);  ///arrotonda per difetto
  Mating::CrossOver(m_chrom, keep, m_conf.NumberOfParameters());
};

void Population::PairAndMate_Beta() {
  // number of chromosomes to keep
  int keep = floor(m_conf.KeepFraction*m_conf.PopulationSize);  ///arrotonda per difetto
  Mating::CrossOver_Beta(m_chrom, keep, m_conf.NumberOfParameters());
};

void Population::Evolve() {
  // number of genes to change
  int mutat = floor(m_conf.MutationRate*m_conf.PopulationSize*m_conf.ParDomain.NumberOfParameters());
  Mutation::Elite(m_chrom, mutat, m_conf);
};

void Population::Improve_Hyb(TH1F* data, TF1* model) {
  Hybrid::Random_Hybrid(m_chrom, m_conf, data, model);
};
