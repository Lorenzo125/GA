#ifndef Mutation_h
#define Mutation_h

#include <vector>
#include <random>

class Mutation {
public:
  static void Elite(std::vector<Chromosome>& m_chrom, int mutat, Config& m_conf) {
    // random generator engine
    std::random_device rd;
    std::mt19937 rng(rd());
    int ngenes = (m_conf.NumberOfParameters())-1;
    int nchrom = (m_chrom.size())-1;

    for (int i=0; i<mutat; i++){
    std::uniform_int_distribution<int> uni_1(0, ngenes);
    int ra1 = uni_1(rng);

    std::uniform_int_distribution<int> uni_2(1, nchrom); //1 per introdurre elitismo
    int ra2 = uni_2(rng);

    m_chrom [ra2][ra1] = gRandom->Uniform(m_conf.ParDomain[ra1].min_val, m_conf.ParDomain[ra1].max_val);
    };
  };
};

#endif
