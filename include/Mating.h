#ifndef Mating_h
#define Mating_h

#include <vector>
#include <random>

class Mating {
public:
  static void CrossOver(std::vector<Chromosome>& m_chrom, int keep, int ngenes) {
    // random generator engine
    std::random_device rd;
    std::mt19937 rng(rd());

    // shuffle the indeces associated to the chromosomes to keep
    std::vector<int> keep_index(keep); //create a vector of int of keep-dimension
    std::iota(keep_index.begin(), keep_index.end(), 0); //fill with crescent values from 0
    std::shuffle(keep_index.begin(), keep_index.end(), rng); //random shuffle of values

    std::size_t const half_size = keep_index.size() / 2;
    std::vector<int> mother(keep_index.begin(), keep_index.begin() + half_size);
    std::vector<int> father(keep_index.begin() + half_size, keep_index.end());

    // mating (crossover)
    std::uniform_int_distribution<int> uni(0, ngenes);
    int cover = uni(rng);

    std::cout << keep << " - " << cover << std::endl;
    std::vector<double> ma, pa, ch(ngenes);
    for (auto i = 0; i < mother.size(); ++i) {
      cover = uni(rng);
      ma = m_chrom[mother[i]].Genes();
      pa = m_chrom[father[i]].Genes();

      // ma pa
      std::copy(ma.begin(), ma.begin()+cover, ch.begin());
      std::copy(pa.begin()+cover, pa.end(), ch.begin()+cover);
      m_chrom[keep] = Chromosome(ch);
      keep++;

      // pa maprop []
      std::copy(pa.begin(), pa.begin()+cover, ch.begin());
      std::copy(ma.begin()+cover, ma.end(), ch.begin()+cover);
      m_chrom[keep] = Chromosome(ch);
      keep++;
    }
  };

  static void CrossOver_Beta(std::vector<Chromosome>& m_chrom, int keep, int ngenes) {
      int nchrom = (m_chrom.size());
      float k=0.0;
      //calcolo le probabilità sequenziali (da migliorare)
      std::vector<double> prob(keep);
      for (int i=1;i<=keep;i++)
      k+=i;
      for (int i=1; i <= keep; i++)
      prob[i-1]=(keep-i+1)/k;
      for (int i=1; i < keep; i++)
      prob[i]+=prob[i-1];

      // genero i figli (uno solo ad accoppiamento, cosi riesco sempre a rispettare keep, anche se mi viene dispari)
      for (int i=0;i<(nchrom-keep);i++){
        int ma=0, pa=0;
        double ra1 = gRandom->Uniform(0, 1);
        for (int u=1;u<keep;u++){ //scelgo la madre
          if (ra1>prob [u-1] && ra1<=prob [u])
          ma=u;
        };
        for (int u=1;u<keep;u++){ //scelgo il padre
          if (ra1>prob [u-1] && ra1<=prob [u])
          pa=u;
        };

        for (int u=0;u<ngenes;u++){ //uniform crossover
          double beta = gRandom->Uniform(0, 1);
          m_chrom [nchrom-1-i][u]= m_chrom [ma][u] - beta*(m_chrom [ma][u] - m_chrom [pa][u]);
        };
    };
  };
};

#endif
