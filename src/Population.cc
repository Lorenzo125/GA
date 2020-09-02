#include "Population.h"

#include "TH1F.h"
#include <random>

Population::Population(const Config &conf)
   : m_chrom(conf.getPopulationSize(), conf.getNumberOfParameters()), m_conf(conf){};

Population::Population(const Population &p) : m_chrom(p.m_chrom), m_conf(p.m_conf){};

void Population::init()
{

   std::random_device rd;
   std::mt19937       rng(rd());

   for (size_t i = 0; i < m_chrom.size(); ++i) {
      for (size_t j = 0; j < m_conf.getNumberOfParameters(); ++j) {
         std::uniform_real_distribution<double> uni(m_conf.getParameterDomain().getParameterMin(j),
                                                    m_conf.getParameterDomain().getParameterMax(j));
         m_chrom[i].setGene(j, uni(rng));
      };
   };
};

size_t Population::size() const
{
   return m_chrom.size();
};

Config Population::getConfig()
{
   return m_conf;
};

void Population::sort()
{
   std::sort(m_chrom.begin(), m_chrom.end());
};

void Population::crossover()
{ // uniform crossover

   std::random_device                     rd;
   std::mt19937                           rng(rd());
   std::uniform_real_distribution<double> uni(0.0, 1.0);

   // generate offspring
   for (int i = 0; i < (m_chrom.size() - m_conf.getKeep()); i = i + 2) {
      m_chrom[m_chrom.size() - 1 - i].setIndicatorDown();
      m_chrom[m_chrom.size() - 2 - i].setIndicatorDown();
      int    ma = 0, pa = 0;
      double ra1 = uni(rng);
      // mother
      for (int u = 1; u < m_conf.getKeep(); u++) {
         if (ra1 > m_conf.getProb(u - 1) && ra1 <= m_conf.getProb(u)) ma = u;
      };
      double ra2 = uni(rng);
      // father
      for (int u = 1; u < m_conf.getKeep(); u++) {
         if (ra2 > m_conf.getProb(u - 1) && ra2 <= m_conf.getProb(u)) pa = u;
      };
      for (int u = 0; u < m_conf.getNumberOfParameters(); u++) {
         double beta = uni(rng);
         m_chrom[m_chrom.size() - 1 - i].setGene(u, m_chrom[ma].getGene(u) -
                                                       beta * (m_chrom[ma].getGene(u) - m_chrom[pa].getGene(u)));
         m_chrom[m_chrom.size() - 2 - i].setGene(u, m_chrom[pa].getGene(u) +
                                                       beta * (m_chrom[ma].getGene(u) - m_chrom[pa].getGene(u)));
      };
   };
};

void Population::mutation()
{
   // random generator engine
   std::random_device rd;
   std::mt19937       rng(rd());

   int mutat = floor(m_conf.getMutationRate() * m_conf.getPopulationSize() * m_conf.getNumberOfParameters());

   for (int i = 0; i < mutat; i++) {
      std::uniform_int_distribution<int> uni_1(0, m_conf.getNumberOfParameters() - 1);
      int                                ra1 = uni_1(rng);

      std::uniform_int_distribution<int> uni_2(1, m_chrom.size() - 1); //elitism
      int                                ra2 = uni_2(rng);

      std::uniform_real_distribution<double> uni_3(m_conf.getParameterDomain().getParameterMin(ra1),
                                                   m_conf.getParameterDomain().getParameterMax(ra1));
      m_chrom[ra2].setGene(ra1, uni_3(rng));
      m_chrom[ra2].setIndicatorDown();
   };
};

void Population::draw(TH1F *t_data, TF1 *t_model)
{
   TApplication *app = new TApplication("app", 0, 0);
   TCanvas *     c1  = new TCanvas("c1", "", 1500, 500);
   m_chrom[0].setModelBest(t_model);
   t_model->Draw("C");
   t_model->SetTitle("Normalized function and histogram");
   t_data->Draw("SAME");
   app->Run();
};
