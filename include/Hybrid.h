#ifndef Hybrid_h
#define Hybrid_h

#include <iostream>
#include <vector>
#include <random>
#include "Population.h"
#include "Config.h"
#include "Chi2Fit.h"

#include "TH1F.h"

class Hybrid {
public:
  static void Random_Hybrid(std::vector<Chromosome>& m_chrom, Config& m_conf, TH1F* data, TF1* model) {

    Config conf_aus = m_conf;
    ParametersDomain domain_aus;

    //definisco un dominio di ricerca centrato nel miglior risultato
    for (int i=0;i < m_conf.NumberOfParameters(); i++){
      double toll = (m_conf.ParDomain.ViewParMax(i)-m_conf.ParDomain.ViewParMin(i))/10; //10 arbitrario ////aggiungere tutte le funzioni
      double domain_inf =  m_chrom[0].ViewGene(i)-toll;
      double domain_sup =  m_chrom[0].ViewGene(i)+toll;
      if (domain_inf < m_conf.ParDomain.ViewParMin(i))
      domain_inf = m_conf.ParDomain.ViewParMin(i);
      if (domain_sup > m_conf.ParDomain.ViewParMax(i))
      domain_sup = m_conf.ParDomain.ViewParMax(i);
      domain_aus.SetParDomain(i,"", domain_inf, domain_sup);
    };

    Population pop_aus(conf_aus);
    pop_aus.Init();

    for (int i=0; i<10;i++){ //meta-genetic algortithm
      pop_aus.PairAndMate_Beta();
      Chi2Fit::ComputeCost(pop_aus, data, model);
      pop_aus.Sort();
    };

  if (m_chrom[0].ViewCost() > pop_aus.m_chrom[0].ViewCost())
  for (int i=0;i<m_chrom.size();i++)
  m_chrom[0].SetGene = pop_aus.m_chrom[0].ViewGene(i);

  };

};

#endif
