#ifndef Hybrid_h
#define Hybrid_h

#include <iostream>
#include <vector>
#include <random>
#include "Config.h"
#include "Chi2Fit.h"
#include "TH2F.h"

class Hybrid {
public:

  static void Improve_Hyb_Migrad(Population& pop, TH2F* data, TF2* model) {

    Config conf_aus = pop.Configuration();

    std::random_device rd;
    std::mt19937 rng(rd());

    //definisco un dominio di ricerca centrato nel miglior risultato
    model -> SetParLimits(0, conf_aus.ParDomain.ViewParMin(0), conf_aus.ParDomain.ViewParMax(0));
    std::uniform_real_distribution<double> uni1 (conf_aus.ParDomain.ViewParMin(0), conf_aus.ParDomain.ViewParMin(0));
    model -> SetParameter(0, uni1(rng));

    for (int i=1;i<conf_aus.NumberOfParameters();i++){
      double toll=(conf_aus.ParDomain.ViewParMax(i)-conf_aus.ParDomain.ViewParMin(i))/10;
      double domain_inf =  pop.AccessChromosome(0).ViewGene(i)-toll;
      double domain_sup =  pop.AccessChromosome(0).ViewGene(i)+toll;

      if (domain_inf < conf_aus.ParDomain.ViewParMin(i))
      domain_inf = conf_aus.ParDomain.ViewParMin(i);

      if (domain_sup > conf_aus.ParDomain.ViewParMax(i))
      domain_sup = conf_aus.ParDomain.ViewParMax(i);

      model -> SetParLimits(i, domain_inf, domain_sup);
      std::uniform_real_distribution<double> uni2 (domain_inf, domain_sup);
      model -> SetParameter(i, uni2(rng));
    };

    data -> Fit ("model","Q");

    if (pop.AccessChromosome(0).ViewCost() > model -> GetChisquare()){
      for (int i=0;i<conf_aus.NumberOfParameters();i++){
        pop.AccessChromosome(0).SetGene(i,model -> GetParameter(i));
      };
      pop.AccessChromosome(0).Cost(model -> GetChisquare());
    };
  };


  static void Improve_Hyb_Random(Population& pop, TH2F* data, TF2* model){  //random search
    Config conf_aus = pop.Configuration();
    ParametersDomain domain_aus(conf_aus.NumberOfParameters());

    //definisco un dominio di ricerca centrato nel miglior risultato

    domain_aus.SetParDomain(0,"", pop.AccessChromosome(0).ViewGene(0)-5, pop.AccessChromosome(0).ViewGene(0)+5);

    for (size_t i=1; i < conf_aus.NumberOfParameters(); i++){ //per altri parametri definisco tolleranza relativa
      double toll = (conf_aus.ParDomain.ViewParMax(i)-conf_aus.ParDomain.ViewParMin(i))/10;
      double domain_inf =  pop.AccessChromosome(0).ViewGene(i)-toll;
      double domain_sup =  pop.AccessChromosome(0).ViewGene(i)+toll;
      if (domain_inf < conf_aus.ParDomain.ViewParMin(i))
      domain_inf = conf_aus.ParDomain.ViewParMin(i);
      if (domain_sup > conf_aus.ParDomain.ViewParMax(i))
      domain_sup = conf_aus.ParDomain.ViewParMax(i);
      domain_aus.SetParDomain(i,"", domain_inf, domain_sup);
    };

    Population pop_aus(conf_aus);
    pop_aus.Init();

    Chi2Fit::ComputeCost(pop_aus, data, model);
    pop_aus.Sort();

    if (pop.AccessChromosome(0).ViewCost() > pop_aus.AccessChromosome(0).ViewCost()) {
      pop.AccessChromosome(0) = pop_aus.AccessChromosome(0);
    };
  };
};

#endif
