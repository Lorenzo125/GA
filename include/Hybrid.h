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

    model -> SetParLimits(0, conf_aus.ParDomain.ViewParMin(0), conf_aus.ParDomain.ViewParMax(0));

    //definisco un dominio di ricerca centrato nel miglior risultato
    for (int i=1;i<conf_aus.NumberOfParameters();i++){
      model -> SetParameter(i, pop.AccessChromosome(0).ViewGene(i));
      double toll=0.2;
      double domain_inf =  pop.AccessChromosome(0).ViewGene(i)-toll;
      double domain_sup =  pop.AccessChromosome(0).ViewGene(i)+toll;
      if (domain_inf < conf_aus.ParDomain.ViewParMin(i))
      domain_inf = conf_aus.ParDomain.ViewParMin(i);
      if (domain_sup > conf_aus.ParDomain.ViewParMax(i))
      domain_sup = conf_aus.ParDomain.ViewParMax(i);
      model -> SetParLimits(i, domain_inf, domain_sup);
    };

    data -> Fit ("model","Q");
    for (int i=0;i<conf_aus.NumberOfParameters();i++){
      pop.AccessChromosome(0).SetGene(i,model -> GetParameter(i));
    };
    pop.AccessChromosome(0).DownIndicator();
  };


  static void Improve_Hyb_Random(Population& pop, TH2F* data, TF2* model){  //random search
    Config conf_aus = pop.Configuration();
    ParametersDomain domain_aus(conf_aus.NumberOfParameters());

    //definisco un dominio di ricerca centrato nel miglior risultato
    for (size_t i=0; i < conf_aus.NumberOfParameters(); i++){
      double toll=0.1;
      //toll = (conf_aus.ParDomain.ViewParMax(i)-conf_aus.ParDomain.ViewParMin(i))/10; //arbitrario
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
