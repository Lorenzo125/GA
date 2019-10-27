#ifndef Hybrid_h
#define Hybrid_h

#include <iostream>
#include <vector>
#include <random>
#include "Config.h"
#include "Chi2Fit.h"

#include "TH1F.h"

class Hybrid {
public:

  static void Improve_Hyb(Population& pop, TH1F* data, TF1* model) {

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

    /*for (int i=0; i<10;i++){ //meta-genetic algortithm
      pop_aus.PairAndMate_Beta();
      Chi2Fit::ComputeCost(pop_aus, data, model);
      pop_aus.Sort();
      pop_aus.Evolve();
      pop_aus.Sort();
    };*/

    Chi2Fit::ComputeCost(pop_aus, data, model); //random search
    pop_aus.Sort();

    if (pop.AccessChromosome(0).ViewCost() > pop_aus.AccessChromosome(0).ViewCost()) {
      pop.AccessChromosome(0) = pop_aus.AccessChromosome(0);
    };
  };

};

#endif
