#include <iostream>
#include "Config.h"
#include "DataGenerator.h"
#include "Chi2Fit.h"
#include "Hybrid.h"
#include <algorithm>

int main(int argc, char *argv[]) {
  std::cout << ">>> FitGA++" << std::endl;

  Config conf;
  conf.PopulationSize = 64;
  conf.MutationRate = 0.25;
  conf.KeepFraction = 0.6;


  DataGenerator dg("exp(-1.5*x)+0.2*(x-0.6)^2", 0, 1);
  TH1F* data = new TH1F("data", "", 800, 0., 1.);
  dg.FillTH1F(data, 10000);

  TF1* model = new TF1("model", "[0]*(exp([1]*x)+[2]*(x-[3])^2)", 0, 1);
  ParametersDomain domain(model->GetNumberFreeParameters());
  domain.SetParDomain(0, "p0", 1, 500);
  domain.SetParDomain(1, "p1", -2.25, -0.75);
  domain.SetParDomain(2, "p2", 0.1, 0.3);
  domain.SetParDomain(3, "p3", 0.3, 0.9);
  conf.ParDomain = domain;

  Population pop(conf);
  pop.Init();

  Chi2Fit::ComputeCost(pop, data, model);
  pop.Sort();

  std::cout << "--- Original population ---" << std::endl;
  std::cout << pop << std::endl;

  for (auto i = 0; i < 200; ++i) {
    pop.PairAndMate_Beta();
    //Chi2Fit::ComputeCost(pop, data, model); //non serve senza elitismo
    //pop.Sort();
    pop.Evolve();
    Chi2Fit::ComputeCost(pop, data, model);
    pop.Sort();
    Hybrid::Improve_Hyb(pop, data, model);
    std::cout << "--- Iteration " << i+1 << std::endl;
    std::cout << pop << std::endl;
  };

  delete data;
  delete model;

  return EXIT_SUCCESS;
}
