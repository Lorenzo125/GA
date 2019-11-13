#include <iostream>
#include "Config.h"
#include "Chi2Fit.h"
#include "Hybrid.h"
#include <algorithm>

int main(int argc, char *argv[]) {
  std::cout << ">>> FitGA++" << std::endl;

  Config conf;
  conf.PopulationSize = 16;
  conf.MutationRate = 0.25;
  conf.KeepFraction = 0.50;


  TF2* f = new TF2("f","exp(2*x)*(3^x)*exp(4*y)", 0., 1., 0., 1.);
  TH2F* data = new TH2F("data", "", 100, 0., 1., 100, 0., 1.);
  data -> FillRandom("f", 10000000);

  TF2* model = new TF2("model","[0]*exp([1]*x)*([2]^x)*exp([3]*y)", 0., 1., 0., 1.);
  ParametersDomain domain(model->GetNumberFreeParameters());
  domain.SetParDomain(0, "p0", 1, 100000);
  domain.SetParDomain(1, "p1", 1.0, 3.0);
  domain.SetParDomain(2, "p2", 2.0, 4.0);
  domain.SetParDomain(3, "p3", 3.0, 5.0);
  conf.ParDomain = domain;

  Population pop(conf);
  pop.Init();

  Chi2Fit::ComputeCost(pop, data, model);
  pop.Sort();

  std::cout << "--- Original population ---" << std::endl;
  std::cout << pop << std::endl;

  for (auto i = 0; i < 15; ++i) {
    pop.PairAndMate_Beta();
    Chi2Fit::ComputeCost(pop, data, model);
    pop.Sort();
    pop.Evolve();
    Chi2Fit::ComputeCost(pop, data, model);
    pop.Sort();
    Hybrid::Improve_Hyb_Migrad(pop, data, model);
    std::cout << "--- Iteration " << i+1 << std::endl;
    std::cout << pop << std::endl;
  };

  delete data;
  delete model;
  return EXIT_SUCCESS;
}
