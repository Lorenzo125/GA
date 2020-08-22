#include <iostream>
#include "Config.h"
#include "Evaluate.h"
#include "TFile.h"
#include "TMath.h"
#include "TApplication.h"
#include <algorithm>

int main(int argc, char *argv[])
{

   Config conf;
   conf.setPopulationSize(200);
   conf.setMutationRate(0.10);
   conf.setKeepFraction(0.40);
   conf.setConfigCrossover();

   // max number of iterations allowed
   size_t iter_max = 10000;
   // max number of iterations with the same best cost (with a tollerance toll)
   size_t stuck_max = 100;
   double toll      = 1e-3;

   TApplication *app = new TApplication("app", 0, 0);

   // import the dataset from a .root  file
   TFile *f   = new TFile("XXX.root", "READ");
   TH1F * aux = (TH1F *)f->Get("data");
   //

   // clone the initial data in order to normalize the histogram and use
   // it without working on the original dataset
   TH1F *data = (TH1F *)aux->Clone();
   if (data->GetSumw2N() == 0) data->Sumw2(kTRUE);
   data->Scale(1 / (data->GetEntries() * data->GetBinWidth(1)));
   //

   // define the model as a parametric function
   TF1 *model = new TF1("model", "XXX", XX, XX);
   //

   ParametersDomain domain(model->GetNumberFreeParameters());

   // for every parameter in the model
  domain.setParameterDomain(a, "",XX, XX));
  domain.setParameterDomain(b, "",XX, XX));
  //

  conf.setParameterDomain(domain);

  double Ndof = (data->GetNbinsX() - model->GetNumberFreeParameters());

  std::vector<double> chi_run;

  Population pop(conf);
  pop.init();

  Evaluate::computeCostFit(pop, data, model, Ndof);
  pop.sort();
  chi_run.push_back(pop[0].getCost());

  size_t iter  = 0;
  size_t stuck = 0;

  do {
     iter++;
     pop.crossover();
     Evaluate::computeCostFit(pop, data, model, Ndof);
     pop.sort();
     pop.mutation();
     Evaluate::computeCostFit(pop, data, model, Ndof);
     pop.sort();
     chi_run.push_back(pop[0].getCost());
     if (abs(chi_run[iter - 1] - chi_run[iter]) < toll) {
        stuck++;
     } else
        stuck = 0;
  } while (stuck < stuck_max & iter < iter_max);

  std::cout << pop << std::endl;

  pop.draw();

  delete data;
  delete model;
  return EXIT_SUCCESS;
}
