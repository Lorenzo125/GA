#include <iostream>
#include "Config.h"
#include "Evaluate.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TMath.h"
#include "TStopwatch.h"
#include "TCanvas.h"
#include "TApplication.h"
#include <algorithm>

int main(int argc, char *argv[])
{

   std::cout << ">>> FitGA++" << std::endl;

   Config conf;
   conf.setPopulationSize(20);
   conf.setMutationRate(0.10);
   conf.setKeepFraction(0.50);
   conf.setConfigCrossover();

   size_t iter_max = 20; // max iterazioni in cui si può avere sempre lo stesso risultato
   double toll     = 1e-10;
   int    num_run  = 100;

   TGraphAsymmErrors *chi_trend  = new TGraphAsymmErrors();
   TGraphAsymmErrors *time_trend = new TGraphAsymmErrors();
   TFile *            file       = TFile::Open("GA_20_0.2_0.6.root", "RECREATE");

   std::vector<double> x = {10, 20, 40, 60, 80, 100}; // dominio centrato su parametro reale con larghezza +-x %

   TApplication *app = new TApplication("app", 0, 0);

   TStopwatch *timer = new TStopwatch();

   TFile *f   = new TFile("data_J_10000.root", "READ");
   TH1F * aux = (TH1F *)f->Get("data");

   // clone the initial data in order to normalize and use it without work on the original source

   TH1F *data = (TH1F *)aux->Clone();
   data->Scale(1 / (data->GetEntries() * data->GetBinWidth(1)));

   TF1 *model = new TF1("model",
                        "[12]*([0]*([1]*ROOT::Math::crystalball_function(x, [6], [5], [7]*[3], [8]+[4]) + "
                        "(1.-[1])*ROOT::Math::crystalball_function(x, [6], [5], [7], [8])) + (1.-[0])*([2]*exp([11]*x) "
                        "+ (1.-[2])*(1.-TMath::Erf((x-[10])/[9]))*0.5))",
                        1.50, 8.70);

   ParametersDomain domain(model->GetNumberFreeParameters());

   // domain of parameters

   std::vector<double> dom(13);

   dom[0]  = 4.3588e-01;
   dom[1]  = 1.3231e-02;
   dom[2]  = 3.5251e-01;
   dom[3]  = 1.0;
   dom[4]  = 5.8000e-01;
   dom[5]  = 1.99999e+01;
   dom[6]  = 1.0464;
   dom[7]  = 6.9119e-02;
   dom[8]  = 3.0980;
   dom[9]  = 6.8136e-01;
   dom[10] = 2.5249;
   dom[11] = -4.9824e-01;
   dom[12] = 1.5598;

   for (size_t coeff = 0; coeff < x.size(); coeff++) {

      std::vector<double> chi_GA; // memorize the best chi for each analysis on a coefficient
      std::vector<double> time_GA;

      for (size_t j = 0; j < dom.size(); j++) {
         domain.setParameterDomain(j, "", dom[j] - abs(dom[j] * x[coeff] / 100), dom[j] + abs(dom[j] * x[coeff] / 100));
      };

      conf.setParameterDomain(domain);

      double Ndof = (data->GetNbinsX() - model->GetNumberFreeParameters());

      for (size_t k = 0; k < num_run; k++) {

         std::vector<double> chi_run; // memorize the chi progression to undesrstand when to stop

         Population pop(conf);
         pop.init();

         Evaluate::computeCostFit(pop, data, model);
         pop.sort();
         chi_run.push_back(pop[0].getCost() / Ndof);

         size_t iter  = 0;
         size_t stack = 0;

         timer->Reset();
         timer->Start();

         do {
            iter++;
            pop.crossover();
            Evaluate::computeCostFit(pop, data, model);
            pop.sort();
            pop.mutation();
            Evaluate::computeCostFit(pop, data, model);
            pop.sort();
            chi_run.push_back(pop[0].getCost() / Ndof);
            if (chi_run[iter] / Ndof - chi_run[iter - 1] / Ndof <
                toll) // if the result doesn't change significally for 10 times, the algorithm stops
               stack++;
            else
               stack = 0;
         } while (stack < 20 & iter < iter_max);

         time_GA.push_back(timer->RealTime());
         chi_GA.push_back(pop[0].getCost() / Ndof);
      };

      double sum_chi   = 0.;
      double sum2_chi  = 0.;
      double sum_time  = 0.;
      double sum2_time = 0.;

      for (int l = 0; l < num_run; l++) {
         sum_chi += chi_GA[l];
         sum2_chi += TMath::Power(chi_GA[l], 2);
         sum_time += time_GA[l];
         sum2_time += TMath::Power(time_GA[l], 2);
      };

      double average_chi  = sum_chi / num_run;
      double average_time = sum_time / num_run;

      double variance_chi  = (sum2_chi - num_run * TMath::Power(average_chi, 2)) / (num_run - 1);
      double variance_time = (sum2_time - num_run * (TMath::Power(average_time, 2))) / (num_run - 1);

      chi_trend->SetPoint(coeff, x[coeff], average_chi);
      chi_trend->SetPointError(coeff, 0, 0, TMath::Power(variance_chi / num_run, 0.5),
                               TMath::Power(variance_chi / num_run, 0.5));

      time_trend->SetPoint(coeff, x[coeff], average_time);
      time_trend->SetPointError(coeff, 0, 0, TMath::Power(variance_time / num_run, 0.5),
                                TMath::Power(variance_time / num_run, 0.5));
   };

   file->cd();
   chi_trend->Write("chi_trend");
   time_trend->Write("time_trend");
   file->Close();

   TCanvas *c1 = new TCanvas("c1", "", 1000, 500);
   c1->Divide(2, 1);
   c1->cd(1);
   chi_trend->GetXaxis()->SetTitle("% domain");
   chi_trend->GetYaxis()->SetTitle("#chi^{2}/NDOF");
   chi_trend->SetFillColor(4);
   chi_trend->SetFillStyle(3005);
   chi_trend->SetTitle("#chi^{2}/NDOF trend versus % domain");
   chi_trend->Draw("a3");

   c1->cd(2);

   time_trend->GetXaxis()->SetTitle("% domain");
   time_trend->GetYaxis()->SetTitle("time");
   time_trend->SetFillColor(2);
   time_trend->SetFillStyle(3005);
   time_trend->SetTitle("time trend versus % domain");
   time_trend->Draw("a3");

   app->Run();

   delete data;
   delete model;
   return EXIT_SUCCESS;
}
