#ifndef Evaluate_h
#define Evaluate_h

#include "Population.h"

class Evaluate {
public:
   /// Compute the chi2/Ndof using genes of chromosome as fitting parameters  <br> <br>
   /// The chi2/Ndof is assigned to chromosomes as cost to be minimized
   static void computeCostFit(Chromosome &t_chr, TH1F *t_data, TF1 *t_model, double t_Ndof)
   {
      double chi2 = 0., x = 0., y = 0., c = 0., e = 0.;
      if (t_chr.getIndicator() == 0) {
         TAxis *xAxis = t_data->GetXaxis();
         t_chr.setModel(t_model);
         chi2 = 0.;
         for (size_t m = 1; m <= t_data->GetNbinsX(); m++) {
            x = xAxis->GetBinCenter(m);
            c = t_data->GetBinContent(m);
            e = t_data->GetBinError(m);
            if (e == 0) continue;
            chi2 += TMath::Power((c - t_model->Eval(x)) / e, 2);
         };
         t_chr.setCost(chi2 / t_Ndof);
         t_chr.setIndicatorUp();
      };
   };

   /// Compute the chi2/Ndof  for each chromosome in the population <br> <br>
   static void computeCostFit(Population &t_pop, TH1F *t_data, TF1 *t_model, double t_Ndof)
   {
      for (size_t i = 0; i < pop.size(); ++i) {
         computeCostFit(pop[i], t_data, t_model, t_Ndof);
      };
   };
};

#endif
