#ifndef Evaluate_h
#define Evaluate_h

#include "Population.h"

class Evaluate {
public:
   static void computeCostFit(Chromosome &chr, TH1F *data, TF1 *model, double Ndof)
   {
      double chi2 = 0., x = 0., y = 0., c = 0., e = 0.;
      if (chr.getIndicator() == 0) {
        TAxis *xAxis = data->GetXaxis();
        chr.updateModel(model);
        chi2 = 0.;
         for (size_t m = 1; m <= data->GetNbinsX(); m++) {
            x = xAxis->GetBinCenter(m);
            c = data->GetBinContent(m);
            e = data->GetBinError(m);
            if (e == 0) continue;
            chi2 += TMath::Power((c - model->Eval(x)) / e, 2);
         };
      chr.setCost(chi2/Ndof);
      chr.setIndicatorUp();
      };
   };

   static void computeCostFit(Population &pop, TH1F *data, TF1 *model, double Ndof)
   {
      for (size_t i = 0; i < pop.size(); ++i) {
         computeCostFit(pop[i], data, model, Ndof);
      };
   };
};

#endif
