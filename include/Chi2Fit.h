#ifndef Chi2Fit_h
#define Chi2Fit_h

#include "Population.h"

class Chi2Fit {
public:
  static void ComputeCost(Chromosome& chr, TH2F* data, TF2* model) {
    double chi2 = 0., x = 0., y = 0., c = 0., e = 0.;
    TAxis* xAxis = data->GetXaxis();
    TAxis* yAxis = data->GetYaxis();
    chr.UpdateModel(model);
    chi2 = 0.;
    if (chr.ViewIndicator()==0){
      for (size_t m = 1; m <= data->GetNbinsX(); m++) {
        x = xAxis->GetBinCenter(m);
        for (size_t n = 1; n <= data->GetNbinsY(); n++){
          y = yAxis->GetBinCenter(n);
          c = data->GetBinContent(m,n);
          e = data->GetBinError(m,n);
          if (e == 0) continue;
          chi2 += TMath::Power((c-model->Eval(x,y))/e, 2);
        };
      };
      chr.Cost(chi2);
      chr.UpIndicator();
    };
  };

  static void ComputeCost(Population& pop, TH2F* data, TF2* model) {
    for (size_t i = 0; i < pop.Size(); ++i) {
      ComputeCost(pop[i], data, model);
    }
  };

};

#endif
