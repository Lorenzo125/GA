#ifndef Chi2Fit_h
#define Chi2Fit_h

#include "Population.h"
#include "TH1F.h"

class Chi2Fit {
public:
  static void ComputeCost(Chromosome& chr, TH1F* data, TF1* model) {
    double chi2 = 0., x = 0., y = 0., e = 0.;
    chr.UpdateModel(model);
    chi2 = 0.;
    for (size_t b = 1; b <= data->GetNbinsX(); ++b) {
      x = data->GetBinCenter(b);
      y = data->GetBinContent(b);
      e = data->GetBinError(b);
      if (e == 0) continue;
      chi2 += TMath::Power((y - model->Eval(x))/e, 2);
    }
    chr.Cost(chi2);
  };

  static void ComputeCost(Population& pop, TH1F* data, TF1* model) {
    double chi2 = 0., x = 0., y = 0., e = 0.;
    for (size_t i = 0; i < pop.Size(); ++i) {
      ComputeCost(pop[i], data, model);
    }
  };

};

#endif
