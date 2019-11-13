#include "DataGenerator.h"
#include "TRandom.h"

DataGenerator::DataGenerator(const std::string& formula, double xmin, double xmax, double ymin, double ymax) {
  m_func = new TF2("f", formula.c_str(), xmin, xmax, ymin, ymax);
};

DataGenerator::~DataGenerator() {
  delete m_func;
};

void DataGenerator::FillTH2F(TH2F* h, size_t nevents) {
  TF1* f1 = new TF1("f1", "exp(2*x)*3^x", 0., 1.);
  TF1* f2 = new TF1("f2", "exp(4*x)", 0., 1.);
  for (size_t i = 0; i < nevents; ++i) {
    h -> Fill(f1->GetRandom(), f2->GetRandom());
  };
  /*  for (size_t i=1;i<=h->GetNbinsX(); i++) {   // aggiungo il rumore
  h->SetBinContent(i, h->GetBinContent(i)+gRandom->Gaus(0, 0.1*(h->GetBinContent(i))));
};*/
};
