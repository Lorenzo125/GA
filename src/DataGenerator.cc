#include "DataGenerator.h"

DataGenerator::DataGenerator(const std::string& formula, double vmin, double vmax) {
  m_func = new TF1("f", formula.c_str(), vmin, vmax);
};

DataGenerator::~DataGenerator() { 
  delete m_func;
};

void DataGenerator::FillTH1F(TH1F* h, size_t nevents) {
  for (size_t i = 0; i < nevents; ++i) {
    h->Fill(m_func->GetRandom());
  };
};
