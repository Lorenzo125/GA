#ifndef DataGenerator_h
#define DataGenerator_h

#include "TF1.h"
#include "TF2.h"
#include "TH2F.h"

class DataGenerator {
public:
  DataGenerator(const std::string&, double, double, double,double);

  ~DataGenerator();

  void FillTH2F(TH2F*, size_t);

private:
  TF2* m_func;
};

#endif
