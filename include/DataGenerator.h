#ifndef DataGenerator_h
#define DataGenerator_h

#include "TF1.h"
#include "TH1F.h"

class DataGenerator {
public:
  DataGenerator(const std::string&, double, double); 

  ~DataGenerator();

  void FillTH1F(TH1F*, size_t);

private:
  TF1* m_func;
};

#endif
