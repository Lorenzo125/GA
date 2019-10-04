// costruisce l'istogramma che traccia il valore medio del chi

void fill_with_mean (int num_iteraz, int num_chrom){
  double sum = 0.;
  for (int i=0; i<num_chrom; i++) {
    sum+=evaluate(&pop[i]);
  }
  hmean.Fill(num_iteraz,sum/num_chrom);
}
