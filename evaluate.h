//valuta il cromosoma

double evaluate (struct chrom *chr){
  double chi = 0.;
  for (int m=0; m<ngene; m++) {
    //assegno al modello di fit i valori dei parametri (geni) del cromosoma
    f.SetParameter (m, chr->gene[m]);
  }
  for (int i=1; i<=BIN; i++) {
    if (h_noise.GetBinError(i) == 0) continue;
    chi+=TMath::Power((f.Eval(h_noise.GetBinCenter(i),0,0,0)-h_noise.GetBinContent(i))/h_noise.GetBinError(i),2)/(BIN-ngene); //valuto la funzione nella coordinata x=h_noise.GetBinCenter(i) con i parametri del cromosoma
  }
  return chi;
}
