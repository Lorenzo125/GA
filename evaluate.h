double evaluate (struct chrom *chr){
  int i,m;
  double chi=0;
  for (m=0;m<ngene;m++)
  f.SetParameter (m, chr->gene[m]);
  for (i=1;i<=BIN;i++)
  chi+=TMath::Power((f.Eval(h_noise.GetBinCenter(i),0,0,0)-h_noise.GetBinContent(i))/h_noise.GetBinError(i),2)/(BIN-ngene);
  return chi;
}
