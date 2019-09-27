//stampa geni miglior cromosoma, funzione di fit e istogramma con rumore

void print_pop (){
  c1.cd(1);
  int m,u;
  for (u=0;u<ngene;u++)
  printf ("%f_",pop [0].gene[u]);
  cout << endl;
  h_noise.Draw("HIST");
  for (m=0;m<ngene;m++)
  f.SetParameter (m, pop[0].gene[m]);
  f.Draw("Same");
}
