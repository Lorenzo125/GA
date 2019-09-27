void arrange(){
  struct chrom pop_aus [N_CHROM_MAX]; // popolazione ausiliaria
  int i,u,tmin; //tmin registra la riga del minimo
  struct chrom chr_aus; //mi serve per lo scambio
  for (i=0;i<nchrom;i++){ //cerco i crom con costo minore tra i seguenti, prendo il minimo e lo metto nella prima riga,
    //quello che era nella prima riga lo metto al posto del crom con costo minimo
    tmin=i;
    for (u=i+1;u<nchrom;u++)
    if (evaluate(&pop[u])<evaluate(&pop[tmin]))
    tmin=u;
    chr_aus=pop[i];
    pop[i]=pop[tmin];
    pop[tmin]=chr_aus;
  };
}
