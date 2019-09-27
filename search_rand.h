// usa metodo random search per minimizzare cromosoma con buone qualità

void search_rand (struct chrom *chr) {
  int u,i,m,l,tmin;
  struct chrom chr_best=*chr;
  struct chrom pop_aus [N_CHROM_MAX]; //popolazione ausiliaria
  double dominio_aus [N_GENE_MAX][2]; //dominio ausiliario
  for (u=0;u<ngene;u++){ //definisco il dominio di ricerca (ausiliario) centrato nel cromosoma da minimizzare in modo da non fuoriuscire dai bordi
    if (chr->gene[u]-TOLL>dominio[u][0])
    dominio_aus [u][0]=(chr->gene[u])-TOLL;
    else dominio_aus [u][0]=dominio [u][0];
    if (chr->gene[u]+TOLL<dominio[u][1])
    dominio_aus [u][1]=(chr->gene[u])+TOLL;
    else dominio_aus [u][1]=dominio [u][1];
  };
  for (m=0;m<ITERAZ_RAND;m++){
    for (l=0;l<ngene;l++){ //genero la popolazione
      for (i=0;i<nchrom;i++){
        pop_aus [i].gene[l]=r()*(dominio_aus[l][1]-dominio_aus[l][0])+dominio_aus[l][0];
      };
    };
    tmin=0; //trovo cromosoma con costo minore
    for (u=0;u<nchrom;u++)
    if (evaluate(&pop_aus[u])<evaluate(&pop_aus[tmin]))
    tmin=u;
    if (evaluate(&pop_aus[tmin])<evaluate(&chr_best))
    chr_best=pop_aus[tmin];
  };
  *chr=chr_best;
}
