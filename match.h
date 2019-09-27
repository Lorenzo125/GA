void match (){ // Xfactor e G pari a 1/2, non è un parametro troppo importante e rischia di essere "pericoloso"
double prob [nchrom/2];
int i,u,ma,pa;
double n,k=0.0f;
float alpha;

for (i=1;i<=(nchrom/2);i++) //calcolo le probabilit� sequenziali
k+=i;
for (i=0;i<(nchrom/2);i++)
prob [i]=((nchrom/2)-i)/k;
for (i=1;i<(nchrom/2);i++)
prob [i]+=prob[i-1];
for (i=0;i<nchrom/2;i+=2){ // genero le COPPIE di figli
  n= r();
  for (u=1;u<nchrom/2;u++){
    if (n<=prob[0])
    ma=0;
    if (n>prob [u-1] && n<=prob [u])
    ma=u;
  };
  n= r();
  for (u=1;u<nchrom/2;u++){
    if (n<=prob[0])
    pa=0;
    if (n>prob [u-1] && n<=prob [u])
    pa=u;
  };
  for (u=0;u<ngene;u++){ //uniform crossover
    alpha=r();
    pop [nchrom/2+i].gene[u]=pop[ma].gene[u]-alpha*(pop[ma].gene[u]-pop[pa].gene[u]);
    pop [nchrom/2+i+1].gene[u]=pop[pa].gene[u]+alpha*(pop[ma].gene[u]-pop[pa].gene[u]);
  };
};
}
