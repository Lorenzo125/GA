// mutazioni

void mutat (){
int ra1,ra2,u,v;
for (u=0;u<num_mut;){
  ra1=rand()%nchrom;
  ra2=rand()%ngene;
  for (v=0;v<ngene;v++){
    if (ra1!=0 && ra2==v){ //escluso il primo numero per l'elitismo
      pop[ra1].gene[v]=r()*(dominio[v][1]-dominio[v][0])+dominio[v][0];
      u++;
    };
  };
};
}
