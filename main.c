// ROULETTE WHEEL CON ELITISMO IBRIDIZZATO CON RANDOM SEARCH

#include <iostream>
#include <cstdlib>
#include <string>
using namespace std;

#define N_CHROM_MAX 200
#define N_GENE_MAX 200
#define ITERAZ 10
#define ITERAZ_RAND 10 //numero iterazioni dentro al random searcher
#define TOLL  0.1 //dimensione intorno spazio di ricerca del random searcher
#define BIN 1200 //numero bin nell'istogramma dei dati
#define EVENTI 100000 //numero di eventi che vanno a definire l'istogramma


double r(); // numero random tra 0 a 1 nei reali
double evaluate (struct chrom *chr); //valuta il cromosoma
void arrange (); //ordina dal minore al maggiore
void match (); // seleziona e accoppia
void mutat (); // mutazioni
void print_pop (); //stampa popolazione
void print_hist_best (); //stampa istogramma andamento del best
void fill_with_mean (int num_iteraz, int num_chrom); // costruisce l'istogramma che traccia il valore medio del chi
void print_hist_mean (); // stampa istogramma andamento chi medio
void search_rand (struct chrom *chr); // usa metodo random search

struct chrom {
  double gene [N_GENE_MAX];
};

struct chrom pop [N_CHROM_MAX]; //popolazione di cromosomi

double dominio [N_GENE_MAX][2]; //tiene traccia degli estremi del dominio, [0][1] contiene x_max, [1][0] y_min

int nchrom, ngene, num_mut; //numero cromosomi, numero geni, numero mutazioni

TF1 f_start("f_start","exp(-1.5*x)+0.2*(x-0.6)^2",0,1); //funzione da cui ricavo i dati

TF1 f("f", "(exp([0]*x)+[1]*(x-[2])^2)*[3]",0,1); //modello funzione di fit

TF1 f_gaus ("f_gaus", "TMath::Gaus(x,0,[0]*0.01)",-10,10); //l'ideale sarebbe mettere in una define lo 0.01 e il dominio tale che si adatti correttamente

TH1F h("h", "Data and fit", BIN, 0, 1); //istogramma con dati

TH1F h_noise("h_noise", "Data and fit with noise", BIN, 0, 1); //istogramma con valori aggiornati con rumore Gaussiano

TH1F hbest ("hbest", "Best #chi^{2}/ndof", ITERAZ, 0, ITERAZ); //registra valore del chi migliore per ogni generazione

TH1F hmean ("hmean", "Mean #chi^{2}/ndof", ITERAZ, 0, ITERAZ); //registra valore del chi medio per ogni generazione

TCanvas c1("c1", "", 1500,500); //ambiente di lavoro

#include "arrange.h"
#include "search_rand.h"
#include "match.h"
#include "mutat.h"
#include "r.h"
#include "print_pop.h"
#include "print_hist_best.h"
#include "fill_with_mean.h"
#include "print_hist_mean.h"
#include "evaluate.h"

int main()
{
  srand(time_t(NULL)); //inizializzazione
  c1.Divide (4,1);
  int i,l,m;
  hbest.GetXaxis()->SetTitle("# generzione");
  c1.cd(1);
  h.Draw();
  print_hist_best ();
  print_hist_mean ();

  for (i=0;i<EVENTI; i++) { //genero istogramma
    double value = f_start.GetRandom();
    h.Fill(value);
  };

  c1.cd(4);
  for (i=1;i<=BIN;i++){ //aggiungo termine di noise
    f_gaus.SetParameter(0,h.GetBinContent(i));
    f_gaus.Draw();
    h_noise.Fill((i-1)*(h_noise.GetBinWidth(i)),h.GetBinContent(i)+f_gaus.GetRandom());
    //h_noise.Fill((i-1)*(h_noise.GetBinWidth(i)),h.GetBinContent(i)); //senza rumore
    c1.Modified();
    c1.Update();
  };
  c1.cd(1);
  h_noise.Draw("HIST");

  nchrom=84;
  ngene=4;
  num_mut=17;
  dominio [0][0]=-2.25;//canc
  dominio [0][1]=-0.75;//canc
  dominio [1][0]=0.1;//canc
  dominio [1][1]=0.3;//canc
  dominio [2][0]=0.3;//canc
  dominio [2][1]=0.9;//canc

  double int_hist, int_func, domin_middle; //stima del valore attorno cui si sviluppa N
  int_hist=h_noise.Integral(0,BIN,"width");
  domin_middle=0;
  do {
    domin_middle+=10;
    dominio[ngene-1][0]=domin_middle;
    for (m=0;m<ngene;m++)
    f.SetParameter (m, dominio[m][0]);
  } while (int_hist/f.Integral(f.GetXmin(),f.GetXmax())>1);
  dominio [3][0]=domin_middle*2/3; //valori di approssimazione
  dominio [3][1]=domin_middle*4/3;

  for (l=0;l<ngene;l++){ //genero la popolazione
    for (i=0;i<nchrom;i++){
      pop [i].gene[l]=r()*(dominio[l][1]-dominio[l][0])+dominio[l][0];
    };
  };

  for (i=0;i<ITERAZ;i++){ //ciclo in cui si accoppiano, mutano e ordinano i cromosomi
    match();
    if(i!=(ITERAZ-1))
    mutat ();
    search_rand (&pop[0]);
    arrange();
    hbest.Fill(i,evaluate(&pop[0]));
    fill_with_mean (i,nchrom);
    print_pop ();
    print_hist_best ();
    print_hist_mean ();
    c1.Modified();
    c1.Update();
  };
  return 0;
}
