#include <iostream>
// #include <fstream>
#include <Eigen/Dense>
#include <math.h>  // sqrt, round
// #include <tuple>
// #include <complex>
#include <iomanip> // setprecision

using namespace std;
using namespace Eigen;
typedef complex<double> dcomp;


/* Invertiert die ersten @l Stellen von @zahl vor dem Komma in Binärdarstellung
 * und gibt diese als integer Zahl wieder zurück.
 */
int flip_bits(int zahl, int l)
{
  string tmp = "";
  // Lese Zahl bitweise in umgekehrter Reihenfolge aus
  for (int i=0; i<l; i++){
    tmp += to_string(((zahl >> i) & 1));
  }
  // Konvertiere in umgekehrter Reihenfolge zurueck in Integer
  return stoi(tmp, nullptr, 2);
}


/* Nummeriert die Einträge von @input um und schreibt diese in @output.
 * Dabei werden die ersten @m Stellen der Binärzahldarstellung des Index
 * invertiert und wieder in eine integer Zahl umgewandelt
 */
void umnummerierung(VectorXcd &input, VectorXcd &output, int m){
  double dim = pow(2, m);
  for(int i=0; i<dim; i++){
    output(i) = input(flip_bits(i, m));
  }
}


/* diskrete Fouriertransformation als nicht-speicheroptimierter Algorithmus aus der
 * Vorlesung mit 2^m Diskretisierungspunkten.
 */
VectorXcd discrete_fft(double m, VectorXcd &f){
  double N = pow(2, m);
  const dcomp compz (0., 1.);  // complexe Zahl i
  VectorXcd start = VectorXcd::Zero(N);
  umnummerierung(f, start, m);

  // Initialisierung
  MatrixXcd neu = MatrixXcd::Zero(N,N);
  MatrixXcd alt = MatrixXcd::Zero(N,N);
  for (int i=0; i<N; i++) {
    alt.col(i) = start;
  }

  // Iterative Berechnung
  for(int k=1; k<=m; k++)
  {
    for(int j=0; j<pow(2,k); j++)
    {
      for(int l=0; l<pow(2,m-k); l++)
      {
        // Einträge j<2^k
        neu(j,l) = alt(j,2*l)+exp(2*M_PI*compz*static_cast<double>(j)/pow(2,k))*alt(j,2*l+1);

        // Einträge j>=2^k
        for(int i=1; i<pow(2,m-k); i++) {
          neu(j+i*pow(2,k),l) = neu(j,l);
        }
      }
    }
    alt = neu;
  }
  return alt.col(0);
}


/* Verschiebung der Elemente des Ergebnisvektors von discrete_fft,
 * um diese in eine geordnete Reihenfolge zu bringen
 */
void sortiere(VectorXcd &input, VectorXcd &output, int N)
{
  for (int i=0; i<N/2; i++)
  {
    output(2*i) = input(i);  // gerade k
    output(2*i+1) = input(i+N/2);  // ungerade k
  }
}


/* Fast Fourier Transformation
 * INPUT:     m       Legt Anzahl Diskretisierungpunkte N=2^m fest
 *            f       Funktionswerte, Dimension N=2^m
 *            TODO: Wertebereich mit übergeben
 * OUTPUT:    F       Fouriertransformierte
 */
VectorXcd fft(double m, VectorXcd &f){
  // VectorXcd F = discrete_fft(m, f);
  int N = pow(2,m);
  VectorXcd F = VectorXcd::Zero(N);
  sortiere(f, F, N);
  // TODO: Daten mit Phasenfaktor multiplizieren

  return F;
}


int main()
{
  cout << "Aufgabe 2: FFT" << endl;
  cout << "\tTeil a)" << endl;

  cout << "\tdirekte Auswertung" << endl;
  double m = 3;
  double N = pow(2, m);
  const dcomp i (0., 1.);
  VectorXcd f =  VectorXd::LinSpaced(N, 1, N).cwiseSqrt();
  VectorXcd F = VectorXcd::Zero(N);
  for(double j=0; j<N; j++){
    for (double l=0; l<N; l++){
      F(j) += exp(2*M_PI*i*l*j/N)*f(l);
    }
  }
  cout << F << endl << endl;

  cout << "\tFFT" << endl;
  VectorXcd Ftest = discrete_fft(m, f);
  cout << Ftest << endl << endl;
  // TODO: Herausfinden, warum fft komisches Ergebnis hat...

  cout << "\tTeil b)" << endl;
  VectorXcd Fvoll = fft(m, f);


      // // Auslesen in eine txt-Datei
      // filename = "build/bild_"+to_string(k[l])+".txt";
      // file.open(filename, ios::trunc);
      // //file << "# Array" << endl;
      // file.close();
  return 0;
}
