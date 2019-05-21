#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <math.h>  // sqrt, round
// #include <tuple>
// #include <complex>
// #include <iomanip> // setprecision

using namespace std;
using namespace Eigen;


/* Berechnet die eingeschlossene Fläche eines 2 dimensionalen Polygonzuges
 * INPUT: Matrix r mit Dimension 2x(#Punkte)
 * OUTPUT: Flächeninhalt der eingeschlossenen Fläche
 */
double flaechePolygonzug(Eigen::MatrixXd &r)
{
  // Berechne Schwerpunkt
  VectorXd r_com = r.rowwise().sum() / r.rows();
  double A = 0;
  // Addiere Flächeninhalte der Dreiecke (r_com, r_i, r_i+1) auf
  for (int i=0; i<r.cols(); i++)
  {
    if(i==r.cols()-1)  // r_0 und r_N sind miteinander verbunden
    {
      A += 0.5*((r(1,0)-r_com(1))*(r(0,i)-r_com(0))+(r_com(0)-r(0,0))*(r(1,i)-r_com(1)));
    } else {
      A += 0.5*((r(1,i+1)-r_com(1))*(r(0,i)-r_com(0))+(r_com(0)-r(0,i+1))*(r(1,i)-r_com(1)));
    }
  }
  return A;
}


/* Initialisiere Polygonzug in 2 Dimensionen
 * Dabei liegen @N Punkte gleichmäßig verteilt auf einem Quadrat der Kantenlänge @l,
 * das im Urpsprung zentriert ist.
 * INPUT    N   # Punkte
 *          l   Seitenlänge des Quadrats
 * OUTPUT   r   2xN Matrix mit den Punkten des Polygonzuges
 */
Eigen::MatrixXd initPoly(int N, double l)
{
  double bor = l/2;  // Grenze des Quadrats in x- und y
  int pps = N/4;  // Punkte pro Seite
  MatrixXd r = MatrixXd::Zero(2, N);

  // Obere Seite des Quadrats
  r.block(0,0,1,pps) = ArrayXd::LinSpaced(pps+1, bor, -bor).head(pps).transpose();
  r.block(1,0,1,pps) = bor*ArrayXd::Ones(pps).transpose();
  // Linke Seite des Quadrats
  r.block(0,pps,1,pps) = -bor*ArrayXd::Ones(pps).transpose();
  r.block(1,pps,1,pps) = ArrayXd::LinSpaced(pps+1, bor, -bor).head(pps).transpose();
  // Untere Seite des Quadrats
  r.block(0,2*pps,1,pps) = ArrayXd::LinSpaced(pps+1, -bor, bor).head(pps).transpose();
  r.block(1,2*pps,1,pps) = -bor*ArrayXd::Ones(pps).transpose();
  // Rechte Seite des Quadrats
  r.block(0,3*pps,1,pps) = bor*ArrayXd::Ones(pps).transpose();
  r.block(1,3*pps,1,pps) = ArrayXd::LinSpaced(pps+1, -bor, bor).head(pps).transpose();
  return r;
}


/* Berechne 2 dimensionalen Gradienten von Funktion f am Ort x mit Schrittweite h
 */
Eigen::VectorXd gradient_2d(double (*f)(double, double), Eigen::VectorXd x, double h)
{
  VectorXd g = VectorXd::Zero(2);
  g(0) = 0.5*(f(x(0)+h,x(1))-f(x(0)-h,x(1)))/h;
  g(1) = 0.5*(f(x(0),x(1)+h)-f(x(0),x(1)-h))/h;
  return g;
}


/* Optimaler Polygonzug mit vorgegebenem Flächeninhalt
 * Augmented Lagrangian Method
 */
Eigen::VectorXd opti_poly(Eigen::MatrixXd r, function<double(VectorXd)> f) // TODO: Funktion schreiben
{
  double mu = 1;
  double lam = 1;

  return VectorXd::Ones(9);
}


/* Berechnet Energie des Polygonzuges
 * INPUT:   r   Punkte des Polygonzuges, Dimension der Matrix 2x(#Punkte)
 */
double energie(Eigen::MatrixXd r)
{
  double en = 0;
  for (int i=0; i<r.cols(); i++)
  {
    if(i==r.cols()-1)  // letzter Punkt, mit erstem vernküpft
    {
      en += abs((r.col(0)-r.col(i)).sum());
    } else {
      en += abs((r.col(i+1)-r.col(i)).sum());
    }
  }
  return en;
}


/* Zu minimierende Funktion der Augmented Lagrangian Method
 * INPUT:   r   Punkte des Polygonzuges, Dimension der Matrix 2x(#Punkte)
 *          lam Parameter, siehe Skript
 *          mu  Parameter, siehe Skript
 *          A0  Optimaler Flächeninhalt des Polygonzuges
 * OUTPUT:  Vektor mit Funktionswerten
 */
double func(Eigen::MatrixXd r, double lam, double A0, double mu)
{
  double A = flaechePolygonzug(r);
  return energie(r) - lam*(A-A0) + 0.5*mu*(A-A0)*(A-A0);
}


double testfunction(double x1, double x2)
{
  return x1*x1 + 2*x2*x2;
}


int main() {
  cout << "Beginn Aufgabe 2!\n" << endl;

  double h = 1e-6;  // Schrittweite des Gradienten
  int N = 80;  // #Punkte
  double l = 1.;  // Länge der Seiten des Quadrates
  double A0 = M_PI;  // Idealer Flächeninhalt
  ofstream stream;

  // Testumgebung, ob Flächenberechnung klappt!
  // Testobjekte: Matrix mit den einzelnen Ortsvektoren in 2D
  // MatrixXd r_test = 0.5*MatrixXd::Ones(2, 5);
  // r_test(0, 1) *= -1;
  // r_test(0, 2) *= -1;
  // r_test(1, 2) *= -1;
  // r_test(1, 3) *= -1;
  // r_test(0, 4) = 0.5;
  // r_test(1, 4) = 0.25;
  // cout << r_test << endl;
  // double A_test = flaechePolygonzug(r_test);
  // cout << "A sollte 1 ergeben: " << A_test << endl;

  // Testumgebung, ob Energieberechnung klappt!
  // MatrixXd r_test = 0.5*MatrixXd::Ones(2, 4);
  // r_test(0, 1) *= -1;
  // r_test(0, 2) *= -1;
  // r_test(1, 2) *= -1;
  // r_test(1, 3) *= -1;
  // // r_test(0, 4) = 0.5;
  // // r_test(1, 4) = 0.25;
  // cout << r_test << endl;
  // double E_test = energie(r_test);
  // cout << "E sollte 4 ergeben: " << E_test << endl;

  // Teste Gradientenmethode
  // VectorXd x_test = VectorXd::Ones(2);
  // cout << gradient_2d(testfunction, x_test, h) << endl;

  // Initialisiere Polygonzug auf Quadrat
  MatrixXd r = initPoly(N, l);
  // cout << r << endl;
  // cout << flaechePolygonzug(r) << endl;

  // Speichere Startkonfiguration
  stream.open("build/aufg2-l1-start.txt");
  stream << "x y" << endl;
  for(int i=0; i<r.cols(); i++)
  {
    for(int j=0; j<r.rows(); j++)
    {
      stream << r(j, i);
      if(j != r.rows()-1){ stream << " "; }
    }
    stream << endl;
  }
  stream.close();

  cout << "\nEnde Aufgabe 2!" << endl;
  return 0;
}
