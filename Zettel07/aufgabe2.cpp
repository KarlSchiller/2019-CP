#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <math.h>  // sqrt, round
// #include <tuple>
// #include <complex>
// #include <iomanip> // setprecision

using namespace std;
using namespace Eigen;

// Von Aufgabe 1
/*
* f           zu minierende Funktion
* x           aktueller Schritt
* grad        Schrittrichtung
* lambda      zu minimierender Parameter
* epsilon_f   Minimierungstoleranz
*/
// TODO: umschreiben auf double func(Eigen::VectorXd r, double lam, double A0, double mu)
void minimize(function<double(VectorXd, double, double, double)> f,
VectorXd x, VectorXd grad, double &lambda, double epsilon_f,
double lam, double A0, double mu){
  //Minimiere mit Newtonverfahren
  double lambda_prev, lambda_next,acc,h;
  acc = epsilon_f;
  h = 1e-4;
  do {
    double f_prime, f_2prime;
    f_prime = (f(x + (lambda + h) * grad, lam, A0, mu) - f(x + (lambda - h) * grad, lam, A0, mu)) / (2 * h);
    f_2prime = f(x + (lambda + h) * grad, lam, A0, mu) - 2 * f(x + lambda * grad, lam, A0, mu) + f(x + (lambda - h) * grad, lam, A0, mu);
    f_2prime = f_2prime / (h * h);
    lambda_prev = lambda;
    lambda_next = lambda - f_prime / f_2prime;
    lambda = lambda_next;
    cout << abs(f(x+lambda*grad, lam, A0, mu) - f(x+lambda_prev*grad, lam, A0, mu)) << endl;
  } while(abs(f(x+lambda*grad, lam, A0, mu) - f(x+lambda_prev*grad, lam, A0, mu)) >= acc);
}

// Von Aufgabe 1 TODO: funktion und gradient anpassen
        // func(r,lam,A0,mu),
        // grad_alm(r,h,func,lam,A0,mu),

/*
implementation of BFGS
* f           minimized function
* g           gradient function
* x_0         Starting point
* c_0         Initial inverse Hesse-Matrix
* epsilon_f   Tolerance for one-dimensional minimization
* epsilon_g   Tolerance for gradientnorm
*/
VectorXd bfgs(function<double(VectorXd, double, double, double)> f,
              function<VectorXd(VectorXd, double, function<double(VectorXd, double, double, double)> f1,
              double, double, double)> g,
              VectorXd x_0,
              MatrixXd c_0, double epsilon_f, double epsilon_g,
            double lam_func, double A0, double mu, double h){
  /* Fahrplan:
  * p bestimmen
  * Liniensuchschritt durchlaufen
  * Updaten
  */

  VectorXd p, b, x_i, s_k, y_k, b_0;
  double lam = 0, rho = 0;
  int iteration = 0;
  MatrixXd c(c_0.rows(), c_0.cols()), I;
  I = MatrixXd::Identity(c_0.rows(), c_0.cols());

  // Erstes b und p erzeugen
  // Dabei ist p die Minimierungsrichtung

  b_0 = g(x_0, h, f, lam_func, A0, mu);
  p = -c_0*b_0;

  // Liniensuchschritt durchlaufen
  cout << "Bis zur Minimierung alles gut!" << endl;
  minimize(f, x_0, p, lam, epsilon_f, lam_func, A0, mu);
  x_i = x_0 + lam*p;

  // Neues b, und dann s_k und y_k
  b = g(x_i, h, f, lam_func, A0, mu);
  s_k = x_i - x_0;
  y_k = b - b_0;
  c = c_0;
  // Bestimmung von C_k+1
  // in Schleife eingebettet
  do{
  rho = 1/y_k.dot(s_k);

  c = (I-rho*s_k*y_k.transpose())*c*(I-rho*y_k*s_k.transpose()) + rho*s_k*s_k.transpose();

  p = -c*b;
  // Liniensuchschritt durchlaufen und dabei x_0 auf den alten wert setzen
  x_0 = x_i;
  minimize(f, x_i, p, lam, epsilon_f, lam_func, A0, mu);
  x_i = x_0 + lam*p;

  // Bestimmung von s_k und y_k
  s_k = x_i - x_0;
  b_0 = b;
  b = g(x_i, h, f, lam_func, A0, mu);
  y_k = b - b_0;
  iteration++;
  cout << b.norm() << endl;
  }while(b.norm() > epsilon_g);
  //cout << "Iterationen: " << iteration << endl;
  return x_i;
}


/* Berechnet die eingeschlossene Fläche eines 2 dimensionalen Polygonzuges
 * INPUT: Vektor r mit Dimension 2x(#Punkte)
 * OUTPUT: Flächeninhalt der eingeschlossenen Fläche
 */
double flaechePolygonzug(Eigen::VectorXd orte)
{
  // Resizen, weil die Funktion für eine 2x(#Punkte)-Matrix konzipiert ist
  MatrixXd r = orte;
  r.resize(2, orte.size()/2);

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
 * OUTPUT   r   2*N Vektor mit den Punkten des Polygonzuges
 */
Eigen::VectorXd initPoly(int N, double l)
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

  // Resizen der Matrix in einen Vektor
  r.resize(r.rows()*r.cols(), 1);
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


/* Gradient für die Augmentierte Lagrangemethode
 */
Eigen::VectorXd grad_alm(
    VectorXd r,
    double h,
    function<double(VectorXd,double,double,double)> func,
    double lam,
    double A0,
    double mu)
{
  VectorXd grad(r.size());
  VectorXd x_plus, x_minus;

  // Gradienten für jede Komponente berechnen
  for (int i=0; i<r.size(); i++)
  {
    x_plus = r;
    x_plus(i) += h;
    x_minus = r;
    x_minus(i) -= h;
    grad(i) = 0.5*((func(x_plus,lam,A0,mu)-func(x_minus,lam,A0,mu))/h);
  }
  return grad;
}


/* Optimaler Polygonzug mit vorgegebenem Flächeninhalt
 * Augmented Lagrangian Method
 */
Eigen::VectorXd opti_poly(Eigen::VectorXd r, double h, double A0,
    function<double(Eigen::VectorXd,double,double,double)> func)
{
  double lam = 1;
  double mu = 1;
  VectorXd r_neu = r;
  MatrixXd c_0 = MatrixXd::Identity(r.size(), r.size());
  cout << "Vor-Schleife" << endl;
  while((flaechePolygonzug(r)-A0)/A0 >= 1e-6)
  {
    cout << "Durchlauf" << endl;
    r_neu = bfgs(
        func,
        grad_alm,
        r,
        c_0,
        1e-6,
        1e-6,
        lam, A0, mu, h
        );  // TODO: bfgs so umschreiben, dass es hier funktioniert
    lam = lam - mu*(flaechePolygonzug(r)-A0);
    mu *= 2;
    r = r_neu;
  }
  lam = lam - mu*(flaechePolygonzug(r)-A0);

  return r;
}


/* Berechnet Energie des Polygonzuges
 * INPUT:   r   Punkte des Polygonzuges, Dimension des Vektors 2x(#Punkte)
 */
double energie(Eigen::VectorXd orte)
{
  // Resizen, weil die Funktion für eine 2x(#Punkte)-Matrix konzipiert ist
  MatrixXd r = orte;
  r.resize(2, r.size()/2);
  //cout << r << endl;

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
 * INPUT:   r   Punkte des Polygonzuges, Dimension des Vektors 2x(#Punkte)
 *          lam Parameter, siehe Skript
 *          mu  Parameter, siehe Skript
 *          A0  Optimaler Flächeninhalt des Polygonzuges
 * OUTPUT:  Vektor mit Funktionswerten
 */
double func(Eigen::VectorXd r, double lam, double A0, double mu)
{
  double A = flaechePolygonzug(r);
  return energie(r) - lam*(A-A0) + 0.5*mu*(A-A0)*(A-A0);
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
  // r_test.resize(r_test.rows()*r_test.cols(), 1);
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
  // r_test.resize(r_test.rows()*r_test.cols(), 1);
  // cout << r_test << endl;
  // double E_test = energie(r_test);
  // cout << "E sollte 4 ergeben: " << E_test << endl;

  // Initialisiere Polygonzug auf Quadrat
  MatrixXd r = initPoly(N, l);
  VectorXd test = opti_poly(r, h, A0, func);
  //cout << test << endl;
  //cout << flaechePolygonzug(r) << endl;

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
