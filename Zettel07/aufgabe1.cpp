#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <math.h>  // sqrt()

using namespace std;
using namespace Eigen;

// Analytisch bestimmte Gradientenfunktion
VectorXd gradient(VectorXd x){
  VectorXd temp(5);
  temp(0) = 2*1*(x(0)-1) + -3*x(1);
  temp(1) = 2*3*(x(1)+1) + -3*x(0) + 2*x(2);
  temp(2) = 2*2*(x(2)-4) + 2*x(1)   -2*x(3);
  temp(3) = 2*2*(x(3)+2) + -2*x(2) + 1*x(4);
  temp(4) = 2*1*(x(4)+3) + 1*x(3);
  return temp;
}

// Funktion f
double f(VectorXd x){
  VectorXd a(5), b(5), d(4);
  double temp = 0;
  a << 1, 3, 2, 2, 1;
  b << 1, -1, 4, -2, -3;
  d << -3, 2, -2, 1;
  for(int i = 0; i < a.size(); i++){
    temp += a(i)*(x(i)-b(i))*(x(i)-b(i));
  }
  for(int i = 0; i < d.size(); i++){
    temp += d(i)*x(i)*x(i+1);
  }
  return temp;
}

/*
* f           zu minierende Funktion
* x           aktueller Schritt
* grad        Schrittrichtung
* lambda      zu minimierender Parameter
* epsilon_f   Minimierungstoleranz
*/
void minimize(function<double(VectorXd)> f,VectorXd x, VectorXd grad, double &lambda, double epsilon_f){
  //Minimiere mit Newtonverfahren
  double lambda_prev, lambda_next,acc,h;
  acc = epsilon_f;
  h = 1e-4;
  do {
    double f_prime, f_2prime;
    f_prime = (f(x + (lambda + h) * grad) - f(x + (lambda - h) * grad)) / (2 * h);
    f_2prime = f(x + (lambda + h) * grad) - 2 * f(x + lambda * grad) + f(x + (lambda - h) * grad);
    f_2prime = f_2prime / (h * h);
    lambda_prev = lambda;
    lambda_next = lambda - f_prime / f_2prime;
    lambda = lambda_next;
  } while(abs(f(x+lambda*grad) - f(x+lambda_prev*grad)) >= acc);
}

/*
implementation of BFGS
* f           minimized function
* g           gradient function
* x_0         Starting point
* c_0         Initial inverse Hesse-Matrix
* epsilon_f   Tolerance for one-dimensional minimization
* epsilon_g   Tolerance for gradientnorm
* fstream     file for iteration steps
*/
VectorXd bfgs(double (*f)(VectorXd), VectorXd (*g)(VectorXd), VectorXd x_0,
              MatrixXd c_0, double epsilon_f, double epsilon_g, int &iteration){
  /* Fahrplan:
  * p bestimmen
  * Liniensuchschritt durchlaufen
  * Updaten
  */

  VectorXd p, b, x_i, s_k, y_k, b_0;
  double lam = 0, rho = 0;
  MatrixXd c(c_0.rows(), c_0.cols()), I;
  I = MatrixXd::Identity(c_0.rows(), c_0.cols());

  // Erstes b und p erzeugen
  // Dabei ist p die Minimierungsrichtung

  b_0 = g(x_0);
  p = -c_0*b_0;

  // Liniensuchschritt durchlaufen
  //cout << "Bis zur Minimierung alles gut!" << endl;
  minimize(f, x_0, p, lam, epsilon_f);
  x_i = x_0 + lam*p;

  // Neues b, und dann s_k und y_k
  b = g(x_i);
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
  minimize(f, x_i, p, lam, epsilon_f);
  x_i = x_0 + lam*p;

  // Bestimmung von s_k und y_k
  s_k = x_i - x_0;
  b_0 = b;
  b = g(x_i);
  y_k = b - b_0;
  iteration++;
  }while(b.norm() > epsilon_g);
  cout << "Iterationen: " << iteration << endl;
  return x_i;
}

int main() {
  cout << "Beginn des Programms!" << endl;
  // Initialisieren der Größen
  VectorXd x_0(5), x_i(5);
  MatrixXd c_0;
  double epsilon_f = 1e-6, epsilon_g = 1e-6;
  int iteration = 0;

  x_0 << 0, 0, 0, 0, 0;
  c_0 = MatrixXd::Identity(5, 5);
  c_0 = f(x_0)*c_0;
  x_i = bfgs(f, gradient, x_0, c_0, epsilon_f, epsilon_g, iteration);
  cout << f(x_0) << " " << f(x_i) << endl;

  cout << "Ende des Programms!" << endl;
  return 0;
}
