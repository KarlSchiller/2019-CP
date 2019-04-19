#include <iostream>
// #include <fstream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <math.h>

using namespace std;
using namespace Eigen;

// init hamilton in position representation
MatrixXd hamilton_pos_repr(float delta, float lam, int dim)
{
  MatrixXd h(dim,dim);
  // non-diagonal
  h.diagonal<1>() = -pow(delta,-2)*ArrayXd::Ones(dim-1);
  h.diagonal<-1>() = -pow(delta,-2)*ArrayXd::Ones(dim-1);
  // diagonal
  ArrayXd tempvector = ArrayXd::Ones(dim);
  for (int i=0; i<dim; i++)
  {
    tempvector(i) = 2*pow(delta,-2)+delta*delta*i*i+lam*pow(delta,4)*pow(i,4);
  }
  h.diagonal() = tempvector;
  return h;
}


int main()
{
  cout << "Aufgabe 2" << endl;

  const int L = 10;
  const float delta = 0.1;  // discretization steps
  const float lam = 0.2;  // lambda, distortion

  cout << "\tTeil a)" << endl;
  const int dim = L/delta;
  cout << "\tDimension " << dim << endl;

  // init hamilton
  MatrixXd h = hamilton_pos_repr(delta, lam, dim);

  // compute eigenvalues
  VectorXd eivals = h.eigenvalues().real();
  std::sort(eivals.data(), eivals.data()+eivals.size());
  cout << "\t10 minor eigenvalues:" << endl
    << "\t" << eivals.head(10).transpose() << endl;

  cout << "\tTeilb)" << endl;

  // Save to file
  // filename = "build/bild_"+to_string(k[l])+".txt";
  // file.open(filename, ios::trunc);
  // //file << "# Array" << endl;
  // file.close();
  return 0;
}
