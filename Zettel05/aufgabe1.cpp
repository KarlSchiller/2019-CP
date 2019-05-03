#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <math.h>  // sqrt()
#include <tuple>
#include "Dateien/service.cpp"

using namespace std;
using namespace Eigen;


// Maximales Nebendiagonalelement finden
tuple<int, int, double> finde_neben(MatrixXd &M){
    // Dimension herausfinden
    int N = M.cols();
    int zeile, spalte;
    double max = 0;
    for(int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            if (i != j){
                if(abs(M(i, j)) > max){
                    max = abs(M(i, j));
                    zeile = i;
                    spalte = j;
                }
            }
        }
    }
    return make_tuple(zeile, spalte, max);
}
// Aufruf der Methode mit:
// tie(zeile_max, spalte_max, max) = finde_neben(A);

// Approximation der Singulärwerte
VectorXd sort_vec(VectorXd &vector, int k){
    VectorXd approx(k);
        for (int i = 0; i < k; i++){
            approx(i) = vector(i);
        }
    return approx;
}

// Approximation der Matrizen
MatrixXd sort_mat(MatrixXd &M, int k) {
    MatrixXd approx(k, k);
    for (int i = 0; i < k; i++){
        for (int j = 0; j < k; j++){
            approx(i, j) = M(i, j);
        }
    }
    return approx;
}

int main()
{
    cout << "Beginn des Programms!" << endl;
    // Einlesen der Daten
    MatrixXd M;
    int k[3] = {10, 20, 50};
    loadData(M, "Dateien/Bild", 512, 512);
    cout << M.rows() << "x" << M.cols() << endl;
    // SVD durchführen
    BDCSVD<MatrixXd> svd(M, ComputeFullU|ComputeFullV);
    VectorXd sing = svd.singularValues();
    MatrixXd U = svd.matrixU();
    MatrixXd V = svd.matrixV();

    // Initialisieren der Variablen
    VectorXd approx;
    MatrixXd approxW, approxU, approxV, approxA;
    ofstream file;
    string filename;

    // Durchführen für die verschiedenen k-Werte
    for (int l=0; l<3; l++){
      // Approximation der Singulärwerte und umschreiben in eine Diagonalmatrix
      approx = sort_vec(sing, k[l]);
      approxW = approx.asDiagonal();

      // Approximation der U- und der V-Matrix
      approxU = sort_mat(U, k[l]);
      approxV = sort_mat(V, k[l]);

      // Transformieren der Matrix
      approxA = approxU*approxW*approxV.transpose();

      // Auslesen in eine txt-Datei
      filename = "build/bild_"+to_string(k[l])+".txt";
      file.open(filename, ios::trunc);
      //file << "# Array" << endl;
      for (int i = 0; i < k[l]; i++){
          file << i << ";";
      }
      file << endl;
      for (int i=0; i < k[l]; i++){
          for (int j = 0; j < k[l]; j++){
              file << approxA(i, j) << "; ";
          }
          file << endl;
      }
      file.close();
    }
    cout << "Ende des Programms!" << endl;
    return 0;
}
