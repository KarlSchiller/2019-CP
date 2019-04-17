#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include "Dateien/service.cpp"

using namespace std;
using namespace Eigen;


int main()
{
    cout << "Aufgabe 2" << endl;
    // Teil a)
    // Einlesen
    MatrixXd TRAIN;
    cout << "\tImportiere Trainingsdaten" << endl;
    int test = loadData(TRAIN, "Dateien/Training", 112*92, 360);
    if(test==0){
      cout << "\tFehler beim Einlesen der Daten!" << endl;
      return 1;
    }
    // cout << TRAIN.rows() << "x" << TRAIN.cols() << endl;

    // SVD durchführen
    cout << "\tFuehre SVD durch" << endl;
    MatrixXd A = TRAIN;
    BDCSVD<MatrixXd> svd(A, ComputeThinU);
    VectorXd sing = svd.singularValues();
    // cout << "U Thin " << svd.matrixU().rows() << "x" << svd.matrixU().cols() << endl;
    // MatrixXd U = svd.matrixU();

    // Speichere Sprektrum zum plotten im build-Ordner
    cout << "\tSpeichern der Eigenwerte" << endl;
    ofstream outfile;
    outfile.open("build/aufg2-eigenvalues.txt", ios::trunc);
    outfile << "evalues" << endl;
    for (int i=0; i<360; i++){
      outfile << sing(i) << endl;
    }
    outfile.close();

    // Transformiere erstes Bild des Trainingsdatensatzes
    cout << "\tTrafo des ersten Bildes des Trainingsdatensatzes" << endl;
    VectorXd training_pic = A.col(0);
    // Nehme die ersten k Eigenwerte
    int k[2] = {200, 300};
    for(int l=0; l<2; l++)
    {
      VectorXd transformed_pic = VectorXd::Zero(10304);
      for(int i=0; i<k[l]; i++){
        transformed_pic = transformed_pic + training_pic.dot(svd.matrixU().col(i))*svd.matrixU().col(i);
      }

      // Speichere das erste Bild und die transformierte Version
      cout << "\tSpeichern des ersten Bildes und der Trafo" << endl;
      string filename = "build/aufg2-k" + to_string(k[l]) + ".txt";
      outfile.open(filename, ios::trunc);
      for (int i=0; i<10304; i++){
        outfile << training_pic(i) << "; "<< transformed_pic(i) << endl;
      }
      outfile.close();
    }


    // Teil b)

    // BDCSVD<MatrixXd> svd(M, ComputeFullU|ComputeFullV);
    // VectorXd sing = svd.singularValues();
    // VectorXd approx = sort_vec(sing, k);
    // MatrixXd approxW = approx.asDiagonal();

    // MatrixXd U = svd.matrixU();
    // MatrixXd V = svd.matrixV();
    // MatrixXd approxU = sort_mat(U, k);
    // MatrixXd approxV = sort_mat(V, k);

    // MatrixXd approxA = approxU*approxW*approxV.transpose();
    return 0;
}
