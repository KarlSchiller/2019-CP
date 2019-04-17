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
    BDCSVD<MatrixXd> svd(A, ComputeFullV|EigenvaluesOnly);
    VectorXd sing = svd.singularValues();
    // MatrixXd V = svd.matrixV();

    // Speichere Sprektrum zum plotten im build-Ordner
    cout << "\tSpeichern der Eigenwerte" << endl;
    ofstream outfile;
    outfile.open("build/aufg2-eigenvalues.txt", ios::trunc);
    outfile << "evalues" << endl;
    for (int i=0; i<360; i++){
      outfile << sing(i) << endl;
    }
    outfile.close();

    // Nehme die ersten k Eigenwerte
    int k = 200;
    // MatrixXd reducedV = svd.matrixV().block(0,0,k,10304);
    // cout << reducedV.rows() << "x" << reducedV.cols() << endl;

    // Transformiere erstes Bild des Trainingsdatensatzes
    cout << "\tTrafo des ersten Bildes des Trainingsdatensatzes" << endl;
    VectorXd training_pic = A.col(0);
    VectorXd transformed_pic = VectorXd::Zero(10304);
    cout << "transformed pic " << transformed_pic.size() << endl;
    cout << "V column 0      " << svd.matrixV().row(0).size() << endl;
    // for(int i=0; i<k; i++){
      // transformed_pic = transformed_pic + training_pic.dot(svd.matrixV().col(i))*svd.matrixV().col(i);
    // }

    // // Speichere das erste Bild und die transformierte Version
    // cout << "\tSpeichern des ersten Bildes und der Trafo" << endl;
    // // ofstream firstpic;
    // outfile.open("build/aufg2-firstpic.txt", ios::trunc);
    // outfile << "input; k200" << endl;
    // for (int i=0; i<10304; i++){
      // outfile << training_pic(i) << transformed_pic(i) << endl;
    // }
    // outfile.close();


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
