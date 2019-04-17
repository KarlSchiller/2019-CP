#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <math.h>  // sqrt()
#include "Dateien/service.cpp"

using namespace std;
using namespace Eigen;


int main()
{
    cout << "Aufgabe 2" << endl;
    // Teil a)
    cout << "\tTeil a)" << endl;
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
    // Verwende nun vollständigen Eigenraum mit 360 Basisvektoren
    cout << "\tTeil b)" << endl;
    cout << "\tBerechne Entwicklungskoeffizienten der Trainingsdaten" << endl;
    MatrixXd transformedTRAIN = svd.matrixU().transpose()*TRAIN;

    // Lese Testdaten ein
    // Einlesen
    MatrixXd TEST;
    cout << "\tImportiere Testdaten" << endl;
    test = loadData(TEST, "Dateien/Testdata", 112*92, 40);
    if(test==0){
      cout << "\tFehler beim Einlesen der Daten!" << endl;
      return 1;
    }
    cout << "\tBerechne Entwicklungskoeffizienten der Testdaten" << endl;
    MatrixXd transformedTEST = svd.matrixU().transpose()*TEST;

    // Berechne Abstaende der Testbilder zu den Trainingsbildern
    // Leider keine tolle Eigen-Funktion dafuer oder fuer die Euklid Distanz gefunden :/
    cout << "\tBerechne Distanzen Training <-> Test" << endl;
    MatrixXd dist = MatrixXd::Zero(360, 40);
    VectorXd difference;
    for (int i=0; i<40; i++){
      for (int j=0; j<360; j++){
        difference = transformedTEST.col(i)-transformedTRAIN.col(j);
        dist(j, i) = sqrt(difference.dot(difference));
      }
    }

    // Suche für jedes Testbild den Index der minimalen Distanz
    // Das ist dann der Index des zugehörigen Trainingsbildes
    cout << "\tSpeichere Distanzen" << endl;
    Eigen::MatrixXd::Index min_index;
    outfile.open("build/aufg2-distanzen.txt", ios::trunc);
    int wrong = 0;  // Anzahl falsch zugeordneter Bilder
    for (int i=0; i<40; i++){
      dist.col(i).minCoeff(&min_index);
      outfile << min_index << endl;
      if (min_index < 9*(i) || min_index >= 9*(i+1)){
        wrong++;
      }
    }
    outfile.close();
    cout << "\t" << wrong << " Bilder falsch zugeordnet" << endl;

    return 0;
}
