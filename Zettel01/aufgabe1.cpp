#include <iostream>
#include <Eigen/Dense>

// using Eigen::Matrix3i;
// using Eigen::PartialPivLU;
using namespace std;
using namespace Eigen;

int main()
{
    cout << "Beginn des Programms" << endl;
    // Aufgabenteil a)
    //  Intitialisieren der benötigten Variablen
    const int N = 4;
    int A[4][4] = {
        {1, 5, -4, 2},
        {-4, -17, 14, -5},
        {2, 16, -10, 7},
        {6, 33, -22, 13}
    };
    int L[N][N]{};
    int U[N][N]{};
    for (int i = 0; i < N; i++){
        L[i][i] = 1;
    }
    for (int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            U[i][j] = A[i][j];
        }
    }
    // Intitialisieren beendet
    double rij = 0;
    for (int i = 0; i < N-1; i++){
        for (int j = i + 1; j < N; j++){
            // Das Element rij bestimmen
            rij = U[j][i]/U[i][i];
            for (int k = 0; k < N; k++){
                // Fall abfangen, dass U[i][i] 0 ist, denn ansonsten ist Pivotisierung notwendig
                if (U[i][i] != 0){
                    U[j][k] = (U[j][k])-(U[i][k]*rij);
                    L[j][i] = rij;
                }
                else{
                    cout << "Pivotisierung notwendig!";
                    break;
                }
            }
        }
    }
    // Nur für die Ausgabe der in a) bearbeiteten Matrizen
    //for (int i = 0; i < N; i++){
    //    for (int j = 0; j < N; j++){
    //        cout << U[i][j] << " ";
    //    }
    //    cout << "\n";
    //}
    //cout << endl;

    // Aufgabenteil b)
    // initialisieren
    int Ab[N][N] = {
        {1, 5, -4, 2},
        {1, 5, -22, 13},
        {-4, 17, 14, 5},
        {2, 16, -10, 7}
    };
    int P[N][N]{};
    for (int i = 0; i < N; i++){
        P[i][i] = 1;
    }
    double Lb[N][N]{};
    double Ub[N][N]{};
    for (int i = 0; i < N; i++){
        Lb[i][i] = 1;
    }
    for (int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            Ub[i][j] = Ab[i][j];
        }
    }
    int max = 0;
    int max_i = 0;
    // Initialisierung beendet
    for (int i = 0; i < N-1; i++){
        // Betragsgrößtes element finden
                max = abs(Ub[i+1][i]);
                max_i = i+1;
                for (int m = i+1; m < N - 1; m++){
                        if(abs(Ub[m+1][i]) > max){
                            max_i = m + 1;
                            max = abs(Ub[m+1][i]);
                        }
                    }
                /* Nachgeprüfen, ob der maximale Wert != 0,
                 * ansonsten hat auch die Pivotisierung nicht geholfen
                 * und es wird durch 0 geteilt */
                if(max != 0){
                    // zeilen tauschen
                    for (int l = 0; l < N; l++){
                        int c = Ub[i][l];
                        Ub[i][l] = Ub[max_i][l];
                        Ub[max_i][l] = c;
                    }

                    for (int l = 0; l < N; l++){
                        int c = P[l][i];
                        P[l][i] = P[l][max_i];
                        P[l][max_i] = c;
                    }
                    for (int l = 0; l < i; l++){
                        int c = Lb[i][l];
                        Lb[i][l] = Lb[max_i][l];
                        Lb[max_i][l] = c;
                    }
                }
                else {
                    cout << "Pivotisierung fehlgeschlagen!" << endl;
                }
                // LU-Zerlegung durchführen
                for (int j = i + 1; j < N; j++){
                    // Das Element rij bestimmen
                    // cout << "Ubji " << Ub[j][i] << " " << "Ubii " << Ub[i][i] << " " << "rij " << Ub[j][i]/Ub[i][i] << endl;
                    rij = Ub[j][i]/Ub[i][i];
                    for (int k = 0; k < N; k++){
                        Ub[j][k] = (Ub[j][k])-(Ub[i][k]*rij);
                        Lb[j][i] = rij;
                       // cout << "Ubik " << (Ub[i][k]) << " " << "abgez. " << (Ub[i][k]*rij) << endl;
                    }
                }
        }
        cout << endl;

    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            cout << P[i][j] << " ";
        }
        cout << "\n";
    }
    cout << endl;

    // Aufgabenteil c)
    // Matrix A initialisieren
    MatrixXd Ac(3,3);
    Ac << 1, 2, 4, -2, 1, 0, 4, 2, -4;
    MatrixXd M = Ac;
    PartialPivLU<Ref<MatrixXd> > lu(M);

    // Definition von U und L, gewonnen aus M
    MatrixXd Uc(3,3);
    Uc << 4, 2, -4, 0, 2, -2, 0, 0, 6.5;
    MatrixXd Lc(3,3);
    Lc << 1, 0, 0, -0.5, 1, 0, 0.25, 0.75, 1;
    // Lc*Uc ergibt das Pivotisierte A
    //cout << Lc*Uc << endl;
    //cout << endl;
    // Daher kenne ich die Pivotisierungsmatrix :)
    //MatrixXd Pc(3, 3);
    //Pc << 0, 0, 1, 0, 1, 0, 1, 0, 0;
    //cout << Pc*Lc*Uc << endl;
    /*Man kann aber auch ohne Vorwissen über die Pivotisierungsmatrix die Gleichung lösen: */
    /* Über Transponieren und umstellen der Gleichung PA = LU die Pivotisierungmatrix
    rechts von der A-Matrix schreiben, um über solve die Lösung zu bekommen.*/
    MatrixXd At = Ac.transpose();
    MatrixXd Ac1 = Lc*Uc;
    MatrixXd Ac2 = Ac1.transpose();
    lu.compute(At);
    MatrixXd Loesung = lu.solve(Ac2);
    // Transponierte Lösung anzeigen, ist egal, ob nochmal transponiert wird in diesem Fall
    cout << Loesung.transpose() << endl;
    //MatrixXd Ax(4, 4);
    //Ax << 1, 5, -4, 2, 1, 5, -22, 13, -4, 17, 14, 5, 2, 16, -10, 7;
    //PartialPivLU<Ref<MatrixXd> > test(Ax);
    //cout << Ax << endl;

    cout << "Ende des Programms" <<endl;
    return 0;
}
