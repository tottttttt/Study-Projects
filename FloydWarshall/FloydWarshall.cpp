// FloydWarshall.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//


#include <bits/stdc++.h>
using namespace std;
#define INF INT_MAX /*implementation defined*/


void printMatrix(int** matrix, int numberOfVert) {
    for (int i = 0; i < numberOfVert; i++) {
        for (int j = 0; j < numberOfVert; j++) {
            if (matrix[i][j] == INF) {
                cout << "INFINITY" << " ";
            }
            else {
                cout << matrix[i][j] << " ";
            }
        }
        cout << endl;
    }
}
double FloydWarshall(int** matrix, int numberOfVert) {
    
    double start_time = clock();
    for (int k = 0; k < numberOfVert; k++) {
        for (int i = 0; i < numberOfVert; i++) {
            for (int j = 0; j < numberOfVert; j++) {
                matrix[i][j] = min(matrix[i][j], matrix[i][k] + matrix[k][j]);
            }
        }
    }
    return (clock() - start_time) / 1000;
}



int intRand(int min, int max) {
    int a = (int)(rand()) / INT_MAX * (max - min) + min;
    return a;
}
double test(int numberOfVert) {
    
    
    int** matrix = new int* [numberOfVert];
    for (int i = 0; i < numberOfVert; i++) {
        matrix[i] = new int[numberOfVert];
    }


    for (int i = 0; i < numberOfVert; i++) {
        for (int j = 0; j < numberOfVert; j++) {
            matrix[i][j] = rand() % 101 + 1;
        }
    }

    for (int i = 0; i < numberOfVert; i++) {
        matrix[i][i] = 0;
    }
    
 
    cout << "Old matrix" << endl;
    printMatrix(matrix, numberOfVert);

    double result = FloydWarshall(matrix, numberOfVert);

    cout << "New matrix" << endl;

    printMatrix(matrix, numberOfVert);

    for (int i = 0; i < numberOfVert; i++) {
        delete[] matrix[i];
    }
    delete[] matrix;



    return result;
}
int main() {
    
   
    int n;
    
    cin >> n;
    test(n);

    return 0;
}

