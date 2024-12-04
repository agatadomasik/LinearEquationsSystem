#include <iostream>
#include <fstream>
#include "Matrix.h"

int main()
{
    // Ex 1
    int N = 977;
    Matrix m = Matrix(N, 10, -1, -1);
    double* b = new double[N];
    double* x = new double[N];
    for (int i = 0; i < N; i++)
        b[i] = sin(i * 4);

    // Ex 2
    int it = 0;
    double duration;
    double* res = new double[N];
    //x = m.solveJacobi(b, 1000, 0.000000001, &it, &duration, res);
    //std::cout << it << " " << duration << std::endl;
    //x = m.solveGaussSeidel(b, 1000, 0.000000001, &it, &duration, res);
    //std::cout << it << " " << duration << std::endl;


    // Ex 3
    Matrix m2 = Matrix(N, 3, -1, -1);
    //x = m2.solveJacobi(b, 100, 0.000000001, &it, &duration, res);
    //x = m2.solveGaussSeidel(b, 100, 0.000000001, &it, &duration, res);

    // Ex 4
    x = m2.solveDirect(b, 100, 0.000000001, &it, &duration, res);
    std::cout << duration << std::endl;

    // Ex 5
    //std::ofstream file("time.txt");

    //Matrix* m3 = new Matrix(100, 10, -1, -1);
    //x = m3->solveJacobi(b, 100, 0.000000001, &it, &duration, res);
    //file << duration << std::endl;
    //x = m3->solveGaussSeidel(b, 100, 0.000000001, &it, &duration, res);
    //file << duration << std::endl;
    //x = m3->solveDirect(b, 100, 0.000000001, &it, &duration, res);
    //file << duration << std::endl << std::endl;

    //std::cout << "1";

    //m3 = new Matrix(500, 10, -1, -1);
    //x = m3->solveJacobi(b, 100, 0.000000001, &it, &duration, res);
    //file << duration << std::endl;
    //x = m3->solveGaussSeidel(b, 100, 0.000000001, &it, &duration, res);
    //file << duration << std::endl;
    //x = m3->solveDirect(b, 100, 0.000000001, &it, &duration, res);
    //file << duration << std::endl << std::endl;

    //std::cout << "2";


    //m3 = new Matrix(1000, 10, -1, -1);
    //x = m3->solveJacobi(b, 100, 0.000000001, &it, &duration, res);
    //file << duration << std::endl;
    //x = m3->solveGaussSeidel(b, 100, 0.000000001, &it, &duration, res);
    //file << duration << std::endl;
    //x = m3->solveDirect(b, 100, 0.000000001, &it, &duration, res);
    //file << duration << std::endl << std::endl;

    //std::cout << "3";


    //m3 = new Matrix(2000, 10, -1, -1);
    //x = m3->solveJacobi(b, 100, 0.000000001, &it, &duration, res);
    //file << duration << std::endl;
    //x = m3->solveGaussSeidel(b, 100, 0.000000001, &it, &duration, res);
    //file << duration << std::endl;
    //x = m3->solveDirect(b, 100, 0.000000001, &it, &duration, res);
    //file << duration << std::endl << std::endl;

    //std::cout << "4";


    //m3 = new Matrix(3000, 10, -1, -1);
    //x = m3->solveJacobi(b, 100, 0.000000001, &it, &duration, res);
    //file << duration << std::endl;
    //x = m3->solveGaussSeidel(b, 100, 0.000000001, &it, &duration, res);
    //file << duration << std::endl;
    //x = m3->solveDirect(b, 100, 0.000000001, &it, &duration, res);
    //file << duration << std::endl << std::endl;

    //file.close();


}

