#pragma once

class Matrix {
	double** A;
	int size;

public:
	Matrix(int N, double a1, double a2, double a3);
	Matrix(const Matrix& m);
	~Matrix();
	void print();
	double* solveJacobi(double* b, int it_max, double err_max, int *it, double* duration, double* res);
	double* solveGaussSeidel(double* b, int it_max, double err_max, int *it, double* duration, double* res);
	double* solveDirect(double* b, int it_max, double err_max, int* it, double* duration, double* res);
	double* operator*(double* v);
};