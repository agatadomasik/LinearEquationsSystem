#pragma once
#include "Matrix.h"
#include <iostream>
#include <chrono>
#include <fstream>



Matrix::Matrix(int N, double a1, double a2, double a3) {
	this->size = N;
	this->A = new double* [N];
	for (int i = 0; i < N; i++)
		this->A[i] = new double[N];

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (j == i) this->A[i][j] = a1;
			else if (j == i + 1 || j == i - 1) this->A[i][j] = a2;
			else if (j == i + 2 || j == i - 2) this->A[i][j] = a3;
			else this->A[i][j] = 0;
		}
	}
}

Matrix::Matrix(const Matrix& m)
{
	this->size = m.size;
	this->A = new double* [this->size];

	for (int i = 0; i < this->size; ++i)
		this->A[i] = new double[this->size];

	for (int i = 0; i < this->size; ++i)
		for (int j = 0; j < this->size; ++j)
			this->A[i][j] = m.A[i][j];
	
}

void Matrix::print() {
	for (int i = 0; i < this->size; i++) {
		for (int j = 0; j < this->size; j++)
			std::cout << this->A[i][j] << " ";
		std::cout << std::endl;
	}
}

double* Matrix::operator*(double* b) {
	double* c = new double[size];
	for (int i = 0; i < size; ++i) {
		double a = 0;
		for (int j = 0; j < size; ++j) {
			a += A[i][j] * b[j];
		}
		c[i] = a;
	}
	return c;
}

Matrix::~Matrix() {
	for (int i = 0; i < size; ++i)
		delete[] A[i];
	delete[] A;
}


double* Matrix::solveJacobi(double* b, int it_max, double err_max, int *it, double* duration, double* res) {
	*it = 1;
	double* x = new double[size];
	double* x_prev = new double[size];
	for (int i = 0; i < size; ++i) x_prev[i] = 1;

	std::ofstream file("jacobi_error.txt");
	auto start = std::chrono::high_resolution_clock::now();

	while ((*it) < it_max) {
		for (int i = 0; i < size; ++i) {
			double s = 0;
			for (int j = 0; j < size; ++j) {
				if (j != i)
					s += A[i][j] * x_prev[j];
			}
			x[i] = (b[i] - s) / A[i][i];
		}
		for (int i = 0; i < size; ++i) x_prev[i] = x[i];

		res = (*this) * x;
		for (int i = 0; i < size; ++i) res[i] -= b[i];


		// norma Frobeniusa
		double err_norm = 0;
		for (int i = 0; i < size; ++i) err_norm += pow(res[i], 2);
		err_norm = sqrt(err_norm);
		file << err_norm << std::endl;

		if (err_norm < err_max) break;

		(*it)++;
	}

	auto end = std::chrono::high_resolution_clock::now();
	auto difference = end - start;
	*duration = std::chrono::duration<double, std::milli>(difference).count();
	file.close();

	return x;
}


double* Matrix::solveGaussSeidel(double* b, int it_max, double err_max, int* it, double* duration, double* res) {

	*it = 1;
	double s;
	double* x = new double[size];
	double* x_prev = new double[size];
	for (int i = 0; i < size; ++i) x_prev[i] = 1;

	std::ofstream file("gauss_error.txt");
	auto start = std::chrono::high_resolution_clock::now();

	while (*it < it_max) {

		for (int i = 0; i < size; ++i)
		{
			s = 0;

			for (int j = 0; j < i; ++j) 
				s += A[i][j] * x[j];

			for (int j = i + 1; j < size; ++j) 
				s += A[i][j] * x_prev[j];

			x[i] = (b[i] - s) / A[i][i];
		}

		for (int i = 0; i < size; ++i) x_prev[i] = x[i];

		// wektor rezydualny
		res = (*this) * x;
		for (int i = 0; i < size; ++i) res[i] -= b[i];

		// norma Frobeniusa
		double err_norm = 0;
		for (int i = 0; i < size; ++i) err_norm += pow(res[i], 2);
		err_norm = sqrt(err_norm);
		file << err_norm << std::endl;

		if (err_norm < err_max) break;

		(*it)++;

	}

	auto end = std::chrono::high_resolution_clock::now();
	auto difference = end - start;
	*duration = std::chrono::duration<double, std::milli>(difference).count();
	file.close();

	return x;
}

double* Matrix::solveDirect(double* b, int it_max, double err_max, int* it, double* duration, double* res) {
	Matrix* U = new Matrix(*this);
	Matrix* L = new Matrix(size, 1, 0, 0);
	double* x = new double[size];


	std::ofstream file("direct_error.txt");
	auto start = std::chrono::high_resolution_clock::now();

	for (int i = 0; i < size - 1; ++i) {
		for (int j = i + 1; j < size; ++j) {
			L->A[j][i] = U->A[j][i] / U->A[i][i];

			for (int k = i; k < size; ++k)
				U->A[j][k] = U->A[j][k] - L->A[j][i] * U->A[i][k];
		}
	}

	double* y = new double[size];
	for (int i = 0; i < size; ++i)
	{
		double s = 0;

		for (int j = 0; j < i; ++j) 
			s += L->A[i][j] * y[j];

		y[i] = (b[i] - s) / L->A[i][i];
	}

	for (int i = size - 1; i >= 0; --i)
	{
		double s = 0;

		for (int j = i + 1; j < size; ++j) 
			s += U->A[i][j] * x[j];

		x[i] = (y[i] - s) / U->A[i][i];
	}

	auto end = std::chrono::high_resolution_clock::now();
	auto difference = end - start;
	*duration = std::chrono::duration<double, std::milli>(difference).count();

	// wektor rezydualny
	res = (*this) * x;
	for (int i = 0; i < size; ++i) res[i] -= b[i];

	// norma Frobeniusa
	double err_norm = 0;
	for (int i = 0; i < size; ++i) err_norm += pow(res[i], 2);
	err_norm = sqrt(err_norm);
	file << err_norm << std::endl;
	file.close();

	return x;
}