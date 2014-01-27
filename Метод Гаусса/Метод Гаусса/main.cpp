#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <time.h>
#include <iostream>

using namespace std;

const int n = 10;
double A[n][n], B[n];
double Res[n];

double drand(double dmin, double dmax){
	double f = (double)rand() / RAND_MAX;
	return dmin + f * (dmax - dmin);
}

int main(){
	srand(time(NULL));

	for (int i = 0; i < n; ++i)
	for (int j = 0; j < n; ++j)
		A[i][j] = drand((double)1, (double)100);
	for (int i = 0; i < n; ++i)
		B[i] = drand((double)1, (double)100);
	for (int i = 0; i < n - 1; ++i)
		for (int j = i + 1; j < n; ++j){
			double koef = A[j][i] / A[i][i];
			A[j][i] = 0;
			for (int k = i + 1; k < n; ++k)
				A[j][k] -= A[i][k] * koef;
			B[j] = B[i] * koef;
		}

	for (int i = 0; i < n; ++i){
		for (int j = 0; j < n; ++j)
		cout << A[i][j] << " ";
		cout << endl;
	}

	for (int i = n - 1; i >= 0; --i){
		Res[i] = B[i];
		for (int j = n - 1; j > i; --j)
			Res[i] -= Res[j] * A[i][j];
		Res[i] /= A[i][i];
	}
	cout << endl;
	for (int i = 0; i < n; ++i)
		cout << Res[i] << " ";

	return 0;
}