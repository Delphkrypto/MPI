#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <time.h>
#include <iostream>
#include <mpi.h>
#include <memory.h>
#pragma comment (lib, "mpi.lib")
#pragma comment (lib, "mpe.lib")
#pragma comment (lib, "fmpich2.lib")
#pragma comment (lib, "fmpich2g.lib")
#pragma comment (lib, "cxx.lib")
#pragma comment (lib, "irlog2rlog.lib")
#pragma comment (lib, "rlog.lib")
#pragma comment (lib, "TraceInput.lib")

using namespace std;

MPI_Status status;
//const int n = 6;
double *A, *C, *B, *D, *At;
double *Res, tmp;
int n, Bt, FLAG;
int ProcNum = -1;
int ProcRank = -1;

int main(int argc, char **argv){
	srand(time(NULL));
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	if(!ProcRank){
		freopen("input.txt", "rt", stdin);
		scanf("%d", &n);
		for (int i = 1; i < ProcNum; ++i)
			MPI_Send(&n, 1, MPI_INT,i, 0, MPI_COMM_WORLD);
	}
	else
		MPI_Recv(&n, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	
	A = (double *)malloc(n * n * sizeof(double));
	C = (double *)malloc(n * n * sizeof(double));
	B = (double *)malloc(n * sizeof(double));
	D = (double *)malloc(n * sizeof(double));
	At = (double *)malloc(n * sizeof(double));
	Res = (double *)malloc(n * sizeof(double));

	if(!ProcRank){
		for (int i = 0; i < n; ++i)
			if(!(i % ProcNum))
				for (int j = 0; j < n; ++j)
					scanf("%lf", &A[(i / ProcNum) * n + j]);
			else{
				for (int j = 0; j < n; ++j)
					scanf("%lf", &At[j]);
				MPI_Send(At, n, MPI_DOUBLE, i % ProcNum, 0, MPI_COMM_WORLD);
			}
		for (int i = 0; i < n; ++i)
			if(!(i % ProcNum))
				scanf("%lf", &B[i / ProcNum]);
			else{
				scanf("%lf", &Bt);
				MPI_Send(&Bt, 1, MPI_DOUBLE, i % ProcNum, 0, MPI_COMM_WORLD);
			}
	}
	else{
		for (int i = 0; i < n / ProcNum; ++i){
			MPI_Recv(At, n*n, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			memcpy(A + (i * n) * sizeof(double), At, n * sizeof(double));
		}
		for (int i = 0; i < n / ProcNum; ++i){
			MPI_Recv(&B[i], n, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	///////////////////
	//if(!ProcRank){
	//for (int i = 0; i < n; ++i){
	//	for (int j = 0; j < n; ++j)
	//	cout << A[i * n + j] << " \t";
	//	cout << endl;
	//}
	//cout << "!" << endl;
	//for (int j = 0; j < n; ++j)
	//	cout << B[j] << " \t";
	//cout << endl;
	//cout << "!!!" << endl;
	//}
	//MPI_Barrier(MPI_COMM_WORLD);
	//if(ProcRank == 1){
	//for (int i = 0; i < n; ++i){
	//	for (int j = 0; j < n; ++j)
	//	cout << C[i * n + j] << " \t";
	//	cout << endl;
	//}
	//cout << endl;
	//cout << "!" << endl;
	//for (int j = 0; j < n; ++j)
	//	cout << D[j] << " \t";
	//cout << endl;
	//cout << "!!!" << endl;
	//}
	//MPI_Barrier(MPI_COMM_WORLD);
	//if(ProcRank == 2){
	//for (int i = 0; i < n; ++i){
	//	for (int j = 0; j < n; ++j)
	//	cout << C[i * n + j] << " \t";
	//	cout << endl;
	//}
	//cout << endl;
	//cout << "!" << endl;
	//for (int j = 0; j < n; ++j)
	//	cout << D[j] << " \t";
	//cout << endl;
	//cout << "!!!" << endl;
	//}
	///////////////////////

	while(FLAG < n)
	{
		if(ProcRank == FLAG % ProcNum){
			for (int i = 0; i < n; ++i)
				At[i] = A[(FLAG / ProcNum) * n + i]; 
			for (int i = 0; i < ProcNum; ++i)
				if(i != ProcRank)
				{
					MPI_Send(At, n, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
					MPI_Send(&B[FLAG / ProcNum], 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
				}
			int koef = 0;
			for (int i = (FLAG / ProcNum) + 1; i < n / ProcNum; ++i){
				koef = A[i * n] / A[(FLAG / ProcNum) * n]; 
				for (int j = 0; j < n; ++j)
					A[i * n + j] -= A[(FLAG / ProcNum) * n + j] * koef;
				B[i] -= B[FLAG / ProcNum] * koef; 
			}
		}
		else{
			MPI_Recv(At, n, MPI_DOUBLE, FLAG % ProcNum, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			MPI_Recv(&Bt, 1, MPI_DOUBLE, FLAG % ProcNum, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			int koef = 0;
			for (int i = (FLAG / ProcNum) + 1; i < n / ProcNum; ++i){
				koef = A[i * n] / At[0]; 
				for (int j = 0; j < n; ++j)
					A[i * n + j] -= A[(FLAG / ProcNum) * n + j] * koef; 
				B[i] -= Bt * koef; 
			}
		}

	}

	//for (int i = 0; i < n - 1; ++i){
	//	if(ProcRank){
	//		for (int j = (i + 1) / (ProcNum - 1); j < n / (ProcNum - 1); ++j)
	//			if ((j * (ProcNum - 1) + (ProcRank - 1)) > i){
	//				int t = j * (ProcNum - 1) + (ProcRank - 1);
	//				double koef = C[t * n + i] / C[i * n + i];
	//				C[t * n + i] = 0;
	//				for (int k = i + 1; k < n; ++k)
	//					C[j * n + k] -= C[i * n + k] * koef;
	//				D[j] -= D[i] * koef;
	//				//MPI_Send(&t, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	//				//MPI_Send(C, n*n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	//				//MPI_Send(D, n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	//				//cout << "send number " << t << " from " << ProcRank << endl;
	//			}
	//		MPI_Send(C, n*n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	//		MPI_Send(D, n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	//		if(ProcRank == 2){
	//			for (int i = 0; i < n; ++i){
	//				for (int j = 0; j < n; ++j)
	//					cout << C[i * n + j] << " \t";
	//				cout << endl;
	//			}
	//			cout << endl << endl;
	//		}
	//		//cout << "send number " << t << " from " << ProcRank << endl;
	//	}
	//	else{
	//		for (int j = 1; j < ProcNum; ++j){
	//			MPI_Recv(C, n * n, MPI_DOUBLE, j, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	//			MPI_Recv(D, n, MPI_DOUBLE, j, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	//			for (int k = 0; k < n / (ProcNum - 1); ++k){
	//				memcpy(A + (((k * (ProcNum - 1) + j - 1) * n) * sizeof(double)) , C + (((k * (ProcNum - 1) + j - 1) * n) * sizeof(double)), n * sizeof(double));
	//				B[k * (ProcNum - 1) + j - 1] = D[k * (ProcNum - 1) + j - 1];
	//			}
	//		}
	//		
	//	}

	//	MPI_Barrier(MPI_COMM_WORLD);

	//	if(!ProcRank)
	//		for (int i = 1; i < ProcNum; ++i){
	//			MPI_Send(A, n*n, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
	//			MPI_Send(B, n, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
	//		}
	//	else{
	//		MPI_Recv(C, n*n, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	//		MPI_Recv(D, n, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	//	}
	//}

	if (!ProcRank){
		for (int i = 0; i < n; ++i){
			for (int j = 0; j < n; ++j)
			cout << A[i * n + j] << " ";
			cout << endl;
		}

		for (int i = n - 1; i >= 0; --i){
			Res[i] = B[i];
			for (int j = n - 1; j > i; --j)
				Res[i] -= Res[j] * A[i * n + j];
			Res[i] /= A[i * n + i];
		}
		cout << endl;
		for (int i = 0; i < n; ++i)
			cout << Res[i] << " ";
	}
	
	MPI_Finalize();
	return 0;
}