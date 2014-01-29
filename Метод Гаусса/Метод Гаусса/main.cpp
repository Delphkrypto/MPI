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
double *A, *C, *B, *D, *At;
double *Res, tmp, Bt, tm;
int n, FLAG;
int ProcNum = -1;
int ProcRank = -1;

int main(int argc, char **argv){
	srand(time(NULL));
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	if(!ProcRank){
		freopen("input.txt", "rt", stdin);
		freopen("output.txt", "wt", stdout);
		scanf("%d", &n);
	}
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

	A = (double *)malloc(n * n * sizeof(double));
	C = (double *)malloc(n * n * sizeof(double));
	B = (double *)malloc(n * sizeof(double));
	D = (double *)malloc(n * sizeof(double));
	At = (double *)malloc(n * sizeof(double));
	Res = (double *)malloc(n * sizeof(double));

	MPI_Barrier(MPI_COMM_WORLD);
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

	}
	else{
		for (int i = 0; i < n / ProcNum; ++i){
			MPI_Recv(At, n, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			memcpy(&A[i * n], At, n * sizeof(double));
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if(!ProcRank)
		for (int i = 0; i < n; ++i)
			if(!(i % ProcNum))
				scanf("%lf", &B[i / ProcNum]);
			else{
				scanf("%lf", &Bt);
				MPI_Send(&Bt, 1, MPI_DOUBLE, i % ProcNum, 0, MPI_COMM_WORLD);
			}
	else
		for (int i = 0; i < n / ProcNum; ++i){
			MPI_Recv(&B[i], n, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		}

	MPI_Barrier(MPI_COMM_WORLD);
	tm = MPI_Wtime();
	FLAG = 0;
	while(FLAG < n)
	{
		if(ProcRank == (FLAG % ProcNum)){
			for (int i = 0; i < n; ++i)
				At[i] = A[(FLAG / ProcNum) * n + i]; 
			Bt = B[FLAG / ProcNum];
			for (int i = 0; i < ProcNum; ++i)
				if(i != ProcRank)
				{
					MPI_Send(At, n, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
					MPI_Send(&Bt, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
				}
			double koef = 0;
			for (int i = (FLAG / ProcNum) + 1; i < n / ProcNum; ++i){
				koef = A[i * n + FLAG] / At[FLAG]; 
				A[i * n + FLAG] = 0;
				for (int j = FLAG + 1; j < n; ++j)
					A[i * n + j] -= At[j] * koef;
				B[i] -= Bt * koef; 
			}
		}
		else{
			MPI_Recv(At, n, MPI_DOUBLE, FLAG % ProcNum, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			MPI_Recv(&Bt, 1, MPI_DOUBLE, FLAG % ProcNum, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			double koef = 0;
			for (int i = (FLAG / ProcNum); i < n / ProcNum; ++i){
				if((i * ProcNum + ProcRank) > FLAG){
					koef = A[i * n + FLAG] / At[FLAG]; 
					A[i * n + FLAG] = 0;
					for (int j = FLAG + 1; j < n; ++j)
						A[i * n + j] -= At[j] * koef; 
					B[i] -= Bt * koef;
				}
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		++FLAG;

	}
	MPI_Barrier(MPI_COMM_WORLD);
	FLAG = n - 1;
	while(FLAG >= 0) {
		if(ProcRank == (FLAG % ProcNum)){
			Res[FLAG] = B[FLAG / ProcNum] / A[(FLAG / ProcNum) * n + FLAG];
		for (int i = 0; i < ProcNum; ++i)
			if(i != ProcRank)
				MPI_Send(&Res[FLAG], 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
		for (int i = FLAG / ProcNum - 1; i >= 0; --i){
				B[i] -= Res[FLAG] * A[n * i + FLAG]; 
			}
		}
		else{
			MPI_Recv(&Res[FLAG], 1, MPI_DOUBLE, FLAG % ProcNum, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			for (int i = FLAG / ProcNum; i >= 0; --i){
				B[i] -= Res[FLAG] * A[n * i + FLAG]; 
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		--FLAG;
	}

	if (!ProcRank)
		for (int i = 0; i < n; ++i)
			cout << Res[i] << " ";

	if(!ProcRank){
		tm = MPI_Wtime() - tm;
		cout << endl << tm;
	}
	MPI_Finalize();
	return 0;
}