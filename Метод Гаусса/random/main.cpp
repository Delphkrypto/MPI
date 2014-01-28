#include <stdlib.h>
#include <stdio.h>
#include <time.h>

double drand(double dmin, double dmax){
	double f = (double)rand() / RAND_MAX;
	return dmin + f * (dmax - dmin);
}
const int n = 6;
int main(){
	srand(time(NULL));
	freopen("input.txt", "wt", stdout);
	printf("%d\n", n);
	for (int i = 0; i < n; ++i){
		for (int j = 0; j < n; ++j)
			printf("%lf ", drand((double)1, (double)100));
		printf("\n");
	}
		for (int i = 0; i < n; ++i)
			printf("%lf ", drand((double)1, (double)100));
	return 0;
}