#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N 5
#define IT_MAX 30

double licz_R(double *a, double *b, int n, double x){
	double R2;
	b[n]=0;
	for (int i=n-1;i>=0;i--){
		b[i]=a[i+1] +x*b[i+1];

	}
	R2 = a[0]+ x*b[0];
	return R2;
}//funkcja liczaca Rj

int main(void)
{
	
	double a[N+1],b[N+1];
	double x0,x1,x2,R0,R1,R2;
	int n;

	a[0]=240.0;
	a[1]=-196.;
	a[2]=-92.;
	a[3]=33.;
	a[4]=14.;
	a[5]=1.;
	//inicjalizacja wektora współczynnników

	for (int i=1;i<=N;i++){
		n=N-i+1;
		x0=0.;
		x1=0.1;
		
		R0 = licz_R(a,b,n,x0);
		R1 = licz_R(a,b,n,x1);
		//printf("%f\t %f\t \n", R0,R1);

		for (int j=1;j<IT_MAX;j++){
			x2=x1-R1*(x1-x0)/(R1-R0);

			R2 = licz_R(a,b,n,x2);
			R0=R1;
			R1=R2;
			x0=x1;
			x1=x2;

			printf("%d\t %d\t %g\t %g\t \n", i,j,x2,R2);	//wypisanie potrzebnych danych, indeks zewnetrznej i wewnetrznej petli, x2 oraz R2

			if(fabs(x1-x0)<1.0e-7) break;
		}
		printf("\n");
		for(int k=0; k<=(n-1); k++)a[k]=b[k];
	}
	


	return 0;

}