#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N 7

void print(double A[][7]){
		for (int i = 0; i < N; i++){
			for(int j=0;j<N;j++)
				printf("%g\t", A[i][j]);
		printf("\n");
	}
}

void multiply(double A[][N], double *b, double *y){
	for(int i=0;i<N;i++){
		y[i]=0;
		for(int j=0;j<N;j++)y[i]+=A[i][j]*b[j];
		}
}

void matrix_multiply(double A[][N], double X[][N], double T[][N]){
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			T[i][j]=0;
			for(int k=0;k<N;k++) T[i][j]+=A[i][k]*X[k][j];
		}
	}
}

void transpose(double Xt[][N], double X[][N]){
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++) Xt[j][i] = X[i][j];
	}
}

double scalar_multiply(double *x, double *y){
	double result=0;
	for(int i=0;i<N;i++) result+= (x[i]*y[i]);
	return result;
}

void add_vect(double*x, double *y, double a, double* pom){
	for(int i=0;i<N;i++) pom[i]= (x[i]+a*y[i]);
	
}


int main(void)
{
	double A[N][N],x[N],x1[N],X[N][N],temp[N], T[N][N], Xt[N][N],D[N][N];
	double a,lam,n;
	
	for (int i=0;i<N;i++){
		for (int j=0;j<N;j++){
			A[i][j] = 1/sqrt(2+abs(i-j));
		}
	}
	// printf("Macierz A: \n");	
	// print(A);

	// printf("Wartości lambda w pierwszej i drugiej iteracji: \n");	
	for (int k=0;k<N;k++){
		for(int i=0;i<N;i++)x[i]=1;
	
		for (int i=1;i<=12;i++){
			multiply(A,x,x1);

			for (int j=0;j<k;j++){
				for(int p=0;p<N;p++)temp[p]=X[p][j];
				a=scalar_multiply(x1,temp);
				add_vect(x1,temp,-a,x1);
			}
			lam = scalar_multiply(x1,x)/scalar_multiply(x,x);
			//printf("%d	%.7f \n",i,lam);				
			n=sqrt(scalar_multiply(x1,x1));

			for (int i=0;i<N;i++) x[i] = x1[i]/n;
			for (int i=0;i<N;i++) X[i][k] = x[i];
		}
		//printf("\n\n");
	}
	// printf("Macierz X zawierająca wektory własne: \n");	
	// print(X); 			
	
	
	transpose(Xt,X);
	matrix_multiply(Xt,A,T);
	matrix_multiply(T,X,D);
	// printf("Macierz D: \n");	
	print(D);

	return 0;

}