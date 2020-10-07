#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "nrutil.h"
#include "nrutil.c"
#include "gaussj.c" 


#define N 1000 // rozmiar macierzy M: NxN
#define M  5// m przyjete w tresci zadania

void print(float **A){
		for (int i = 1; i <= N; ++i){
			for(int j=1;j<=N;j++)
				printf("%.4f\t", A[i][j]);
		printf("\n");
	}
}

int abs(int x){
	if (x>0) return x;
	else return -x;
}

int max(int x,int y){
	if (x>y) return x;
	else return y;
}
int min(int x,int y){
	if (x>y) return y;
	else return x;
}

void multiply1(float **A, float **b, float **y){
	for(int i=1;i<=N;i++){
		int jmin=max(1,i-M);
		int jmax=min(i+M,N);
		y[i][1]=0;
		for(int j=jmin;j<=jmax;j++)y[i][1]+=A[i][j]*b[j][1];
		}
}

float scalar_multiply(float **x, float **y){
	float result=0;
	for(int i=1;i<=N;i++) result+= (x[i][1]*y[i][1]);
	return result;
}

void add_vect(float **x, float **y, float a, float ** pom){

	for(int i=1;i<=N;i++) pom[i][1]= (x[i][1]+a*y[i][1]);
	
}


int main(void)
{
	float **A, **b, **x, **r, **v, **temp, **Av;
	float alfa, beta, s;

	//	Alokacja macierzy i wektorów
	A = matrix(1, N, 1, N);
	b = matrix(1,N,1,1);
    x = matrix(1,N,1,1);
    r = matrix(1,N,1,1);
	v = matrix(1,N,1,1);
	temp = matrix(1,N,1,1);
	Av = matrix(1,N,1,1);

	// 	Wypelnienie macierzy A
	for (int i = 1; i <= N; ++i)
	{
		for (int j = 1; j <= N; ++j){
			if (abs(i-j)<=M) A[i][j] = (1.0/(1 + abs(i-j)));
			else  A[i][j] = 0;
		}
	}


//printf("Macierz A: \n"); print(A);printf("\n");

//ustawienie wektorów wyrazów wolnych
	for (int i = 1; i <= N; ++i) {
		b[i][1]=(float)i;
		x[i][1]=0;
		//temp[i][1]=0;
		//Av[i][1]=0;	
    }

	clock_t t;
	t = clock();

	gaussj(A,N,b,1);

	t = clock() - t;
	printf("czas wykonywania metoda eliminacji zupelnej: %li, w sekundach: %f \n",t, (float)t/CLOCKS_PER_SEC);

	for (int i = 1; i <= N; ++i)
	{
		for (int j = 1; j <= N; ++j){
			if (abs(i-j)<=M) A[i][j] = (1.0/(1 + abs(i-j)));
			else  A[i][j] = 0;
		}
	}


//printf("Macierz A: \n"); print(A);printf("\n");

//ustawienie wektorów wyrazów wolnych
	for (int i = 1; i <= N; ++i) {
		b[i][1]=(float)i;
		x[i][1]=0;
		//temp[i][1]=0;
		//Av[i][1]=0;	
    }

	//for (int i = 1; i <= N; ++i) printf("%.2f \n", b[i][1]);
	//for (int i = 1; i <= N; ++i) printf("%.2f \n", x[i][1]);

 	multiply1(A,x,temp);
	//inicjalizacja wektorów r i v
	for (int i = 1; i <= N; ++i) {
		r[i][1]=b[i][1]-temp[i][1];
		v[i][1]=r[i][1];
		//printf("%.2f \n", temp[i][1]);
		//printf("%.2f \n", r[i][1]);
	}
	//for (int i = 1; i <= N; ++i) printf("%.2f \n", r[i][1]);
	//for (int i = 1; i <= N; ++i) printf("%.2f \n", v[i][1]);
	int k=1;
	//printf("%f \n", sqrt(scalar_multiply(r,r))); 
	
	//
	
	t = clock();

	while(sqrt(scalar_multiply(r,r)) > 1e-6){
		s=scalar_multiply(r,r);
		multiply1(A,v,Av);
		alfa = s/scalar_multiply(v,Av);

		add_vect(x,v,alfa,temp);
		for (int i = 1; i <= N; ++i) x[i][1] = temp[i][1];
		
		add_vect(r,Av,-alfa,temp);
		for (int i = 1; i <= N; ++i) r[i][1] = temp[i][1];

		beta = scalar_multiply(r,r)/s;

		add_vect(r,v,beta,temp);
		for (int i = 1; i <= N; ++i) v[i][1] = temp[i][1];
		
		//printf("%d\t%f\t%f\t%f\t%f\n",k,sqrt(scalar_multiply(r,r)),alfa,beta,sqrt(scalar_multiply(x,x)));
		k++;
	}

	t = clock() - t;
	printf("czas wykonywania metoda wstegowa: %li, w sekundach: %f",t, (float)t/CLOCKS_PER_SEC);

	//	Zwolnienie pamieci
	free_matrix(A, 1, N, 1, N);
	free_matrix(b, 1, N, 1, 1);
	free_matrix(r, 1, N, 1, 1);
	free_matrix(v, 1, N, 1, 1);
	free_matrix(x,1,N,1,1);
	free_matrix(temp, 1, N, 1, 1);
	free_matrix(Av, 1, N, 1, 1);


	return 0;

}
