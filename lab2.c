#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "nrutil.h"
#include "nrutil.c" 
#include "lubksb.c"
#include "ludcmp.c"


#define N 3 // rozmiar macierzy M: NxN

float **lu(float **A, int *indx, float *d){
	ludcmp(A,N,indx,d);
	return A;
}

void print(float **A){
		for (int i = 1; i <= N; ++i){
			for(int j=1;j<=N;j++)
				printf("%g\t", A[i][j]);
		printf("\n");
	}
}

float max(float** A){
	float result=fabs(A[1][1]);
		for (int i = 1; i <= N; ++i){
			for(int j=1;j<=N;j++){
				if(result<=fabs(A[i][j]))result= fabs(A[i][j]);
			}
	}
	return result;
}

int main(void)
{
	float **A, **B, **A1, **B1;
	int * indxA, *indxB;
	//	Alokacja macierzy
	A = matrix(1, N, 1, N);
	A1 = matrix(1, N, 1, N);
	B = matrix(1, N, 1, N);
	B1 = matrix(1, N, 1, N);
	indxA = ivector(1, N);
	indxB = ivector(1, N);

	float w=1,dA=1,dB=1;

	// 	Wypelnienie macierzy A,A1,B,B1
	for (int i = 1; i <= N; ++i)
	{
		for (int j = 1; j <= N; ++j){
			A[i][j] = w;
			A1[i][j] =w;
			B[i][j] =w;
			B1[i][j]=w;
			w+=1;
		}
	}

	B[1][1]=B1[1][1]=1.1;

	A = lu(A,indxA,&dA);
	B= lu(B,indxB,&dB);

printf("Macierz LU z A: \n");print(A);printf("\n");
printf("Macierz A: \n");print(A1);printf("\n");
printf("Macierz LU z B: \n") ;print(B);printf("\n");
printf("Macierz B: \n");print(B1);printf("\n");

//inicjalizacja wektorw wyrazow wolnych
float *b1, *b2, *b3;
//ustalenie macierzy odwrotnych od A(LU) i B(LU)
float **AA, **BB;

	b1 = vector(1, N);
	b2 = vector(1, N);
	b3 = vector(1, N);

AA = matrix(1, N, 1, N);
BB = matrix(1, N, 1, N);

//ustawienie wektorów wyrazów wolnych
	for (int i = 1; i <= N; ++i) {
		b1[i]=0;
		b2[i]=0;
		b3[i]=0;
	}

	b1[1]=1;
	b2[2]=1;
	b3[3]=1;

//	Rozwiazanie ukladu rownan Mx=b - wywolanie procedury:
	lubksb(A, N, indxA, b1);
	lubksb(A, N, indxA, b2);
	lubksb(A, N, indxA, b3);

	for (int i = 1; i <= N; ++i)AA[i][1]=b1[i];
	for (int i = 1; i <= N; ++i)AA[i][2]=b2[i];
	for (int i = 1; i <= N; ++i)AA[i][3]=b3[i];
//ponowne ustawienie wektorów wyrazów wolnych
for (int i = 1; i <= N; ++i) {
		b1[i]=0;
		b2[i]=0;
		b3[i]=0;
	}

	b1[1]=1;
	b2[2]=1;
	b3[3]=1;

	lubksb(B, N, indxB, b1);
	lubksb(B, N, indxB, b2);
	lubksb(B, N, indxB, b3);
	for (int i = 1; i <= N; ++i)BB[i][1]=b1[i];
	for (int i = 1; i <= N; ++i)BB[i][2]=b2[i];
	for (int i = 1; i <= N; ++i)BB[i][3]=b3[i];


printf("Macierz odwrotna z A: \n");print (AA);printf("\n");
printf("Macierz odwrotna z B: \n");print(BB);printf("\n");

float wA, wAA, wB, wBB;

	wA = max(A1);
	wAA =  max(AA);
	wB = max(B1);
	wBB =  max(BB);

//wspolczynniki normy 
printf("Wsp normy z A: \n");
printf("%.2f \n",wA);

float kA = wA*wAA;

printf("%.2f \n",wAA);
printf("%.2f \n",kA);
float kB = wB*wBB;
printf("Wsp normy z B: \n");
printf("%.2f \n",wB);
printf("%.2f \n",wBB);
printf("%.2f \n",kB);


//mnozenie maciezy
float **C, **CC;
	C = matrix(1, N, 1, N);
	CC = matrix(1, N, 1, N);

	for (int i = 1; i <= N; ++i){
			for(int j=1;j<=N;j++){
				for(int k=1;k<=N;k++) {
					C[i][j] += A1[i][k] * AA[k][j];
					CC[i][j] += B1[i][k] * BB[k][j];
				}
			}
	}

	printf("Iloczyn dla  A: \n");print(C);
	printf("Iloczyn dla B: \n");print(CC);

	//	Zwolnienie pamieci
	free_matrix(A, 1, N, 1, N);
	free_matrix(A1, 1, N, 1, N);
	free_matrix(B, 1, N, 1, N);
	free_matrix(B1, 1, N, 1, N);
//	free_matrix(indxA, 1, N, 1, N);
//	free_matrix(indxB, 1, N, 1, 1);

	return 0;

}
