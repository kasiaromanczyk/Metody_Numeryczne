#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "nrutil.h"
#include "nrutil.c"
#include "gaussj.c" 

#define N 5 //liczba wezlów

float f1(float x){
    return 1/(1+x*x);
}

float pochodna(float x){

    float dx=0.01,wynik;
    return (f1(x+dx) - f1(x-dx))/(2*dx);

}

float sklej(int i, float x, float h, float* xw){

    float d;

    if(x>=xw[i-2] && x<xw[i-1])  d= (x-xw[i-2])*(x-xw[i-2])*(x-xw[i-2]);
    else if (x>=xw[i-1] && x<xw[i]) d= h*h*h +3*h*h*(x - xw[i-1]) + 3*h*(x - xw[i-1])*(x - xw[i-1]) - 3*(x - xw[i-1])*(x - xw[i-1])*(x - xw[i-1]);
    else if (x>=xw[i] && x<xw[i+1]) d= h*h*h +3*h*h*(xw[i+1] - x) + 3*h*(xw[i+1] - x)*(xw[i+1] - x) - 3*(xw[i+1] - x)*(xw[i+1] - x)*(xw[i+1] - x);
    else if (x>=xw[i+1] && x<xw[i+2]) d=(xw[i+2] - x)*(xw[i+2] - x)*(xw[i+2] - x);
    else return 0;

    return d/(h*h*h);

}

float interp(float x, float h, float* c, float *xw){
    float result =0;
    for(int i=0;i<=N+1;i++)result+=c[i]*sklej(i,x,h,xw);

    return result;
}



int main(){


    float *xw, *yw, **b, **A,*c;
    float xmin=-5, xmax = 5;

    float h = (xmax - xmin) / (N-1);
	

	//	Alokacja macierzy i wektorów
	xw = vector(-2, N+3);
	yw = vector(1,N);
    b = matrix(1, N,1,1);
	A = matrix(1,N,1,N);
    c = vector(0,N+1);

    printf("Wektor xw:\n");
    for(int i=-2;i<=N+3;i++) {
        xw[i]=xmin+h*(i-1);
        //printf("%.2f\n", xw[i]);
    }
    //printf("Wektor yw:\n");
    for(int i=1;i<=N;i++) {
        yw[i]=f1(xw[i]);
        //printf("%.5f\n", yw[i]);

    }
    //printf("Macierz a\n");
    for(int i=1;i<=N;i++) {
        for(int j=1;j<=N;j++){
            A[i][j]=0;
            if (i==j) A[i][j]=4;
            if(abs(i-j)==1) A[i][j]=1;
            A[1][2]=2;
            A[N][N-1]=2;
            printf("%.0f ", A[i][j]);
        }
        printf("\n");
    }

    for(int i=1;i<=N;i++) {
        b[i][1] = yw[i];
        
    }

    float alfa=pochodna(xmin);
    float beta=pochodna(xmax);

    printf("alfa = %f\n", alfa);
    printf("beta = %f\n", beta);

    b[1][1]=yw[1]+h*alfa/3;
    b[N][1]=yw[N]-h*beta/3;

    printf("Wektor b:\n");
    for(int i=1;i<=N;i++) printf("%f\n", b[i][1]);

    gaussj(A,N,b,1);


    printf("\n\nWektor c:\n");

    for(int i=1;i<=N;i++)c[i]=b[i][1];
    c[0]=c[2] - h*alfa/3;
    c[N+1]=c[N-1] + h*beta/3;

    for(int i=0;i<=N+1;i++) printf("%f\n", c[i]);

    for(float i=-5;i<=5;i+=0.01) {
        //printf("%f\t%f\n",i,interp(i,h,c,xw));
    }
    //wypisanie do pliku x oraz s(x)

	free_vector(xw, -2, N+3);
    free_vector(yw, 1, N);
    free_vector(c, 0, N+1);
    free_matrix(b, 1, N,1,N);
    free_matrix(A, 1, N,1,N);


    return 0;
}