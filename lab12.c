#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 7 //wymiar tablicy

double f(double x){
    if(fabs(x)<=10e-7)x+=10e-7;
    return sin(x)/x;
}

double f2(double x){
    if(fabs(x)<=10e-7)x+=10e-7;
    return (cos(x)-exp(x))/sin(x);
}

double f3(double x){
    if(fabs(x)<=10e-7)x+=10e-7;
    return 1./(x*exp(1/x));
}

int main(){

    double a=0,b=1;     //przedział calkowania


    double D[N][N];
    double h,s=0;       //h- krok całkowania; s - do sumowania w pierwszej pętli
    for(int i=0;i<N;i++)
        for(int j=0;j<N;j++) D[i][j]=0; //początkowe wypełnienie tablicy zerami

    D[0][0]=(1./2.)*(b-a)*(f3(a)+f3(b));
    //ustalenie pierwszego elementu tablicy
    for(int i=1;i<N;i++){
        h=(b-a)/(pow(2,i));
        s=0;
        for(int j=1;j<=pow(2,i-1);j++) s+=f3(a+(2*j-1)*h);

        D[i][0]=(1./2.)*D[i-1][0]+h*s;
    }//wypełnienie pierwszej kolumny tabeli

    for(int i=1;i<N;i++){
        for(int j=1;j<=i;j++)
            D[i][j]=(pow(4,j)*D[i][j-1]-D[i-1][j-1])/(pow(4,j)-1);
    }//wypełnienie koijnych kolumn

    for(int i=0;i<N;i++){
        for(int j=0;j<=i;j++) printf("%g  ",D[i][j]);
        printf("\n");
    }//wypisanie interesujących nas wyników

    return 0;
}