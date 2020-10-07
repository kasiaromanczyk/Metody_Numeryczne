#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#include"gammln.c"
#include"nrutil.h"
#include"nrutil.c"
#include"gammp.c"
#include"gcf.c"
#include"gser.c"


#define N 10000

unsigned long RJ(unsigned long x, int a, int c, unsigned long m){
    return (a*x+c)%m;
}

double gestosc(double x, double mi, double sig){

    return (1/(sig*sqrt(2*M_PI)))*exp(-(pow(x-mi,2)/(2*pow(sig,2))));
}

double F(double x, double mi, double sig){

    return (1-erf((x-mi)/(sqrt(2)*sig)))/2;
}

int main(void){

    int a=123, c=1;
    unsigned long Xp=10, Xi, m=pow(2,15);
    double xp, xi, sum=0,mi;
    double tab[N];

    tab[0]=(RJ(10, a, c, m))/(m+1.);
    sum=tab[0];
    for(int i=1;i<N;i++) {
        Xi=RJ(Xp, a, c, m);
        //tab[i-1]=Xp/(m+1.);
        tab[i]=Xi/(m+1.) ;
        //printf("%f\t%f\n",tab[i-1],tab[i]);       //wypisanie wartości xi-1,xi
        sum+=tab[i];
        Xp=Xi;
    }

    mi=(double)sum/(double)N; 
    sum=0;

    for(int i=0;i<N;i++) {
        sum+=pow(tab[i]-mi,2);
    }

    double sig=sqrt(sum/N);

    printf("\n%f\t%f\n\n",mi,sig);
    //wypisanie wart, średniej i odchylenia standardowego

    int k=12;
    double xmin=0,xmax=1;
    double d=(xmax-xmin)/(double)k;

    int traf[k],j;

     for(int i=0;i<k;i++)traf[i]=0;
     for(int i=0;i<N;i++) {
        j=(int)((tab[i]-xmin)/d);
        traf[j]=traf[j]+1;
    }
    
    //for(int i=0;i<k;i++)printf("%f\t%f\n",xmin+d*(i+0.5), (double)traf[i]/N);

////////////////////////////////////////////////////////////////////////////////////////
//rozkład normalny 

    int dd=1;           //ograniczenie od góry
    double mi0=0.2, sig0=0.5;       //teorteycznie wartości mi,sig
    xmin=mi0-3*sig0;                //wyznaczenie krańców przedziału
    xmax=mi0+3*sig0;
    double wid=xmax-xmin;           
    double U1,u1,U2,u2;
    double Up1=10,Up2=10;
    double p[k];
    sum=0;
    int n=0;

    while(n<N){
        U1=RJ(Up1, a, c, m);
        u1=(U1/(m+1.))*wid+xmin;
        U2=RJ(U1, a, c, m);
        u2=(U2/(m+1.));
        if(u2<=gestosc(u1,mi0,sig0)){
            tab[n]=u1;
            sum+=tab[n];
            n++;}
        
        // printf("%f\t%f\n",tab[i-1],tab[i]);
        Up1=U2;
    }

    mi=(double)sum/(double)N; 
    sum=0;

    for(int i=0;i<N;i++) {
        //printf("%f\n",tab[i]);
        sum+=pow(tab[i]-mi,2);
    }

    double war=(sum/N);
    sig=sqrt(war);

    //printf("\n%f\t%f\t%f\n\n",mi,war,sig);
    //wymisanie wart. średniej , wariancji i odchylenia standardowego   


    //do histogramu:

    d=(xmax-xmin)/(double)k;

    for(int i=0;i<k;i++)traf[i]=0;
     for(int i=0;i<N;i++) {
        j=(int)((tab[i]-xmin)/d);
        traf[j]=traf[j]+1;
    }
    
    double S=0;

    for(int i=0;i<k;i++){
        p[i]=F(xmin+i*d,mi0,sig0)-F(xmin+(i+1)*d,mi0,sig0);
        //printf("%f\t%f\n",xmin+d*(i+0.5), (double)traf[i]/N);
        S+=pow(traf[i]-N*p[i],2)/(N*p[i]);

    }

    printf("%f\n",S);
    double poz=gammp((k-2-1)/2,S/2);
    double alfa = 1 - poz;

    //printf("%f\t%f\t%f\n",S,poz,alfa);
    return 0;
}
