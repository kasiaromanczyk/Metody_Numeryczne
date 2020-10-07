#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N 100 //liczba wezlów
#define Ms 5
#define Mc 5
#define Ms2 30
#define Mc2 30

 double f1( double x,  double alfa){
    return 2*sin(x) + sin(2*x) + 2*sin(3*x) +alfa;
}

 double f2( double x,  double alfa){
    return 2*sin(x) + sin(2*x) + 2*cos(x) +cos(2*x);
}

 double f3( double x,  double alfa){
    return 2*sin(1.1*x) + sin(2.1*x) + 2*sin(3.1*x) +alfa;
}

 double f4( double x,  double alfa){
    return 2*sin(x) + sin(2*x) + 2*sin(3*x) +alfa;
}

double approx(double x, double* ak, double* bj){
    double sum1=0,sum2=0;
    for(int k=1;k<Ms+1;k++) {
            sum1+=ak[k]*sin(k*x);
        }
    
    for(int k=1;k<Mc+1;k++) {
            sum1+=bj[k]*cos(k*x);
        }
    
    return sum1+sum2+bj[0]/2;
}

double approx2(double x, double* ak, double* bj){
    double sum1=0,sum2=0;
    for(int k=0;k<Ms2+1;k++) {
            sum1+=ak[k]*sin(k*x);
        }
    
    for(int k=0;k<Mc2+1;k++) {
            sum1+=bj[k]*cos(k*x);
        }
    
    return sum1+sum2+bj[0]/2;
}

int main(){


   double wx[N], wy[N], ak[Ms+1], bj[Mc+1],ak2[Ms2+1], bj2[Mc2+1];
    double xmin=0, xmax = 2*M_PI;
    double sum=0;
     double h = (xmax - xmin) / (N-1);
	
    //double alpha=rand()/(RAND_MAX+1.)-0.5;

   
    for(int i=0;i<N;i++) {
        wx[i]=xmin+h*i;
        //if (i<5 || i>95)printf("%g\n", wx[i]);
    }
    //printf("Wektory wx, wy:\n");
    for(int i=0;i<N;i++) {
        wy[i]=f1(wx[i],0);
        //printf("%.2f\t%.2f\n",wx[i],wy[i]);
    }
    
    for(int i=0;i<Ms+1;i++) {
        sum=0;
        for(int j=0;j<N;j++) {
            sum+=wy[j]*sin(i*wx[j]);
        }
        ak[i]=2*sum/N;
        //if (i<5 || i>95)printf("%g\n",ak[i]);
    }

    for(int i=0;i<Mc+1;i++) {
        sum=0;
        for(int j=0;j<N;j++) {
            sum+=wy[j]*cos(i*wx[j]);
        }
        bj[i]=2*sum/N;
        //if (i<5 || i>95)printf("%g\n",bj[i]);
    }
    
    
    //printf("Współczynniki a,b (5):\n");
    //for(int i=0;i<Mc+1;i++) printf("%d\t%.3f\t%.3f\n",i, ak[i],bj[i]);
    //printf("Aproksymacja: (5)\n");
   for(double i=0;i<2*M_PI;i+=0.01){
       printf("%f\t%.4f\n",i,approx(i, ak, bj));
   }

   //Cześć dla Ms=Mc=30 i f4 (dla tej samej alfy)
   printf("\n\n");

       for(int i=0;i<Ms2+1;i++) {
        sum=0;
        for(int j=0;j<N;j++) {
            sum+=wy[j]*sin(i*wx[j]);
        }
        ak2[i]=2*sum/N;
        //if (i<5 || i>95)printf("%g\n",ak[i]);
    }

    for(int i=0;i<Mc2+1;i++) {
        sum=0;
        for(int j=0;j<N;j++) {
            sum+=wy[j]*cos(i*wx[j]);
        }
        bj2[i]=2*sum/N;
        //if (i<5 || i>95)printf("%g\n",bj[i]);
    }
    
    printf("Współczynniki a,b (30):\n");
    for(int i=0;i<Mc2+1;i++) printf("%d\t%.3f\t%.3f\n",i, ak2[i],bj2[i]);
    printf("Aproksymacja (30):\n");
   for(double i=0;i<2*M_PI;i+=0.01){
       printf("%f\t%.4f\n",i,approx2(i, ak2, bj2));
   }


    return 0;
}
