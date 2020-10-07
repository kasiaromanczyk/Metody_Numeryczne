#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define eps1 1.e-2
#define eps2 1.e-3
#define delta 1.e-4 

double f( double x,  double y){
    return (5./2.)*pow((x*x - y),2) +pow((1-x),2);
}

double dx(double x, double y){

    return (f(x+delta,y)-f(x-delta,y))/(2*delta);
}


double dy(double x, double y){

    return (f(x,y+delta)-f(x,y-delta))/(2*delta);
}


int main(){


   double x0,y0,x1,y1;
   double h=0.1;


	
    x0=-0.75;
    y0=1.75;

    for(int i=1;i<=1000;i++) {
        x1=x0-h*dx(x0,y0);
        y1=y0-h*dy(x0,y0);
        //if(i==1)printf("%g,\t%g\n",dx(x0,y0),dy(x0,y0));

        printf("%d\t%g\t%g\n",i,x1,y1);

        if(sqrt(pow(x1-x0,2)+pow(y1-y0,2))<eps2) break;

        x0=x1;
        y0=y1;
        h-=0.0001;

    }


    return 0;
}