#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N 5
#define IT_MAX 30

double f(double x){
	return 1/(1+x*x);
}

double iloraz_rozn(double *xm, double* ym, int j){
	double sum=0, il=1;

	for (int i=0;i<=j;i++){
		il=1;
		for(int k=0;k<=j;k++){
			if(k!=i) il*=1/(xm[i] - xm[k]);
		}
		sum+=ym[i]*il;
	}

	return sum;
}

double wielomian(double x, int n, double* xm, double* fm){
	double sum=0,il=1;

	for(int j=0;j<=n;j++){
		il=1;
		for(int i=0;i<j;i++) il*=(x-xm[i]);
		sum+=fm[j]*il;
	}

	return sum;
}

int main(void)
{
	int n=6;
	double xm[n+1], ym[n+1], fm[n+1];
	double xmin = -5, xmax = 5;
	double h = (xmax - xmin)/n;



	xm[0]=-5;
	xm[1]=-2;
	xm[2]=-0.5;
	xm[3]=0;
	xm[4]=0.5;
	xm[5]=2;
	xm[6]=5;

	printf("Ilorazy różnicowe dla węzłów nieregularnych: \n");
	for(int i=0;i<n+1;i++) {
		ym[i] = f(xm[i]);
		fm[i]=iloraz_rozn(xm,ym,i);
		printf("%f\n", fm[i]);
		
	}

	//for(float d=-5;d<5;d+=0.01) printf("%f \t  %f \n",d,wielomian(d,n,xm,fm));
	//współczynniki
	double X[n+1],s=0,k=0,l=0,m=0,nn=0;
	X[n] = fm[n];

	for(int i=0;i<=5;i++)s+=xm[i];
	X[n-1]= fm[n-1] - fm[n]*s;

	s=0;
	for(int i=0;i<=4;i++){
		for(int j=i+1;j<5;j++)k+=xm[j];
		s+= fm[n-1] - fm[n]*k;
		k=0;
	}
	X[n-2]= fm[n-2] - fm[n-1]*s;


	s=0;k=0;l=0;
	for(int i=0;i<=3;i++){
		for(int j=i+1;j<=4;j++){
			for(int p=j+1;p<=5;p++)l+=xm[p];
			k+=xm[j]*(fm[n-1]-fm[n]*l);
			l=0;
		}

		s+=xm[i] * (fm[n-2]-k);
		k=0;
	}
	X[n-3]= fm[n-3] -s;


	s=0;k=0;l=0,m=0;
	for(int i=0;i<=2;i++){
		for(int j=i+1;j<=3;j++){
			for(int p=j+1;p<=4;p++){
				for(int r=p+1;p<=5;p++)m+=xm[r];
				l+=xm[p]*(fm[n-1]-fm[n]*m);
				m=0;
			}
			k+=xm[j]*(fm[n-2]-fm[n-1]*l);
			l=0;
		}
		s+=xm[i] * (fm[n-3]-k);
		k=0;
	}
	X[n-4]= fm[n-4]-s;





	s=0;k=0;l=0,m=0,nn=0;
	for(int i=0;i<=1;i++){
		for(int j=i+1;j<=2;j++){
			for(int p=j+1;p<=3;p++){
				for(int r=p+1;p<=4;p++){
					for(int t=r+1;t<=5;t++)nn+=xm[p];
					m+=xm[p]*(fm[n-1]-fm[n]*nn);
					nn=0;
				}
				l+=xm[p]*(fm[n-2]-fm[n-1]*m);
				m=0;
			}
			k+=xm[j]*(fm[n-3]-fm[n-3]*l);
			l=0;
		}
		s+=xm[i] * (fm[n-4]-k);
		k=0;
	}
	X[n-5]= fm[n-5]-s;

	X[0] = fm[0] - xm[0]*( fm[1] - xm[1]*(fm[2] - xm[2]*(fm[3] - xm[3]*(fm[4] - xm[4]*(fm[5] - xm[5]*(fm[6]))))));

	//for(int i=0;i<=n;i++)printf("%f \t", X[i]);

	printf("Ilorazy różnicowe dla węzłów regularnych: \n");
	for(int i=0;i<n+1;i++) {
		xm[i] = xmin + h*i;
		ym[i] = f(xm[i]);
		fm[i]=iloraz_rozn(xm,ym,i);
		printf("%f\n", fm[i]);
		

	}

	//for(float d=-5;d<5;d+=0.01) printf("%f \t  %f \n",d,wielomian(d,n,xm,fm));
	for(int i=0;i<=n;i++)X[i]=0;
	s=0;k=0;l=0;m=0;nn=0;
	X[n] = fm[n];

	printf("\n\n");

	for(int i=0;i<=5;i++)s+=xm[i];
	X[n-1]= fm[n-1] - fm[n]*s;

	s=0;
	for(int i=0;i<=4;i++){
		for(int j=i+1;j<5;j++)k+=xm[j];
		s+= fm[n-1] - fm[n]*k;
		k=0;
	}
	X[n-2]= fm[n-2] - fm[n-1]*s;


	s=0;k=0;l=0;
	for(int i=0;i<=3;i++){
		for(int j=i+1;j<=4;j++){
			for(int p=j+1;p<=5;p++)l+=xm[p];
			k+=xm[j]*(fm[n-1]-fm[n]*l);
			l=0;
		}

		s+=xm[i] * (fm[n-2]-k);
		k=0;
	}
	X[n-3]= fm[n-3] -s;


	s=0;k=0;l=0,m=0;
	for(int i=0;i<=2;i++){
		for(int j=i+1;j<=3;j++){
			for(int p=j+1;p<=4;p++){
				for(int r=p+1;p<=5;p++)m+=xm[r];
				l+=xm[p]*(fm[n-1]-fm[n]*m);
				m=0;
			}
			k+=xm[j]*(fm[n-2]-fm[n-1]*l);
			l=0;
		}
		s+=xm[i] * (fm[n-3]-k);
		k=0;
	}
	X[n-4]= fm[n-4]-s;





	s=0;k=0;l=0,m=0,nn=0;
	for(int i=0;i<=1;i++){
		for(int j=i+1;j<=2;j++){
			for(int p=j+1;p<=3;p++){
				for(int r=p+1;p<=4;p++){
					for(int t=r+1;t<=5;t++)nn+=xm[p];
					m+=xm[p]*(fm[n-1]-fm[n]*nn);
					nn=0;
				}
				l+=xm[p]*(fm[n-2]-fm[n-1]*m);
				m=0;
			}
			k+=xm[j]*(fm[n-3]-fm[n-3]*l);
			l=0;
		}
		s+=xm[i] * (fm[n-4]-k);
		k=0;
	}
	X[n-5]= fm[n-5]-s;

	X[0] = fm[0] - xm[0]*( fm[1] - xm[1]*(fm[2] - xm[2]*(fm[3] - xm[3]*(fm[4] - xm[4]*(fm[5] - xm[5]*(fm[6]))))));

	//for(int i=0;i<=n;i++)printf("%f \t", X[i]);

	return 0;

}