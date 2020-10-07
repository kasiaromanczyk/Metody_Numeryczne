#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_eigen.h>

#define N 1
#define L 10
#define n 200


int deltaK(int i, int j){
	if (i==j) return 1;
	else return 0;
}

double ro(double x, int alfa){

	return 1 + 4*alfa*x*x;
}

void print(gsl_matrix *A){
		for (int i = 0; i < n; i++){
			for(int j=0; j< n; j++)
				if (i>190 && j>190)printf("%.2f\t", gsl_matrix_get(A, i, j));
		if (i>190)printf("\n");
	}
}


void function(int alfa,gsl_matrix *A, gsl_matrix *B){
	double dX =(double) L/(n+1);

	for (int i=0;i<n;i++){
		for (int j=0;j<n;j++){
			gsl_matrix_set(A, i, j, (-deltaK(i,j+1) + 2*deltaK(i,j) -deltaK(i,j-1))/(dX*dX)); 
			gsl_matrix_set(B, i, j, (ro((-L/2) + dX*(i+1),alfa)*deltaK(i,j)/N)); 
		}
	}
	gsl_vector *eval = gsl_vector_calloc(n);
	gsl_matrix *evec = gsl_matrix_calloc(n,n);
	gsl_eigen_gensymmv_workspace *w = gsl_eigen_gensymmv_alloc(n);

	gsl_eigen_gensymmv(A, B, eval, evec, w);
	gsl_eigen_gensymmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);


	if(alfa==0 ||alfa==100) {
		for(int j=0;j<n;j++){
			printf("%f ",(float) (-L/2) + dX*(j+1));
			for(int i=0;i<6;i++) printf("%f ", gsl_matrix_get(evec,j,i));
			printf("\n");
		}
		printf("\n\n");
	}

	gsl_matrix_free(evec);
	gsl_vector_free(eval);
	gsl_eigen_gensymmv_free(w);
}

int main(void)
{

	gsl_matrix *A = gsl_matrix_calloc(n, n);
	gsl_matrix *B = gsl_matrix_calloc(n, n);


	for(int alfa=0;alfa<=100;alfa+=2){
		function(alfa, A, B);		
	}


	gsl_matrix_free(A);
	gsl_matrix_free(B);
	
	
	return 0;

}
