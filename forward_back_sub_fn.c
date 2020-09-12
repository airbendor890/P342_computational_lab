#include<stdio.h>
#include"utility.h"

void forward_back_sub(int n,double a[n][n],int *P,double * b){
	/* first make the permutation matrix from P*/
	int perm[n][n];
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			if(j==P[i])
				perm[i][j]=1;
			else
				perm[i][j]=0;

		}
	}

	/* now permute the vector b*/
	double c[n];
	for(int i=0;i<n;i++){
		double sum=0.0;
		for(int j=0;j<n;j++)
			sum+=perm[i][j]*b[j];
		c[i]=sum;
	}
	for(int i=0;i<n;i++)b[i]=c[i];

	
	double sum=0.0;
	/*forward substitution Ly=b */
	for(int i=0;i<n;i++){
		sum=b[i];
		for(int j=0;j<i;j++){
			sum=sum-a[i][j]*b[j];
		}
		b[i]=sum;
	}

	/* back subtitution Ux=y */

	for(int i=n-1;i>=0;i--){
		sum=b[i];
		for(int j=i+1;j<n;j++){
			sum=sum-a[i][j]*b[j];
		}

		b[i]=sum/a[i][i];
	}

}
