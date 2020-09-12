#include<stdio.h>
#include<stdlib.h>
#include"utility.h"

double lu_dcmpsn(int n,double a[n][n],int *P){
	/* P will hold the permutation infromation that will occur during 
	partial pivoting*/

	double sum=0.0;
	double det=1.0;
	for(int j=0;j<n;j++){	/*loop over columns of crouts method*/
		for(int i=0;i<=j;i++){
			sum=a[i][j];
			for(int k=0;k<i;k++) sum=sum-a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		
		if(abs(a[j][j])<=1.0e-12) lu_partial_pivot(n,a,P,j);

		for(int i=j+1;i<n;i++){
			sum=a[i][j];
			for(int k=0;k<j;k++) sum=sum-a[i][k]*a[k][j];
			//still if diag element is 0 its a singular matrix
			//return 0  
			if(a[j][j]==0.0){	
				printf("singular matrix !!");
				return 0;
			}
			else
			a[i][j]=sum/a[j][j];
		}

	}

	//its not singular ..return its determinant
	for(int i=0;i<n;i++)det=det*a[i][i];

	return det;

}	
