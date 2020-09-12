#include<stdio.h>
#include<stdlib.h>
#include"utility.h"



//-----------------------------------------------------------------------

void lu_partial_pivot(int n,double arr[n][n],int* P,int r){
		int r1=r+1;
		while(abs(arr[r][r])<1.0e-12 && r1<n){
			if(abs(arr[r1][r])>abs(arr[r][r])){
				double temp;
				for(int c=0;c<n;c++){
					temp=arr[r1][c];
					arr[r1][c]=arr[r][c];
					arr[r][c]=temp;
				}

				int dum;
				dum=P[r1];
				P[r1]=P[r];
				P[r]=dum;
			}
			r1++;
		}



}
//-----------------------------------------------------------------------

