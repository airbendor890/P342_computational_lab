#include<stdio.h>
#include<stdlib.h>
//gauss jordan method

//partial pivot
void partial_pivot(int n,double [n][n],double* b,int r);
//(pointer to matrix A,pointer to b,r,n)


//gauss jordan method

void gauss_jordan(int n,double [n][n],double* b);
//(pointer to matrix A,pointer to b,n)

void mat_mul(int n,double [n][n],double [n][n]);

int main(){
	printf("-------------------------GAUSS JORDAN METHOD--------------------\n");
// Q1-------------------------------------------------------------------------	

	double A_q1[3][3];double b_q1[3];
	//A_q1 is matrix A of Q1 and b_q1 is b of Q1,similarly for Q2	
	FILE *file1;FILE *file2;
	file1=fopen("A_q1.txt","r");
	file2=fopen("b_q1.txt","r");
	//reading file
	for(int i=0;i<3;i++){
		fscanf(file2,"%lf",&b_q1[i]);
		for(int j=0;j<3;j++)
		fscanf(file1,"%lf",&A_q1[i][j]);
	}
	fclose(file1);fclose(file2);
	//augmented matrix(initial)
	printf("\nInitial augmented matrix for Q1\n");
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++)
		printf("%lf  ",A_q1[i][j]);
		printf("%lf\n\n",b_q1[i]);
	}
	
	gauss_jordan(3,A_q1,b_q1);
	
	//augmented matrix reduced
	printf("\nreduced augmented matrix for Q1\n");
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++)
		printf("%lf  ",A_q1[i][j]);
		printf("%lf\n\n",b_q1[i]);
	}

//Q2-------------------------------------------------------------------------

	double A_q2[3][3];double b_q2[3];
	file1=fopen("A_q2.txt","r");
	file2=fopen("b_q2.txt","r");
	//reading file
	for(int i=0;i<3;i++){
		fscanf(file2,"%lf",&b_q2[i]);
		for(int j=0;j<3;j++)
		fscanf(file1,"%lf",&A_q2[i][j]);
	}
	fclose(file1);fclose(file2);
	//augmented matrix(initial)
	printf("\nInitial augmented matrix for Q2\n");
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++)
		printf("%lf  ",A_q2[i][j]);
		printf("%lf\n\n",b_q2[i]);
	}
	
	gauss_jordan(3,A_q2,b_q2);
	
	//augmented matrix reduced
	printf("\nreduced augmented matrix for Q2\n");
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++)
		printf("%lf  ",A_q2[i][j]);
		printf("%lf\n\n",b_q2[i]);
	}

//Q3---------------------------------------------------------------------------
	double A_q3[3][3];
	double b1_q3[3]={1,0,0};double b2_q3[3]={0,1,0};double b3_q3[3]={0,0,1};
	double A_inv_q3[3][3];//will store A inverse

	file1=fopen("A_q3.txt","r");
	for(int i=0;i<3;i++)
	for(int j=0;j<3;j++)
	fscanf(file1,"%lf",&A_q3[i][j]);
	fclose(file1);
	gauss_jordan(3,A_q3,b1_q3);
	file1=fopen("A_q3.txt","r");
	for(int i=0;i<3;i++)
	for(int j=0;j<3;j++)
	fscanf(file1,"%lf",&A_q3[i][j]);
	fclose(file1);
	gauss_jordan(3,A_q3,b2_q3);
	file1=fopen("A_q3.txt","r");
	for(int i=0;i<3;i++)
	for(int j=0;j<3;j++)
	fscanf(file1,"%lf",&A_q3[i][j]);
	fclose(file1);
	gauss_jordan(3,A_q3,b3_q3);
	file1=fopen("A_q3.txt","r");
	for(int i=0;i<3;i++)
	for(int j=0;j<3;j++)
	fscanf(file1,"%lf",&A_q3[i][j]);
	fclose(file1);

	
	printf("\n\nfinding invers\n\n");

	//make A inverse from b1,b2,b3
	for(int i=0;i<3;i++)
	for(int j=0;j<3;j++){
		if(i==0)
		A_inv_q3[j][i]=b1_q3[j];
		if(i==1)
		A_inv_q3[j][i]=b2_q3[j];
		if(i==2)
		A_inv_q3[j][i]=b3_q3[j];
	}

	//print
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++)
		printf("%lf  ",A_q3[i][j]);
		printf("|  ");
		for(int k=0;k<3;k++)
		printf("%lf  ",A_inv_q3[i][k]);
		printf("\n\n");
	}

	//check
	printf("\n check with matrix multiplication\n");
	mat_mul(3,A_q3,A_inv_q3);

	return 0;
}

//-----------------------------------------------------------------------

void partial_pivot(int n,double arr[n][n],double* b,int r){
	if(arr[r][r]==0){
		int r1=r+1;
		while((arr[r][r]==0 && r1<n)){
			if(abs(arr[r1][r])>abs(arr[r][r])){
				double temp;
				for(int c=0;c<n;c++){
					temp=arr[r1][c];
					arr[r1][c]=arr[r][c];
					arr[r][c]=temp;
				}
				temp=b[r1];
				b[r1]=b[r];
				b[r]=temp;
			}
			r1++;
		}
	}



}
//-----------------------------------------------------------------------

void gauss_jordan(int n,double arr[n][n],double* b){
	for(int r=0;r<n;r++){
		partial_pivot(n,arr,b,r);
		double pivot=arr[r][r];
		for(int c=r;c<n;c++)
		arr[r][c]=arr[r][c]/pivot;
		b[r]=b[r]/pivot;

		for(int r1=0;r1<n;r1++){
			if(r1==r || arr[r1][r]==0)continue;
			else{
				double factor=arr[r1][r];
				for(int c=r;c<n;c++)
				arr[r1][c]=arr[r1][c]-factor*arr[r][c];
				b[r1]=b[r1]-factor*b[r];
													
			}
		}

	}
}
//--------------------------------------------------------------------

void mat_mul(int n,double arr1[n][n],double arr2[n][n]){

	double mul[n][n];
	for(int i=0;i<n;i++)
	for(int j=0;j<n;j++)mul[i][j]=0;

	for(int i=0;i<n;i++)
	for(int j=0;j<n;j++)
	for(int k=0;k<n;k++)
	mul[i][j]+=arr1[i][k]*arr2[k][j];

	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++)
		printf("%lf  ",mul[i][j]);
		printf("\n\n");
	}
}


