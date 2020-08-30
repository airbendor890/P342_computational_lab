#include<stdio.h>
#include<stdlib.h>

void partial_pivot(int n,double arr[n][n],double* b,int r);
void gauss_jordan(int n,double arr[n][n],double* b);

int main(){
	double arr[4][4];
	double b[4];
	FILE *file;
	file=fopen("test_matrix.txt","r");
	for(int i=0;i<4;i++)
	for(int j=0;j<4;j++){
		fscanf(file,"%lf",&arr[i][j]);
	}

	fclose(file);

	file=fopen("test_vector.txt","r");
	for(int i=0;i<4;i++)
	fscanf(file,"%lf",&b[i]);

	gauss_jordan(4,arr,b);

	for(int i=0;i<4;i++){
		for(int j=0;j<4;j++){
			printf("%lf ",arr[i][j]);
		}
		printf("%lf\n\n",b[i]);
	}

	return 0;

}


//-----------------------------------------------------------------------

void partial_pivot(int n,double arr[n][n],double* b,int r){
	if(arr[r][r]==0){
		int r1=r+1;
		while((arr[r][r]==0 && r1<n))
		//for(int r1=r+1;r1<n;r1++)				
					{
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

