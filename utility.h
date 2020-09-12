void mat_mul(int n,double arr1[n][n],double arr2[n][n]);
double lu_dcmpsn(int n,double a[n][n],int * P);
void lu_partial_pivot(int n,double arr[n][n],int * P,int r);
void forward_back_sub(int n,double a[n][n],int *P,double *b);
