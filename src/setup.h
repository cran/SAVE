float *vector();
double *dvector();
unsigned long int *ulivector();
int *ivector();
int ***iarray3();
float ***array3();
double ***darray3();
int **iarray2();
float **array2();
double **darray2(int n1, int n2);
void free_vector();
void free_dvector();
void free_array2(float **a, int n1, int n2);
void free_darray2(double **a, int n1, int n2);
void free_ivector();
void nrerror();