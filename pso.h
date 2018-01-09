
typedef double (*fptr)(double*,int);
typedef double(*fptr_1d)(double);
fptr functionLookup(char* fName);

void update_best(double * g, double * x, int nd);

void initialize_pso_domain(char* bounds_file, int* type, double* lb, double* ub ,double * gbest, int nd);

int lower_bound_check(double g, double l, double u, double epsilon);
int upper_bound_check(double g, double l, double u, double epsilon);



void adjust_domain(int* type, double* lb, double* ub, double * gbest,int nd, double alpha);

double fminbr(double a, double b, fptr_1d f, double tol);


void pso_optimization(fptr fun, int np, int nd, int ni, double* lb, double* ub,
                      double* value, double* gbest, int ni_stop, double tol);
