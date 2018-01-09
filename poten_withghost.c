#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "randomlib.h"
#include "matmul2.h"
#include "pso.h"
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))


//#define MAX  5000
/* global variable */

#define numofgeom 282 // change this number if number of data points are different


/* Global Variables*/
double V[numofgeom],rH1I[numofgeom],rH2I[numofgeom],rOI[numofgeom],rGI[numofgeom];
double xH1[numofgeom],xH2[numofgeom],xO[numofgeom],xI[numofgeom],xG[numofgeom];
double yH1[numofgeom],yH2[numofgeom],yO[numofgeom],yI[numofgeom],yG[numofgeom];
double zH1[numofgeom],zH2[numofgeom],zO[numofgeom],zI[numofgeom],zG[numofgeom];
double approxV[numofgeom];
//struct matrix* data; /*if data is needed in application*/
double * gbest;
int nd;
int ig;
fptr fun ;
fptr_1d f;

double fun_1D(double);
double poten_fitting(double , int);
double poten_energy(double *, double, double, double, double);
void read_coord_data(void);
void local_optimizer(int, int , double*, double*, double*, int , double);


int  main (int argc, char *argv[]){
int i,j,k;
int seed1,seed2;
srand((unsigned)time(NULL));
seed1= rand()%30000;
seed2= rand()%30000;
//seed1= 17795;
//seed2= 1795;

//seed1 = 17373;
//seed2 = 17496;
RandomInitialise(seed1,seed2);
/* Reading command line inputs */
char func[80];
strcpy(func, argv[1]);
fun  = functionLookup(func);
f = &fun_1D;

int np     =  atoi(argv[2]);
   nd      =  atoi(argv[3]);
int ni     = atoi(argv[4]);
int maxRun = atoi(argv[5]);
int option = atoi(argv[6]);

read_coord_data();


/* If inter-atomic  distances are provided, we can use the following.
char file[20];
strcpy(file,"datafile.txt");
data = read_matrix(file,row,col);
*/

/* Alocationg memory spaces */
int* type = (int*)malloc(sizeof(int)*nd);
double* lb = (double*)malloc(sizeof(double)*nd);
double* ub = (double*)malloc(sizeof(double)*nd);
gbest = (double*)malloc(sizeof(double)*nd);/*gbest is global variable */

char bounds_file[40];
strcpy(bounds_file,"bounds_withghost.txt");
initialize_pso_domain(bounds_file,type,lb,ub,gbest,nd);
/*for( i = 0; i < nd; i ++){

   printf("%d\t %f\t %f\t %f\n",type[i],lb[i], ub[i], gbest[i]);
}
*/
double start,finish;
double tol = 0.000001;
double value = 9999999;
int ni_stop = 100;
int ntries = 10;// number of times the FORWARD and BACKWARD is repeated in local optimizer
/*Start Optimization */

FILE* fp;
fp=fopen("output_withghost.txt","w");

start=time(NULL);
if(option==0){
    fprintf(fp,"\n**================== Standalone PSO Is Used ==================**\n");
        pso_optimization(fun,np,nd,ni,lb,ub,&value,gbest,ni_stop,tol );
        }


else if(option==1){

         fprintf(fp,"\n**================== Multi-Run PSO Is Used ==================**\n");
          double alpha = 1;

         fprintf(fp,"\n Run\tRMSE\n");
          for ( k = 0 ; k < maxRun; k++){

          pso_optimization(fun, np , nd, ni, lb, ub, &value, gbest, ni_stop,tol);
          fprintf(fp, "%d\t%lf\n", k+1, value);
          adjust_domain(type, lb , ub, gbest,nd,alpha);
          alpha = alpha/1.1;
    }
}


 else if(option==2){

         fprintf(fp,"\n**===== Multi-Run PSO with Local Optimizer Is Used ========**\n");
          double alpha = 1;
          double value1, value2 ;
         fprintf(fp,"\n Run\t RMSE1\t\t RMSE2\n");
          for ( k = 0 ; k < maxRun; k++){

          pso_optimization(fun, np , nd, ni, lb, ub, &value, gbest, ni_stop,tol);

          value1 = value;

          local_optimizer(nd, ntries, lb, ub, &value, ni_stop, tol);

          value2 = value;

          fprintf(fp, "%d\t%lf\t%lf\n", k+1,value1,value2);
          adjust_domain(type, lb , ub, gbest,nd,alpha);
          alpha = alpha/1.1;


    }


}

else
    {fprintf(fp,"\n No option is found !!!\n");}

 finish=time(NULL);

double time=difftime(finish, start);

//double approxV[MAX];
for(i= 0;i< numofgeom;i++){
approxV[i]= poten_energy(gbest,rH1I[i],rH2I[i],rOI[i],rGI[i]);
}
//fprintf(fp,"%lf\t%lf\n",V[i],approxV[i]);}



gbest[7]=round(gbest[6]+gbest[7]),  gbest[6]=  round(gbest[6]);
gbest[15]=round(gbest[14]+gbest[15]); gbest[14]= round(gbest[14]);
gbest[23]=round(gbest[22]+gbest[23]);gbest[22]= round(gbest[22]);



/* Write outputs */

fprintf(fp,"\n**================== PSO Parameters Used ==================**\n");
fprintf(fp,"seed1:%d\n",seed1);
fprintf(fp,"seed2:%d\n",seed2);
fprintf(fp,"dimension : %d\n",nd);
fprintf(fp,"num_particles: %d\n",np);
fprintf(fp,"max_num_iterations per run allowed: %d\n",ni);
fprintf(fp,"maxRun for multi-run PSO: %d\n",maxRun);
fprintf(fp,"total wall-clock time: %lf\n",time);
fprintf(fp,"best rmse : %lf\n",value);

fprintf(fp,"\n**================= Optimal Model Parameters =================**\n");
for(j=0;j<nd;j++){fprintf(fp,"%lf\n",gbest[j]);}

fprintf(fp,"\n**================= DFT Energy and Fitted Energy =================**\n");
//double approxV[MAX];
for(i= 0;i< numofgeom;i++){

 fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",rH1I[i],rH2I[i],rOI[i],rGI[i],V[i],approxV[i]);
}
fclose(fp);
free(lb);free(ub);free(gbest);free(type);
//free_matrix(data);

return 0;
}


void read_coord_data(void) {
 int i;
 FILE* input,*inputV;

    input  = fopen("geom.txt","r");
    inputV = fopen("energy.txt","r");
    char tempblank[1];

    for( i=0;i<numofgeom;i++)
    {
        fscanf(input,"%lf\t %lf\t %lf",&xO[i],&yO[i],&zO[i]);
        fscanf(input,"%lf\t%lf\t%lf",&xH1[i],&yH1[i],&zH1[i]);
        fscanf(input,"%lf\t%lf\t%lf",&xH2[i],&yH2[i],&zH2[i]);
        fscanf(input,"%lf\t%lf\t%lf",&xI[i],&yI[i],&zI[i]);
        fscanf(input,"%c",tempblank);
        fscanf(inputV,"%lf",&V[i]);


  /* calculating distances: */
        xG[i] = 0.00000; yG[i]= 0.00000; zG[i]= 0.00000;
        rH1I[i]= sqrt(pow((xH1[i]-xI[i]),2) +pow((yH1[i]-yI[i]),2)+pow((zH1[i]-zI[i]),2));
        rH2I[i]= sqrt(pow((xH2[i]-xI[i]),2)+pow((yH2[i]-yI[i]),2)+pow((zH2[i]-zI[i]),2));
        rOI[i] = sqrt(pow((xO[i] - xI[i]),2) +pow((yO[i]-yI[i]),2)+pow((zO[i]-zI[i]),2));
        rGI[i]= 0.0;

    }
  fclose(input);
  fclose(inputV);


}

double poten_fitting(double * x, int nd){/* nd is not necessary in this function but kept to be consistent for PSO input */
  int i;
     double rmse = 0;
     double sse = 0;

     for( i= 0; i< numofgeom; i++) {

          xG[i] =0.00000; yG[i]  = x[24]; zG[i] = 0.0000; // x[24] is alpha parameter

          rGI[i] = sqrt(pow((xG[i] - xI[i]),2) + pow((yG[i]-yI[i]),2)+pow((zG[i]-zI[i]),2));

          approxV[i] = poten_energy(x,rH1I[i],rH2I[i],rOI[i],rGI[i]);
          sse = sse + pow(approxV[i] - V[i],2);
       }

      rmse = sqrt(sse/numofgeom);
      return rmse;

}



double poten_energy(double * x, double r1, double r2, double r3, double r4){

     /* x -> parameter vector to be optimized
        r -> vector of inter-atomic distance */

    double  m1,m2,m3; /* integer parameters */
    double  n1,n2,n3; /* integer parameters */

    double pe = 0;
  //Calculating Nitrogen-Carbon potential energy

    m1  = round(x[6]);    n1 = round(x[6] + x[7]);

    m2  = round(x[14]);   n2 = round(x[14] + x[15]);

    m3  = round(x[22]);   n3 = round(x[22] + x[23]);


									// I-H1
            pe = pe + x[0] * exp(-x[1]*r1) + x[2]/ pow((r1+ x[4]), m1) + x[3]/ pow((r1+x[5]),n1);
                                    // I-H2
            pe = pe + x[0] * exp(-x[1]*r2) + x[2]/ pow((r2+ x[4]), m1) + x[3]/ pow((r2+x[5]),n1);
                                    // I-O
            pe = pe + x[8] * exp(-x[9]*r3) + x[10]/ pow((r3+ x[12]), m2) + x[11]/ pow((r3+x[13]),n2);
                                   // I-G
            pe = pe + x[16] * exp(-x[17]*r4) + x[18]/ pow((r4+ x[20]), m3) + x[19]/ pow((r4+x[21]),n3);

    return pe;

}



double fun_1D(double t) {
	double sum;
	gbest[ig] = t;
	sum = fun(gbest, nd);
	return sum;
}


void  local_optimizer(int nd, int ntries, double* lb, double* ub, double* value, int ni_stop, double tol){
		 /* this is extern global variable*/
		int  k;
		double l, u;
		double temp_t, old_value;
		double eps;

		double best = fun(gbest, nd);
		double alpha = 0.1;

		old_value = best;
		double val;
		//printf("\nbefore local optimizer %f\n",val_old);
		for (k = 0; k< ntries; k++) {

			for (ig = 0; ig < nd; ig ++) {

				eps = alpha*(ub[ig] - lb[ig ]);

				l = MAX(gbest[ig ] - eps, lb[ig]);

				u = MIN(gbest[ig] + eps, ub[ig]);

				temp_t = gbest[ig];

				gbest[ig] = fminbr(l, u, f, tol);

				val = fun(gbest, nd);

				if (val > best) { gbest[ig] = temp_t; }
				else { best = val; }


			}

			for ( ig= nd - 1; ig >= 0; ig--) {

                eps = alpha*(ub[ig] - lb[ig ]);

				l = MAX(gbest[ig ] - eps, lb[ig]);

				u = MIN(gbest[ig] + eps, ub[ig]);

				temp_t = gbest[ig];

				gbest[ig] = fminbr(l, u, f, tol);

				val = fun(gbest, nd);

				if (val > best) { gbest[ig] = temp_t; }
				else { best = val; }



			}

			for ( ig= nd - 1; ig >= 0; ig--) {

                eps = alpha*(ub[ig] - lb[ig ]);

				l = MAX(gbest[ig ] - eps, lb[ig]);

				u = MIN(gbest[ig] + eps, ub[ig]);

				temp_t = gbest[ig];

				gbest[ig] = fminbr(l, u, f, tol);

				val = fun(gbest, nd);

				if (val > best) { gbest[ig] = temp_t; }
				else { best = val; }



			}

			for (ig = 0; ig < nd; ig ++) {

				eps = alpha*(ub[ig] - lb[ig ]);

				l = MAX(gbest[ig ] - eps, lb[ig]);

				u = MIN(gbest[ig] + eps, ub[ig]);

				temp_t = gbest[ig];

				gbest[ig] = fminbr(l, u, f, tol);

				val = fun(gbest, nd);

				if (val > best) { gbest[ig] = temp_t; }
				else { best = val; }


			}
			best = fun(gbest, nd);
			// printf("\n val_old  %f\t  val_best%f\n",val_old,val_best);
			if ((k > 0) && (k % ni_stop == 0)) {
				double estimate = fabs((best - old_value) / (best + tol));
				if (estimate < tol) break;
			}
			old_value = best;
		}
		best = fun(gbest, nd);

		(*value) = best;
	}



fptr functionLookup(char* fName){
/*if (strcmp(fName,"rosenbrock")==0){return &rosenbrock;}*/
if (strcmp(fName,"poten_fitting")==0){return &poten_fitting;}
else {printf("Function Name %s not recognized\n",fName);return NULL;}
}

