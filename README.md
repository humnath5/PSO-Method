# PSO-Method
This optimization software can be implemented for non-linear least square fitting of analytic potential energy function. 


### How to compile potential fitting without ghost atom?
g++  poten_noghost.c  pso.c  matmul2.c  randomlib.c  -o  poten_noghost

### How to compile potential fitting with ghost atom?
g++  poten_withghost.c  pso.c  matmul2.c  randomlib.c -o  poten_withghost


### How to run potential fitting without ghost atom ?

./poten_noghost  poten_fitting  50  16  500  10  0

./poten_noghost function np nd ni maxRun option

### How to run potential fitting with ghost atom ?


./poten_withghost  poten_fitting  50  25  500  10  0

./ poten_withghost function np nd ni maxRun option

